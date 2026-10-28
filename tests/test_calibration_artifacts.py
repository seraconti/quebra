"""The calibration artifacts: what their dependence axes measure and which rows they carry.

A dependence axis must read 0 with no dependence and mean the same at every event count, or
curves for different n are drawn against different offsets. C5's unstudentized row is a
bench-only comparison and must reach no artifact.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import quebra.analyzers.checks.c5_rank_autocorr as c5
from quebra.analyzers.calibration_summary import (
    ARM_E,
    power_vs_dependence,
    read_dependence,
    size_vs_n,
)

RESULTS = Path(__file__).resolve().parents[1] / "jobs" / "bench" / "results"


def _arm_e_power_table(*, rhos, ns, lag1):
    """Arm E power rows for C5 (both variants) and C6, with a chosen `mean_induced_lag1`."""
    keys = [
        ("c5_rank_autocorr", "permutation", c5.VARIANT_STUDENTIZED),
        ("c5_rank_autocorr", "permutation", c5.VARIANT_RAW),
        ("c6_exchangeability", "permutation", np.nan),
    ]
    return pd.DataFrame(
        [
            {
                "arm": ARM_E,
                "n_target": n,
                "rho": rho,
                "check": check,
                "calibration": calibration,
                "variant": variant,
                "rejection_rate": 0.05 + rho,
                "mean_induced_lag1": lag1,
            }
            for n in ns
            for rho in rhos
            for check, calibration, variant in keys
        ]
    )


@pytest.mark.unit
def test_the_power_axis_is_the_copula_rank_correlation_and_zero_without_dependence():
    """Oracle: Spearman's correlation under a Gaussian copula, `(6/pi) arcsin(rho/2)`.

    The table's measured lag-1 is set to -0.05 everywhere, which an axis reading it would
    show; the axis must instead be exactly 0 at `rho = 0` and one list for every n.
    """
    artifact = power_vs_dependence(
        _arm_e_power_table(rhos=(0.0, 0.3), ns=(20, 1000), lag1=-0.05)
    )
    assert artifact.rho == [0.0, 0.3]
    assert artifact.neighbour_rank_corr[0] == 0.0
    assert artifact.neighbour_rank_corr[1] == pytest.approx(
        6.0 / np.pi * np.arcsin(0.15), rel=1e-15
    )
    assert all(len(per_n) == 2 for per_n in artifact.power_by_row_and_n.values())


@pytest.mark.statistical
def test_the_measured_lag1_plus_one_over_n_matches_the_copula_rank_correlation_at_n_1000():
    """Oracle: `(6/pi) arcsin(rho/2)`, the rank correlation Arm E's generator induces.

    On one segment `mean_induced_lag1` is C5's uncentred `r_1`, whose exchangeable mean is
    `-1/N`. With `1/mean_n_events` added back, every Arm E power cell at the largest event
    count must sit within 3 Monte Carlo SE of the identity. The SE uses `r_1`'s null SD,
    `1/sqrt(N)`, over `n_used` replicates; an AR(1) dependence shrinks that SD rather than
    inflating it, so the band is if anything wide.
    """
    power = pd.read_csv(RESULTS / "power_table.csv")
    arm_e = power[power["arm"] == ARM_E]
    cells = arm_e[arm_e["n_target"] == arm_e["n_target"].max()].drop_duplicates(
        ["shape", "rho"]
    )
    assert len(cells) >= 2 and (cells["mean_n_segments"] == 1.0).all()
    centred = cells["mean_induced_lag1"] + 1.0 / cells["mean_n_events"]
    truth = 6.0 / np.pi * np.arcsin(cells["rho"] / 2.0)
    se = 1.0 / np.sqrt(cells["mean_n_events"] * cells["n_used"])
    z = (centred - truth) / se
    assert (z.abs() <= 3.0).all(), z.round(2).tolist()


@pytest.mark.unit
def test_no_calibration_artifact_carries_the_unstudentized_row():
    """Contract: C5's unstudentized row is in the bench tables only.

    The committed tables carry it, so a filter that let it through would show here; C5's
    studentized row must still be there, and a row in no known list must raise.
    """
    size = pd.read_csv(RESULTS / "size_table.csv")
    power = pd.read_csv(RESULTS / "power_table.csv")
    assert (size["variant"] == c5.VARIANT_RAW).any()
    assert (power["variant"] == c5.VARIANT_RAW).any()
    sizes, powers = size_vs_n(size), power_vs_dependence(power)
    labels = [
        *sizes.rate_by_row,
        *sizes.se_by_row,
        *sizes.rate_lo_by_row,
        *sizes.rate_hi_by_row,
        *powers.power_by_row_and_n,
        *powers.rows_drawn,
    ]
    assert not [label for label in labels if c5.VARIANT_RAW in label]
    assert any(label.endswith(f"/{c5.VARIANT_STUDENTIZED}]") for label in labels)

    stranger = power.head(1).assign(check="c9_unknown")
    with pytest.raises(ValueError, match="BENCH_ROW_KEYS"):
        power_vs_dependence(pd.concat([power, stranger]))


@pytest.mark.unit
def test_read_dependence_adds_back_one_over_n_per_cell_and_refuses_several_segments():
    """Oracle: on one segment the exchangeable mean of `r_1` is `-1/N` (Dufour and Roy 1985,
    Eq. (2.2), p. 260).

    Cells whose measured lag-1 is exactly `-1/N`, at three event counts and two rows each,
    must centre to 0 at every `rho`; a cell of two segments must raise.
    """
    table = pd.DataFrame(
        [
            {
                "arm": "C_ar1_reads_carved",
                "clock": "in_spec",
                "n_target": int(n_events),
                "rho": rho,
                "check": check,
                "mean_induced_lag1": -1.0 / n_events,
                "mean_n_events": n_events,
                "mean_n_segments": 1.0,
            }
            for n_events in (20.0, 50.0, 400.0)
            for rho in (0.1, 0.5)
            for check in ("c5_rank_autocorr", "c6_exchangeability")
        ]
    )
    artifact = read_dependence(table)
    assert artifact.read_rho == [0.1, 0.5]
    assert artifact.centred_lag1 == pytest.approx([0.0, 0.0], abs=1e-15)
    with pytest.raises(ValueError, match="more than one segment"):
        read_dependence(table.assign(mean_n_segments=2.0))
