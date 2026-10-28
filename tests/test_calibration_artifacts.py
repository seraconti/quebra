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
    ARM_D,
    ARM_E,
    bonferroni_z_crit,
    null_se,
    power_vs_dependence,
    power_vs_trend,
    read_dependence,
    size_vs_n,
)

RESULTS = Path(__file__).resolve().parents[1] / "jobs" / "bench" / "results"


def _arm_e_power_table(*, rhos, ns, lag1):
    """Arm E power rows for C5 (both variants) and C6 at two Weibull shapes, with a chosen
    `mean_induced_lag1`. Shape 1.5 rejects 0.02 more often, so pooling shows as spread."""
    keys = [
        ("c5_rank_autocorr", "permutation", c5.VARIANT_STUDENTIZED),
        ("c5_rank_autocorr", "permutation", c5.VARIANT_RAW),
        ("c6_exchangeability", "permutation", np.nan),
    ]
    return pd.DataFrame(
        [
            {
                "arm": ARM_E,
                "clock": "in_spec",
                "n_target": n,
                "shape": shape,
                "quantised": False,
                "censoring_target": 0.0,
                "rho": rho,
                "check": check,
                "calibration": calibration,
                "variant": variant,
                "rejection_rate": 0.05 + rho + (0.02 if shape == 1.5 else 0.0),
                "mean_induced_lag1": lag1,
            }
            for n in ns
            for shape in (0.75, 1.5)
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


@pytest.mark.unit
def test_power_vs_dependence_names_its_pooling_and_refuses_what_it_cannot_draw():
    """Contract: a power value is a mean over shapes, so the artifact names the pooling,
    carries its spread and refuses a point that pools fewer cells than that. Another arm, a
    curve missing a `rho`, a drawn cell whose every replicate failed and an unknown check
    raise; a failed cell the artifact does not draw is ignored.
    """
    table = _arm_e_power_table(rhos=(0.0, 0.3), ns=(20, 1000), lag1=-0.05)
    artifact = power_vs_dependence(table)
    row = "c6_exchangeability [permutation]"
    assert artifact.pooled_over == ["shape"]
    assert artifact.power_lo_by_row_and_n[row][20] == pytest.approx([0.05, 0.35])
    assert artifact.power_hi_by_row_and_n[row][20] == pytest.approx([0.07, 0.37])
    assert artifact.power_by_row_and_n[row][20] == pytest.approx([0.06, 0.36])

    with pytest.raises(ValueError, match="is not Arm E"):
        power_vs_dependence(table, arm="D_trp_power_law")
    gap = (table["check"] == "c6_exchangeability") & (table["n_target"] == 20)
    with pytest.raises(ValueError, match="covers rho"):
        power_vs_dependence(table[~(gap & (table["rho"] == 0.3))])

    # What `run_cell` writes when every replicate of a cell fails: the cell's rows replaced
    # by one row naming no check.
    drawn = (table["n_target"] == 20) & (table["shape"] == 1.5) & (table["rho"] == 0.3)
    failed = (
        table[drawn]
        .head(1)
        .assign(
            check=None,
            calibration=None,
            variant=None,
            rejection_rate=np.nan,
            fail_reason="ValueError: probe",
        )
    )
    with pytest.raises(ValueError, match="failed every replicate"):
        power_vs_dependence(pd.concat([table[~drawn], failed], ignore_index=True))
    elsewhere = pd.concat(
        [table, failed.assign(arm="D_trp_power_law")], ignore_index=True
    )
    assert (
        power_vs_dependence(elsewhere).power_by_row_and_n == artifact.power_by_row_and_n
    )
    with pytest.raises(ValueError, match="pools 1 of 2"):
        power_vs_dependence(table[~(drawn & (table["check"] == "c6_exchangeability"))])

    # Two unknown keys that differ only where one is missing must still raise a ValueError.
    strangers = table.head(2).assign(
        check="c9_unknown", calibration=["permutation", None]
    )
    with pytest.raises(ValueError, match="BENCH_ROW_KEYS"):
        power_vs_dependence(pd.concat([table, strangers], ignore_index=True))


@pytest.mark.statistical
def test_the_measured_lag1_plus_one_over_n_matches_the_copula_rank_correlation_at_n_1000():
    """Oracle: `(6/pi) arcsin(rho/2)`, the rank correlation Arm E's generator induces.

    On one segment `mean_induced_lag1` is C5's uncentred `r_1`, whose exchangeable mean is
    `-1/N`. With `1/mean_n_events` added back, every Arm E power cell at the largest event
    count must sit within 3 Monte Carlo SE of the identity. The SE uses `r_1`'s null SD,
    `1/sqrt(N)`, over `n_used` replicates; an AR(1) dependence shrinks that SD rather than
    inflating it, so the band is if anything wide.

    The 3 SE is per cell and not corrected for the 12 cells it reads (two shapes, six rho).
    The cells also sit a little low on average: beyond its `-1/N` offset, `r_1` is biased
    toward 0 under dependence by an amount that grows with `rho` and shrinks with N. Both
    are tolerated here, not modelled. The test reads the committed table and calls no source
    function: it is evidence that the identity describes the generator, and the guard on
    the code that draws it is the axis test above.
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
    studentized row must still be there, and an unknown check in a drawn cell must raise.
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

    stranger = power[power["arm"] == ARM_E].head(1).assign(check="c9_unknown")
    with pytest.raises(ValueError, match="BENCH_ROW_KEYS"):
        power_vs_dependence(pd.concat([power, stranger]))


@pytest.mark.unit
def test_size_vs_n_refuses_a_point_that_pools_fewer_cells_than_it_says():
    """Contract: each size point is a mean over the two Weibull shapes. A drawn cell whose
    every replicate failed, or a row missing from one shape, must raise rather than leave a
    one-shape value labelled as the mean over both.
    """
    size = pd.read_csv(RESULTS / "size_table.csv")
    drawn = (
        (size["arm"] == "A_iid_weibull")
        & (~size["quantised"].astype(bool))
        & (size["censoring_target"] == 0.0)
        & (size["n_target"] == 20)
        & (size["shape"] == 1.5)
    )
    assert size_vs_n(size).pooled_over == ["shape"]
    failed = (
        size[drawn]
        .head(1)
        .assign(
            check=None, calibration=None, variant=None, fail_reason="ValueError: probe"
        )
    )
    with pytest.raises(ValueError, match="failed every replicate"):
        size_vs_n(pd.concat([size[~drawn], failed], ignore_index=True))
    with pytest.raises(ValueError, match="pools 1 of 2"):
        size_vs_n(size[~(drawn & (size["check"] == "c6_exchangeability"))])


@pytest.mark.unit
def test_read_dependence_adds_back_one_over_n_per_cell_and_refuses_several_segments():
    """Oracle: on one segment the exchangeable mean of `r_1` is `-1/N` (Dufour and Roy 1985,
    Eq. (2.2), p. 260).

    Three cells per `rho`, each written as two rows, measure `-1/N + d` with `d` = +0.03, 0
    and -0.03. Centred once per cell they average 0, with an SE across the three cells of
    `0.03/sqrt(3)`; counting each cell once per row would shrink that SE. A cell of two
    segments, a cell whose rows disagree and a cell whose every replicate failed must raise.
    """
    table = pd.DataFrame(
        [
            {
                "arm": "C_ar1_reads_carved",
                "clock": "in_spec",
                "n_target": int(n_events),
                "rho": rho,
                "check": check,
                "mean_induced_lag1": -1.0 / n_events + offset,
                "mean_n_events": n_events,
                "mean_n_segments": 1.0,
            }
            for n_events, offset in ((20.0, 0.03), (50.0, 0.0), (400.0, -0.03))
            for rho in (0.1, 0.5)
            for check in ("c5_rank_autocorr", "c6_exchangeability")
        ]
    )
    artifact = read_dependence(table)
    assert artifact.read_rho == [0.1, 0.5]
    assert artifact.centred_lag1 == pytest.approx([0.0, 0.0], abs=1e-15)
    assert artifact.centred_lag1_se == pytest.approx(
        [0.03 / np.sqrt(3.0)] * 2, rel=1e-12
    )
    with pytest.raises(ValueError, match="more than one segment"):
        read_dependence(table.assign(mean_n_segments=2.0))
    split = table.copy()
    split.loc[0, "mean_induced_lag1"] += 0.01
    with pytest.raises(ValueError, match="differs between rows"):
        read_dependence(split)
    first_cell = (table["n_target"] == 20) & (table["rho"] == 0.1)
    failed = (
        table[first_cell]
        .head(1)
        .assign(check=None, mean_induced_lag1=np.nan, fail_reason="ValueError: probe")
    )
    with pytest.raises(ValueError, match="failed every replicate"):
        read_dependence(pd.concat([table[~first_cell], failed], ignore_index=True))


@pytest.mark.unit
def test_power_vs_trend_reads_b_and_takes_b_1_from_the_uncensored_size_cells():
    """Contract: the x-axis is Arm D's trend exponent `b`. Its `b = 1` point, no trend, is
    Arm D's uncensored size cells; a censored size cell (rate 0.9 here) must not enter it,
    and without an uncensored `b = 1` cell the artifact raises.
    """

    def cells(b, censoring, rate):
        return [
            {
                "arm": ARM_D,
                "clock": "in_spec",
                "n_target": n,
                "shape": shape,
                "quantised": False,
                "censoring_target": censoring,
                "b": b,
                "rho": np.nan,
                "check": check,
                "calibration": "permutation",
                "variant": np.nan,
                "rejection_rate": rate,
            }
            for n in (20, 100)
            for shape in (0.75, 1.5)
            for check in ("c1_lewis_robinson", "c6_exchangeability")
        ]

    power = pd.DataFrame(cells(0.7, 0.0, 0.4) + cells(1.5, 0.0, 0.6))
    size = pd.DataFrame(cells(1.0, 0.0, 0.05) + cells(1.0, 0.25, 0.9))
    artifact = power_vs_trend(power, size)
    assert artifact.b == [0.7, 1.0, 1.5]
    assert artifact.pooled_over == ["shape"]
    assert artifact.power_by_row_and_n["c1_lewis_robinson [permutation]"][20] == (
        pytest.approx([0.4, 0.05, 0.6])
    )
    with pytest.raises(ValueError, match="b = 1"):
        power_vs_trend(power, size[size["censoring_target"] > 0.0])


@pytest.mark.statistical
def test_every_drawn_check_is_at_its_size_at_b_1_on_the_committed_tables():
    """Oracle: a permutation-calibrated check is exact, so at `b = 1` (no trend) each drawn
    row's rejection rate is binomial around alpha.

    Each point is the mean of the two shape cells, so its null SE is
    `sqrt(se_1^2 + se_2^2) / 2` from their `n_used`. The band is the two-sided Bonferroni
    threshold over every drawn (row, n) point. The test reads the committed tables: it is
    evidence about them, and the guard on the code is the test above.
    """
    size = pd.read_csv(RESULTS / "size_table.csv")
    power = pd.read_csv(RESULTS / "power_table.csv")
    artifact = power_vs_trend(power, size)
    at_null = artifact.b.index(1.0)
    null = size[
        (size["arm"] == ARM_D) & (size["b"] == 1.0) & (size["censoring_target"] == 0.0)
    ].assign(
        row=lambda d: [
            f"{c} [{k}" + (f"/{v}]" if isinstance(v, str) and v else "]")
            for c, k, v in zip(d["check"], d["calibration"], d["variant"])
        ]
    )
    z = []
    for row in artifact.rows_drawn:
        for n, rates in artifact.power_by_row_and_n[row].items():
            used = null[(null["row"] == row) & (null["n_target"] == n)]["n_used"]
            assert len(used) == 2
            se = np.sqrt(sum(null_se(artifact.alpha, u) ** 2 for u in used)) / 2.0
            z.append((rates[at_null] - artifact.alpha) / se)
    assert len(z) >= 2 * len(artifact.rows_drawn)
    assert max(abs(v) for v in z) <= bonferroni_z_crit(len(z)), [round(v, 2) for v in z]
