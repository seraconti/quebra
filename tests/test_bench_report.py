"""The bench report's two dependence numbers: the Arm E grid's column and the Arm C finding.

Both must read 0 with no dependence and mean the same at every event count. The grid prints
the generator's own rank correlation; the Arm C finding adds back `r_1`'s exchangeable mean.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from jobs.bench import report
from jobs.bench.arms import ARM_C, ARM_E


@pytest.mark.unit
def test_the_arm_e_grid_prints_the_generators_rank_correlation():
    """Oracle: `(6/pi) arcsin(rho/2)`, the rank correlation of neighbouring durations Arm E's
    generator sets, rounded as the table rounds (3 decimals). The measured
    `mean_induced_lag1` is set to -0.05 here, which a column reading it would show.
    """
    power = pd.DataFrame(
        [
            {
                "arm": ARM_E,
                "check": check,
                "calibration": "permutation",
                "variant": variant,
                "rho": rho,
                "n_target": n,
                "rejection_rate": 0.1,
                "mean_induced_lag1": -0.05,
            }
            for check, variant in (
                ("c5_rank_autocorr", "studentized"),
                ("c6_exchangeability", np.nan),
            )
            for rho in (0.1, 0.5)
            for n in (20, 100)
        ]
    )
    table, _strong = report.dependence_grid_table(report._label_column(power))
    expected = [round(6.0 / np.pi * np.arcsin(rho / 2.0), 3) for rho in table["rho"]]
    assert list(table["neighbour_rank_corr"]) == expected


@pytest.mark.unit
def test_the_arm_c_finding_adds_back_one_over_n_per_cell_and_refuses_what_it_cannot_centre():
    """Oracle: on one segment the exchangeable mean of `r_1` is `-1/N` (Dufour and Roy 1985,
    Eq. (2.2), p. 260).

    Arm C cells measuring `-1/N + d`, `rho = 0` from the size table (where it is blank) and
    `rho > 0` from the power table, two rows each, must come back with `centred_lag1 = d`. A
    cell of two segments, a cell whose rows disagree and a cell whose every replicate failed
    must raise.
    """

    def cell(rho, n_events, d):
        return [
            {
                "arm": ARM_C,
                "clock": "in_spec",
                "n_target": int(n_events),
                "rho": rho,
                "check": check,
                "mean_induced_lag1": -1.0 / n_events + d,
                "mean_n_events": n_events,
                "mean_n_segments": 1.0,
                "fail_reason": "",
            }
            for check in ("c5_rank_autocorr", "c6_exchangeability")
        ]

    size = pd.DataFrame(cell(np.nan, 20.0, 0.01) + cell(np.nan, 400.0, -0.002))
    power = pd.DataFrame(cell(0.5, 20.0, 0.03) + cell(0.5, 400.0, 0.0))
    centred = report.centred_arm_c_cells(size, power).set_index(["n_target", "rho"])
    assert centred["centred_lag1"].to_dict() == pytest.approx(
        {(20, 0.0): 0.01, (400, 0.0): -0.002, (20, 0.5): 0.03, (400, 0.5): 0.0},
        abs=1e-15,
    )

    with pytest.raises(ValueError, match="more than one segment"):
        report.centred_arm_c_cells(size, power.assign(mean_n_segments=2.0))
    split = power.copy()
    split.loc[0, "mean_induced_lag1"] += 0.01
    with pytest.raises(ValueError, match="disagree"):
        report.centred_arm_c_cells(size, split)
    failed = power.head(1).assign(
        check=None, mean_induced_lag1=np.nan, fail_reason="ValueError: probe"
    )
    with pytest.raises(ValueError, match="failed every replicate"):
        report.centred_arm_c_cells(size, pd.concat([power.iloc[2:], failed]))
