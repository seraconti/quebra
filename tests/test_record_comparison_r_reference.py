"""The comparisons across records against R's `survival` and `metafor`.

Oracle, unless a test names another: R, run by `jobs/rscripts/comparison_reference.R` on the
inputs committed beside its values (versions in the fixture's meta rows). The suite never
invokes R. Each record is built as an event table and pushed through the shipped steps in
`analyzers/record_comparison.py`.
"""

from __future__ import annotations

import subprocess

import numpy as np
import pandas as pd
import pytest

from quebra.analyzers import record_comparison
from quebra.analyzers.event_table import EventTable
from quebra.analyzers.windows import SIDE_IN_SPEC
from tests.fixtures import R_COMPARISON_INPUTS, R_COMPARISON_VALUES

pytestmark = pytest.mark.statistical

_INPUTS = pd.read_csv(R_COMPARISON_INPUTS)
_VALUES = pd.read_csv(R_COMPARISON_VALUES)
CASES = sorted(_INPUTS["case"].unique())
PAIR_CASES = [c for c in CASES if c.startswith("pair_")]
POOL_CASES = [c for c in CASES if c.startswith("pool_")]
REL = 1e-9
RULE = record_comparison.TAU_RULE_MIN_LARGEST_AGE


def _tables(case: str) -> list[EventTable]:
    """One event table per record, in the order R wrote them."""
    rows = _INPUTS[_INPUTS["case"] == case]
    tables = []
    for record in pd.unique(rows["record"]):
        mine = rows[rows["record"] == record].sort_values("i")
        tables.append(
            EventTable(
                dataset_id=str(record),
                threshold_label="fixture",
                side=SIDE_IN_SPEC,
                age_s=mine["age"].to_numpy(dtype=float),
                event=mine["event"].to_numpy(dtype=int) == 1,
            )
        )
    return tables


def _r(case: str, quantity: str) -> np.ndarray:
    rows = _VALUES[(_VALUES["case"] == case) & (_VALUES["quantity"] == quantity)]
    if rows.empty:
        raise KeyError(f"no R value {quantity!r} for {case!r}")
    return rows.sort_values("i")["value"].to_numpy(dtype=float)


def _r1(case: str, quantity: str) -> float:
    (value,) = _r(case, quantity)
    return float(value)


def _version(package: str) -> tuple[int, ...]:
    return tuple(
        int(_r1("meta", f"{package}_{k}")) for k in ("major", "minor", "patch")
    )


def test_the_fixture_was_written_by_the_versions_the_docstrings_cite() -> None:
    """Oracle: the versions `analyzers/record_comparison.py` cites."""
    assert _version("survival") == (3, 8, 6)
    assert _version("metafor") == (5, 2, 1)


@pytest.mark.r
def test_the_fixture_versions_are_the_installed_r_s(requires_rscript: str) -> None:
    """Oracle: the installed R, asked directly. A fixture it cannot re-derive is stale."""
    for package in ("survival", "metafor"):
        out = subprocess.run(
            [requires_rscript, "-e", f'cat(unlist(packageVersion("{package}")))'],
            capture_output=True,
            text=True,
            timeout=120,
            check=True,
        ).stdout.split()
        assert tuple(int(v) for v in out) == _version(package), package


@pytest.mark.parametrize("case", CASES)
def test_the_comparison_tau_and_the_record_that_set_it_match_r(case: str) -> None:
    """Oracle: R's min(max(age)) over the case's records, and which.min for the setter."""
    tables = _tables(case)
    tau = record_comparison.comparison_tau(tables, tau_rule=RULE)
    assert tau.rmst_tau_s == pytest.approx(_r1(case, "tau"), rel=REL)
    setter = int(_r1(case, "tau_set_by_index")) - 1
    assert tau.tau_set_by == tables[setter].dataset_id
    assert tau.excluded == ()


@pytest.mark.parametrize("case", PAIR_CASES)
def test_the_log_rank_statistic_and_p_match_r_survdiff(case: str) -> None:
    """Oracle: R `survdiff(rho = 0)`; the sign of z is the sign of O - E for record a."""
    a, b = _tables(case)
    pair = record_comparison.compare_pair(a, b, tau_rule=RULE)
    assert pair.logrank_not_computed_reason is None
    assert pair.logrank_chi2 == pytest.approx(_r1(case, "logrank_chisq"), rel=REL)
    assert pair.logrank_p == pytest.approx(_r1(case, "logrank_p"), rel=REL)
    assert np.sign(pair.logrank_z) == np.sign(_r1(case, "logrank_o_minus_e_a"))


@pytest.mark.parametrize("case", PAIR_CASES)
def test_the_restricted_mean_difference_and_its_se_match_r(case: str) -> None:
    """Oracle: the variance of a difference of independent estimates, on R's per-record values.

    Each record's restricted mean and SE are R `summary(survfit, rmean = tau)` at the
    comparison tau; R combines them as a - b and sqrt(se_a^2 + se_b^2).
    """
    a, b = _tables(case)
    pair = record_comparison.compare_pair(a, b, tau_rule=RULE)
    assert pair.summary_a is not None and pair.summary_b is not None
    assert pair.summary_a.rmst_s == pytest.approx(_r1(case, "rmean_a"), rel=REL)
    assert pair.summary_a.rmst_se_s == pytest.approx(
        _r1(case, "rmean_se_a"), rel=REL, abs=1e-12
    )
    assert pair.summary_b.rmst_s == pytest.approx(_r1(case, "rmean_b"), rel=REL)
    assert pair.summary_b.rmst_se_s == pytest.approx(
        _r1(case, "rmean_se_b"), rel=REL, abs=1e-12
    )
    assert pair.rmst_difference_s == pytest.approx(
        _r1(case, "rmean_diff"), rel=REL, abs=1e-12
    )
    assert pair.rmst_difference_se_s == pytest.approx(
        _r1(case, "rmean_diff_se"), rel=REL
    )


def _summary(case: str) -> record_comparison.RandomEffectsSummary:
    summary = record_comparison.random_effects_summary(_tables(case), tau_rule=RULE)
    assert summary.not_pooled_reason is None
    assert all(row.pooled for row in summary.rows)
    return summary


def _zero_or_rel(ours: float | None, theirs: float) -> None:
    """Equal to R at REL, and exactly 0 where R reports 0 (a truncated tau^2 or bound)."""
    assert ours is not None
    if theirs == 0.0:
        assert ours == 0.0
    else:
        assert ours == pytest.approx(theirs, rel=REL)


@pytest.mark.parametrize("case", POOL_CASES)
def test_the_dersimonian_laird_summary_matches_metafor_rma(case: str) -> None:
    """Oracle: metafor `rma(method = "DL")` on R's own per-record restricted means."""
    summary = _summary(case)
    rows = summary.rows
    np.testing.assert_allclose(
        [row.summary.rmst_s for row in rows], _r(case, "yi"), rtol=REL
    )
    np.testing.assert_allclose(
        [row.summary.rmst_se_s**2 for row in rows], _r(case, "vi"), rtol=REL
    )
    assert summary.k_pooled == len(rows) and summary.q_df == len(rows) - 1
    assert summary.pooled_s == pytest.approx(_r1(case, "dl_mu"), rel=REL)
    assert summary.pooled_se_s == pytest.approx(_r1(case, "dl_se"), rel=REL)
    assert summary.pooled_lower_s == pytest.approx(_r1(case, "dl_ci_lb"), rel=REL)
    assert summary.pooled_upper_s == pytest.approx(_r1(case, "dl_ci_ub"), rel=REL)
    _zero_or_rel(summary.between_var_s2, _r1(case, "dl_tau2"))
    assert summary.q_statistic == pytest.approx(_r1(case, "dl_q"), rel=REL, abs=1e-12)
    _zero_or_rel(summary.i2_fraction, _r1(case, "dl_i2_percent") / 100.0)


@pytest.mark.parametrize("case", POOL_CASES)
def test_the_q_profile_interval_matches_metafor_confint(case: str) -> None:
    """Oracle: metafor `confint` (Q-profile, tol 1e-12, search widened) on the DL fit."""
    summary = _summary(case)
    _zero_or_rel(summary.between_var_lower_s2, _r1(case, "qp_tau2_lb"))
    _zero_or_rel(summary.between_var_upper_s2, _r1(case, "qp_tau2_ub"))


@pytest.mark.parametrize("case", POOL_CASES)
def test_the_prediction_interval_matches_metafor_predict_riley(case: str) -> None:
    """Oracle: metafor `predict(predtype = "Riley")`, t with k - 2 degrees of freedom."""
    summary = _summary(case)
    assert summary.prediction_lower_s == pytest.approx(_r1(case, "pi_lb"), rel=REL)
    assert summary.prediction_upper_s == pytest.approx(_r1(case, "pi_ub"), rel=REL)
