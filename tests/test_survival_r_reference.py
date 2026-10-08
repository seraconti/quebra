"""Kaplan-Meier, Nelson-Aalen and the curve summaries against R's `survival` package.

Oracle, unless a test names another: R `survival` (version pinned in the fixture's meta
rows), run by `jobs/rscripts/survival_reference.R` on the inputs committed beside its
values. The suite never invokes R. Each case is built as an event table and pushed through
the shipped path: `kaplan_meier.run` via `make_inputs_from_event_table`, and the estimators
in `analyzers/survival.py`.
"""

from __future__ import annotations

import subprocess
from statistics import NormalDist

import numpy as np
import pandas as pd
import pytest

from quebra.analyzers import kaplan_meier, survival
from quebra.analyzers.event_table import EventTable
from quebra.analyzers.windows import SIDE_IN_SPEC
from tests.fixtures import R_SURVIVAL_INPUTS, R_SURVIVAL_VALUES, R_TURNBULL_INPUTS

pytestmark = pytest.mark.statistical

_INPUTS = pd.read_csv(R_SURVIVAL_INPUTS)
_VALUES = pd.read_csv(R_SURVIVAL_VALUES)
_TB_INPUTS = pd.read_csv(R_TURNBULL_INPUTS)
CASES = sorted(_INPUTS["case"].unique())
REL = 1e-9


def _table(case: str) -> tuple[EventTable, float]:
    rows = _INPUTS[_INPUTS["case"] == case].sort_values("i")
    table = EventTable(
        dataset_id="r",
        threshold_label=case,
        side=SIDE_IN_SPEC,
        age_s=rows["age"].to_numpy(dtype=float),
        event=rows["event"].to_numpy(dtype=int) == 1,
    )
    return table, float(rows["tau"].iloc[0])


def _r(case: str, quantity: str) -> np.ndarray:
    rows = _VALUES[(_VALUES["case"] == case) & (_VALUES["quantity"] == quantity)]
    return rows.sort_values("i")["value"].to_numpy(dtype=float)


def _same(ours, theirs) -> None:
    """Array fields: NaN exactly where R says NA, the rest equal to R."""
    assert not any(v is None for v in ours), ours
    ours = np.asarray(ours, dtype=float)
    np.testing.assert_array_equal(np.isnan(ours), np.isnan(theirs))
    keep = ~np.isnan(theirs)
    np.testing.assert_allclose(ours[keep], theirs[keep], rtol=REL, atol=1e-12)


def _same_quantiles(ours: tuple[float | None, ...], theirs: np.ndarray) -> None:
    """Quantile fields: None exactly where R says NA (never NaN), the rest equal to R."""
    assert [v is None for v in ours] == np.isnan(theirs).tolist(), (ours, theirs)
    _same([np.nan if v is None else v for v in ours], theirs)


def _fixture_version() -> tuple[int, ...]:
    return tuple(
        int(_r("meta", f"survival_{k}")[0]) for k in ("major", "minor", "patch")
    )


def test_the_fixture_was_written_by_the_survival_version_the_docstrings_cite() -> None:
    """Oracle: the version `analyzers/survival.py` cites, survival 3.8.6."""
    assert _fixture_version() == (3, 8, 6)


@pytest.mark.r
def test_the_fixture_version_is_the_installed_r_s(requires_rscript: str) -> None:
    """Oracle: the installed R, asked directly. A fixture it cannot re-derive is stale."""
    out = subprocess.run(
        [requires_rscript, "-e", 'cat(unlist(packageVersion("survival")))'],
        capture_output=True,
        text=True,
        timeout=120,
        check=True,
    ).stdout.split()
    assert tuple(int(v) for v in out) == _fixture_version()


@pytest.mark.parametrize("case", CASES)
def test_kaplan_meier_and_its_log_log_band_match_r(case: str) -> None:
    table, _ = _table(case)
    curve = kaplan_meier.run(
        kaplan_meier.make_inputs_from_event_table(table, label=case)
    )
    steps = slice(1, None)  # R reports death times only; ours start at (0, 1)
    _same(curve.time_min[steps] * 60.0, _r(case, "km_time"))
    _same(curve.survival[steps], _r(case, "km_surv"))
    _same(curve.band_lower[steps], _r(case, "km_lower"))
    _same(curve.band_upper[steps], _r(case, "km_upper"))


@pytest.mark.parametrize("case", CASES)
def test_nelson_aalen_and_its_aalen_variance_match_r(case: str) -> None:
    """Oracle: R at every death age; analytic H(0) = V(0) = 0 at the origin R omits."""
    table, _ = _table(case)
    na = survival.nelson_aalen(table)
    assert na.time_s[0] == na.cumulative_hazard[0] == na.variance[0] == 0.0
    _same(na.time_s[1:], _r(case, "km_time"))
    _same(na.cumulative_hazard[1:], _r(case, "na_cumhaz"))
    _same(np.sqrt(na.variance[1:]), _r(case, "na_std"))


def test_the_nelson_aalen_log_band_matches_a_hand_computation() -> None:
    """Oracle: a hand computation of H exp(+/- z sqrt(V) / H). R reports no band on H.

    Deaths at 1, 2 and 3 with 4, 3 and 1 at risk (a censoring at 2.5). H is 1/4, 7/12,
    19/12; V is 1/16, 25/144, 169/144; so sqrt(V) / H is 1, 5/7, 13/19. The band is
    undefined at the origin, where H = 0. z comes from the standard library, not scipy.
    """
    table = EventTable(
        dataset_id="",
        threshold_label="h",
        side=SIDE_IN_SPEC,
        age_s=np.array([1.0, 2.0, 2.5, 3.0]),
        event=np.array([True, True, False, True]),
    )
    na = survival.nelson_aalen(table, conf_level=0.9)
    z = NormalDist().inv_cdf(0.95)
    hazard = np.array([1 / 4, 7 / 12, 19 / 12])
    spread = z * np.array([1.0, 5 / 7, 13 / 19])
    np.testing.assert_array_equal(na.time_s, [0.0, 1.0, 2.0, 3.0])
    np.testing.assert_allclose(na.cumulative_hazard, [0.0, *hazard], rtol=1e-12)
    np.testing.assert_allclose(
        na.variance, [0, 1 / 16, 25 / 144, 169 / 144], rtol=1e-12
    )
    assert np.isnan(na.band_lower[0]) and np.isnan(na.band_upper[0])
    np.testing.assert_allclose(na.band_lower[1:], hazard * np.exp(-spread), rtol=1e-12)
    np.testing.assert_allclose(na.band_upper[1:], hazard * np.exp(spread), rtol=1e-12)


@pytest.mark.parametrize("case", CASES)
def test_quantiles_their_intervals_and_the_restricted_mean_match_r(case: str) -> None:
    """Oracle: R; `max_observed_age_s` against the largest input age.

    `tau_past_follow_up` restricts to 100 s a record watched to 2 s: R extends S flat.
    """
    table, tau = _table(case)
    probs = tuple(_r(case, "quantile_p").tolist())
    summary = survival.curve_summaries(table, quantiles=probs, rmst_tau_s=tau)
    assert summary.quantile_p == probs
    assert summary.max_observed_age_s == _INPUTS[_INPUTS["case"] == case]["age"].max()
    _same_quantiles(summary.quantile_s, _r(case, "quantile"))
    _same_quantiles(summary.quantile_lower_s, _r(case, "quantile_lower"))
    _same_quantiles(summary.quantile_upper_s, _r(case, "quantile_upper"))
    _same([summary.rmst_s], _r(case, "rmean"))
    _same([summary.rmst_se_s], _r(case, "rmean_se"))


# The cases R ran at tau = the largest observed age, the default restriction.
DEFAULT_TAU_CASES = [c for c in CASES if _table(c)[1] == np.max(_table(c)[0].age_s)]


@pytest.mark.parametrize("case", DEFAULT_TAU_CASES)
def test_the_default_restriction_is_the_largest_age_and_matches_r(case: str) -> None:
    """Oracle: R, at its own probabilities and tau, which are the defaults."""
    table, tau = _table(case)
    summary = survival.curve_summaries(table)
    assert summary.quantile_p == tuple(_r(case, "quantile_p").tolist())
    assert summary.rmst_tau_s == tau
    _same([summary.rmst_s], _r(case, "rmean"))
    _same([summary.rmst_se_s], _r(case, "rmean_se"))


# ---------------------------------------------------------------------------
# Turnbull. R's EM stops short of the maximum on a flat likelihood and exposes no
# tolerance, so R pins where the mass may sit and a likelihood floor; the masses
# themselves are pinned only where they have a closed form.
# ---------------------------------------------------------------------------

TB_CASES = sorted(_TB_INPUTS["case"].unique())


def _tb_table(case: str) -> EventTable:
    rows = _TB_INPUTS[_TB_INPUTS["case"] == case].sort_values("i")
    lo = rows["lo"].to_numpy(dtype=float)
    hi = rows["hi"].to_numpy(dtype=float)
    event = np.isfinite(hi)
    age = np.where(event, (lo + np.where(event, hi, 0.0)) / 2.0, lo)
    return EventTable(
        dataset_id="r",
        threshold_label=case,
        side=SIDE_IN_SPEC,
        age_s=age,
        event=event,
        age_lo_s=lo,
        age_hi_s=hi,
    )


def _positions(result: survival.Turnbull) -> np.ndarray:
    mid = (result.innermost_lo_s + result.innermost_hi_s) / 2.0
    return np.where(result.innermost_point, result.innermost_lo_s, mid)


def _loglik(table: EventTable, result: survival.Turnbull, mass: np.ndarray) -> float:
    q, p = result.innermost_lo_s, result.innermost_hi_s
    lo, hi = table.age_lo_s[:, None], table.age_hi_s[:, None]
    exact = (table.age_lo_s == table.age_hi_s)[:, None]
    point = result.innermost_point[None, :]
    inside = np.where(point, (q > lo) & (q <= hi), (q >= lo) & (p <= hi))
    inside = np.where(exact, point & (q == lo), inside)
    return float(np.sum(np.log(inside @ mass)))


@pytest.mark.parametrize("case", TB_CASES)
def test_turnbull_puts_mass_only_where_r_does_and_reaches_at_least_r_s_likelihood(
    case: str,
) -> None:
    """Oracle: R `survfit` on `interval2`. Each puts mass where the other does."""
    table = _tb_table(case)
    result = survival.turnbull(table)
    assert result.converged
    ours = _positions(result)
    r_time, r_mass = _r(case, "tb_time"), _r(case, "tb_mass")
    weighted = r_time[r_mass > 1e-6]
    assert all(np.any(np.isclose(ours, t, atol=1e-12)) for t in weighted)
    ours_weighted = ours[result.mass > 1e-6]
    assert all(np.any(np.isclose(weighted, x, atol=1e-12)) for x in ours_weighted)
    r_on_ours = np.array(
        [r_mass[np.isclose(r_time, x, atol=1e-12)].sum() for x in ours]
    )
    assert r_on_ours.sum() == pytest.approx(1.0, abs=1e-9)
    assert result.log_likelihood >= _loglik(table, result, r_on_ours) - 1e-9


def test_turnbull_on_disjoint_intervals_is_the_count_in_each_over_n() -> None:
    """Oracle: analytic. Each window's interval holds exactly one innermost interval."""
    result = survival.turnbull(_tb_table("tb_disjoint"))
    np.testing.assert_allclose(result.mass, [2 / 6, 3 / 6, 1 / 6], atol=1e-8)
    np.testing.assert_allclose(_positions(result), _r("tb_disjoint", "tb_time"))


def test_turnbull_with_zero_tol_stops_once_the_likelihood_stops_rising() -> None:
    """Oracle: analytic. On disjoint intervals one update reaches count/n exactly, so with
    `tol = 0` the next iteration raises the log-likelihood by nothing and the EM stops."""
    result = survival.turnbull(_tb_table("tb_disjoint"), tol=0.0, max_iter=1000)
    assert result.converged
    assert result.last_change <= 0.0
    assert result.n_iter < result.max_iter
    np.testing.assert_allclose(result.mass, [2 / 6, 3 / 6, 1 / 6], atol=1e-12)


def test_turnbull_out_of_iterations_reports_the_last_iterate_as_unconverged() -> None:
    """Oracle: specification; running out of iterations is reported, never raised.

    The returned log-likelihood is the returned masses' own, and lies below a converged
    run's, since EM never lowers the likelihood.
    """
    table = _tb_table("tb_mixed")
    short = survival.turnbull(table, max_iter=5)
    assert not short.converged
    assert short.n_iter == short.max_iter == 5
    assert short.last_change > short.tol
    assert short.mass.sum() == pytest.approx(1.0, abs=1e-12)
    assert short.log_likelihood == pytest.approx(_loglik(table, short, short.mass))
    assert short.log_likelihood < survival.turnbull(table).log_likelihood


def test_turnbull_on_an_empty_table_is_a_typed_empty_result() -> None:
    """Oracle: specification. An empty side is data, not an error."""
    result = survival.turnbull(_deaths(0))
    assert isinstance(result, survival.Turnbull) and result.n_windows == 0
    for name in ("innermost_lo_s", "innermost_hi_s", "innermost_point", "mass"):
        assert len(getattr(result, name)) == 0, name
    assert len(result.survival_after) == 0
    assert result.converged and result.n_iter == 0
    assert np.isnan(result.log_likelihood)


@pytest.mark.parametrize("case", ["ties_and_censoring", "mixed_forty"])
def test_turnbull_on_exact_deaths_and_right_censoring_is_kaplan_meier(
    case: str,
) -> None:
    """Oracle: Kaplan-Meier, itself pinned to R above. With every death exact and every
    censoring right-open, the self-consistent estimate is the product-limit one."""
    exact, _ = _table(case)
    table = EventTable(
        dataset_id="r",
        threshold_label=case,
        side=SIDE_IN_SPEC,
        age_s=exact.age_s,
        event=exact.event,
        age_lo_s=exact.age_s,
        age_hi_s=np.where(exact.event, exact.age_s, np.inf),
    )
    result = survival.turnbull(table, tol=1e-14)
    curve = kaplan_meier.run(
        kaplan_meier.make_inputs_from_event_table(table, label=case)
    )
    drops = -np.diff(curve.survival)
    at = curve.time_min[1:] * 60.0
    ours = np.array(
        [
            result.mass[
                result.innermost_point & np.isclose(result.innermost_lo_s, a)
            ].sum()
            for a in at
        ]
    )
    np.testing.assert_allclose(ours, drops, atol=1e-6)


# ---------------------------------------------------------------------------
# Arguments. Oracle: specification; an empty table refuses what a full one does.
# ---------------------------------------------------------------------------


def _deaths(n_windows: int) -> EventTable:
    return EventTable(
        dataset_id="",
        threshold_label="a",
        side=SIDE_IN_SPEC,
        age_s=np.arange(1.0, n_windows + 1.0),
        event=np.ones(n_windows, dtype=bool),
    )


@pytest.mark.parametrize("n_windows", [0, 2])
@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"rmst_tau_s": -1.0}, "rmst_tau_s must be finite"),
        ({"rmst_tau_s": np.inf}, "rmst_tau_s must be finite"),
        ({"rmst_tau_s": np.nan}, "rmst_tau_s must be finite"),
        ({"conf_level": 0.0}, "conf_level must lie in"),
        ({"conf_level": 1.0}, "conf_level must lie in"),
        ({"conf_level": np.nan}, "conf_level must lie in"),
        ({"quantiles": (0.0,)}, "quantile probability must lie in"),
        ({"quantiles": (0.5, 1.0)}, "quantile probability must lie in"),
        ({"quantiles": (-0.2,)}, "quantile probability must lie in"),
        ({"quantiles": (1.5,)}, "quantile probability must lie in"),
        ({"quantiles": (np.nan,)}, "quantile probability must lie in"),
    ],
    ids=[
        "tau-negative",
        "tau-inf",
        "tau-nan",
        "conf-0",
        "conf-1",
        "conf-nan",
        "p-0",
        "p-1",
        "p-negative",
        "p-above-1",
        "p-nan",
    ],
)
def test_curve_summaries_refuse_an_invalid_argument_whether_or_not_the_table_is_empty(
    n_windows: int, kwargs: dict, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        survival.curve_summaries(_deaths(n_windows), **kwargs)


@pytest.mark.parametrize("n_windows", [0, 2])
@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"tol": -1e-10}, "tol must be finite"),
        ({"tol": np.nan}, "tol must be finite"),
        ({"tol": np.inf}, "tol must be finite"),
        ({"max_iter": 0}, "max_iter must be at least 1"),
        ({"max_iter": -1}, "max_iter must be at least 1"),
    ],
    ids=["tol-negative", "tol-nan", "tol-inf", "max-iter-0", "max-iter-negative"],
)
def test_turnbull_refuses_an_invalid_argument_whether_or_not_the_table_is_empty(
    n_windows: int, kwargs: dict, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        survival.turnbull(_deaths(n_windows), **kwargs)
