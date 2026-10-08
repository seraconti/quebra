"""Kaplan-Meier: the estimator the within-calibration tier depends on and nothing validated.

Measured before this file existed: `kaplan_meier` appeared in exactly two test modules, and
the strongest claim either made about the estimator's output was `len(curve.time_min) > 1`.
The module was 76% line-covered at the time, which is the point - coverage measures
execution, not checking.

Indexed by subject (`kaplan_meier.run`) per AGENTS.md section 7. The oracles differ across
the file and each test names its own in its first docstring line.

The censored case in `test_a_hand_computed_censored_case` was derived independently of this
implementation, from the textbook definition, before being compared against it. An oracle
computed by reading the code it checks is not an oracle.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from quebra.analyzers import event_table, windows
from quebra.analyzers import kaplan_meier as km

# Marked PER TEST, not at module level. R6.2 requires that of a file whose tests answer more
# than one kind of question, and this one does: eleven assert an estimator value against a
# named oracle, one is the carve-to-estimator composition, and four are algebraic or
# contract assertions with no external oracle. A module-level `statistical` mark counted all
# thirty as statistical and inflated the reported tier census.


def _curve(durations, observed, **kw):
    return km.run(
        km.KaplanMeierInputs(
            duration_min=np.asarray(durations, dtype=float),
            death_observed=np.asarray(observed, dtype=bool),
            **kw,
        )
    )


# --------------------------------------------------------------- analytic special cases


@pytest.mark.statistical
@pytest.mark.parametrize("n", [2, 5, 40])
def test_with_no_censoring_it_is_one_minus_the_empirical_cdf(n) -> None:
    """Oracle: the empirical survival function, which KM reduces to when nothing is censored.

    With every observation an observed death the product-limit collapses term by term to
    `(n - i)/n`, so `S(t)` is exactly the fraction of durations strictly greater than `t`.
    Computed here from the durations, not read off the curve.
    """
    durations = np.arange(1, n + 1, dtype=float)
    curve = _curve(durations, np.ones(n, dtype=bool))

    for t, s in zip(curve.time_min, curve.survival):
        expected = float(np.count_nonzero(durations > t)) / n
        assert s == pytest.approx(expected, abs=1e-12), f"S({t})"
    assert curve.survival[0] == 1.0
    assert curve.survival[-1] == pytest.approx(0.0, abs=1e-12)


@pytest.mark.statistical
def test_a_hand_computed_censored_case() -> None:
    """Oracle: a seven-window case worked out by hand from the definition.

    Durations and status, in minutes:

        2 death | 3 death | 3 CENSORED | 5 CENSORED | 7 death | 7 death | 9 CENSORED

    Risk sets, counting an observation as at risk at time `t` when its duration is >= t, so
    the censoring at 3 is still in the risk set for the death at 3:

        t=2  n=7  d=1   S = 6/7                    = 0.857142857
        t=3  n=6  d=1   S = 6/7 * 5/6  = 5/7       = 0.714285714
        t=7  n=3  d=2   S = 5/7 * 1/3  = 5/21      = 0.238095238

    t=5 and t=9 are censoring-only and emit no point. This is the case that catches a
    mis-stepped risk set, which the ECDF reduction above cannot see: with no censoring
    every `n_j` is `n - i` regardless of how the loop counts.

    Had the censoring at 3 been removed BEFORE the death at 3, the risk set there would be
    5 and the curve would read `S(3) = 24/35`, `S(7) = 8/35`. Both are asserted against
    below. The MEDIAN is 7 under either convention, so it cannot guard this and is not
    used to.
    """
    curve = _curve([2, 3, 3, 5, 7, 7, 9], [True, True, False, False, True, True, False])

    assert curve.time_min.tolist() == [0.0, 2.0, 3.0, 7.0]
    assert curve.n_at_risk.tolist() == [7, 7, 6, 3]
    assert curve.survival == pytest.approx([1.0, 6 / 7, 5 / 7, 5 / 21], abs=1e-12)
    assert curve.n_deaths == 4
    assert curve.n_censored == 3

    # Greenwood, derived by hand from the same risk sets: the cumulative sum of
    # d/(n(n-d)) is 1/42, then +1/30 = 2/35, then +2/3 = 76/105. The band follows from
    # S**exp(+/- z sqrt(V)/|log S|), which at t=3 gives (0.258153665, 0.919797456).
    # Without this the variance estimator's only oracle is scipy, which is tier-4
    # evidence; a hand derivation is what makes the tier-1 claim in the instrument
    # report true.
    assert curve.band_lower[2] == pytest.approx(0.25815366545879487, abs=1e-12)
    assert curve.band_upper[2] == pytest.approx(0.9197974560448583, abs=1e-12)

    # The censoring-first convention, named so a silent switch cannot pass.
    assert curve.n_at_risk[2] != 5
    assert curve.survival[2] != pytest.approx(24 / 35, abs=1e-12)
    assert curve.survival[3] != pytest.approx(8 / 35, abs=1e-12)


@pytest.mark.statistical
def test_the_risk_set_at_a_tie_still_contains_the_censored_window() -> None:
    """Oracle: the death-before-censoring tie convention, stated in `kaplan_meier.run`.

    A death and a censoring recorded at the same time: the censored window has not yet
    left, so it counts in that time's risk set. With 2 observations tied at t=1, one death
    and one censoring, the risk set is 2 and `S = 1/2`. Were the censoring removed first
    the risk set would be 1 and `S` would be 0, so the two conventions are distinguishable
    and this asserts which one ships.
    """
    curve = _curve([1, 1], [True, False])

    assert curve.n_at_risk.tolist() == [2, 2]
    assert curve.survival == pytest.approx([1.0, 0.5], abs=1e-12)
    assert curve.survival[-1] != 0.0, (
        "S = 0 would mean the censored window was dropped from its own tie's risk set"
    )


@pytest.mark.unit
@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_the_estimate_is_invariant_to_input_order(seed) -> None:
    """Oracle: KM is a function of the multiset of (duration, status) pairs.

    This is the property the `lexsort` at the head of `run` exists to provide, and the one
    worth asserting: the sort key itself is measurably a no-op for the OUTPUT, because the
    tie block recomputes `n_j` at the block's first index and counts deaths across the
    whole block, so intra-tie ordering cannot reach the result. Pinning the key directly
    would be a control that cannot fail; pinning the property it is there for cannot.
    """
    durations = np.array([2, 3, 3, 5, 7, 7, 9], dtype=float)
    observed = np.array([True, True, False, False, True, True, False])
    reference = _curve(durations, observed)

    perm = np.random.default_rng(seed).permutation(len(durations))
    shuffled = _curve(durations[perm], observed[perm])

    assert shuffled.time_min.tolist() == reference.time_min.tolist()
    assert shuffled.survival == pytest.approx(reference.survival, abs=1e-15)
    assert shuffled.n_at_risk.tolist() == reference.n_at_risk.tolist()


@pytest.mark.statistical
@pytest.mark.parametrize(
    ("durations", "observed", "r_median"),
    [
        ([1, 2], [True, True], 1.5),
        ([10, 20, 30, 40], [True, True, True, True], 25.0),
        ([10, 20, 30, 40], [True, True, False, False], 30.0),
        (list(range(1, 25)), [True] * 24, 12.5),
        (list(range(1, 35)), [True] * 34, 17.5),
    ],
    ids=[
        "flat_to_next_death",
        "sample_median",
        "flat_to_last_follow_up",
        "one_ulp_above_one_half",
        "one_ulp_below_one_half",
    ],
)
def test_the_median_on_a_flat_stretch_at_one_half_is_the_stretch_midpoint(
    durations, observed, r_median
) -> None:
    """Oracle: R survival 3.8.6 `quantile(survfit(Surv(durations, observed) ~ 1), 0.5)`.

    Run in R 4.5.3 on these exact inputs; the last two are `Surv(1:24, rep(1, 24))` and
    `Surv(1:34, rep(1, 34))`. The first three land S exactly on 0.5. The first-crossing
    rule would answer 1, 20 and 20; R takes the midpoint of the flat stretch, which ends at
    the next death or, with none, at the largest observed time. Without censoring that is
    the ordinary sample median.

    The last two land one ulp off 0.5 (R prints S = 0.50000000000000011 at 12 and
    0.49999999999999989 at 17, and so does `run`). Only `QUANTILE_TOL` makes those the flat
    stretch R sees: without it the answers are 13.0 and 17.0.
    """
    curve = _curve(durations, observed)

    assert curve.median_survival_min == pytest.approx(r_median, abs=1e-12)


@pytest.mark.unit
@pytest.mark.parametrize("p", [0.0, 1.0, -0.2, 1.5, float("nan")])
def test_a_quantile_outside_the_open_unit_interval_is_refused(p) -> None:
    """Oracle: specification. The quantile rule is stated for 0 < p < 1 only.

    Unchecked, p = 1.5 returned None (read as "never reached") and p = -0.2 the first step
    time; at p = 0 R's `survival:::findq` special-cases the answer to the first time, which
    this rule does not reproduce.
    """
    with pytest.raises(ValueError, match=r"must lie in \(0, 1\)"):
        km.survival_quantile(
            np.array([0.0, 1.0, 2.0]), np.array([1.0, 0.5, 0.0]), p, 2.0
        )


@pytest.mark.statistical
def test_a_curve_that_never_reaches_one_half_has_no_median() -> None:
    """Negative control for the median: without it, returning a constant would pass above."""
    curve = _curve([1, 2, 3, 4], [True, False, False, False])

    assert curve.survival[-1] == pytest.approx(0.75, abs=1e-12)
    assert curve.median_survival_min is None


@pytest.mark.unit
def test_the_band_is_undefined_at_exactly_the_endpoints_and_defined_between() -> None:
    """Oracle: `loglog_band`'s stated exclusion set, `{S = 1} u {S = 0}`.

    Named for what it pins rather than for the infinite Greenwood term, because that term
    is NOT separately observable. `denom == 0` holds only when `deaths == n_j`, which sets
    `s = 0` and ends the loop, so the point is already excluded by the `survival > 0.0`
    mask before the infinity is ever read. Measured: replacing `np.inf` with `0.0` in the
    accumulator leaves this whole file green. The `np.isfinite(greenwood)` clause in the
    mask is therefore defence in depth, not a third exclusion condition, and this test does
    not pretend to exercise it.

    Asserted as the DISCRIMINATING pair rather than as "the band is NaN there", which holds
    for every excluded point regardless and so could not fail: NaN at exactly the two
    endpoints, finite between. A change making the band NaN everywhere fails this.
    """
    curve = _curve([1, 2, 2], [True, True, True])

    assert curve.survival == pytest.approx([1.0, 2 / 3, 0.0], abs=1e-12)
    nan_at = set(np.flatnonzero(np.isnan(curve.band_lower)).tolist())
    assert nan_at == {0, 2}, (
        "band must be undefined at S=1 and S=0, and defined between"
    )
    assert np.isfinite(curve.band_lower[1]) and np.isfinite(curve.band_upper[1])
    assert curve.band_lower[1] <= curve.survival[1] <= curve.band_upper[1]


# ------------------------------------------ cross-implementation, no new dependency, no R

CASES = {
    "no censoring": ([2, 3, 5, 7, 11], [True] * 5),
    "interleaved": (
        [2, 3, 3, 5, 7, 7, 9],
        [True, True, False, False, True, True, False],
    ),
    "censoring first": ([1, 2, 3, 4], [False, True, False, True]),
    "heavily censored": ([1, 2, 3, 4, 5, 6], [True, False, False, False, False, False]),
    "tied deaths": ([4, 4, 4, 9], [True, True, True, False]),
}


def _scipy_sf(durations, observed):
    from scipy.stats import CensoredData, ecdf

    d = np.asarray(durations, dtype=float)
    o = np.asarray(observed, dtype=bool)
    return ecdf(CensoredData(uncensored=d[o], right=d[~o])).sf


def _common_grid(curve, sf):
    """Union of both step functions' jump points, with KM's synthetic t=0 row dropped.

    `run` prepends a `t = 0`, `S = 1` row that `ecdf` has no counterpart for, and appends a
    point only where `deaths > 0` while `ecdf` steps at every distinct observation time. So
    the arrays are different lengths and an elementwise comparison misaligns every point.
    Both are step functions, so evaluating each on the union grid is the comparison that
    means anything.
    """
    return np.unique(np.concatenate([curve.time_min[1:], sf.quantiles]))


@pytest.mark.statistical
@pytest.mark.parametrize("case", list(CASES))
def test_the_product_limit_estimate_agrees_with_scipys_ecdf(case) -> None:
    """Oracle: `scipy.stats.ecdf` on `CensoredData`, an independent implementation.

    scipy is already a runtime dependency, so this is a reference cross-check that costs
    no new package and does not need R. It is the strongest single piece of evidence in
    this file: an independent implementation of the same estimator, arrived at by a
    different route, agreeing to machine precision.
    """
    durations, observed = CASES[case]
    curve = _curve(durations, observed)
    sf = _scipy_sf(durations, observed)

    grid = _common_grid(curve, sf)
    assert km._step_eval(curve.time_min, curve.survival, grid) == pytest.approx(
        sf.evaluate(grid), abs=1e-12
    )


@pytest.mark.statistical
@pytest.mark.parametrize("case", list(CASES))
def test_the_log_log_band_agrees_with_scipys_log_log_interval(case) -> None:
    """Oracle: `ecdf(...).sf.confidence_interval(method="log-log")`, the same construction.

    Compared only where BOTH are finite. `loglog_band` returns NaN where `S = 1`, where
    `S = 0` and where Greenwood is infinite, and scipy does not use the same exclusion set,
    so the excluded points are a difference in reporting policy rather than in arithmetic.
    Stating the rule here is what stops it becoming a reconciliation later.
    """
    durations, observed = CASES[case]
    curve = _curve(durations, observed)
    lo, hi = _scipy_sf(durations, observed).confidence_interval(method="log-log")

    at = curve.time_min[1:]
    ours_lo, ours_hi = curve.band_lower[1:], curve.band_upper[1:]
    theirs_lo, theirs_hi = lo.evaluate(at), hi.evaluate(at)

    both = np.isfinite(ours_lo) & np.isfinite(theirs_lo)
    assert both.any(), f"{case}: nothing comparable, the test would pass vacuously"
    assert ours_lo[both] == pytest.approx(theirs_lo[both], abs=1e-12)
    assert ours_hi[both] == pytest.approx(theirs_hi[both], abs=1e-12)


@pytest.mark.unit
def test_a_curve_with_no_windows_has_no_distance_to_any_curve() -> None:
    """Oracle: specification. An empty side's S = 1 is a placeholder, not an estimate.

    Its support ends at 0, so unrefused it scores 0 against every curve, including one
    that falls to 0, and `compare` would rank the two as identical.
    """
    empty_table = event_table.EventTable(
        dataset_id="",
        threshold_label="3 µs",
        side=windows.SIDE_OUT_OF_SPEC,
        age_s=np.array([], dtype=float),
        event=np.array([], dtype=bool),
    )
    empty = km.curve_from_event_table(empty_table, label="empty")
    falls = _curve([1.0, 2.0], [True, True], label="falls", threshold_label="3 µs")
    assert empty.n_windows == 0 and falls.survival[-1] == 0.0

    for a, b in ((empty, falls), (falls, empty)):
        with pytest.raises(ValueError, match="holds no windows"):
            km.log_time_separation(a, b)
    for curves in ([falls, empty], [empty]):
        with pytest.raises(ValueError, match="holds no windows"):
            km.compare(curves, "3 µs")


@pytest.mark.unit
def test_run_and_compare_print_nothing(capsys) -> None:
    """Oracle: specification. A step is pure compute: its result is the artifact.

    The widest pair is a field of the comparison (`pair`), not a log line.
    """
    capsys.readouterr()
    near = _curve([1.0, 2.0], [True, True], label="near", threshold_label="t")
    far = _curve([100.0, 200.0], [True, True], label="far", threshold_label="t")
    comparison = km.compare([near, far], "t")

    assert comparison.pair == ("near", "far")
    printed = capsys.readouterr()
    assert printed.out == "" and printed.err == "", printed


@pytest.mark.statistical
def test_the_scipy_comparison_rejects_an_estimator_that_drops_censored_windows() -> (
    None
):
    """Positive control: the comparison above must be able to FAIL.

    AGENTS.md section 4 records a carve control that compared `reference` to `reference`
    and passed with a deliberately broken carve, so a cross-check with no demonstrated
    failure mode is not evidence. The wrong estimator used here is a real one: the crude
    survival curve that DISCARDS right-censored windows, which the within-calibration panel
    drew before Kaplan-Meier replaced it. On a heavily censored record the two must not
    agree.
    """
    durations, observed = CASES["heavily censored"]
    d, o = np.asarray(durations, dtype=float), np.asarray(observed, dtype=bool)
    sf = _scipy_sf(durations, observed)

    grid = np.unique(d)
    deaths_only = d[o]
    crude = np.array([np.mean(deaths_only > g) for g in grid])

    assert not np.allclose(crude, sf.evaluate(grid), atol=1e-6), (
        "dropping censored windows must change the curve, or this comparison proves nothing"
    )


# ------------------------------------------- simulation truth and the carve-to-KM handoff


def _carve_with_every_window_kind():
    """A synthetic carve producing all five window kinds, so the handoff has something to do.

    Reads at 60 s, one long gap. Yields, in order: a `scan_start` endurance bag, an
    observed `up_crossing`/`down_crossing` window, a single-read `up_crossing` window dying
    at `gap_start` with duration 0, a `gap_resume` endurance bag, and an `up_crossing`
    window censored at `scan_end`.
    """
    from quebra.analyzers import windows

    t = np.array(
        [0, 60, 120, 180, 240, 300, 360, 420, 480, 1800, 1860, 1920, 1980, 2040],
        dtype=float,
    )
    v = np.array([5, 5, 5, 1, 5, 5, 5, 1, 5, 5, 5, 1, 5, 5], dtype=float)
    return windows.run(
        windows.WindowsInputs(
            t_rel_s=t,
            values=v,
            thresholds=[("3 µs", 3.0, True)],
            dataset_id="unit",
            gap_mult=3.0,
        )
    )


@pytest.mark.integration
def test_the_handoff_drops_exactly_the_unobserved_births_and_counts_them() -> None:
    """Oracle: the carve's own `birth_type` column, on a trace built to contain both kinds.

    `make_inputs_from_windows` keeps only `up_crossing` births, because an endurance bag
    was already in spec when observation began or resumed and its recorded duration is a
    residual lifetime, not a lifetime. This is the composition that the validation-jobs
    tier of `quebraplan.md` 5.1 would test through a whole job; it is the one place the
    wiring and the estimator can jointly be wrong, and it costs one test here.
    """
    carved = _carve_with_every_window_kind()
    at = carved.windows_in_spec
    assert set(at["birth_type"]) == {"scan_start", "up_crossing", "gap_resume"}

    inputs = km.make_inputs_from_windows(at, threshold_label="3 µs", label="unit")

    kept = at[at["birth_type"] == "up_crossing"]
    assert inputs.n_windows_carved == len(at) == 5
    assert inputs.n_unobserved_birth_dropped == 2
    assert len(inputs.duration_min) == len(kept) == 3
    # durations are seconds on the table and MINUTES on the inputs
    assert inputs.duration_min == pytest.approx(
        kept["duration_s"].to_numpy() / 60.0, abs=1e-12
    )
    assert (
        inputs.death_observed.tolist()
        == (kept["death_type"] == "down_crossing").tolist()
    )


@pytest.mark.integration
def test_an_out_of_spec_table_estimates_that_side_and_a_mixed_table_is_refused() -> (
    None
):
    """Oracle: the carve's own out-of-spec table, and specification for the refusals.

    The side estimated is the table's own `side` column. Both sides concatenated are two
    lifetime laws in one sample, so the handoff raises rather than picking one; a table
    short of a column the event table reads fails the handoff's own column check.
    """
    carved = _carve_with_every_window_kind()
    out = carved.windows_out_of_spec
    inputs = km.make_inputs_from_windows(out, threshold_label="3 µs", label="out")

    assert inputs.side == windows.SIDE_OUT_OF_SPEC
    assert inputs.duration_min == pytest.approx(
        out["duration_s"].to_numpy() / 60.0, abs=1e-12
    )
    assert inputs.death_observed.tolist() == (~out["censored"]).tolist()
    assert km.run(inputs).side == windows.SIDE_OUT_OF_SPEC

    both = pd.concat([carved.windows_in_spec, out], ignore_index=True)
    with pytest.raises(ValueError, match="must hold one side"):
        km.make_inputs_from_windows(both, threshold_label="3 µs", label="mixed")
    with pytest.raises(
        KeyError, match=r"Kaplan-Meier requires columns \['t_birth_s'\]"
    ):
        km.make_inputs_from_windows(
            out.drop(columns="t_birth_s"), threshold_label="3 µs", label="out"
        )


@pytest.mark.unit
@pytest.mark.parametrize(
    ("flag", "code"), [("birth_observed", "birth_type"), ("censored", "death_type")]
)
def test_a_flag_that_disagrees_with_its_code_is_refused(flag, code) -> None:
    """Oracle: the carve's rule. A birth is observed iff it is a crossing; a window is
    censored iff it ended at a gap or at the end of the scan.

    Selection reads the flags. Flipped on an up-crossing window, as here, a flag would
    drop a lifetime or enter a censored window as a death, and only this check notices.
    The other direction, a scan-start or gap-resume birth flagged observed, is refused
    earlier by the event table: its length interval needs the read before birth, which
    that window does not have.
    """
    at = _carve_with_every_window_kind().windows_in_spec.copy()
    assert at.loc[1, "birth_type"] == "up_crossing"
    at.loc[1, flag] = not at.loc[1, flag]

    with pytest.raises(
        ValueError, match=f"{flag} flag that disagrees with their {code}"
    ):
        km.make_inputs_from_windows(at, threshold_label="3 µs", label="u")


@pytest.mark.unit
def test_every_curve_in_a_set_names_the_side_it_estimates() -> None:
    """Oracle: the (label, side) key each curve is stored under, empty sides included.

    Both sides' curves carry the threshold label alone as `label`, so `side` is what tells
    a curve opened alone which side it estimates.
    """
    always_in = windows.run(
        windows.WindowsInputs(
            t_rel_s=np.arange(6.0),
            values=np.full(6, 9.0),
            thresholds=[("3 µs", 3.0, True)],
        )
    )
    for carved, n_empty in ((_carve_with_every_window_kind(), 0), (always_in, 2)):
        kset = km.kaplan_meier_set(event_table.event_tables_from_carve(carved))
        assert {key[1] for key in kset.curves} == set(event_table.SIDES)
        assert sum(c.n_windows == 0 for c in kset.curves.values()) == n_empty
        for (_label, side), curve in kset.curves.items():
            assert curve.side == side


@pytest.mark.statistical
def test_a_zero_duration_window_is_censored_and_leaves_survival_at_one() -> None:
    """Oracle: the carve's death convention, which makes a zero-duration DEATH unreachable.

    `specvalidity08.md` measured 2 of 418 windows with `duration_s == 0.0` on a real
    record, both dying at a gap or at the scan end. That is not an accident of that record:
    `windows.py` sets `t_death_s = t[e]` for a down-crossing and `t[e-1]` otherwise, with
    `e > s` always, so an OBSERVED death always has a positive duration on a strictly
    increasing read clock. A zero-duration window can therefore only be censored, and
    `run` appends no point for it, so `S(0)` stays 1 rather than dropping.
    """
    carved = _carve_with_every_window_kind()
    inputs = km.make_inputs_from_windows(
        carved.windows_in_spec, threshold_label="3 µs", label="u"
    )
    zero = inputs.duration_min == 0.0
    assert zero.sum() == 1
    assert not inputs.death_observed[zero][0], "a zero-duration window must be censored"

    curve = km.run(inputs)
    assert curve.n_zero_duration == 1
    assert curve.time_min[0] == 0.0
    assert curve.survival[0] == 1.0, "a censored window at t=0 must not drop the curve"
    assert 0.0 in curve.censor_time_min.tolist()
    assert curve.censor_survival[curve.censor_time_min == 0.0][0] == 1.0


@pytest.mark.statistical
def test_the_band_covers_a_known_exponential_survival_at_nominal() -> None:
    """Oracle: simulation truth, `S(t) = exp(-t)`, with the durations drawn from it.

    The measurement lives in `instrument_validation.measure_band_coverage` and is IMPORTED,
    not repeated: a tier-3 verdict must quote a number the artifact produced, and that
    module records having been burned once by a measurement that lived only in a test.

    Right-censored at a fixed administrative time, so censoring is non-informative by
    construction and the estimator is entitled to the data.

    **The tolerance is built at the NOMINAL level, not at the observed coverage.** A plug-in
    `sqrt(c(1-c)/N)` shrinks as `c` approaches 1, which makes the window asymmetric and
    lenient in exactly the direction that makes a band unsafe - it would tolerate more
    under-coverage than over-coverage. `sqrt(0.95*0.05/N)` does not depend on the number
    being judged.

    **What this does and does not certify.** It certifies the band's SCALE at t = 0.5:
    measured against a 90% band (Z = 1.6449) coverage falls to 0.8914 and this fails. It
    does NOT certify the Greenwood arithmetic - mutating `n_j*(n_j - d)` to `n_j**2` leaves
    coverage unchanged at these settings, because deaths are usually 1 and the two agree to
    first order. That mutation is caught by the scipy band comparison instead, which is the
    right division of labour and is stated here so neither test is credited with the other's
    reach.

    The claim is pointwise at ONE evaluation time, fixed a priori. Coverage elsewhere on
    the curve is not measured and is not claimed.
    """
    from quebra.analyzers.instrument_validation import (
        BAND_COVERAGE_EVAL_AT,
        BAND_COVERAGE_REPLICATES,
        measure_band_coverage,
    )

    coverage = measure_band_coverage()
    se_at_nominal = float(np.sqrt(0.95 * 0.05 / BAND_COVERAGE_REPLICATES))

    assert abs(coverage - 0.95) <= 3 * se_at_nominal, (
        f"coverage {coverage:.4f} at t={BAND_COVERAGE_EVAL_AT:g} over "
        f"{BAND_COVERAGE_REPLICATES} replicates; nominal 0.95, SE at nominal "
        f"{se_at_nominal:.5f}, tolerance 3 SE = {3 * se_at_nominal:.5f}"
    )


def _crude_survival(windows_min: list[float]) -> list[tuple[float, float]]:
    """The crude estimator the panel drew before Kaplan-Meier: P(W >= t), uncensored only."""
    w = np.asarray(windows_min, dtype=float)
    w = w[np.isfinite(w)]
    return [(float(x), float(np.sum(w >= x) / len(w))) for x in np.unique(w)]


@pytest.mark.statistical
def test_the_crude_estimator_differs_from_kaplan_meier_on_three_named_axes() -> None:
    """Oracle: the two implementations, compared axis by axis on inputs built to isolate each.

    `kaplan_meier.py`'s module docstring names TWO differences from the crude estimator the
    panel drew before Kaplan-Meier (`_crude_survival` below). There is a third, and it is
    the one that is present even when the other two are switched off.

    A. CENSORING. The crude estimator was fed only `~censored` windows; KM keeps them as
       censored observations.
    B. BIRTH TYPE. The crude path did not filter on `birth_type`, so endurance bags entered
       it; KM drops them and counts the drop.
    C. RIGHT-CONTINUITY. `_crude_survival` returns `P(W >= t)`; `run` returns `P(T > t)`.
       They differ by one step on IDENTICAL uncensored inputs, and the crude curve never
       reaches 0.

    The panel draws Kaplan-Meier and no crude estimator remains in `src/`, so
    `_crude_survival` is a test-local reference: the two asserts on `crude` describe it and
    cannot fail on a production change. Every other assert reads `run`, the carve or
    `make_inputs_from_windows`.
    """
    # Axis C in isolation: same durations, no censoring, no windows with unseen births.
    durations = [1.0, 2.0, 3.0]
    crude = dict(_crude_survival(durations))
    curve = _curve(durations, [True, True, True])
    km_at = dict(zip(curve.time_min.tolist(), curve.survival.tolist()))

    assert crude[1.0] == pytest.approx(1.0)  # P(W >= 1) = 3/3
    assert km_at[1.0] == pytest.approx(2 / 3)  # P(T  > 1) = 2/3
    assert min(crude.values()) == pytest.approx(1 / 3), "crude never reaches 0"
    assert curve.survival[-1] == pytest.approx(0.0)

    # Axes A and B on the real carve: what each estimator is even given.
    carved = _carve_with_every_window_kind()
    at = carved.windows_in_spec
    inputs = km.make_inputs_from_windows(at, threshold_label="3 µs", label="u")

    crude_input = at[~at["censored"]]["duration_s"].to_numpy() / 60.0
    assert len(crude_input) == 3  # axis A: 2 censored windows never reach it
    assert (at[~at["censored"]]["birth_type"] != "up_crossing").sum() == 2  # axis B
    assert inputs.n_unobserved_birth_dropped == 2
    assert sorted(crude_input.tolist()) != sorted(inputs.duration_min.tolist()), (
        "the two estimators are not even given the same windows"
    )


# ------------------------------------------------------------------------------- MTBF
#
# Two tests, not four tiers. MTBF is a MEAN, not an estimator of a distribution: it makes
# no distributional claim, so applying the four-tier instrument scheme to it would be a
# tier scheme on a subject that does not have one.


def _events(times_unix_s):
    from quebra.core.types import CalibrationEvent

    return [CalibrationEvent(t_event_unix_s=float(t)) for t in times_unix_s]


@pytest.mark.unit
def test_mtbf_intervals_are_the_sorted_first_difference() -> None:
    """Oracle: `np.diff` on the sorted times, on a deliberately UNSORTED input.

    `run` sorts before differencing, so a pre-sorted fixture would leave that sort
    unexercised and the test would hold for an implementation that omitted it. The
    hand-written expectation is the interval set of the SORTED series.
    """
    from quebra.analyzers import mtbf

    unsorted = [500.0, 100.0, 400.0, 250.0]
    result = mtbf.run(mtbf.MtbfInputs(events=_events(unsorted)))

    assert result.intervals_s.tolist() == [150.0, 150.0, 100.0]
    assert result.event_times_unix_s.tolist() == [250.0, 400.0, 500.0]
    assert result.stats["count"] == 3
    assert result.stats["mean_s"] == pytest.approx(400 / 3)
    assert result.stats["min_s"] == 100.0 and result.stats["max_s"] == 150.0
    # `std_s` is np.std with ddof=0, a POPULATION standard deviation. Recorded rather than
    # changed: whether the spread of a sample of intervals should be the sample SD is a
    # domain decision, and nothing else in the repo depends on which it is today.
    assert result.stats["std_s"] == pytest.approx(float(np.std([150.0, 150.0, 100.0])))


@pytest.mark.unit
def test_mtbf_refuses_fewer_than_two_events() -> None:
    """Oracle: the stated precondition. One event yields no interval, so there is no mean.

    Raised rather than returned empty, per the repo's rule that a silent fallback yields a
    wrong-but-plausible result.
    """
    from quebra.analyzers import mtbf

    with pytest.raises(ValueError, match="at least 2 events"):
        mtbf.run(mtbf.MtbfInputs(events=_events([100.0])))
