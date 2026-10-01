"""Unit tests for the six checks: identities, oracles, and the guards.

The two tests that matter most are the ones that pin a statistic against something
INDEPENDENT of the implementation:

- `test_eq7_is_the_classical_anderson_darling` compares eq (7) at `gamma_hat = 1` with the
  textbook `-N - (1/N) sum (2i-1)[ln u_i + ln(1-u_{N+1-i})]`. Two unrelated algebraic
  forms agreeing to float precision is strong evidence the transcription from the PDF is
  right.
- `test_eq7_at_gamma_one_matches_the_limiting_ad_null` checks the same statistic against
  Marsaglia's limiting AD distribution over 6000 replicates.

Both force `gamma = 1`, so both pin the TRANSCRIPTION and neither pins the shipped path,
which divides by an estimated `gamma_hat`. The distinction matters: reading the gate as
proof of the whole implementation hides that the shipped asymptotic path is measurably
oversized at n = 20. `test_shipped_c2_asymptotic_size_matches_its_measurement` holds the
shipped size to a coarse band around its measured value, and
`test_the_divisor_scales_each_statistic_by_its_own_power` pins which divisor ships.

The third group at the bottom of the file covers `statistic_batch` and
`segments_from_windows`. The second is where a defect merging renewal segments across
unobserved read gaps is invisible to inspection.
"""

from __future__ import annotations

import numpy as np
import pytest

import quebra.analyzers.checks.c1_lewis_robinson as c1
import quebra.analyzers.checks.c2_anderson_darling as c2
import quebra.analyzers.checks.c5_rank_autocorr as c5
import quebra.analyzers.checks.c6_exchangeability as c6
import quebra.analyzers.checks.cvm_cramer_von_mises as cvm
from quebra.analyzers.checks._multiprocess import (
    GAMMA_COMPLETE,
    GAMMA_COMPLETE_SAMPLE,
    GAMMA_DEFAULT,
    GAMMA_ESTIMATORS,
    GAMMA_TRUNCATED,
    gamma_hat,
    gamma_hat_batch,
)
from quebra.analyzers import windows
from quebra.analyzers.checks._multiprocess import segments_from_windows
from quebra.analyzers.checks._permutation import (
    PermutationSet,
    block_permutations,
    permutation_p_value,
    resolve_perm,
)
from quebra.analyzers.checks.battery import ROW_KEYS, row_key, run_battery
from quebra.analyzers.checks.result import (
    CALIB_ASYMPTOTIC,
    CALIB_PERMUTATION,
    CLOCK_CALENDAR,
    CLOCK_IN_SPEC,
    Segment,
)


def exponential_segment(n: int, rng: np.random.Generator, rate: float = 1.0) -> Segment:
    """A genuinely time-truncated exponential segment: tau fixed, event count random."""
    tau = float(n) / rate
    gaps = rng.exponential(1.0 / rate, size=int(n + 10 * np.sqrt(n) + 50))
    kept = gaps[: int(np.searchsorted(np.cumsum(gaps), tau, side="left"))]
    return Segment(x=kept, tau=tau, n_censored_dropped=1)


# --------------------------------------------------------------------------- C2


def _classical_ad(x: np.ndarray, tau: float) -> float:
    u = np.cumsum(x) / tau
    n = len(u)
    i = np.arange(1, n + 1)
    return float(-n - np.sum((2 * i - 1) * (np.log(u) + np.log(1 - u[::-1]))) / n)


@pytest.mark.statistical
@pytest.mark.parametrize("n", [3, 5, 20, 60, 200])
def test_eq7_is_the_classical_anderson_darling(n):
    rng = np.random.default_rng(n)
    gaps = rng.exponential(size=n + 1)
    tau = float(gaps.sum())
    x = gaps[:n]
    assert c2._eq7(x, tau, 1.0) == pytest.approx(_classical_ad(x, tau), rel=1e-9)


@pytest.mark.statistical
def test_ad_limiting_cdf_reproduces_published_critical_values():
    """Oracle: Marsaglia & Marsaglia (2004), JSS 9(2), p. 2, 20-place limiting percentiles.

    They give 1.9329578327, 2.4923671600 and 3.8781250216 for the 90, 95 and 99 percent
    points; their digits agree with the exact limit, the series of Anderson and Darling
    (1954) eq (8), p. 768, as `scripts/verify_gold_standard.py` computes it. The paper
    quotes `adinf` to |error| < 2e-6 (z < 2) and < 8e-7 (z >= 2), but against the exact
    limit its error at these three points is -1.1e-5, +8.1e-6 and -2.6e-6 (and 1.95e-5
    near z = 1), so the tolerance is 2e-5. The
    previous band of 5e-4 also accepted 3.857 (upper tail 0.0102), the 1% point in Stephens
    (1974, JASA 69, Table 1A, p. 732), which Marsaglia & Marsaglia correct to 3.878125.
    """
    for statistic, alpha in [
        (1.9329578327, 0.10),
        (2.4923671600, 0.05),
        (3.8781250216, 0.01),
    ]:
        assert 1.0 - c2.ad_limiting_cdf(statistic) == pytest.approx(alpha, abs=2e-5)


@pytest.mark.statistical
@pytest.mark.parametrize("n", [20, 50])
def test_eq7_at_gamma_one_matches_the_limiting_ad_null(n):
    """Transcription pin for eq (7) - and narrower than it looks, deliberately.

    This forces `gamma = 1`. Conditional on N, exponential gaps under a pre-chosen `tau`
    make `T_i/tau` exactly uniform order statistics, so the statistic's null is the
    classical Anderson-Darling one; `ad_limiting_cdf` is Marsaglia's LIMITING function, so
    the agreement below is excellent-at-this-N rather than exact. Together with
    `test_eq7_is_the_classical_anderson_darling` it pins the transcription.

    It does NOT pin the shipped path, which divides by an estimated `gamma_hat` - see
    `test_shipped_c2_asymptotic_size_matches_its_measurement` for that, and the bench's size
    table for what it costs.
    """
    rng = np.random.default_rng(90210 + n)
    reps = 6000
    statistics = np.array(
        [
            c2._eq7(s.x, s.tau, 1.0)
            for s in (exponential_segment(n, rng) for _ in range(reps))
        ]
    )
    p = 1.0 - c2.ad_limiting_cdf(statistics)
    for alpha in (0.10, 0.05, 0.01):
        rate = float(np.mean(p < alpha))
        se = float(np.sqrt(alpha * (1 - alpha) / reps))
        assert abs(rate - alpha) < 4 * se, (
            f"n={n} alpha={alpha}: rejection {rate:.4f}, expected {alpha} +/- {4 * se:.4f}"
        )


@pytest.mark.statistical
@pytest.mark.parametrize("n, expected", [(20, 0.0557), (50, 0.0502)])
def test_shipped_c2_asymptotic_size_matches_its_measurement(n, expected):
    """Oracle: simulation truth, measured under `GAMMA_DEFAULT` at 40,000 replicates.

    Dividing by an ESTIMATED `gamma_hat` fattens the upper tail: at n = 20 the shipped
    asymptotic path rejects at 0.0557 against a nominal 0.05, and by n = 50 it is at
    nominal, 0.0502. Both figures come from this generator with seed `20260930 + n`,
    independent of the stream below (MC SE 0.0011).

    This is a COARSE bound, and it cannot show the oversize itself: at 4000 replicates the
    4 SE band around 0.0557 also contains the nominal 0.05, the `gamma = 1` size and the
    population-divisor size. It fails only if the shipped size moves by more than about
    0.014. `test_the_divisor_scales_each_statistic_by_its_own_power` pins which divisor
    ships. The oversize at n = 20 is a measured fact of the bench and of the figures above,
    not something this test detects.
    """
    rng = np.random.default_rng(5150 + n)
    reps = 4000
    p = np.array(
        [
            c2.run([exponential_segment(n, rng)], calibration=CALIB_ASYMPTOTIC).p_value
            for _ in range(reps)
        ]
    )
    rate = float(np.mean(p < 0.05))
    se = float(np.sqrt(0.05 * 0.95 / reps))
    assert abs(rate - expected) < 4 * se, (
        f"shipped C2 asymptotic size at n={n} is {rate:.4f}, expected ~{expected}"
    )


@pytest.mark.statistical
@pytest.mark.parametrize(
    "module, power", [(c1, 0.5), (c2, 1.0), (cvm, 1.0)], ids=["c1", "c2", "cvm"]
)
def test_the_divisor_scales_each_statistic_by_its_own_power(module, power):
    """Oracle: analytic. The sample and population `gamma_hat` differ by sqrt(N/(N-1)).

    C1 divides by `gamma_hat`, C2 and CvM by `gamma_hat**2`, so the shipped statistic is
    the 1/N one times `((N-1)/N) ** power`, with power 1/2 for C1 and 1 for the other two.
    A run that silently fell back to the 1/N form gives a ratio of exactly 1 and fails.
    """
    segment = exponential_segment(20, np.random.default_rng(21))
    n = segment.n_events
    shipped = module.run([segment], calibration=CALIB_ASYMPTOTIC).statistic
    population = module.run(
        [segment], calibration=CALIB_ASYMPTOTIC, gamma_estimator=GAMMA_COMPLETE
    ).statistic
    assert shipped / population == pytest.approx(((n - 1) / n) ** power, rel=1e-12)


@pytest.mark.unit
def test_c2_refuses_a_segment_whose_last_event_lands_on_tau():
    x = np.array([1.0, 2.0, 3.0])
    with pytest.raises(ValueError, match="failure censoring"):
        c2.statistic([Segment(x=x, tau=float(x.sum()))])


@pytest.mark.unit
def test_c2_refuses_an_asymptotic_calibration_for_multiple_segments():
    rng = np.random.default_rng(0)
    segments = [exponential_segment(20, rng) for _ in range(3)]
    with pytest.raises(ValueError, match="no asymptotic calibration"):
        c2.run(segments, calibration=CALIB_ASYMPTOTIC)


# --------------------------------------------------------------------------- C1


@pytest.mark.statistical
def test_eq16_reduces_to_eq4_for_a_single_segment():
    """Oracle: Technometrics eq (4), p. 103, by hand; both versions' eq (16), both divisors."""
    rng = np.random.default_rng(7)
    segment = exponential_segment(40, rng)
    x, tau = segment.x, segment.tau
    for estimator in (GAMMA_COMPLETE, GAMMA_COMPLETE_SAMPLE):
        gamma = gamma_hat(x, tau, estimator)
        eq4 = (
            (1.0 / gamma)
            * np.sqrt(12.0)
            / (tau * np.sqrt(len(x)))
            * (np.cumsum(x).sum() - len(x) * tau / 2.0)
        )
        for weights in c1.C1_WEIGHTS:
            assert c1.statistic(
                [segment], gamma_estimator=estimator, weights=weights
            ) == pytest.approx(eq4, rel=1e-12)


@pytest.mark.statistical
def test_c1_detects_the_trend_it_is_built_for():
    """Evidence about `a1_renewal_durations`: C1 detects the trend that assumption forbids.

    Events crowded early give a strongly negative statistic, crowded late a positive one.

    The two segments share one multiset of gaps and differ only in ORDER, so `gamma_hat`,
    `tau` and `N` are identical and the sign difference can only come from the trend term.
    (Equally-spaced gaps would be the sharper trend but have zero dispersion, which makes
    `gamma_hat` zero and the statistic undefined - the guard fires before the test can.)
    """
    rng = np.random.default_rng(3)
    gaps = np.sort(rng.exponential(size=25))
    tau = float(gaps.sum() + 5.0)
    early = c1.statistic([Segment(x=gaps, tau=tau)])  # small gaps first -> events early
    late = c1.statistic([Segment(x=gaps[::-1], tau=tau)])  # large gaps first -> late
    assert early < -2.0, f"expected a strong negative trend statistic, got {early}"
    assert late > 2.0, f"expected a strong positive trend statistic, got {late}"


@pytest.mark.statistical
def test_gamma_hat_eq10_goes_negative_where_the_complete_form_cannot():
    """Pins why the default is a complete-gap form - see `_multiprocess.gamma_hat`.

    Both estimators get the SAME input, which the previous version of this test did not:
    it fed `GAMMA_COMPLETE` a perturbed vector and `GAMMA_TRUNCATED` a constant one, so it
    demonstrated nothing about the difference between them.
    """
    x = np.array([0.9, 1.0, 1.1, 1.0, 1.0])
    tau = 5.5
    # The complete-gap form is a genuine variance and cannot go negative.
    assert gamma_hat(x, tau, GAMMA_DEFAULT) > 0.0
    # Eq (10) on the identical vector is a difference of two large terms and does.
    with pytest.raises(ValueError, match="variance is negative"):
        gamma_hat(x, tau, GAMMA_TRUNCATED)


@pytest.mark.statistical
def test_gamma_hat_distinguishes_a_negative_variance_from_a_constant_vector():
    """The two failures have different causes and must not share a message.

    Clamping a materially negative eq (10) variance to zero would report it as "every
    gap is identical", which is a different and false diagnosis.
    """
    with pytest.raises(ValueError, match="every gap is identical"):
        gamma_hat(np.ones(5), 5.0, GAMMA_COMPLETE)
    with pytest.raises(ValueError, match="variance is negative"):
        gamma_hat(np.array([0.9, 1.0, 1.1, 1.0, 1.0]), 5.5, GAMMA_TRUNCATED)


@pytest.mark.unit
@pytest.mark.parametrize("estimator", GAMMA_ESTIMATORS)
def test_the_batch_estimator_refuses_one_gap_as_its_scalar_twin_does(estimator):
    """At N = 1 the sample divisor is 1/0: the batch path must raise, not return NaN."""
    with pytest.raises(ValueError, match="at least 2 gaps"):
        gamma_hat(np.array([1.0]), 2.0, estimator)
    with pytest.raises(ValueError, match="at least 2 gaps"):
        gamma_hat_batch(np.ones((3, 1)), 2.0, estimator)


# --------------------------------------------------------------------------- shared


@pytest.mark.statistical
def test_permutation_p_value_is_tie_corrected():
    """All-tied null must give p = 1, not p = 0: a tie supports the null."""
    assert permutation_p_value(1.0, np.ones(99)) == pytest.approx(1.0)
    # Strictly smaller null draws: the minimum attainable p is 1/(1+B), never 0.
    assert permutation_p_value(5.0, np.zeros(99)) == pytest.approx(1.0 / 100.0)


@pytest.mark.unit
def test_permutations_stay_inside_their_segment():
    perm = block_permutations([3, 4, 2], n_perm=200, rng=np.random.default_rng(1))
    for lo, hi in perm.blocks():
        block = perm.indices[:, lo:hi]
        assert block.min() >= lo and block.max() < hi
        assert np.all(np.sort(block, axis=1) == np.arange(lo, hi))


@pytest.mark.statistical
@pytest.mark.parametrize(
    "check_module, kwargs",
    [
        (c5, {"variant": c5.VARIANT_STUDENTIZED}),
        (c5, {"variant": c5.VARIANT_RAW}),
        (c6, {}),
    ],
)
def test_rank_checks_are_calibrated_under_exchangeability(check_module, kwargs):
    """Size at alpha = 0.10 on iid gaps, which is these checks' exact null."""
    rng = np.random.default_rng(31415)
    reps = 600
    p = np.array(
        [
            check_module.run(
                [exponential_segment(30, rng)], n_perm=199, rng=rng, **kwargs
            ).p_value
            for _ in range(reps)
        ]
    )
    rate = float(np.mean(p < 0.10))
    se = float(np.sqrt(0.10 * 0.90 / reps))
    assert abs(rate - 0.10) < 4 * se, f"rejection {rate:.4f} at nominal 0.10"


@pytest.mark.statistical
def test_c5_and_c6_find_a_strongly_ordered_sequence():
    """A monotonically increasing duration sequence is maximally non-exchangeable."""
    x = np.linspace(1.0, 20.0, 40)
    segment = Segment(x=x, tau=float(x.sum() + 5.0))
    perm = block_permutations([len(x)], 999, np.random.default_rng(2))
    assert c5.run([segment], perm=perm).p_value <= 0.01
    assert c6.run([segment], perm=perm).p_value <= 0.01


@pytest.mark.unit
def test_battery_returns_every_row_and_shares_one_permutation_set():
    rng = np.random.default_rng(11)
    segment = exponential_segment(40, rng)
    perm = block_permutations([segment.n_events], 199, rng)
    results = run_battery([segment], perm=perm)
    assert [row_key(r) for r in results] == list(ROW_KEYS)
    assert all(r.p_value is not None for r in results)
    # Standalone calls with the same permutation set must reproduce the battery exactly.
    standalone = c1.run([segment], calibration=CALIB_PERMUTATION, perm=perm)
    battery_c1 = next(
        r for r in results if row_key(r) == ("c1_lewis_robinson", CALIB_PERMUTATION, "")
    )
    assert standalone.p_value == pytest.approx(battery_c1.p_value)


@pytest.mark.unit
def test_the_battery_runs_the_shipped_gamma_estimator():
    """Oracle: each check's own `run` under `GAMMA_DEFAULT`, which differs from 1/N here.

    The battery has its own `gamma_estimator` default and is the path every ledger and
    bench p-value takes, so a regression of that default alone must fail here. The
    permutation comparison above cannot see it: at m = 1 gamma is permutation invariant.
    """
    rng = np.random.default_rng(13)
    segment = exponential_segment(20, rng)
    perm = block_permutations([segment.n_events], 99, rng)
    rows = {row_key(r): r for r in run_battery([segment], perm=perm)}
    for module in (c1, c2, cvm):
        shipped, population = (
            module.run(
                [segment], calibration=CALIB_ASYMPTOTIC, gamma_estimator=estimator
            ).p_value
            for estimator in (GAMMA_DEFAULT, GAMMA_COMPLETE)
        )
        assert shipped != pytest.approx(population, rel=1e-6)
        row = rows[(module.CHECK_NAME, CALIB_ASYMPTOTIC, "")]
        assert row.p_value == pytest.approx(shipped, rel=1e-12)
    # And an estimator the caller names must reach every check, not only the default.
    asked = run_battery([segment], perm=perm, gamma_estimator=GAMMA_COMPLETE)
    rows_asked = {row_key(r): r for r in asked}
    for module in (c1, c2, cvm):
        expected = module.run(
            [segment], calibration=CALIB_ASYMPTOTIC, gamma_estimator=GAMMA_COMPLETE
        ).p_value
        row = rows_asked[(module.CHECK_NAME, CALIB_ASYMPTOTIC, "")]
        assert row.p_value == pytest.approx(expected, rel=1e-12)


@pytest.mark.unit
def test_battery_drops_only_the_tau_checks_when_asked():
    rng = np.random.default_rng(12)
    segment = exponential_segment(30, rng)
    perm = block_permutations([segment.n_events], 99, rng)
    keys = [
        row_key(r) for r in run_battery([segment], perm=perm, include_tau_checks=False)
    ]
    # C5, C6 and CvM-by-permutation survive; C1, C2 and CvM-ASYMPTOTIC do not.
    # CvM is the one entitled to this case: its
    # statistic has no `1/(s(1-s))` weight so it stays finite when `tau == T_N`. Its
    # asymptotic row is still dropped here, because a limiting null needs a truncation
    # time chosen independently of the events and this one is not.
    assert all(key[0].startswith(("c5", "c6", "cvm")) for key in keys)
    assert ("cvm_cramer_von_mises", CALIB_PERMUTATION, "") in keys
    assert ("cvm_cramer_von_mises", CALIB_ASYMPTOTIC, "") not in keys
    assert not any(key[0].startswith(("c1", "c2")) for key in keys)
    assert len(keys) == 4


@pytest.mark.unit
@pytest.mark.parametrize(
    "segment, match",
    [
        (Segment(x=np.array([1.0]), tau=3.0), "at least 2"),
        (Segment(x=np.array([1.0, 0.0, 2.0]), tau=9.0), "strictly positive"),
        (Segment(x=np.array([1.0, 2.0]), tau=1.0), "below the last event time"),
        (Segment(x=np.array([1.0, np.nan]), tau=9.0), "finite"),
    ],
)
def test_validate_segment_raises_rather_than_returning_nonsense(segment, match):
    with pytest.raises(ValueError, match=match):
        c1.statistic([segment])


# ------------------------------------------------------- the batch/segment paths


@pytest.mark.statistical
@pytest.mark.parametrize("sizes", [[6], [4, 5], [3, 3, 4]])
def test_statistic_batch_reproduces_the_statistic_under_the_identity_permutation(sizes):
    """A misalignment between `segments` and `perm.blocks()` would compute the observed
    statistic from one assembly and the null from another - silently, and it would corrupt
    every permutation p-value the bench produced. Pinned with an identity permutation, for
    which the batch answer must equal the scalar one exactly."""
    rng = np.random.default_rng(808)
    segments = [exponential_segment(size * 6, rng) for size in sizes]
    total = sum(s.n_events for s in segments)
    identity = PermutationSet(
        sizes=tuple(s.n_events for s in segments),
        indices=np.tile(np.arange(total), (3, 1)),
    )
    for module in (c1, c2):
        batch = module.statistic_batch(segments, identity)
        scalar = module.statistic(segments)
        assert np.allclose(batch, scalar, rtol=0, atol=1e-12), (
            f"{module.CHECK_NAME}: identity-permuted batch {batch[0]} != {scalar}"
        )


def _carve_windows(t, in_spec):
    result = windows.run(
        windows.WindowsInputs(
            t_rel_s=np.asarray(t, dtype=float),
            values=np.where(np.asarray(in_spec, dtype=bool), 1.0, -1.0),
            thresholds=[("thr", 0.0, True)],
        )
    )
    return result.windows, result.diagnostics.get("gap_spans_s")


def _gapped_record():
    """A record whose gap is flanked by OUT-of-spec reads on both sides.

    The carve cannot mark this: `gap_resume` is only emitted next to an in-spec read, so
    the window table contains no trace of the gap at all. That is the whole point.
    """
    pre = [1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0]
    post = [0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0]
    t = np.concatenate(
        [np.arange(len(pre), dtype=float), 500.0 + np.arange(len(post), dtype=float)]
    )
    return _carve_windows(t, pre + post)


@pytest.mark.unit
def test_the_gapped_record_really_hides_its_gap_from_the_window_table():
    """Positive control for the two tests below: without this, they could pass because the
    carve marked the gap and `gap_spans_s` changed nothing."""
    frame, gaps = _gapped_record()
    assert gaps, "the carve must have detected a read gap"
    assert "gap_resume" not in set(frame["birth_type"]), (
        "this record is only interesting because the carve leaves no birth-type trace"
    )


@pytest.mark.unit
@pytest.mark.parametrize("clock", [CLOCK_IN_SPEC, CLOCK_CALENDAR])
def test_segments_split_on_gap_spans_the_birth_taxonomy_cannot_see(clock):
    frame, gaps = _gapped_record()
    merged, _ = segments_from_windows(frame, clock=clock)
    split, _ = segments_from_windows(frame, clock=clock, gap_spans_s=gaps)
    assert len(merged) == 1, "birth types alone should merge this record"
    assert len(split) == 2, f"gap_spans_s should split it; got {len(split)}"


@pytest.mark.unit
def test_the_merged_calendar_segment_contains_a_fabricated_inter_event_gap():
    """Names the damage: on the calendar clock the merged segment reports the 488 s of
    UNOBSERVED time as if it were an observed inter-event interval, and feeds it to
    eqs (4) and (7) as a renewal gap."""
    frame, gaps = _gapped_record()
    merged, _ = segments_from_windows(frame, clock=CLOCK_CALENDAR)
    split, _ = segments_from_windows(frame, clock=CLOCK_CALENDAR, gap_spans_s=gaps)
    gap_before, gap_after = gaps[0]
    unobserved_s = gap_after - gap_before
    assert merged[0].x.max() > unobserved_s, "the merged segment should span the gap"
    for segment in split:
        assert segment.x.max() < unobserved_s, (
            "no split segment may contain an interval longer than the read gap"
        )


@pytest.mark.unit
def test_interior_segments_are_truncated_at_the_gap_start_not_at_an_event():
    """`tau` for an interior segment is when watching STOPPED. Taking it from the last
    window's death is the event-determined boundary eq (7) forbids."""
    frame, gaps = _gapped_record()
    split, _ = segments_from_windows(frame, clock=CLOCK_CALENDAR, gap_spans_s=gaps)
    gap_start_s = gaps[0][0]
    first_births = frame["t_birth_s"].to_numpy(dtype=float)
    assert split[0].tau == pytest.approx(gap_start_s - first_births[0]), (
        "the interior segment must be truncated at the gap start"
    )


@pytest.mark.unit
def test_omitting_gap_spans_keeps_the_old_behaviour():
    """The parameter is optional, and its absence must degrade predictably rather than
    change results for callers that never had a gap."""
    rng = np.random.default_rng(3)
    t = np.arange(40, dtype=float)
    frame, gaps = _carve_windows(t, rng.random(40) > 0.4)
    assert not gaps, "a uniformly spaced record has no gaps"
    without, _ = segments_from_windows(frame, clock=CLOCK_IN_SPEC)
    with_empty, _ = segments_from_windows(frame, clock=CLOCK_IN_SPEC, gap_spans_s=gaps)
    assert len(without) == len(with_empty)
    for a, b in zip(without, with_empty):
        assert np.array_equal(a.x, b.x) and a.tau == b.tau


# -------------------- calendar-clock truncation (R8.1)


# One read per second, so the median spacing is 1.0 and "one read interval" is 1.0.
DT_S = 1.0


def _carve_unit_series(in_spec: list[int]):
    """Carve a unit-spaced read series. Returns (window table, diagnostics)."""
    t = np.arange(len(in_spec), dtype=float) * DT_S
    result = windows.run(
        windows.WindowsInputs(
            t_rel_s=t,
            values=np.where(np.asarray(in_spec, dtype=bool), 1.0, -1.0),
            thresholds=[("thr", 0.0, True)],
        )
    )
    return result.windows, result.diagnostics


def _record_ending_on_one_in_spec_read() -> list[int]:
    """Enough complete excursions to be a segment, then a single in-spec final read.

    The final read is flanked by an out-of-spec read before it and nothing after, so the
    carve births a window at it and kills it at `scan_end` in the same instant.
    """
    body = [1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0]
    return body + [0, 1]


@pytest.mark.unit
def test_the_fixture_really_ends_in_a_zero_duration_window():
    """Positive control. Without it every test below could pass on a record whose final
    window has an ordinary positive duration, asserting nothing about the degeneracy."""
    frame, _ = _carve_unit_series(_record_ending_on_one_in_spec_read())
    last = frame.iloc[-1]
    assert last["duration_s"] == 0.0, (
        f"the fixture must end in a zero-duration window; got {last['duration_s']}"
    )
    assert last["death_type"] == "scan_end"
    assert float(last["t_birth_s"]) == float(last["t_death_s"])


@pytest.mark.unit
def test_calendar_tau_equals_T_N_when_the_final_window_has_zero_duration():
    """Arithmetic: tau - T_N = t_death[-1] - t_birth[-1], which is 0 for a zero-duration
    final window. Asserted as an exact equality, not a tolerance, because it is an identity
    between two subtractions of the same float."""
    frame, _ = _carve_unit_series(_record_ending_on_one_in_spec_read())
    segments, _dropped = segments_from_windows(
        frame, clock=CLOCK_CALENDAR, min_events=2
    )
    assert segments, "the fixture must produce at least one calendar segment"

    final = segments[-1]
    t_n = float(np.cumsum(final.x)[-1])
    assert final.tau == t_n, (
        f"expected tau == T_N exactly; got tau={final.tau!r} T_N={t_n!r}"
    )


@pytest.mark.unit
def test_an_observation_end_at_the_last_read_does_not_clear_the_degeneracy():
    """The correction this measurement forced on SPEC 0008's own R8.1 prescription.

    R8.1 first said the observation end is "the last `t_read_s`". That is exactly where the
    zero-duration window is born, so it leaves tau == T_N untouched. `jobs/bench/arms.py`
    ends "one full read interval PAST the last read" and that is what clears it. Both
    branches are asserted, so the wrong prescription cannot come back unnoticed.
    """
    frame, _ = _carve_unit_series(_record_ending_on_one_in_spec_read())
    last_read_s = float(frame["t_death_s"].to_numpy(dtype=float).max())

    at_last_read, _ = segments_from_windows(
        frame,
        clock=CLOCK_CALENDAR,
        min_events=2,
        observation_end_s=last_read_s,
    )
    still_degenerate = at_last_read[-1]
    assert still_degenerate.tau == float(np.cumsum(still_degenerate.x)[-1]), (
        "an observation end at the last read must leave tau == T_N; if this fails the "
        "degeneracy has moved and R8.1's prescription needs re-measuring"
    )

    one_interval_past, _ = segments_from_windows(
        frame,
        clock=CLOCK_CALENDAR,
        min_events=2,
        observation_end_s=last_read_s + DT_S,
    )
    cleared = one_interval_past[-1]
    slack = cleared.tau - float(np.cumsum(cleared.x)[-1])
    assert slack == pytest.approx(DT_S), (
        f"one read interval past the last read must give exactly that much slack; got {slack}"
    )
    # And the statistic C2 refused above is now computable.
    assert np.isfinite(c2.statistic(one_interval_past))


@pytest.mark.unit
def test_an_interior_block_ending_in_a_zero_duration_window_is_not_reachable_by_observation_end():
    """Why the truncation rule cannot be the whole fix: `observation_end_s` applies to the
    FINAL block only (`_multiprocess.py:222-229`), so an interior block ending in a
    zero-duration window at a gap start stays degenerate under every value of it.

    Measured on the real record that motivated this: 070723_6D2S_qubit4 carries two
    zero-duration windows, one at `gap_start` and one at `scan_end`.
    """
    # In spec at the last read before the gap, so the pre-gap block ends zero-duration.
    pre = [1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1]
    post = [1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0]
    t = np.concatenate(
        [
            np.arange(len(pre), dtype=float),
            500.0 + np.arange(len(post), dtype=float),
        ]
    )
    result = windows.run(
        windows.WindowsInputs(
            t_rel_s=t,
            values=np.where(np.asarray(pre + post, dtype=bool), 1.0, -1.0),
            thresholds=[("thr", 0.0, True)],
        )
    )
    frame = result.windows
    gaps = result.diagnostics.get("gap_spans_s")
    assert gaps, "the fixture must contain a detected read gap"

    zero_duration = frame[frame["duration_s"] == 0.0]
    assert len(zero_duration), "the fixture must produce a zero-duration window"
    assert "gap_start" in set(zero_duration["death_type"]), (
        "the zero-duration window must be the one that ends the INTERIOR block"
    )

    far_past_the_end = float(t.max()) + 1000.0
    segments, _ = segments_from_windows(
        frame,
        clock=CLOCK_CALENDAR,
        min_events=2,
        gap_spans_s=gaps,
        observation_end_s=far_past_the_end,
    )
    interior = segments[0]
    assert interior.tau == float(np.cumsum(interior.x)[-1]), (
        "the interior block must stay at tau == T_N however far the observation end is "
        "pushed; if this fails, observation_end_s has started reaching interior blocks"
    )


# -------------------- the shared permutation-set guard (R8.5a)


SIZES = [3, 4, 2]
N_PERM = 50


def _blocked_perm(sizes=SIZES, seed=1) -> PermutationSet:
    return block_permutations(sizes, N_PERM, np.random.default_rng(seed))


@pytest.mark.unit
def test_the_set_it_builds_matches_calling_block_permutations_directly():
    """The extraction must not have changed WHICH permutations are drawn: the same seed has
    to give the same matrix, or every permutation p-value in the package moves."""
    direct = block_permutations(SIZES, N_PERM, np.random.default_rng(7))
    viahelper = resolve_perm(SIZES, None, N_PERM, np.random.default_rng(7))
    assert np.array_equal(direct.indices, viahelper.indices)


@pytest.mark.unit
def test_a_matching_set_is_returned_unchanged_and_not_rebuilt():
    """Identity, not equality: rebuilding would consume the rng and draw different
    permutations while looking correct."""
    supplied = _blocked_perm()
    assert resolve_perm(SIZES, supplied, N_PERM, np.random.default_rng(99)) is supplied


@pytest.mark.unit
def test_a_mismatched_set_raises_rather_than_being_rebuilt():
    """The load-bearing half. A set blocked for different segment sizes belongs to a
    different record; silently rebuilding it would yield a plausible, wrong p-value."""
    wrong = _blocked_perm(sizes=[5, 5])
    with pytest.raises(ValueError, match="permutation set is blocked as"):
        resolve_perm(SIZES, wrong, N_PERM, np.random.default_rng(1))


@pytest.mark.unit
@pytest.mark.parametrize("sizes", [[3, 4, 2], (3, 4, 2), np.array([3, 4, 2])])
def test_both_call_spellings_agree(sizes):
    """The six sites differed only in whether `sizes` arrived as a list from
    `segment_sizes(segments)` or as a precomputed sequence. Both must normalise the same, or
    a tuple-vs-list comparison would reject a set that matches."""
    supplied = _blocked_perm()
    assert resolve_perm(sizes, supplied, N_PERM, np.random.default_rng(1)) is supplied


@pytest.mark.unit
def test_an_absent_rng_still_raises_block_permutations_own_message():
    """Not duplicated in the helper, so there is one message for that failure rather than two
    that can drift. `AGENTS.md` records why an implicit rng is fatal here."""
    with pytest.raises(ValueError, match="requires an explicit rng"):
        resolve_perm(SIZES, None, N_PERM, None)


@pytest.mark.unit
def test_an_rng_is_not_required_when_a_usable_set_is_supplied():
    """The construct-or-validate asymmetry, asserted: only the building branch needs a
    generator. This is why the block cannot be expressed as a precondition, and why R8.5b
    evaluates a contract library against a different invariant."""
    supplied = _blocked_perm()
    assert resolve_perm(SIZES, supplied, N_PERM, None) is supplied


@pytest.mark.policy
def test_no_call_site_still_carries_its_own_copy():
    """The guard against a seventh copy reappearing, and against this extraction being
    reverted in one file only - which is the `AGENTS.md` section 4 failure exactly."""
    from pathlib import Path

    checks = (
        Path(__file__).resolve().parents[1] / "src" / "quebra" / "analyzers" / "checks"
    )
    carriers = [
        path.name
        for path in sorted(checks.glob("*.py"))
        if "permutation set is blocked as" in path.read_text(encoding="utf-8")
    ]
    assert carriers == ["_permutation.py"], (
        f"the guard should live only in _permutation.py; also found in {carriers}"
    )
