"""Tier 3: does each instrument hold its nominal level under a null it is entitled to?

Tier 2 (`test_checks_published_values.py`) pins the ARITHMETIC against a published worked
example. It cannot tell you whether the p-value means what it says. That is this file.

**Two different claims, held to two different standards, and the difference is the point.**

- C5 and C6 are calibrated by PERMUTATION, which is exact under exchangeability. Their level
  is not an approximation that improves with n, it is a combinatorial identity, so a
  deviation is a defect and these are ASSERTED.
- C1, C2 and CvM are calibrated by an ASYMPTOTIC limit, divided by an estimated `gamma_hat`.
  Their finite-n level misses nominal, and BY HOW MUCH AND IN WHICH DIRECTION DEPENDS ON THE
  GAP DISTRIBUTION, which is why no direction is asserted here.

  Under WEIBULL gaps the bench measures, at n = 20, alpha = 0.05, in-spec clock, arm
  A_iid_weibull, unquantised, uncensored, asymptotic calibration - four cells, not two,
  because shape varies:
      C1  shape 0.75 -> 0.0580 (mc_se 0.0052),  shape 1.50 -> 0.0635 (0.0055)
      C2  shape 0.75 -> 0.0575 (0.0052),        shape 1.50 -> 0.0790 (0.0060)
  Above nominal in all four, two of them by less than 1.5 MC SE. An earlier draft of this docstring quoted "0.069" and
  "0.0757", which are the POOLED MEANS over the two shapes and are in no cell of the table -
  the first instance CLAUDE.md's claims-discipline section records, reintroduced here.

  Under EXPONENTIAL gaps, which is what `measure_asymptotic_size` generates, tau = 20 gives
  C1 0.0650, C2 0.0575, CvM 0.0575, the tier-3 cells of
  `jobs/bench/results/instrument_report.md`: C1 above nominal, C2 and CvM about one MC SE
  (0.0067) above it. Both are legitimate iid nulls (these check trend against renewal, not
  exponentiality), so the size of an asymptotic check here still depends on the gap
  distribution through the estimated gamma_hat, and the assertion below stays a wide bound
  rather than a direction.

  Those are at the shipped defaults (seed 777, 1200 replicates) and come from a
  TIME-TRUNCATED null: tau fixed in advance, event count random. A tau derived from the draw
  instead pins the fraction of the window past the last event at n/(n+1), which removes the
  leftover window's variance - eq (7)'s tail term, for C2 - and turns the measured rate
  conservative. That is the size of a scheme the ASYMPTOTIC theory does not cover; the
  permutation checks below are indifferent to it, and `_iid_segments` still builds it for
  them.

**Why not a KS test on the p-values.** At B = 199 a permutation p-value lives on a 200-point
grid, so its distribution differs from the continuous uniform by up to 1/(2B) in sup norm by
construction. KS would be testing the grid, not the check. The exact statement for a discrete
uniform is `P(p <= k/(B+1)) = k/(B+1)` at the grid points, and that is what is tested - at
alphas chosen to LAND on the grid: 0.01 = 2/200, 0.05 = 10/200, 0.10 = 20/200.

**Aggregation and multiplicity, stated before any verdict is read off them.** The exactness
test evaluates 2 checks x 3 layouts x 3 alphas = 18 comparisons against one shared band. The
band is `bonferroni_z_crit(18)` standard errors of the NULL rate, reusing
`calibration_summary.null_se` and `bonferroni_z_crit` so this file and the bench share ONE
definition of calibrated. Without that correction the largest of 18 deviations is about 2.7
SE by chance and a correct check fails as often as a broken one - the failure this repo has
already made once, recorded in CLAUDE.md.

The resolution is honest rather than flattering: at 2000 replicates the band is +/-0.020 at
alpha = 0.10 and +/-0.0067 at alpha = 0.01, so the 0.01 cell can only catch a gross error.
That is stated rather than hidden behind a passing assertion.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

import quebra.analyzers.checks.c1_lewis_robinson as c1
import quebra.analyzers.checks.c2_anderson_darling as c2
import quebra.analyzers.checks.c5_rank_autocorr as c5
import quebra.analyzers.checks.c6_exchangeability as c6
import quebra.analyzers.checks.cvm_cramer_von_mises as cvm
from quebra.analyzers.calibration_summary import bonferroni_z_crit, null_se
from quebra.analyzers.checks._permutation import block_permutations, permutation_p_value
from quebra.analyzers.checks.result import CALIB_ASYMPTOTIC, Segment
from quebra.analyzers.instrument_validation import (
    measure_asymptotic_size,
)

pytestmark = pytest.mark.statistical

# B is chosen so every alpha below is exactly a grid point: p = k/(B+1) = k/200.
N_PERM = 199
ALPHAS = (0.01, 0.05, 0.10)
REPLICATES = 2000

# Three layouts, because the permutation is BLOCKED by segment and a bug in the blocking
# shows up only when there is more than one block. "many short segments" is the shape a
# heavily gapped record actually produces here.
LAYOUTS = {
    "ungapped n=30": [30],
    "three segments of 10": [10, 10, 10],
    "six segments of 5": [5, 5, 5, 5, 5, 5],
}

# 2 checks x 3 layouts x 3 alphas.
N_COMPARISONS = 2 * len(LAYOUTS) * len(ALPHAS)

# One fixed seed per cell, written out rather than derived, so a rerun reproduces exactly.
SEEDS = {
    (check, layout): 90_000 + 10 * i + j
    for i, check in enumerate(("c5", "c6"))
    for j, layout in enumerate(LAYOUTS)
}


def _iid_segments(sizes: list[int], rng: np.random.Generator) -> list[Segment]:
    """Exponential gaps: iid, hence exchangeable, which is C5's and C6's exact null."""
    segments = []
    for n in sizes:
        g = rng.exponential(size=n)
        segments.append(Segment(x=g, tau=float(g.sum() + g.mean())))
    return segments


def _null_p_values(check_module, sizes: list[int], seed: int, **kwargs) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = np.empty(REPLICATES, dtype=float)
    for i in range(REPLICATES):
        segments = _iid_segments(sizes, rng)
        perm = block_permutations([s.n_events for s in segments], N_PERM, rng)
        out[i] = check_module.run(segments, perm=perm, **kwargs).p_value
    return out


@pytest.mark.parametrize("layout_name", list(LAYOUTS))
@pytest.mark.parametrize("check_name", ["c5", "c6"])
def test_the_permutation_p_value_is_exact_under_exchangeability(
    check_name, layout_name
):
    """`P(p <= k/(B+1)) = k/(B+1)`, the identity a permutation test is supposed to satisfy.

    Exact means exact: under exchangeability the observed statistic is equally likely to
    occupy any rank among itself and the B permuted values, so this holds combinatorially
    and does not improve with n. A miss here is a defect in the calibration, not a small-n
    cost - which is why this is asserted and the asymptotic checks below are not.
    """
    module = {"c5": c5, "c6": c6}[check_name]
    kwargs = {"variant": c5.VARIANT_STUDENTIZED} if check_name == "c5" else {}
    # EXPLICIT seed table, not `hash(...)`. Python randomises string hashing per process,
    # so `hash((check_name, layout_name))` gave a different seed on every run: measured
    # 90936 / 90550 / 90044 across three interpreters. That made an EXACTNESS test - one
    # whose whole claim is that a rate lands inside a band - silently flaky, and it is the
    # same defect `block_permutations` raises over when a caller omits an rng.
    seed = SEEDS[(check_name, layout_name)]
    p = _null_p_values(module, LAYOUTS[layout_name], seed, **kwargs)

    z_crit = bonferroni_z_crit(N_COMPARISONS)
    for alpha in ALPHAS:
        # `p <= alpha`, not `<`: alpha is a grid point, and the identity is about the mass
        # AT OR BELOW it. Using `<` would silently test the next grid point down.
        rate = float(np.mean(p <= alpha))
        band = z_crit * null_se(alpha, REPLICATES)
        assert abs(rate - alpha) <= band, (
            f"{check_name} {layout_name}: rejection {rate:.4f} at nominal {alpha} "
            f"(band +/-{band:.4f} = {z_crit:.2f} null SE, Bonferroni over "
            f"{N_COMPARISONS} comparisons)"
        )


def test_the_permutation_p_value_floor_is_exactly_one_over_b_plus_one():
    """The `+1` in `(1 + #{null >= obs})/(1 + B)`, tested where it is actually visible.

    MUTATION-DRIVEN. The exactness test above cannot see this term: dropping it gives
    `k/B <= alpha`, and at the grid alphas `k/199 <= 0.05` and `(1+k)/200 <= 0.05` both
    reduce to `k <= 9`. The rejection rates are IDENTICAL at 0.01, 0.05 and 0.10 - an
    arithmetic coincidence, not a weak test - so a mutant that removes the `+1` passes every
    assertion in this file's first half.

    Where it is visible is the floor. The observed statistic is one of the `B + 1` values
    being ranked, so it can never be strictly beaten by all of them: the smallest attainable
    p is `1/(B+1)`, and a rule without the `+1` can return exactly 0, which is not a
    p-value. Over 2000 null replicates the `k = 0` case occurs about ten times, so this sees
    the floor rather than assuming it.
    """
    p = _null_p_values(c6, [30], seed=4242)
    floor = 1.0 / (N_PERM + 1)
    assert p.min() == pytest.approx(floor, rel=1e-12), (
        f"smallest observed p {p.min():.6f} != the exact floor {floor:.6f}"
    )
    assert not np.any(p <= 0.0), (
        "a permutation p-value of 0 is arithmetic, not evidence"
    )


def test_a_permuted_statistic_that_TIES_the_observed_counts_against_rejection():
    """The `>=` in the same expression, which continuous data can never exercise.

    MUTATION-DRIVEN, the second half of the same finding. On exponential gaps no permuted
    statistic ever exactly equals the observed one, so replacing `>=` with `>` is a no-op
    across every replicate above and the mutant survives. It stops being a no-op the moment
    the durations are quantised - which is the case this project actually has, since the
    T2* ladder reads a discretised metric.

    Counting a tie AGAINST rejection is the conservative and correct choice: a permuted
    arrangement that reproduces the observed statistic is not evidence that the observed one
    is extreme.
    """
    observed = 2.0
    null = np.array([2.0, 2.0, 2.0, 1.0, 0.0])
    # 3 ties + the observed itself = 4 of 6; the strict form would report 1/6.
    assert permutation_p_value(observed, null) == pytest.approx(4.0 / 6.0)


# --------------------------------------------------------------------------------------
# The asymptotic instruments: MEASURED, direction asserted, magnitude reported.
# --------------------------------------------------------------------------------------


# The measurement lives in `analyzers/instrument_validation.py` and is IMPORTED here, not
# defined in `tests/`. The pipeline cannot import `tests/`, so a measurement defined there
# would leave the report's tier-3 verdicts resting on no number at all. Importing it in this
# direction is what makes "the figure and the test measure the same thing" true rather than
# asserted.


@pytest.mark.parametrize("check_name", ["c1", "c2", "cvm"])
def test_the_asymptotic_size_is_in_the_documented_range(check_name):
    """A wide bound, deliberately, and NO direction claim. See the module docstring.

    The bound catches a BROKEN limit - a mis-transcribed CDF, a wrong scaling - not a
    small-n cost. A tight band would encode one gap distribution's miss as correct, and the
    size depends on the gap distribution (the module docstring gives both), so a direction
    measured on one generator is not a property of the check.

    The magnitudes are not asserted here. They are measured into the artifact by
    `analyzers.instrument_validation.measure_all_asymptotic_sizes` and drawn by the report,
    which is the same function this file imports.
    """
    module = {"c1": c1, "c2": c2, "cvm": cvm}[check_name]
    size = measure_asymptotic_size(module, tau=30.0, seed=777)
    assert 0.005 <= size <= 0.20, f"{check_name} asymptotic size at tau=30: {size:.4f}"


class _SegmentRecorder:
    """A stand-in check that records the segments handed to it and rejects nothing.

    `measure_asymptotic_size` returns only a rate, so the null it simulates is otherwise
    unobservable. Passing this in is what puts the segments under assertion.
    """

    CHECK_NAME = "recorder"

    def __init__(self) -> None:
        self.segments: list[Segment] = []

    def run(self, segments, **kwargs):
        self.segments.extend(segments)
        return SimpleNamespace(p_value=None)


def test_the_size_generator_truncates_on_TIME_not_on_events():
    """The size generator truncates on time: tau fixed in advance, event count random.

    `checks/result` requires tau to be chosen without reference to the events, and nothing
    downstream detects a violation. Both observable consequences are asserted: the event
    count varies, and the fraction of the window used up by the last event varies. A tau
    derived from the draw pins the second at n/(n+1) and the first at n.

    Asserted on what `measure_asymptotic_size` actually builds, not on `_exponential_segment`
    directly: the helper is not what this commit changed, so a test of the helper stays green
    through a revert of the caller.
    """
    recorder = _SegmentRecorder()
    measure_asymptotic_size(recorder, tau=20.0, seed=777, replicates=40)

    assert len(recorder.segments) == 40
    counts = [len(s.x) for s in recorder.segments]
    used = [float(s.x.sum() / s.tau) for s in recorder.segments]

    assert len(set(counts)) > 1, (
        f"event count is fixed at {counts[0]}: tau is being derived from the events"
    )
    assert np.std(used) > 1e-6, (
        f"the leftover window is pinned at {used[0]:.9f}, so it carries no variance and "
        f"the measured rate is not a size"
    )
    assert all(0.0 < u <= 1.0 for u in used)
    assert all(s.tau == 20.0 for s in recorder.segments), (
        "tau must be the value asked for, identically on every replicate"
    )


def test_the_asymptotic_checks_reject_a_gross_trend_at_every_n_they_ship_at():
    """Power, not size: a check that never rejects would pass every test above."""
    for n in (20, 60, 200):
        x = np.linspace(1.0, 20.0, n)
        segment = Segment(x=x, tau=float(x.sum() + 5.0))
        for module in (c1, c2, cvm):
            result = module.run([segment], calibration=CALIB_ASYMPTOTIC)
            assert result.p_value is not None and result.p_value < 0.05, (
                f"{module.CHECK_NAME} missed a monotone ramp at n={n}: p={result.p_value}"
            )
