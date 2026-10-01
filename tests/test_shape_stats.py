"""Shape statistics: xi, Spearman, distance correlation, tie fraction.

Written BEFORE the tie-corrected form (Chatterjee 2021) went into `chatterjee_xi`,
deliberately. Nothing else would
have caught that change at the time. `StaleArtifactGuard` compares field NAMES, so every
existing pickle in `output/` keeps loading through a change of this kind, and nothing else
asserts these values. `Job.job_code_hash` now covers `quebra.analyzers.shape_stats` - it is in
the import closure of every job that reaches these statistics - so an edit here does move the
run identity; that makes cached artifacts stale, which is a different signal from a wrong
number being caught.

The xi tests come in two groups, and the distinction is the point:

- **Tie-free**: the tie-corrected form reduces algebraically to `1 - 3*sum|dr|/(n^2-1)`, so both the old
  and the new estimator must agree, and both must match hand-computed values. Measured on
  the real T2* ladder, EVERY window is tie-free on both axes (0 of 279 at 3 us), so these
  are the tests that pin current behaviour.
- **Tied**: the two forms diverge, and only the tie-corrected form is correct. These are the
  tests that would fail against the tie-free reduction.

`eq8` in the test names below denotes the tie-corrected form.
"""

from __future__ import annotations

import numpy as np
import pytest

from quebra.analyzers.permutation import paired_permutation_test
from quebra.analyzers.shape_stats import (
    DCOR_CONSTANT,
    DCOR_OK,
    DCOR_TOO_SHORT,
    XI_METHOD_ASYMPTOTIC,
    XI_METHOD_NONE,
    XI_METHOD_PERMUTATION,
    chatterjee_xi,
    dcor,
    dcor_with_reason,
    tie_fraction,
    xi_p_value,
    xi_p_value_asymptotic,
)

pytestmark = pytest.mark.statistical


# ------------------------------------------------------------------ xi, tie-free


def _xi_tie_free_reference(x: np.ndarray, y: np.ndarray) -> float:
    """`1 - 3*sum|r_{i+1}-r_i|/(n^2-1)`, Chatterjee's tie-free reduction.

    Independent of the implementation: written from the formula, not by calling the code.
    The tie-corrected form must equal this whenever there are no ties in y.
    """
    from scipy.stats import rankdata

    n = len(x)
    ranks = rankdata(
        np.asarray(y)[np.argsort(np.asarray(x), kind="stable")], method="ordinal"
    )
    return float(1.0 - 3.0 * np.abs(np.diff(ranks)).sum() / (n * n - 1.0))


@pytest.mark.parametrize("n", [10, 20, 100, 1000])
def test_a_perfect_function_hits_the_exact_finite_n_maximum(n):
    """`xi -> 1` for a functional relation only ASYMPTOTICALLY. At finite n the maximum
    attainable value is `1 - 3/(n+1)`, because a perfectly ordered rank vector still has
    `sum|r_{i+1} - r_i| = n - 1` and `1 - 3(n-1)/(n^2-1) = 1 - 3/(n+1)`.

    Pinning the exact ceiling rather than "close to 1" is what makes this a real test: at
    n = 20 the ceiling is 0.8571, and an implementation that genuinely reached 1.0 there
    would be wrong. This project's windows run from 3 to a few dozen reads, so the ceiling
    is well below 1 in practice and any reading of xi has to allow for it.

    TIE-FREE ONLY. `y = 2x + 1` has all-distinct values by construction. With ties in y the
    ceiling is HIGHER, because the tie-corrected denominator shrinks - measured, n = 12 with two
    distinct y values reaches 0.833 against a tie-free ceiling of 0.769.
    """
    x = np.arange(1.0, n + 1.0)
    assert chatterjee_xi(x, 2.0 * x + 1.0) == pytest.approx(
        1.0 - 3.0 / (n + 1.0), rel=1e-12
    )


def test_xi_is_near_one_for_a_non_monotone_function_where_spearman_is_not():
    """The whole reason xi ships beside Spearman: a symmetric hump has rho ~ 0 by
    construction, whatever the predictability, while xi still sees the structure."""
    from scipy.stats import spearmanr

    x = np.linspace(-1.0, 1.0, 40)
    y = x**2
    assert abs(float(spearmanr(x, y).statistic)) < 0.1
    assert chatterjee_xi(x, y) > 0.5


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_xi_matches_the_tie_free_reduction_when_there_are_no_ties(seed):
    """Pins current behaviour on the regime this project's data actually occupies."""
    rng = np.random.default_rng(seed)
    x = rng.normal(size=30)
    y = rng.normal(size=30)
    assert len(np.unique(y)) == len(y), "fixture must be tie-free for this comparison"
    assert chatterjee_xi(x, y) == pytest.approx(_xi_tie_free_reference(x, y), rel=1e-12)


def test_xi_is_nan_on_degenerate_input():
    assert np.isnan(chatterjee_xi(np.array([1.0, 2.0]), np.array([1.0, 2.0])))
    assert np.isnan(chatterjee_xi(np.ones(10), np.arange(10.0)))


# ---------------------------------------------------------------------- xi, tied


def _xi_eq8_reference(x: np.ndarray, y: np.ndarray) -> float:
    """Chatterjee (2021), the tie-corrected form, written from the formula rather than the code.

        xi = 1 - n * sum_i |r_{i+1} - r_i| / (2 * sum_i l_i (n - l_i))

    with `r_i = #{j : y_j <= y_(i)}` and `l_i = #{j : y_j >= y_(i)}` after sorting the
    pairs by x. This is the tie-corrected form; it reduces to the expression in
    `_xi_tie_free_reference` when y has no ties.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    order = np.argsort(x, kind="stable")
    ys = y[order]
    r = np.array([np.sum(y <= v) for v in ys], dtype=float)
    ell = np.array([np.sum(y >= v) for v in ys], dtype=float)
    return float(1.0 - n * np.abs(np.diff(r)).sum() / (2.0 * np.sum(ell * (n - ell))))


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_eq8_and_the_tie_free_reduction_agree_when_y_has_no_ties(seed):
    """The algebraic identity the whole tie-free/tie-corrected distinction rests on.

    If this fails, the two reference implementations in this file disagree about a case
    where they are provably the same, and nothing else in the file can be trusted.
    """
    rng = np.random.default_rng(100 + seed)
    x, y = rng.normal(size=25), rng.normal(size=25)
    assert _xi_eq8_reference(x, y) == pytest.approx(
        _xi_tie_free_reference(x, y), rel=1e-12
    )


@pytest.mark.parametrize("seed", [0, 1, 2, 3])
def test_the_implementation_matches_eq8_on_TIED_data(seed):
    """The one behaviour this increment changed, pinned against the formula.

    Everything else here compares the implementation to the tie-FREE form, which it also
    satisfied by the tie-free reduction. Without this test the tied path - the entire point of
    correcting the estimator - is unpinned, and an analyzer edit would not move the run
    identity to signal it.
    """
    rng = np.random.default_rng(200 + seed)
    x = rng.normal(size=35)
    y = rng.integers(0, 4, size=35).astype(float)  # heavy ties in y
    assert tie_fraction(y) > 0.8
    assert chatterjee_xi(x, y) == pytest.approx(_xi_eq8_reference(x, y), rel=1e-12)


def test_the_two_forms_diverge_once_y_is_tied():
    """The failure the tie-corrected form exists to remove.

    A heavily tied y is what a quantised metric produces - a fidelity ladder where most
    windows are one read long. The T2* ladder has no ties at all, which is why this had to
    be constructed rather than taken from the data.
    """
    rng = np.random.default_rng(7)
    x = rng.normal(size=40)
    y = np.repeat(np.arange(4.0), 10)  # four distinct values, 36 tied observations
    rng.shuffle(y)
    assert tie_fraction(y) == pytest.approx(1.0)
    assert _xi_eq8_reference(x, y) != pytest.approx(
        _xi_tie_free_reference(x, y), rel=1e-6
    )


# ------------------------------------------------------------------ tie_fraction


def test_tie_fraction_counts_observations_not_values():
    assert tie_fraction(np.array([1.0, 1.0, 2.0, 3.0])) == pytest.approx(0.5)
    assert tie_fraction(np.arange(10.0)) == pytest.approx(0.0)
    assert tie_fraction(np.ones(5)) == pytest.approx(1.0)
    assert np.isnan(tie_fraction(np.array([1.0])))


# --------------------------------------------------------------- distance correlation


def test_dcor_is_one_for_an_exact_linear_relation():
    x = np.arange(1.0, 21.0)
    assert dcor(x, 3.0 * x - 2.0) == pytest.approx(1.0, abs=1e-12)


def test_dcor_sees_a_non_monotone_relation_that_spearman_misses():
    from scipy.stats import spearmanr

    x = np.linspace(-1.0, 1.0, 40)
    y = x**2
    assert abs(float(spearmanr(x, y).statistic)) < 0.1
    assert dcor(x, y) > 0.4


def test_dcor_is_symmetric_and_bounded():
    rng = np.random.default_rng(3)
    x, y = rng.normal(size=30), rng.normal(size=30)
    value = dcor(x, y)
    assert dcor(y, x) == pytest.approx(value, rel=1e-12)
    assert 0.0 <= value <= 1.0


def test_dcor_reason_codes_distinguish_their_causes():
    """A NaN attributed to the wrong cause is worse than no reason at all - the reason is
    part of a typed contract the panel reads."""
    rng = np.random.default_rng(5)
    assert dcor_with_reason(rng.normal(size=20), rng.normal(size=20))[1] == DCOR_OK
    assert dcor_with_reason(np.arange(3.0), np.arange(3.0))[1] == DCOR_TOO_SHORT
    assert dcor_with_reason(np.ones(10), np.arange(10.0))[1] == DCOR_CONSTANT
    # Non-finite input must be reported as too_short after masking, not as constant.
    x = np.array([1.0, 2.0, np.nan, np.nan, np.nan])
    assert dcor_with_reason(x, np.arange(5.0))[1] == DCOR_TOO_SHORT


# ------------------------------------------------------------- the calibration branch


def test_a_tie_free_window_uses_the_closed_form():
    """The path every T2* window takes - measured, 0 of 279 windows tied."""
    rng = np.random.default_rng(11)
    x, y = rng.normal(size=25), rng.normal(size=25)
    p, method = xi_p_value(x, y, seed=1)
    assert method == XI_METHOD_ASYMPTOTIC
    assert p == pytest.approx(xi_p_value_asymptotic(chatterjee_xi(x, y), 25), rel=1e-12)


def test_a_tied_window_switches_to_permutation():
    """Design rule, `.claude/qre_checks_reference.tex` Section 10.3 caveat 1: where ties
    dominate, the closed form must not be used."""
    rng = np.random.default_rng(12)
    x = rng.normal(size=30)
    y = rng.integers(0, 3, size=30).astype(float)
    p, method = xi_p_value(x, y, seed=1)
    assert method == XI_METHOD_PERMUTATION
    assert 0.0 < p <= 1.0


def test_the_reported_method_is_never_a_lie():
    """Whatever the branch, the label must name the calibration that produced the p."""
    rng = np.random.default_rng(13)
    for y in (rng.normal(size=20), rng.integers(0, 3, size=20).astype(float)):
        x = rng.normal(size=20)
        p, method = xi_p_value(x, y, seed=3)
        if method == XI_METHOD_ASYMPTOTIC:
            assert p == pytest.approx(
                xi_p_value_asymptotic(chatterjee_xi(x, y), 20), rel=1e-12
            )
        else:
            assert method == XI_METHOD_PERMUTATION
            assert tie_fraction(x) > 0 or tie_fraction(y) > 0


def test_a_degenerate_window_reports_no_method_rather_than_a_number():
    p, method = xi_p_value(np.ones(10), np.arange(10.0), seed=1)
    assert np.isnan(p) and method == XI_METHOD_NONE


def test_the_constant_y_guard_fires():
    """New with the tie-corrected form: a constant y makes the denominator exactly zero."""
    assert np.isnan(chatterjee_xi(np.arange(10.0), np.ones(10)))


def test_the_xi_permutation_is_reproducible_and_seed_dependent():
    rng = np.random.default_rng(14)
    x, y = rng.normal(size=30), rng.integers(0, 3, size=30).astype(float)
    a, _ = xi_p_value(x, y, seed=5)
    b, _ = xi_p_value(x, y, seed=5)
    c, _ = xi_p_value(x, y, seed=6)
    assert a == b, "same seed must give the same p-value"
    assert a != c, "a different seed must give a different draw"


# ------------------------------------------------------- the shared permutation step


def test_paired_permutation_refuses_a_statistic_that_goes_non_finite():
    """Substituting a sentinel would bias the p-value; -inf is anti-conservative."""
    rng = np.random.default_rng(15)
    with pytest.raises(ValueError, match="non-finite value on a permuted resample"):
        paired_permutation_test(
            rng.normal(size=10),
            rng.normal(size=10),
            lambda a, b: float("nan"),
            seed=1,
        )


@pytest.mark.parametrize("n, expect_exact", [(3, True), (4, True), (6, False)])
def test_small_windows_get_an_enumerated_exact_null(n, expect_exact):
    """scipy enumerates when `(n!)^2 <= n_resamples`, so n <= 4 is EXACT at the default B.

    A three-read window is the modal case here, so this is the common path. The p-value
    floor differs between the two regimes - `1/B` enumerated, `1/(B+1)` sampled - and
    `resolution` has to report the one that applied.
    """
    rng = np.random.default_rng(16)
    r = paired_permutation_test(
        np.arange(float(n)), rng.normal(size=n), chatterjee_xi, seed=1
    )
    assert r.exact is expect_exact
    expected_floor = 1.0 / r.n_resamples if expect_exact else 1.0 / (r.n_resamples + 1)
    assert r.resolution == pytest.approx(expected_floor)


def test_paired_permutation_detects_a_real_dependence():
    x = np.linspace(0.0, 1.0, 40)
    r = paired_permutation_test(x, x**2, chatterjee_xi, seed=1)
    assert r.p_value <= r.resolution + 1e-12


# --------------------------------------------------- the closed-form null constant


def test_the_asymptotic_null_variance_is_two_fifths():
    """Chatterjee (2021) Thm 2.2: under independence `sqrt(n) xi -> N(0, 2/5)`.

    ADDED AFTER MUTATION TESTING. Replacing `2.0 / 5.0` with `1.0 / 5.0` left all 35 tests
    of this module passing - the constant that scales EVERY p-value the tool actually emits
    was unpinned. It is the asymptotic branch that runs in production: all 279 windows of
    the real T2* ladder are tie-free, so `xi_p_method` reads `asymptotic` on every one of
    them, and `xi_p_value_asymptotic` is where their significance is decided.

    Simulated rather than restated. Asserting `p == norm.sf(sqrt(n) xi / sqrt(2/5))` would
    be the same arithmetic twice and would survive the mutation; measuring the variance of
    the statistic under its own null is an independent route to the constant.
    """
    n, replicates = 200, 1500
    rng = np.random.default_rng(11)
    xi_values = np.array(
        [
            chatterjee_xi(rng.normal(size=n), rng.normal(size=n))
            for _ in range(replicates)
        ]
    )
    # SE of a variance estimate over 1500 draws is about 0.015, so this band is ~3 SE and
    # the 1/5 mutant sits 13 SE outside it.
    assert float(np.var(np.sqrt(n) * xi_values)) == pytest.approx(0.4, abs=0.05)


def test_the_closed_form_holds_its_nominal_level_under_independence():
    """The operational form of the same statement: does the p-value mean what it says?

    `.claude/qre_checks_reference.tex` Section 10.3 caveat 2 calls n of 35 to 355, the
    range this project runs in, disputed for the normal approximation (no source located
    for the dispute), so this is measured rather than assumed. Type-I comes out
    at or slightly below nominal here - the safe direction - which is also what the
    ledger's OPEN entry on the missing n-floor relies on.
    """
    n, replicates = 60, 1500
    rng = np.random.default_rng(12)
    p_values = np.array(
        [
            xi_p_value_asymptotic(
                chatterjee_xi(rng.normal(size=n), rng.normal(size=n)), n
            )
            for _ in range(replicates)
        ]
    )
    for alpha in (0.05, 0.10):
        # 3 binomial SE, one-sided upward: being conservative is not a failure here.
        se = float(np.sqrt(alpha * (1 - alpha) / replicates))
        assert float(np.mean(p_values <= alpha)) <= alpha + 3 * se
