"""Cramer-von Mises: the fourth functional of the same bridge as C1 and C2.

Two independent pins, which is what makes this more than a second transcription:

- The STATISTIC at `gamma = 1` must equal the classical one-sample Cramer-von Mises `W^2`
  on `u_i = T_i/tau`, computed by `scipy.stats.cramervonmises`. That is an external
  implementation of the same mathematical object, arrived at by a different route.
- The limiting CDF must reproduce the critical values Anderson and Darling (1952) print in
  their Table 1, p. 203.

Kvaloy and Lindqvist (Technometrics 62(1):101-115, 2020) print eq (6) with a `/tau` where
the derivation gives `/tau^2`, in both Technometrics (p. 103) and arXiv:1802.08339v1 (p. 5);
`test_the_bracket_is_the_bridge_integral` is what settles that, by integrating the
definition numerically instead of trusting either version.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import cramervonmises

import quebra.analyzers.checks.c2_anderson_darling as c2
import quebra.analyzers.checks.cvm_cramer_von_mises as cvm
from quebra.analyzers.checks._permutation import PermutationSet, block_permutations
from quebra.analyzers.checks.result import CALIB_ASYMPTOTIC, CALIB_PERMUTATION, Segment

pytestmark = pytest.mark.statistical


def _segment(n: int, seed: int, slack: float = 1.0) -> Segment:
    g = np.random.default_rng(seed).exponential(size=n)
    return Segment(x=g, tau=float(np.cumsum(g)[-1] + slack * g.mean()))


@pytest.mark.parametrize("n, seed", [(13, 3), (36, 7), (120, 11)])
def test_at_gamma_one_it_is_the_classical_cramer_von_mises(n, seed):
    """The external cross-check. scipy computes W^2 by a different route entirely."""
    seg = _segment(n, seed)
    ours = cvm._cvm(seg.x, seg.tau, 1.0)
    theirs = float(cramervonmises(np.cumsum(seg.x) / seg.tau, "uniform").statistic)
    assert ours == pytest.approx(theirs, rel=1e-12)


@pytest.mark.parametrize("n, seed", [(9, 1), (25, 4)])
def test_the_bracket_is_the_bridge_integral(n, seed):
    """Settles the paper's printed `/tau` versus `/tau^2` by integrating the definition.

    `integral_0^1 (N(s tau) - sN)^2 ds` computed numerically on the step function, against
    the closed form the implementation uses. A `/tau` middle term does not reproduce this.
    """
    seg = _segment(n, seed)
    u = np.cumsum(seg.x) / seg.tau
    integrand = lambda s: (float(np.sum(u <= s)) - s * n) ** 2  # noqa: E731
    # Integrate between the jumps: the counting process is a step function, and handing
    # quad its discontinuities is what makes this an accurate check rather than a noisy one.
    knots = np.concatenate([[0.0], u, [1.0]])
    numeric = sum(
        quad(integrand, lo, hi)[0] for lo, hi in zip(knots[:-1], knots[1:]) if hi > lo
    )
    assert cvm._cvm(seg.x, seg.tau, 1.0) == pytest.approx(numeric / n, rel=1e-9)


def test_the_limiting_cdf_reproduces_the_published_critical_values():
    """Oracle: Anderson and Darling (1952), Annals of Math. Stat. 23, Table 1, p. 203.

    The four points asserted are that table's rows a1(z) = .90, .95, .99 and .999; it has no
    .975 row. The 2.5%
    point is therefore not asserted here as a published value. The series matches the
    printed points to within 3.1e-6 in tail probability (the table's own 5-decimal
    rounding of z), so the tolerance is set at 1e-5 rather than a looser band that would
    also accept a wrong series.
    """
    for statistic, alpha in [
        (0.34730, 0.10),
        (0.46136, 0.05),
        (0.74346, 0.01),
        (1.16786, 0.001),
    ]:
        assert 1.0 - cvm.cvm_limiting_cdf(statistic) == pytest.approx(alpha, abs=1e-5)


def test_the_limiting_cdf_is_monotone_and_bounded():
    z = np.linspace(0.02, 3.0, 200)
    values = cvm.cvm_limiting_cdf(z)
    assert np.all(np.diff(values) >= -1e-12)
    assert values[0] >= 0.0 and values[-1] <= 1.0


def test_the_limiting_cdf_refuses_a_non_finite_statistic():
    with pytest.raises(ValueError, match="must be finite"):
        cvm.cvm_limiting_cdf(np.array([np.nan]))


@pytest.mark.parametrize("sizes", [[6], [4, 5], [3, 3, 4]])
def test_statistic_batch_reproduces_the_statistic_under_the_identity(sizes):
    """Guards the `segments`/`perm.blocks()` zip, as the C1 and C2 tests do."""
    segments = [_segment(size * 5, 40 + size) for size in sizes]
    total = sum(s.n_events for s in segments)
    identity = PermutationSet(
        sizes=tuple(s.n_events for s in segments),
        indices=np.tile(np.arange(total), (3, 1)),
    )
    assert cvm.statistic_batch(segments, identity) == pytest.approx(
        cvm.statistic(segments), rel=0, abs=1e-12
    )


def test_it_survives_a_record_that_c2_declines():
    """CvM is well posed on a record that C2 declines.

    This is not the source paper's reason for preferring CvM at m > 1. Kvaloy and Lindqvist
    (Technometrics 2020, Section 4.2, p. 106) give level properties: the normal
    approximation works "less well for the Anderson-Darling test due to the very skew
    distribution".

    The AD integrand carries a `1/(s(1-s))` weight, so the `i = N-1` term of eq (7),
    `ln((tau - T_{N-1})/(tau - T_N))`, is `+inf` when the last event lands on the
    truncation time. CvM's integrand has no such weight, so the same
    record is perfectly well posed for it. On the in-spec clock of a carved record that
    case is the common one, not a corner.
    """
    g = np.random.default_rng(5).exponential(size=20)
    tight = Segment(x=g, tau=float(np.cumsum(g)[-1]))
    with pytest.raises(ValueError):
        c2.statistic([tight])
    assert np.isfinite(cvm.statistic([tight]))


def test_it_detects_a_strong_trend():
    """Events crowded early are a gross departure from uniformity on [0, tau]."""
    early = np.diff(
        np.concatenate([[0.0], np.sort(np.random.default_rng(2).uniform(0, 10, 30))])
    )
    seg = Segment(x=early, tau=100.0)
    assert cvm.run([seg], calibration=CALIB_ASYMPTOTIC).p_value < 0.01


def test_the_permutation_calibration_runs_and_is_seeded():
    seg = _segment(30, 9)
    perm = block_permutations([seg.n_events], 199, np.random.default_rng(1))
    a = cvm.run([seg], calibration=CALIB_PERMUTATION, perm=perm)
    b = cvm.run([seg], calibration=CALIB_PERMUTATION, perm=perm)
    assert a.p_value == b.p_value
    assert 0.0 < a.p_value <= 1.0


def test_the_asymptotic_calibration_is_refused_for_multiple_segments():
    segments = [_segment(15, 20 + i) for i in range(3)]
    with pytest.raises(ValueError, match="no asymptotic calibration"):
        cvm.run(segments, calibration=CALIB_ASYMPTOTIC)


def test_cvm_is_registered_in_the_bench_row_schema():
    """CvM is promoted; this test pins the promoted behaviour and is inverted, not
    deleted, so the change of decision is visible in the history rather than silent.

    Registration alone is not the promotion. `ROW_KEYS` is the schema of the BENCH TABLES,
    so a registered check with no bench rows makes `bench_acceptance_at_n` return None and
    every CvM ledger cell reads `underpowered / no bench cell`. The second assertion is
    therefore the one that matters: the bench must actually carry CvM size rows.
    """
    from pathlib import Path

    import pandas as pd

    from quebra.analyzers.checks.battery import ROW_KEYS

    assert any(key[0] == cvm.CHECK_NAME for key in ROW_KEYS)

    table = (
        Path(__file__).resolve().parents[1]
        / "jobs"
        / "bench"
        / "results"
        / "size_table.csv"
    )
    rows = pd.read_csv(table)
    assert (rows["check"] == cvm.CHECK_NAME).any(), (
        "CvM is registered in ROW_KEYS but jobs/bench/results/size_table.csv has no CvM rows, "
        "so every CvM verdict will read `no bench cell`. Re-run `python jobs/bench/runner.py`."
    )


# ----------------------------------------- the tail, where the series stops converging


def test_the_p_value_is_monotone_across_the_whole_range_including_the_tail():
    """Found by a cold reviewer, and invisible to every other test in this file.

    The four pinned critical values all sit in [0.34, 0.75], so nothing here looked above
    `z = 1`. Out there the 12-term series stops converging - `w = (4j+1)^2/(16z)` goes to
    zero, `K_{1/4}` diverges, and the summed value turns around. Before the fix the p-value
    bottomed out at `z = 7.66` and CLIMBED BACK: p(100) = 3.8e-3, p(500) = 0.100,
    p(1e5) = 0.705. A p-value that rises with the statistic is not a precision problem, it
    is a verdict flip - the reviewer built the case: 682 near-constant durations with a step
    trend give a statistic of 794.4, which returned p = 0.156 and FAILED TO REJECT.
    """
    z = np.concatenate([np.linspace(0.02, 5.0, 400), np.geomspace(5.0, 1e6, 200)])
    p = 1.0 - cvm.cvm_limiting_cdf(z)
    assert np.all(np.diff(p) <= 1e-12), (
        "the p-value must never increase with the statistic"
    )


def test_a_saturated_p_value_says_so_in_its_notes():
    """Oracle: the module's own saturation policy at `_SERIES_Z_MAX`.

    Above the validated range the series is not evaluated at all: `cvm_limiting_cdf`
    replaces `z` with 1.0 and returns a hardcoded 1.0, so `p = 0.0` arrives by GUARD
    CLAUSE, not by computation. The verdict is right and the magnitude is not a
    measurement, so the result must say which it is rather than let a reader quote the
    zero as a number.

    This replaces a test that asserted only `1 - cvm_limiting_cdf(794.4) == 0.0`. That
    assertion is implied by the guard clause for every input above the cutoff, it names a
    reviewer rather than an oracle, and no mutation was found that it caught and
    `test_the_p_value_is_monotone_across_the_whole_range_including_the_tail` did not.
    What nothing asserted, and this does, is that the note is attached.
    """
    x = np.full(682, 1.0)
    x[341:] = 3.0  # a step trend: large statistic, comfortably past the cutoff
    seg = Segment(x=x, tau=float(x.sum() + x.mean()))
    result = cvm.run([seg], calibration=CALIB_ASYMPTOTIC)

    assert result.statistic > cvm._SERIES_Z_MAX
    assert result.p_value == 0.0
    assert f"p_saturated(z>{cvm._SERIES_Z_MAX:g})" in result.notes


def test_an_unsaturated_p_value_carries_no_saturation_note():
    """Negative control for the note above: below the cutoff it must be absent.

    Without this, a change that attached the note unconditionally would pass the test
    above while making the annotation meaningless.
    """
    rng = np.random.default_rng(20260906)
    x = rng.exponential(size=200)
    seg = Segment(x=x, tau=float(x.sum() + x.mean()))
    result = cvm.run([seg], calibration=CALIB_ASYMPTOTIC)

    assert result.statistic < cvm._SERIES_Z_MAX
    assert "p_saturated" not in result.notes


def test_the_series_is_still_exact_where_it_was_pinned():
    """The saturation must not have moved anything inside the validated range."""
    for statistic, alpha in [(0.34730, 0.10), (0.46136, 0.05), (0.74346, 0.01)]:
        assert 1.0 - cvm.cvm_limiting_cdf(statistic) == pytest.approx(alpha, abs=5e-4)


def test_the_cutoff_agrees_with_an_independent_representation():
    """Karhunen-Loeve: `W^2 = sum_k Z_k^2/(k pi)^2` for iid standard normal `Z_k`.

    A genuinely different construction of the same limit - eigenvalues of the CvM kernel
    rather than the Anderson-Darling series - so agreement is evidence rather than a
    restatement. Kept cheap; the tolerance is Monte Carlo error at this draw count.
    """
    rng = np.random.default_rng(0)
    lam = 1.0 / ((np.arange(1, 61) * np.pi) ** 2)
    w2 = (rng.standard_normal((200_000, 60)) ** 2) @ lam
    for z, tol in [(0.35, 0.003), (0.46, 0.002), (0.74, 0.001)]:
        assert float(np.mean(w2 > z)) == pytest.approx(
            1.0 - float(cvm.cvm_limiting_cdf(z)), abs=tol
        )
