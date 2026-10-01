"""C1 - Lewis-Robinson test for trend, time-censored and multi-process.

Kvaloy & Lindqvist, "A class of tests for trend in time censored recurrent event data",
Technometrics 62(1):101-115 (2020). Single process, their eq (4), p. 103:

    LR = (1/gamma_hat) * sqrt(12)/(tau*sqrt(N)) * [ sum_i T_i - (N/2)*tau ]

Multi-process, their eq (13), p. 105, is `LR_m = sum_j A_j LR_j` for weights with
`sum_j a_j^2 = 1`. Two weightings exist and they are DIFFERENT statistics for m > 1:

    WEIGHTS_TECHNOMETRICS  journal eq (14)-(16), p. 105: A_j proportional to
        sqrt(N_j)/gamma_hat_j, giving
        LR_m = sqrt(12) / sqrt( sum_k N_k / gamma_hat_k^2 )
               * sum_j (1/gamma_hat_j^2) [ sum_i T_ij/tau_j - N_j/2 ]
    WEIGHTS_ARXIV_V1  preprint arXiv:1802.08339v1 eq (14)-(16), p. 8: A_j proportional to
        gamma_hat_j * tau_j * sqrt(N_j), giving
        LR_m = sqrt(12) / sqrt( sum_k gamma_hat_k^2 * tau_k^2 * N_k )
               * sum_j [ sum_i T_ij - (N_j/2) * tau_j ]

The journal weights are the published ones; the paper calls them optimal for power-law
trend alternatives `t^b` as `b -> 1` (p. 105, using Section 7.1.2), and they are the ones
behind its m = 6 example (Section 8.2, Table 4, LR p = 0.019). `tests/test_checks_published_values.py` reproduces
that p-value with the journal weights and shows the preprint weights do not. Both are
asymptotically N(0, 1) under the null for any weights converging to constants `a_j` with
`sum a_j^2 = 1` (eq (12), p. 105), so the preprint form is not wrong; it is a different,
superseded test.

Written once as `sqrt(12) sum_j c_j U_j / sqrt(sum_j c_j^2 gamma_j^2 tau_j^2 N_j)` with
`U_j = sum_i T_ij - N_j tau_j/2`, `c_j = 1` (preprint) or `c_j = 1/(tau_j gamma_j^2)`
(journal). Setting m = 1 gives eq (4) exactly for either `c_1`, so a separate
single-process path would be a second copy of the same formula, free to drift.
`tests/test_checks_statistics.py::test_eq16_reduces_to_eq4_for_a_single_segment` asserts the
identity numerically rather than trusting the algebra above.

What the statistic is: `sum_i T_i` is the total of the event times, and `(N/2)*tau` is its
expectation when events are uniform on `[0, tau]` (the null of no trend). So LR is a
centred, scaled "are the events early or late" contrast - negative for a decreasing
intensity, positive for an increasing one - and it is two-sided here. The `1/gamma_hat`
is what makes it a test of TREND rather than a test of the Poisson assumption: dividing by
the estimated coefficient of variation removes the renewal distribution's dispersion,
which is the whole difference between Lewis-Robinson and the plain Laplace test (Section
3.1, p. 103).

Both calibrations ship. Asymptotic is N(0,1) per the paper, and the paper's asymptotics are
for every tau_j growing (eq (12), p. 105), not for many short segments. Permutation is EXACTLY valid here,
not merely asymptotically: `tau_j`, `N_j` and `gamma_hat_j`, hence both weightings, are all
invariant to reordering gaps within a segment, so the only thing a permutation moves is
`sum_i T_ij`. That makes the asymptotic-vs-permutation gap in the bench a clean measurement
of asymptotic error rather than a comparison of two approximations.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from quebra.analyzers.checks._multiprocess import (
    GAMMA_DEFAULT,
    gamma_hat,
    gamma_hat_batch,
)
from quebra.analyzers.checks.result import (
    CALIB_ASYMPTOTIC,
    CALIB_PERMUTATION,
    CLOCK_IN_SPEC,
    CheckResult,
    Segment,
    concatenated_gaps,
    segment_sizes,
    validate_segment,
)
from quebra.analyzers.checks._permutation import (
    resolve_perm,
    DEFAULT_N_PERM,
    PermutationSet,
    check_permuted,
    two_sided_p_value,
)

CHECK_NAME = "c1_lewis_robinson"

_SQRT12 = float(np.sqrt(12.0))

# Multi-process weightings; see the module docstring. Recorded in `CheckResult.notes`.
WEIGHTS_TECHNOMETRICS = "technometrics_eq14"
WEIGHTS_ARXIV_V1 = "arxiv_v1_eq14"
C1_WEIGHTS = (WEIGHTS_TECHNOMETRICS, WEIGHTS_ARXIV_V1)
WEIGHTS_DEFAULT = WEIGHTS_TECHNOMETRICS


def _numerator_term(x: np.ndarray, tau: float) -> float:
    """`sum_i T_i - (N/2)*tau` for one segment."""
    return float(np.sum(np.cumsum(x)) - 0.5 * len(x) * tau)


def _numerator_term_batch(x_matrix: np.ndarray, tau: float) -> np.ndarray:
    n = x_matrix.shape[1]
    return np.cumsum(x_matrix, axis=1).sum(axis=1) - 0.5 * n * tau


def _coefficient(gamma, tau: float, weights: str):
    """`c_j` of the module docstring, scalar or per permutation row."""
    if weights == WEIGHTS_TECHNOMETRICS:
        return 1.0 / (tau * gamma**2)
    if weights == WEIGHTS_ARXIV_V1:
        return np.ones_like(gamma) if isinstance(gamma, np.ndarray) else 1.0
    raise ValueError(f"unknown C1 weights {weights!r}; known: {list(C1_WEIGHTS)}")


def statistic(
    segments: list[Segment],
    *,
    gamma_estimator: str = GAMMA_DEFAULT,
    weights: str = WEIGHTS_DEFAULT,
) -> float:
    """Eq (13) with the `weights` form of eq (16). Raises via `validate_segment`."""
    if not segments:
        raise ValueError("C1 needs at least one segment")
    numerator = 0.0
    denominator = 0.0
    for segment in segments:
        validate_segment(segment, require_strict_tau=False)
        x = np.asarray(segment.x, dtype=float)
        gamma = gamma_hat(x, segment.tau, gamma_estimator)
        c = _coefficient(gamma, segment.tau, weights)
        numerator += c * _numerator_term(x, segment.tau)
        denominator += c**2 * gamma**2 * segment.tau**2 * len(x)
    if denominator <= 0.0:
        raise ValueError("C1 denominator is not positive; every segment is degenerate")
    return _SQRT12 * numerator / float(np.sqrt(denominator))


def statistic_batch(
    segments: list[Segment],
    perm: PermutationSet,
    gamma_estimator: str = GAMMA_DEFAULT,
    permuted: np.ndarray | None = None,
    weights: str = WEIGHTS_DEFAULT,
) -> np.ndarray:
    """`statistic` for every permutation, vectorised over the `(B, total)` index matrix.

    `permuted` lets a caller supply the `(B, total)` gathered matrix it already built.
    At n = 355 that gather is a 355k-element copy and every permutation row wants the same one,
    so `checks/battery.py` builds it once and hands it round; passing None rebuilds it.
    """
    check_permuted(perm, permuted, "c1.statistic_batch")
    if permuted is None:
        permuted = perm.apply(concatenated_gaps(segments))
    numerator = np.zeros(perm.n_perm, dtype=float)
    denominator = np.zeros(perm.n_perm, dtype=float)
    for segment, (lo, hi) in zip(segments, perm.blocks()):
        block = permuted[:, lo:hi]
        gamma = gamma_hat_batch(block, segment.tau, gamma_estimator)
        c = _coefficient(gamma, segment.tau, weights)
        numerator += c * _numerator_term_batch(block, segment.tau)
        denominator += c**2 * gamma**2 * segment.tau**2 * (hi - lo)
    return _SQRT12 * numerator / np.sqrt(denominator)


def run(
    segments: list[Segment],
    *,
    calibration: str = CALIB_ASYMPTOTIC,
    clock: str = CLOCK_IN_SPEC,
    gamma_estimator: str = GAMMA_DEFAULT,
    perm: PermutationSet | None = None,
    n_perm: int = DEFAULT_N_PERM,
    rng: np.random.Generator | None = None,
    permuted: np.ndarray | None = None,
    weights: str = WEIGHTS_DEFAULT,
) -> CheckResult:
    observed = statistic(segments, gamma_estimator=gamma_estimator, weights=weights)
    n_events = int(sum(len(s.x) for s in segments))
    n_censored = int(sum(s.n_censored_dropped for s in segments))
    notes = f"gamma={gamma_estimator}"
    # The weighting only exists for m > 1; at m = 1 both forms are eq (4).
    if len(segments) > 1:
        notes += f" weights={weights}"

    if calibration == CALIB_ASYMPTOTIC:
        p_value = float(2.0 * stats.norm.sf(abs(observed)))
    elif calibration == CALIB_PERMUTATION:
        perm = resolve_perm(segment_sizes(segments), perm, n_perm, rng)
        null = statistic_batch(segments, perm, gamma_estimator, permuted, weights)
        p_value = two_sided_p_value(observed, null)
        notes += f" B={perm.n_perm}"
    else:
        raise ValueError(
            f"C1 supports {CALIB_ASYMPTOTIC!r} and {CALIB_PERMUTATION!r}; "
            f"got {calibration!r}"
        )

    return CheckResult(
        check=CHECK_NAME,
        statistic=observed,
        p_value=p_value,
        calibration=calibration,
        clock=clock,
        n_events=n_events,
        n_segments=len(segments),
        n_censored_dropped=n_censored,
        notes=notes,
    )
