"""C2 - Kvaloy-Lindqvist Anderson-Darling test for renewal, time-censored.

Kvaloy & Lindqvist, "A class of tests for trend in time censored recurrent event data",
Technometrics 62(1):101-115 (2020), eq (7), p. 104. Transcribed from the journal PDF:

    AD = (1/gamma_hat^2) * (1/N) * {
            sum_{i=1}^{N-1} [ (N-i)^2 * ln((tau - T_i)/(tau - T_{i+1}))
                              + i^2 * ln(T_{i+1}/T_i) ]
            + N^2 * [ ln(tau/(tau - T_1)) + ln(tau/T_N) - 1 ]
         }

What it is: with `u_i = T_i/tau`, eq (7) at `gamma_hat = 1` is exactly the classical
Anderson-Darling statistic `A^2 = N * integral (F_N(u) - u)^2 / (u(1-u)) du` testing the
event times for uniformity on `[0, tau]`.
`tests/test_checks_statistics.py::test_eq7_is_the_classical_anderson_darling` asserts that
against the `-N - (1/N) sum (2i-1)[ln u_i + ln(1 - u_{N+1-i})]` form of Anderson and
Darling (1954), JASA 49, eq (2), p. 765. That test
enforces `rel=1e-9` and measures 1e-16 to 1e-13 across its sizes, so the two are the same
statistic and the paper's form is shipped because it is the citation. The `1/gamma_hat^2` generalises it from Poisson to renewal.

**Two things about the calibration are decided by the paper, not by preference.**

1. *Asymptotic, m = 1.* The limit is the classical AD distribution (Section 3.4,
   pp. 103-104), evaluated with
   Marsaglia & Marsaglia (2004) `adinf` - their LIMITING function, deliberately not their
   finite-N correction `errfix`. That correction is derived for the AD statistic of N iid
   uniforms; here the statistic is additionally divided by an ESTIMATED `gamma_hat^2`, so
   the finite-N null is not the finite-N AD null and applying a correction fitted to the
   latter would be a precision claim the derivation does not support. Measuring that
   finite-N gap is what the bench is for.

2. *Multi-process, m > 1: permutation only.* The paper's own Section 4.2 (p. 106) says the
   normal approximation to a weighted sum works "less well for the Anderson-Darling test
   due to the very skew distribution", and its simulation study (Section 6, p. 106) leaves
   AD out for m > 1 in favour of Cramer-von Mises. The paper's own
   m = 6 example (Section 8.2, Table 4) does report AD, as a sum weighted by `tau_j` with a
   null simulated from the limiting law. Rather than ship the normal approximation its
   source finds works "less well", the summed statistic here is calibrated by within-segment permutation,
   which needs no asymptotic distribution at all and is exact under iid gaps. Asking for
   `asymptotic` with m > 1 raises.

The multi-process statistic is the unweighted SUM of per-segment eq (7) values. There is no
eq (16) analogue to borrow: eq (14)'s optimal weights are derived for the trend
statistic's normal limit, and the AD statistics being summed are not even asymptotically
normal, which is the paper's point. An unweighted sum makes no distributional claim and
the permutation calibration handles whatever distribution it has.

**A segment that ends at an event** (`tau == T_N`, `Segment.ends_at_event`) is failure
truncated, where eq (7) on all `N` events is `+inf`. It contributes eq (7) on its first
`n - 1` event times with `tau = T_n`; every log term is then finite, since
`tau - T_{n-1}` is the last gap. `gamma_hat` stays on all `n` gaps. Given `T_n` under a
constant rate the `n - 1` times are uniform order statistics on `(0, T_n)` (Rigdon and
Basu 2000, Sect. 4.3), so at `gamma_hat = 1` this is the classical A^2 of those `n - 1`
uniforms.

Validity assumptions

- Assumption: a segment that ends at an event stopped there by a rule that depends on its
  gaps only through their number and sum: a read gap or the record's end on the side
  clock, a birth reaching the stop on the calendar clock, the record's end for calibration
  cycles. Given `T_n`, its first `n - 1` event times are then what is read, and a
  permutation of all `n` gaps leaves `T_n` and `gamma_hat` unchanged, so the permutation
  p-value stays exact under iid gaps.
  Diagnostic: none from the data; the stopping rules are those of
  `_multiprocess.segments_from_windows` and `cycles.battery_segments`.
  Consequence of violation: a stop that depends on the order of the gaps makes the
  reshuffled records the wrong reference set, and the p-value loses its stated size.
  Reference: Rigdon and Basu (2000), Statistical Methods for the Reliability of Repairable
  Systems, Sect. 4.3, the failure-truncated case.
"""

from __future__ import annotations

import numpy as np

from quebra.analyzers.checks._multiprocess import (
    GAMMA_DEFAULT,
    gamma_hat,
    gamma_hat_batch,
)
from quebra.analyzers.checks._permutation import (
    resolve_perm,
    DEFAULT_N_PERM,
    PermutationSet,
    check_permuted,
    permutation_p_value,
)
from quebra.analyzers.checks.result import (
    CALIB_ASYMPTOTIC,
    CALIB_PERMUTATION,
    CLOCK_IN_SPEC,
    CheckResult,
    Segment,
    TAU_MARGIN,
    concatenated_gaps,
    events_read,
    segment_sizes,
    validate_segment,
)

CHECK_NAME = "c2_anderson_darling"


def ad_limiting_cdf(z: float | np.ndarray) -> np.ndarray:
    """`P(A^2 <= z)` for the LIMITING Anderson-Darling distribution.

    Marsaglia & Marsaglia, "Evaluating the Anderson-Darling Distribution", Journal of
    Statistical Software 9(2):1-5, 2004 - their `adinf` (p. 3). The paper quotes it to
    |error| < 2e-6 for z < 2 and < 8e-7 above. Measured against the exact limit
    (numerical inversion of its characteristic function,
    `scripts/verify_gold_standard.py`), the largest error is 1.95e-5, at z = 0.97, below
    z = 2, and 9.1e-6, at z = 2.59, above it: immaterial at alpha = 0.05, but the
    figures to quote. The 13 coefficients have not been compared with the paper's code
    archive, which is not in hand. Against the tabled critical values (Stephens 1974,
    JASA 69:730-737, Table 1A part 1.0, p. 732) it returns upper tails 0.10001 at z =
    1.933 and 0.05001 at z = 2.492. The same table's 3.070 and 3.857 are NOT the 2.5%
    and 1% points of the limit: this function gives 0.0252 and 0.0102 there, and
    Marsaglia & Marsaglia (p. 2) state that the published 3.857 "is actually 3.878125".
    The exact limiting points, evaluated from the series in Anderson & Darling (1954),
    JASA 49, eq (8), p. 768, are 3.0775 at 2.5% (no published source located) and 3.8781
    at 1% (Marsaglia & Marsaglia, p. 2).
    """
    z = np.asarray(z, dtype=float)
    # Finiteness FIRST: `nan <= 0.0` is False, so a nan statistic would slip past the
    # positivity guard and emerge as `p = 1 - nan = nan` on a CheckResult.
    if not np.all(np.isfinite(z)):
        n_bad = int((~np.isfinite(z)).sum())
        raise ValueError(
            f"{n_bad} of {z.size} Anderson-Darling statistics are not finite. This is a "
            "singular eq (7) and must be refused upstream, not turned into a p-value."
        )
    if np.any(z <= 0.0):
        raise ValueError("the Anderson-Darling statistic is positive; got z <= 0")
    small = z < 2.0
    out = np.empty_like(z)
    zs = np.where(small, z, 1.0)  # keep the unused branch away from the exp overflow
    out_small = (
        np.exp(-1.2337141 / zs)
        / np.sqrt(zs)
        * (
            2.00012
            + (
                0.247105
                - (0.0649821 - (0.0347962 - (0.011672 - 0.00168691 * zs) * zs) * zs)
                * zs
            )
            * zs
        )
    )
    zl = np.where(small, 2.0, z)
    out_large = np.exp(
        -np.exp(
            1.0776
            - (
                2.30695
                - (0.43424 - (0.082433 - (0.008056 - 0.0003146 * zl) * zl) * zl) * zl
            )
            * zl
        )
    )
    out = np.where(small, out_small, out_large)
    return out


def _eq7(x: np.ndarray, tau: float, gamma: float) -> float:
    T = np.cumsum(x)
    n = len(x)
    i = np.arange(1, n, dtype=float)  # i = 1 .. N-1
    head = np.sum(
        (n - i) ** 2 * np.log((tau - T[:-1]) / (tau - T[1:]))
        + i**2 * np.log(T[1:] / T[:-1])
    )
    tail = n**2 * (np.log(tau / (tau - T[0])) + np.log(tau / T[-1]) - 1.0)
    return float((head + tail) / (gamma**2 * n))


def _eq7_batch(x_matrix: np.ndarray, tau: float, gamma: np.ndarray) -> np.ndarray:
    T = np.cumsum(x_matrix, axis=1)
    n = x_matrix.shape[1]
    # Each row is cumsummed in its own permuted order, and float addition is not
    # associative, so `T[:, -1]` varies by ulps across rows even though the multiset is
    # identical. `TAU_MARGIN` is ~4500 ulps and therefore covers this, but the guard is
    # explicit: a row that reached `tau - T <= 0` would emit inf/nan into the null and
    # `permutation_p_value` would report it as "a bug in the statistic" without saying why.
    if np.any(tau - T[:, -1] <= 0.0):
        raise ValueError(
            "a permuted row reached tau <= T_N. The multiset is permutation-invariant, so "
            "this is float summation order, not a property of the data - it means tau "
            f"clears T_N by less than the {TAU_MARGIN:.0e} relative margin the guards use."
        )
    i = np.arange(1, n, dtype=float)[None, :]
    head = (
        (n - i) ** 2 * np.log((tau - T[:, :-1]) / (tau - T[:, 1:]))
        + i**2 * np.log(T[:, 1:] / T[:, :-1])
    ).sum(axis=1)
    tail = n**2 * (np.log(tau / (tau - T[:, 0])) + np.log(tau / T[:, -1]) - 1.0)
    return (head + tail) / (gamma**2 * n)


def statistic(
    segments: list[Segment], *, gamma_estimator: str = GAMMA_DEFAULT
) -> float:
    """Eq (7), summed unweighted over segments.

    Eq (7)'s `ln((tau - T_{N-1})/(tau - T_N))` term, with `T_0 = 0`, is `+inf` when the last
    event lands exactly on the truncation time. Such a segment is failure truncated and
    contributes eq (7) on its first `n - 1` event times (`events_read`), where the term is
    finite; `gamma_hat` reads all `n` gaps.
    """
    if not segments:
        raise ValueError("C2 needs at least one segment")
    total = 0.0
    for segment in segments:
        validate_segment(segment, require_strict_tau=False)
        x = np.asarray(segment.x, dtype=float)
        read = events_read(x, segment.ends_at_event)
        total += _eq7(read, segment.tau, gamma_hat(x, segment.tau, gamma_estimator))
    return total


def statistic_batch(
    segments: list[Segment],
    perm: PermutationSet,
    gamma_estimator: str = GAMMA_DEFAULT,
    permuted: np.ndarray | None = None,
) -> np.ndarray:
    """Eq (7) summed over segments, for every permutation.

    `permuted` is the shared `(B, total)` gathered matrix - see `c1.statistic_batch`.
    """
    check_permuted(perm, permuted, "c2.statistic_batch")
    if permuted is None:
        permuted = perm.apply(concatenated_gaps(segments))
    total = np.zeros(perm.n_perm, dtype=float)
    for segment, (lo, hi) in zip(segments, perm.blocks()):
        block = permuted[:, lo:hi]
        gamma = gamma_hat_batch(block, segment.tau, gamma_estimator)
        read = events_read(block, segment.ends_at_event)
        total += _eq7_batch(read, segment.tau, gamma)
    return total


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
) -> CheckResult:
    observed = statistic(segments, gamma_estimator=gamma_estimator)
    n_events = int(sum(s.n_events_read for s in segments))
    n_censored = int(sum(s.n_censored_dropped for s in segments))
    notes = f"gamma={gamma_estimator}"

    if calibration == CALIB_ASYMPTOTIC:
        if len(segments) > 1:
            raise ValueError(
                "C2 has no asymptotic calibration for m > 1 segments. Kvaloy & "
                "Lindqvist (Technometrics 2020) Section 4.2 find the normal "
                "approximation works 'less well' for the summed Anderson-Darling "
                "statistic ('very skew distribution'), and Section 6 drops AD for "
                "m > 1. Use calibration='permutation'."
            )
        p_value = float(1.0 - ad_limiting_cdf(observed))
        notes += " limiting_AD"
    elif calibration == CALIB_PERMUTATION:
        perm = resolve_perm(segment_sizes(segments), perm, n_perm, rng)
        null = statistic_batch(segments, perm, gamma_estimator, permuted)
        p_value = permutation_p_value(observed, null)
        notes += f" B={perm.n_perm}"
    else:
        raise ValueError(
            f"C2 supports {CALIB_ASYMPTOTIC!r} and {CALIB_PERMUTATION!r}; "
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
