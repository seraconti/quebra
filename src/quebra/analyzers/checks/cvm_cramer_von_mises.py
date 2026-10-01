"""CvM - Cramer-von Mises type test for renewal, time censored.

Another functional of the same Brownian bridge that gives C1 and C2. Kvaloy and
Lindqvist, "A class of tests for trend in time censored recurrent event data",
Technometrics 62(1):101-115 (2020), Section 3, derive a CLASS of tests by applying
different functionals to the tied-down normalised counting process;
`integral W0(s)^2 ds` is this one (Section 3.3, eq (6)), `integral W0(s)^2/(s(1-s)) ds`
is C2, `integral W0(s) ds` is C1.

**Why it is here at all.** For m > 1 the paper REVERSES its single-process preference: its
simulation study (Section 6, p. 106) leaves
Anderson-Darling out for several processes "as the Cramer-von Mises test had better level
properties in this case". Gapped records are routine here and m > 1 is the gap case, so CvM
is the functional the paper itself keeps for m > 1 (Section 6, p. 106), the regime this
project actually operates in. The
paper's m > 1 CvM and AD statistics are sums weighted by `tau_j` with a simulated null
(Section 8.2, Table 4); ours are unweighted sums with a permutation null. The paper
names permutation as an alternative for very small samples in its conclusion (Section 9,
p. 113).

**Derived, not transcribed.** The paper prints its eq (6) middle term as
`- i N (T^2_{i+1} - T^2_i)/tau`, which is dimensionally inconsistent with the first term
`i^2 X_{i+1}/tau`. The typo is in the source, in both versions (Technometrics eq (6),
p. 103; arXiv v1 eq (6), p. 5); the authors' R code (`CvMtestobs` in
github.com/jtkgithub/trendtests) divides by `tau^2`. Integrating from the definition gives
`/tau^2`:

    integral_0^1 (N(s tau) - sN)^2 ds
        = sum_{i=0}^{N-1} [ i^2 (u_{i+1} - u_i) - i N (u^2_{i+1} - u^2_i) ]
          + N^2 [ u_N^2 - u_N + 1/3 ],        u_i = T_i / tau, u_0 = 0

which reproduces the paper's own tail term exactly. Verified numerically against
`scipy.integrate.quad` on the step function and against `scipy.stats.cramervonmises`.

**The scaling.** The normalised bridge is `V0(s) = (N(s tau) - sN) / (gamma_hat sqrt(N))`,
so `integral V0^2 = (1/(gamma_hat^2 N)) integral (N(s tau) - sN)^2`, which is the
`(1/gamma_hat^2)(1/N)` prefactor of eq (6).

At `gamma_hat = 1` the result is EXACTLY the classical one-sample Cramer-von Mises
statistic `W^2` on `u_i = T_i/tau` - agreement to 12 decimals against
`scipy.stats.cramervonmises`, which is the external cross-check and is stronger evidence
than a second transcription of the same formula would be.

**Registered in `battery.ROW_KEYS`.** That tuple is the schema of the bench
tables, so registration alone would have made `bench_acceptance_at_n` return None and every
CvM ledger row read `underpowered / no bench cell`. `jobs/bench/results/size_table.csv`
carries 219 CvM size rows and the power table 348. Over the twelve arm=A_iid_weibull /
in_spec / unquantised / uncensored asymptotic cells there, mean |size - 0.05| is 0.0055 for
CvM, against 0.0051 for C1 and 0.0058 for C2: margins smaller than the Monte Carlo error of
a twelve-cell mean, about 0.0014, so this aggregate does not rank the three. Single cells and
the envelope's worst cells do differ.

Its PERMUTATION row runs even when `include_tau_checks` is False, which is the reason it was
promoted: `tau == T_N` on the in-spec clock of a carved record silences C1 and C2, and this
integrand carries no `1/(s(1-s))` weight. Its ASYMPTOTIC row stays inside that gate, because
a limiting null still needs a truncation time chosen independently of the events.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import gamma as gamma_fn
from scipy.special import kv

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
    concatenated_gaps,
    segment_sizes,
    validate_segment,
)

CHECK_NAME = "cvm_cramer_von_mises"

# Terms of the limiting series. Convergence is fast NEAR THE CRITICAL VALUES and slower at
# the top of the retained range: the j = 3 term is 8.4e-28 at z = 0.35, 9.2e-11 at z = 1, but
# 1.9e-03 at z = 5. (An earlier comment here claimed "below 1e-12 across the range that
# matters", which is nine orders out at the cutoff.) 12 is still generous - the summed value
# is identical to a 40- and an 80-term sum to every printed digit throughout z <= 5 - but the
# margin comes from the later terms, not from the third.
_CDF_TERMS = 12

# Above this the series must NOT be evaluated. Each term carries `exp(-w) K_{1/4}(w)` with
# `w = (4j+1)^2/(16 z)`, so `w -> 0` as `z` grows, `K_{1/4}` diverges, and the truncated sum
# stops converging: measured, `p` bottoms out at `z = 7.66` and then CLIMBS BACK, giving
# p(100) = 3.8e-3, p(500) = 0.100, p(1e5) = 0.705. That is a p-value that fails to reject on
# an overwhelming statistic - a verdict flip, not a precision loss. `scipy`'s own
# `_cdf_cvm_inf` breaks the same way (in scipy 1.17.1 its p rises from about z = 4.75), so
# delegating does not help.
#
# 5.0 is where the series is still monotone and the true survival is already negligible.
# Evidence, all measured rather than assumed: the four pinned critical values all lie in
# [0.34, 1.17] and reproduce to 3.1e-06; an INDEPENDENT Karhunen-Loeve Monte Carlo (the
# `sum_k Z_k^2/(k pi)^2` representation, 4e6 draws) agrees to 4e-05 absolute through z = 2
# and produced no draw above 2.904; and the dominant eigenvalue `1/pi^2` bounds the tail by
# `exp(-pi^2 z / 2)`, which is 2e-11 at z = 5.
_SERIES_Z_MAX = 5.0


def cvm_limiting_cdf(z: float | np.ndarray) -> np.ndarray:
    """`P(W^2 <= z)` for the LIMITING Cramer-von Mises distribution.

    Anderson and Darling (1952), Annals of Mathematical Statistics 23:193-212, eq (4.35),
    p. 202:

        F(z) = 1/(pi sqrt(z)) * sum_j [Gamma(j+1/2)/(Gamma(1/2) j!)] sqrt(4j+1)
               * exp(-(4j+1)^2/(16 z)) * K_{1/4}((4j+1)^2/(16 z))

    with `K` the modified Bessel function of the second kind. Pinned by
    `tests/test_checks_cvm.py` at four critical values from Anderson and Darling (1952)
    Table 1, p. 203: 0.34730 (alpha .10), 0.46136 (.05), 0.74346 (.01) and 1.16786 (.001).
    No 2.5% point is pinned. There the series gives 0.5806147, which agrees with the 0.581
    of Stephens (1974), JASA 69:730-737, Table 1A part 1.0, to that table's three
    decimals.

    VALID ONLY FOR `z <= _SERIES_Z_MAX`; above it the routine saturates at 1.0 rather than
    evaluating a series that has stopped converging. See that constant for the measurement.

    As with C2's `adinf`, no finite-N correction is applied: the statistic here is divided
    by an ESTIMATED `gamma_hat^2`, so the finite-N null is not the finite-N CvM null, and a
    correction fitted to the latter would be a precision claim the derivation cannot carry.
    """
    z = np.asarray(z, dtype=float)
    if not np.all(np.isfinite(z)):
        raise ValueError(
            "the Cramer-von Mises statistic must be finite; a non-finite value means a "
            "singular segment and must be refused upstream, not turned into a p-value."
        )
    if np.any(z <= 0.0):
        raise ValueError("the Cramer-von Mises statistic is positive; got z <= 0")
    # Evaluate only where the series converges; saturate above it. Returning 1.0 (p = 0)
    # beyond `_SERIES_Z_MAX` reports that the true survival is below what this routine can
    # represent, NOT that it was computed - and it is the only choice that keeps the CDF
    # monotone, which is what a p-value has to be to mean anything.
    inside = z <= _SERIES_Z_MAX
    safe = np.where(inside, z, 1.0)
    total = np.zeros_like(safe)
    for j in range(_CDF_TERMS):
        coefficient = gamma_fn(j + 0.5) / (gamma_fn(0.5) * float(math.factorial(j)))
        w = (4.0 * j + 1.0) ** 2 / (16.0 * safe)
        total = total + coefficient * np.sqrt(4.0 * j + 1.0) * np.exp(-w) * kv(0.25, w)
    value = np.clip(total / (np.pi * np.sqrt(safe)), 0.0, 1.0)
    return np.where(inside, value, 1.0)


def _cvm(x: np.ndarray, tau: float, gamma: float) -> float:
    """One segment's statistic. See the module docstring for the derivation."""
    u = np.cumsum(x) / tau
    n = len(x)
    ue = np.concatenate([[0.0], u])
    i = np.arange(n, dtype=float)
    brace = float(
        np.sum(i**2 * np.diff(ue) - i * n * np.diff(ue**2))
        + n**2 * (u[-1] ** 2 - u[-1] + 1.0 / 3.0)
    )
    return float(brace / (gamma**2 * n))


def _cvm_batch(x_matrix: np.ndarray, tau: float, gamma: np.ndarray) -> np.ndarray:
    u = np.cumsum(x_matrix, axis=1) / tau
    n = x_matrix.shape[1]
    zeros = np.zeros((u.shape[0], 1), dtype=float)
    ue = np.hstack([zeros, u])
    i = np.arange(n, dtype=float)[None, :]
    brace = (i**2 * np.diff(ue, axis=1) - i * n * np.diff(ue**2, axis=1)).sum(axis=1)
    brace = brace + n**2 * (u[:, -1] ** 2 - u[:, -1] + 1.0 / 3.0)
    return brace / (gamma**2 * n)


def statistic(
    segments: list[Segment], *, gamma_estimator: str = GAMMA_DEFAULT
) -> float:
    """Summed unweighted over segments, exactly as C2's multi-process form is.

    `require_strict_tau=False`: unlike eq (7), this integrand has no `1/(s(1-s))` weight,
    so nothing here is singular when the last event lands on the truncation time. That is
    this project's reason to carry CvM at m > 1: it survives records C2 declines. The paper's
    reason is level properties (Section 6, p. 106).
    """
    if not segments:
        raise ValueError("CvM needs at least one segment")
    total = 0.0
    for segment in segments:
        validate_segment(segment, require_strict_tau=False)
        x = np.asarray(segment.x, dtype=float)
        total += _cvm(x, segment.tau, gamma_hat(x, segment.tau, gamma_estimator))
    return total


def statistic_batch(
    segments: list[Segment],
    perm: PermutationSet,
    gamma_estimator: str = GAMMA_DEFAULT,
    permuted: np.ndarray | None = None,
) -> np.ndarray:
    check_permuted(perm, permuted, "cvm.statistic_batch")
    if permuted is None:
        permuted = perm.apply(concatenated_gaps(segments))
    total = np.zeros(perm.n_perm, dtype=float)
    for segment, (lo, hi) in zip(segments, perm.blocks()):
        block = permuted[:, lo:hi]
        gamma = gamma_hat_batch(block, segment.tau, gamma_estimator)
        total += _cvm_batch(block, segment.tau, gamma)
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
    notes = f"gamma={gamma_estimator}"

    if calibration == CALIB_ASYMPTOTIC:
        if len(segments) > 1:
            raise ValueError(
                "CvM has no asymptotic calibration for m > 1 segments here. The paper "
                "(Section 4.2) gives a normal approximation to the weighted sum, but this "
                "implementation sums unweighted; use the permutation calibration, which "
                "is exact under iid gaps and needs no limit."
            )
        p_value = float(1.0 - cvm_limiting_cdf(observed))
        notes += " limiting_CvM"
        if observed > _SERIES_Z_MAX:
            # p is 0.0 by SATURATION, not by computation - the series is not evaluated out
            # here. The verdict (reject) is right and the magnitude is not a measurement, so
            # the result has to say which it is rather than let a reader quote the zero.
            notes += f" p_saturated(z>{_SERIES_Z_MAX:g})"
    elif calibration == CALIB_PERMUTATION:
        perm = resolve_perm(segment_sizes(segments), perm, n_perm, rng)
        null = statistic_batch(segments, perm, gamma_estimator, permuted)
        p_value = permutation_p_value(observed, null)
        notes += f" B={perm.n_perm}"
    else:
        raise ValueError(
            f"CvM supports {CALIB_ASYMPTOTIC!r} and {CALIB_PERMUTATION!r}; "
            f"got {calibration!r}"
        )

    return CheckResult(
        check=CHECK_NAME,
        statistic=observed,
        p_value=p_value,
        calibration=calibration,
        clock=clock,
        n_events=int(sum(len(s.x) for s in segments)),
        n_segments=len(segments),
        n_censored_dropped=int(sum(s.n_censored_dropped for s in segments)),
        notes=notes,
    )
