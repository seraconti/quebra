"""C6 - permutation exchangeability test with a portmanteau rank statistic.

    Q = sum_{h=1..H} z_h^2,    z_h = (r_h - mu_h) / sd_h

where `r_h` is the within-segment rank autocorrelation at lag `h` and `(mu_h, sd_h)` are that
lag's mean and standard deviation across the augmented permutation ensemble, the observed row
included. `z_h` is `_rank_serial.studentize`, the same values whose absolute maximum C5's
studentized variant takes. Calibrated by within-segment permutation, upper tail.

Q is a Ljung-Box-type sum, each lag centred at and scaled by its own permutation moments.
Ljung and Box (1978, Eq. (1.1), p. 297) sum squared autocorrelations that are not centred.
Under exchangeability an autocorrelation does not average zero: Dufour and Roy (1985) define
it on one unbroken series as `autocorrelations` does, one fully centred denominator shared
by every lag (Eq. (2.1), p. 260), give its exact mean (Eq. (2.2), p. 260), note that the
result depends only on exchangeability (p. 261), carry it to rank autocorrelations
(Section 3), note that the Ljung-Box standard error is not exact when the observations are
centred (p. 268), and recommend the exact means and variances for tests of randomness
(p. 271). Across a gap only each segment is exchangeable, so their closed forms do not
cover a gapped record. The permutation moments take their place and hold with gaps and ties.

The centring matters most on a gapped record. When segments sit at different levels, every
within-segment product leans positive whatever the order inside each segment, so the null
mean of `r_h` moves well away from zero. A sum of uncentred squares is not symmetric about
that offset: dependence of one sign inflates it and the other sign deflates it. Permutation
calibration keeps its size exact either way, so what the offset distorts is power, in a
direction set by the sign of the dependence.

The null is exchangeability of the durations within each segment, which is what a
permutation test tests directly: no distributional assumption enters, and the p-value is
exact at finite B. C5's studentized variant reads the same `z_h` and is exact in the same
sense.

Against C5's studentized variant: both read the same `z_h`. C5's max sees dependence
concentrated at one lag at full strength; C6's sum gains where dependence is spread thinly
over several lags and loses where it sits entirely in one. That is the trade that makes both
ship rather than one standing in for the other.

Q is symmetric in the sign of each `z_h`, so it is one-sided by construction and says nothing
about the DIRECTION of any dependence it finds. The direction is read from the sign of `z_1`
in `extra["z"]`, which carries each lag's `z_h` with its sign. C5's `extra["r_lag1"]` is the
uncentred `r_1`, whose null mean is not zero, so its sign is not the direction.

Validity assumptions

- Assumption: under the null the observed row and its B permutations are exchangeable.
  Each row is standardised with the moments of all B + 1 rows together, so the
  standardisation is a symmetric function of the augmented sample and the permutation
  p-value stays exact at finite B.
  Diagnostic: `tests/test_tier3_calibration.py`, the exactness tests, one segment and
  several, with and without a level difference between segments.
  Consequence of violation: the p-value loses its stated size.
  Reference: `_rank_serial` module docstring.
- Assumption: a difference in level between segments is part of the null. The copies
  reshuffle durations only inside a segment, so no permutation-based check in the battery
  (C1, C2 and CvM by permutation, C5, C6) can see such a difference; here the centring
  absorbs it by construction.
  Diagnostic: none among the permutation-based checks.
  Consequence of violation: a record whose segments differ only in level is not rejected,
  at any size of shift.
  Reference: `_permutation` module docstring, on why permutation stays within a segment.
- Assumption: every lag has spread across the augmented ensemble.
  Diagnostic: `studentize` raises a ValueError when a lag's range across the B + 1 rows is
  at most `_rank_serial.ZERO_SPREAD_TOL`, as when every segment holds two events, or
  every segment that pairs at that lag is constant inside. C5's studentized variant raises
  on the same arrays. `run_battery` then raises, so the ledger writes the whole cell as one
  `not computed` row with the reason and reports no other check for it, C3 included, and
  the bench drops the replicate from every row.
  Consequence of violation: `z_h` is undefined, and there is no fallback.
  Reference: `_rank_serial.studentize`.
"""

from __future__ import annotations

import numpy as np

from quebra.analyzers.checks._permutation import (
    resolve_perm,
    DEFAULT_N_PERM,
    PermutationSet,
    permutation_p_value,
)
from quebra.analyzers.checks._rank_serial import (
    MAX_LAG_CAP,
    autocorrelations,
    global_ranks,
    lag_layout,
    studentize,
)
from quebra.analyzers.checks.result import (
    CALIB_PERMUTATION,
    CLOCK_IN_SPEC,
    CheckResult,
    Segment,
    concatenated_gaps,
    segment_sizes,
    validate_segment,
)

CHECK_NAME = "c6_exchangeability"


def run(
    segments: list[Segment],
    *,
    clock: str = CLOCK_IN_SPEC,
    perm: PermutationSet | None = None,
    n_perm: int = DEFAULT_N_PERM,
    max_lag: int = MAX_LAG_CAP,
    rng: np.random.Generator | None = None,
    autocorr: tuple[np.ndarray, np.ndarray, object] | None = None,
) -> CheckResult:
    """`autocorr` is `(observed_r, null_r, layout)` shared with C5 - see `c5.run`."""
    if not segments:
        raise ValueError("C6 needs at least one segment")
    for segment in segments:
        validate_segment(segment, require_strict_tau=False)

    sizes = segment_sizes(segments)

    if autocorr is None:
        layout = lag_layout(sizes, max_lag)
        ranks = global_ranks(concatenated_gaps(segments))
        perm = resolve_perm(sizes, perm, n_perm, rng)
        observed_r = autocorrelations(ranks, layout)[0]
        null_r = autocorrelations(perm.apply(ranks), layout)
    else:
        observed_r, null_r, layout = autocorr
    n_perm_used = int(null_r.shape[0])

    z_observed, z_null = studentize(observed_r, null_r)
    observed = float(np.sum(z_observed**2))
    null = (z_null**2).sum(axis=1)
    p_value = permutation_p_value(observed, null)

    return CheckResult(
        check=CHECK_NAME,
        statistic=observed,
        p_value=p_value,
        calibration=CALIB_PERMUTATION,
        clock=clock,
        n_events=int(sum(sizes)),
        n_segments=len(segments),
        n_censored_dropped=int(sum(s.n_censored_dropped for s in segments)),
        notes=f"lags={list(layout.lags)} B={n_perm_used}",
        extra={
            "lags": list(layout.lags),
            "n_pairs_per_lag": list(layout.n_pairs),
            "z": [float(z) for z in z_observed],
        },
    )
