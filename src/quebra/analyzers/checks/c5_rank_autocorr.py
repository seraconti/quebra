"""C5 - studentized max-lag rank autocorrelation, permutation calibrated.

    T = max_{h=1..H} | (r_h - mu_h) / sd_h |

where `r_h` is the within-segment rank autocorrelation at lag `h` and `(mu_h, sd_h)` are
that lag's moments across the augmented permutation ensemble (see `_rank_serial`).

Directed at LOCALIZED serial dependence: if the durations carry a lag-2 relation and
nothing else, the max picks it up at full strength while a portmanteau (C6) dilutes it
across every lag it also looked at. The two are shipped as a pair for that reason, not
as redundancy.

**The unstudentized variant runs beside it in the bench, on purpose - and measurement
changed what is claimed for it.** The expectation was that the SCALE of `max_h |r_h|` would
misbehave generally. It does not. Because every lag shares one denominator, the per-lag null
spreads are nearly identical on an ungapped record (SD ratio 1.04 across lags 1-5 at
N = 60), so there studentizing barely changes any lag's scale. It still moves each lag's
centre, and that alone separates the two variants on one segment: the unstudentized one is
uncentred (below). In `jobs/bench/results/power_table.csv` it rejects less often than the
studentized one in 84 of the 88 Arm E cells where the two differ, and in 63 of 65 on Arm D,
every one of them a single segment.

Where they genuinely separate is short, numerous segments - the censored and gapped case.
At 10 segments of 5 the lag pair counts run 40 down to 10, the null SD ratio reaches 1.64,
and the unstudentized max lands on lag 1 in 43% of replicates against 21-29% once each lag
is on its own scale. So the studentized form is a safeguard for the multi-process regime,
not a general-purpose improvement, and running both in the bench is what lets it say which
regime is which instead of taking the claim on trust.

Both variants read the SAME permutation ensemble, so their difference in any bench cell is
a paired difference and its Monte Carlo error is smaller than that of two independent runs.

**The unstudentized variant is also uncentred, which is why it runs only in the bench.**
Under exchangeability `r_h` does not average zero. On one segment its mean is
`-(n - h)/(n(n - 1))`, which is `-1/n` at lag 1 (Dufour and Roy 1985, Eq. (2.2), p. 260).
On segments at different levels every within-segment product leans positive, so the mean
moves up by an amount that grows with the difference in level; the closed form is
pinned by `tests/test_checks_statistics.py::test_the_gapped_null_mean_is_the_closed_form`.
`max_h |r_h|` then leans as an uncentred sum does: dependence of one sign is inflated and
the other deflated. Permutation calibration keeps its size exact, so what the offset
distorts is power, in a direction set by the sign of the dependence. It is the bench's
measure of what studentizing buys, and nothing else: `battery.run_battery` computes it
only with `include_c5_raw=True`, which the bench passes, so no ledger, survey or band
carries it.

`extra["r_lag1"]` is the same uncentred `r_1` of the observed record, with the same
offsets: about `-1/n` on one segment, and the level effect across gaps. Its sign is not the
direction of the dependence. Read that from the sign of `z_1`, which C6 reports in
`extra["z"]`.
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

CHECK_NAME = "c5_rank_autocorr"

VARIANT_STUDENTIZED = "studentized"
VARIANT_RAW = "unstudentized"
VARIANTS = (VARIANT_STUDENTIZED, VARIANT_RAW)


def run(
    segments: list[Segment],
    *,
    variant: str = VARIANT_STUDENTIZED,
    clock: str = CLOCK_IN_SPEC,
    perm: PermutationSet | None = None,
    n_perm: int = DEFAULT_N_PERM,
    max_lag: int = MAX_LAG_CAP,
    rng: np.random.Generator | None = None,
    autocorr: tuple[np.ndarray, np.ndarray, object] | None = None,
) -> CheckResult:
    """`autocorr` is `(observed_r, null_r, layout)` when a caller already computed them.

    The studentized variant, the unstudentized variant and C6 are three aggregations of
    ONE lag-autocorrelation matrix, so `checks/battery.py` computes it once and passes it
    to all three. Left at None, this recomputes it.
    """
    if variant not in VARIANTS:
        raise ValueError(f"unknown C5 variant {variant!r}; known: {list(VARIANTS)}")
    if not segments:
        raise ValueError("C5 needs at least one segment")
    for segment in segments:
        # No strict-tau demand: C5 never touches tau. It is a test of the ORDER of the
        # observed durations, so a segment whose last event lands on the truncation time
        # is perfectly testable here even though C2 declines it.
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

    if variant == VARIANT_STUDENTIZED:
        observed_row, null_rows = studentize(observed_r, null_r)
    else:
        observed_row, null_rows = observed_r, null_r

    observed = float(np.max(np.abs(observed_row)))
    null = np.max(np.abs(null_rows), axis=1)
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
        notes=f"variant={variant} lags={list(layout.lags)} B={n_perm_used}",
        extra={
            "lags": list(layout.lags),
            "n_pairs_per_lag": list(layout.n_pairs),
            # The uncentred lag-1 autocorrelation, which the bench averages into
            # `mean_induced_lag1`. It carries the null offset in the module docstring, so
            # its sign is not the direction of the dependence.
            "r_lag1": float(observed_r[0]),
        },
    )
