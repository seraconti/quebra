"""The shared vocabulary of the six checks: what a record looks like and what a
check returns.

A check answers one question - "is this duration sequence consistent with a renewal
process / with exchangeable draws?" - and the answer is only interpretable alongside how
it was calibrated and how many events it rested on. So `CheckResult` carries the
calibration and the counts, not just a p-value.

Two structural facts drive the design:

- **A record is a LIST of segments, not one process.** Kvaloy & Lindqvist's eqs (4)/(7)
  (Technometrics 62(1):101-115, 2020) admit exactly one incomplete gap, the terminal one. A record with interior read gaps
  has an incomplete gap at every gap, so it is m independent time-censored processes and
  needs the multi-process extension (eqs 13-16). A gap-free record is the m = 1 case of
  the same code, so there is no separate single-process path to keep in sync.
- **Time censoring is not event censoring.** `tau` is a TRUNCATION time, chosen (or
  imposed by a gap) without reference to the events. The trailing censored window's
  duration is not an event; it is the residual `tau - T_N`. Feeding it in as an event
  would treat "we stopped watching" as "it failed".

All times are SI seconds. Nothing here prints - the bench calls these hundreds of
thousands of times.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


# Clock labels. Which one a result used is part of the result: the two ask different
# questions of the same record and can legitimately disagree.
CLOCK_IN_SPEC = "in_spec"
CLOCK_CALENDAR = "calendar"
CLOCKS = (CLOCK_IN_SPEC, CLOCK_CALENDAR)

# Calibration labels.
CALIB_ASYMPTOTIC = "asymptotic"
CALIB_PERMUTATION = "permutation"
CALIB_R_COPULA = "r_copula"

# How far `tau` must clear `T_N` in relative terms before eq (7) is considered well posed.
#
# This exists because `np.sum(x)` and `np.cumsum(x)[-1]` are DIFFERENT numbers: numpy's sum
# is pairwise, cumsum is sequential, and they disagree by ulps. `validate_segment` would
# check `tau > np.sum(x)` while `_eq7` divides by `tau - np.cumsum(x)[-1]`, so a segment
# could pass validation and then produce `inf` (asymptotic p = 0.0, a false REJECTION) or
# `nan`. Reproduced on trials 6 and 7 of a 40,000-segment sweep at N = 64.
#
# Both halves are needed and they do different jobs. `last_event_time` pins the CONSUMERS
# (`validate_segment`, `gamma_hat`, `Segment.residual`, `_eq7`) to eq (7)'s own summation
# order, so they agree with each other. The margin covers the PRODUCERS, which are outside
# this module and cannot be pinned: `segments_from_windows` builds `tau` with
# `durations.sum()`, and pairwise-vs-sequential summation disagrees in either direction on
# 23% of random vectors.
#
# 1e-12 is ~4500 ulps at double precision - far above that disagreement (~1e-16 relative)
# and far below any residual a real record produces here, where the trailing gap is either
# exactly zero or at least one read interval (~1e2 s) against a T_N of at most ~1e8 s, six
# orders of headroom. A record with a genuinely sub-picosecond relative residual would be
# refused; these timestamps cannot produce one.
TAU_MARGIN = 1e-12


def last_event_time(x: np.ndarray) -> float:
    """`T_N`, computed the way `c2._eq7` computes it.

    Deliberately `cumsum(x)[-1]` and not `sum(x)`: `c2._eq7` forms `T = cumsum(x)` and divides
    by `tau - T[-1]`, so this is the quantity every guard must be written against. See
    `TAU_MARGIN`.
    """
    return float(np.cumsum(np.asarray(x, dtype=float))[-1])


@dataclass(frozen=True)
class Segment:
    """One time-censored renewal process: complete gaps `x` observed over `[0, tau]`.

    `x` are the COMPLETE inter-event gaps in observation order, so `T = cumsum(x)` are
    the event times and `tau - T[-1]` is the incomplete trailing gap. `tau >= sum(x)`
    always; `tau > sum(x)` strictly whenever the segment ended by truncation rather than
    by an event landing exactly on the boundary.
    """

    x: np.ndarray
    tau: float
    # Windows folded into `tau` instead of used as events (0 or 1 per segment - only the
    # LAST window of a segment can be censored, since a gap_start/scan_end death is what
    # ends a segment).
    n_censored_dropped: int = 0

    @property
    def n_events(self) -> int:
        return int(len(self.x))

    @property
    def event_times(self) -> np.ndarray:
        return np.cumsum(self.x)

    @property
    def residual(self) -> float:
        """The incomplete trailing gap, `tau - T_N`. Uses `last_event_time`, not `sum`."""
        return float(self.tau - last_event_time(self.x))


@dataclass(frozen=True)
class CheckResult:
    check: str
    statistic: float
    p_value: float | None
    calibration: str
    clock: str
    n_events: int
    n_segments: int
    # Windows that were folded into a `tau` rather than counted as events. Not a
    # failure - it is the whole point of time censoring - but it bounds how much the
    # record actually told us.
    n_censored_dropped: int = 0
    notes: str = ""
    extra: dict[str, object] = field(default_factory=dict)


def validate_segment(
    segment: Segment, *, require_strict_tau: bool, min_events: int = 2
) -> None:
    """Raise unless `segment` is a well-posed time-censored renewal observation.

    Raising rather than returning a sentinel is the repo convention and it matters
    here: every failure below yields a finite, plausible, WRONG statistic if waved
    through. A zero first gap makes eq (7)'s `ln(T_2/T_1)` term `+inf`, and a
    permutation can put any gap first; a `tau` below `T_N` makes
    `ln(tau - T_i)` complex; a negative gap flows onto a survival curve as a lifetime.

    `require_strict_tau` is C2's extra demand. Eq (7) contains `ln((tau - T_{N-1})/(tau - T_N))` (with `T_0 = 0`),
    which is `+inf` at `tau == T_N` - verified numerically, not assumed. That case is
    failure (type II) censoring, not the time censoring the equation is derived for, so
    C2 declines it while C1 (eq 4, no such term) accepts it.
    """
    x = np.asarray(segment.x, dtype=float)
    if x.ndim != 1:
        raise ValueError(f"segment gaps must be 1-D; got shape {x.shape}")
    if len(x) < min_events:
        raise ValueError(
            f"this check needs at least {min_events} complete gaps in a segment; "
            f"got {len(x)}"
        )
    if not np.all(np.isfinite(x)):
        raise ValueError(
            f"segment gaps must all be finite; {int((~np.isfinite(x)).sum())} of "
            f"{len(x)} are not"
        )
    if not np.all(x > 0.0):
        n_bad = int((x <= 0.0).sum())
        raise ValueError(
            f"segment gaps must be strictly positive; {n_bad} of {len(x)} are <= 0. "
            "A zero-duration window is producible by the real carve (a single-read "
            "scan_end window has duration_s exactly 0.0) and makes eq (7) singular "
            "once it is the first gap, which a permutation can make it, "
            "so it has to be excluded deliberately upstream rather than here."
        )
    if not np.isfinite(segment.tau):
        raise ValueError(f"tau must be finite; got {segment.tau}")
    total = last_event_time(x)
    # The margin applies to BOTH comparisons. Producers build `tau` by other routes -
    # `segments_from_windows` uses `durations.sum()` (pairwise) while `total` here is
    # `cumsum(x)[-1]` (sequential) - and those disagree by ulps in either direction:
    # measured, `sum < cumsum[-1]` on 23% of random vectors. Without the margin an
    # all-down_crossing block, where `tau == T_N` mathematically, was hard-refused as
    # "tau is below the last event time" on a coin flip of summation order.
    if segment.tau < total * (1.0 - TAU_MARGIN):
        raise ValueError(
            f"tau ({segment.tau:.6g}) is below the last event time T_N ({total:.6g}) by "
            f"more than the {TAU_MARGIN:.0e} relative margin. tau is a truncation time "
            "and must cover every event it censors."
        )
    if require_strict_tau and segment.tau <= total * (1.0 + TAU_MARGIN):
        raise ValueError(
            f"tau ({segment.tau:.6g}) does not clear T_N ({total:.6g}) by the required "
            f"relative margin {TAU_MARGIN:.1e}, so the trailing gap is empty or within "
            "float noise of it. Eq (7) contains ln((tau - T_{N-1})/(tau - T_N)), which is +inf at "
            "equality: this is failure censoring, not the time censoring the statistic "
            "is derived for. Use the calendar clock, or extend tau to the real end of "
            "observation."
        )


def total_events(segments: list[Segment]) -> int:
    return int(sum(s.n_events for s in segments))


def segment_sizes(segments: list[Segment]) -> list[int]:
    return [s.n_events for s in segments]


def concatenated_gaps(segments: list[Segment]) -> np.ndarray:
    """Every segment's gaps end to end - the vector the permutation harness shuffles.

    Order matters and is preserved: the permutation matrix is block-structured over
    exactly this layout, so a check may only reassemble segments by slicing it in the
    same order.
    """
    if not segments:
        return np.zeros(0, dtype=float)
    return np.concatenate([np.asarray(s.x, dtype=float) for s in segments])
