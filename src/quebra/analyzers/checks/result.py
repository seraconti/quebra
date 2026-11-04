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
- **Time censoring is not event censoring.** `tau` is where observation stopped. When it
  is a TRUNCATION time, chosen (or imposed by a gap) without reference to the events, the
  trailing censored window's duration is not an event; it is the residual `tau - T_N`.
  Feeding it in as an event would treat "we stopped watching" as "it failed". When
  observation stopped at an event, `tau == T_N` and the segment is failure truncated
  (`Segment.ends_at_event`): C1, C2 and CvM read its first `n - 1` event times
  (`events_read`).

All times are SI seconds. Nothing here prints - the bench calls these hundreds of
thousands of times.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


# Clock labels. Which one a result used is part of the result: the two ask different
# questions of the same record and can legitimately disagree. Each runs on either bag
# (`_multiprocess.segments_from_windows`): on the side clock an event is a window ended by
# the crossing that leaves its side, on the calendar clock a window birth.
CLOCK_SIDE = "side"
CLOCK_CALENDAR = "calendar"
CLOCKS = (CLOCK_SIDE, CLOCK_CALENDAR)
# The side clock's name where it was calibrated and seeded as the in-spec bag's clock: the
# bench tables (`jobs/bench/results/*.csv`), the bench's own cells, and the permutation
# streams (`check_ledger.stream_for`). `BENCH_CLOCK` maps each clock to the name those
# lookups use, so the in-spec bag keeps its bench cells and its p-values. The checks'
# `clock=` defaults are bench names too, because the bench calls them without one.
CLOCK_IN_SPEC = "in_spec"
BENCH_CLOCK = {CLOCK_SIDE: CLOCK_IN_SPEC, CLOCK_CALENDAR: CLOCK_CALENDAR}

# Calibration labels.
CALIB_ASYMPTOTIC = "asymptotic"
CALIB_PERMUTATION = "permutation"
CALIB_R_COPULA = "r_copula"

# How far `tau` must clear `T_N` in relative terms before eq (7) on all `N` events is
# considered well posed. Within it, a segment ends at its last event (`at_last_event`).
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


def at_last_event(tau: float, t_n: float) -> bool:
    """True when `tau` is the last event time `t_n` within `TAU_MARGIN`: observation stopped
    at an event (failure truncation), not at a time chosen apart from the events.

    The one test of it. `validate_segment`, `Segment.ends_at_event` and the published-record
    tier (`instrument_validation`) all read this inequality, so they cannot disagree.
    """
    return bool(tau <= t_n * (1.0 + TAU_MARGIN))


def events_read(x: np.ndarray, at_event: bool) -> np.ndarray:
    """The gaps whose running sums are the event times C1, C2 and CvM read.

    Every gap of a time-truncated segment; every gap but the last of one that ends at an
    event, whose end `T_n` is then the truncation time. Works on one segment's gaps or on a
    `(B, n)` block of permuted rows. Only the event times shrink: validation, `gamma_hat`
    and the permutation blocks keep all `n` gaps, which are all complete either way.
    """
    return x[..., :-1] if at_event else x


@dataclass(frozen=True)
class Segment:
    """One renewal process: complete gaps `x` observed over `[0, tau]`.

    `x` are the COMPLETE inter-event gaps in observation order, so `T = cumsum(x)` are
    the event times and `tau - T[-1]` is the incomplete trailing gap. `tau >= sum(x)`
    always; `tau > sum(x)` strictly when the segment is time truncated, and `tau == T_N`
    within `TAU_MARGIN` when observation stopped at its last event (`ends_at_event`).
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

    @property
    def ends_at_event(self) -> bool:
        """The segment stopped at its last event: `tau == T_N` within `TAU_MARGIN`.

        C1, C2 and CvM then read its first `n - 1` event times on `(0, T_n)`, with
        `tau = T_n` (`events_read`): given `T_n`, under a constant rate those are uniform
        order statistics (Rigdon and Basu 2000, Statistical Methods for the Reliability of
        Repairable Systems, Sect. 4.3, the failure-truncated case). A segment ends this way
        when the record leaves the bag's side, when a cycles segment ends with a complete
        cycle, and when the final window is one read long (zero duration), censored or not.
        """
        return at_last_event(self.tau, last_event_time(self.x))

    @property
    def n_events_read(self) -> int:
        """The events C1, C2 and CvM read: `n - 1` when the segment ends at an event."""
        return self.n_events - 1 if self.ends_at_event else self.n_events


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

    `require_strict_tau` refuses a segment that ends at an event (`at_last_event`). Eq (7)
    on all `N` events contains `ln((tau - T_{N-1})/(tau - T_N))` (with `T_0 = 0`), which is
    `+inf` at `tau == T_N` - verified numerically, not assumed. That case is failure (type
    II) truncation, not the time censoring the equation is derived for. No check asks for
    it: C1, C2 and CvM read such a segment's first `N - 1` events (`events_read`), where
    the term is finite.
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
    if require_strict_tau and at_last_event(segment.tau, total):
        raise ValueError(
            f"tau ({segment.tau:.6g}) does not clear T_N ({total:.6g}) by the required "
            f"relative margin {TAU_MARGIN:.1e}, so the trailing gap is empty or within "
            "float noise of it. Eq (7) contains ln((tau - T_{N-1})/(tau - T_N)), which is +inf at "
            "equality: this is failure truncation, not the time censoring the statistic "
            "is derived for. Read the first N - 1 events with tau = T_N instead "
            "(`events_read`)."
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
