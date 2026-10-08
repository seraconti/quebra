"""Turning a carved record into segments, and the coefficient-of-variation estimators.

`analyzers/windows.py` produces a window table; the checks need time-censored renewal
segments. This module is the only place that mapping lives, and it makes two choices that
the results depend on, so both are stated rather than buried:

**Segment boundaries are read gaps**, taken from `diagnostics["gap_spans_s"]` and from the
birth taxonomy together. Neither alone is enough. The taxonomy is authoritative when the
carve recorded it, but `analyzers/windows.carve` only emits `gap_resume` next to an in-spec
read - a gap flanked by out-of-spec reads leaves NO trace in the window table, and reading
births alone silently merged two renewal processes separated by unobserved hours and folded
an interior censored window into `tau`. Neither source re-derives the gap policy: the gap
threshold, the strict `>` test and gap-wins-ties all stay in `analyzers/windows.py`.

With both sources in play the intended consequence holds again: only the LAST window of a
segment can be censored, because a `gap_start` or `scan_end` death is precisely what ends
one. Callers that omit `gap_spans_s` fall back to the birth-only behaviour and lose that
guarantee.

**Two clocks, because the mapping from windows to a renewal process is not unique.**
See `CLOCK_IN_SPEC` / `CLOCK_CALENDAR` below. Reporting one alone would hide that the
answer depends on the choice.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from quebra.analyzers.windows import (
    BIRTH_GAP_RESUME,
    BIRTH_SCAN_START,
    SIDE_IN_SPEC,
)
from quebra.analyzers.checks.result import (
    CLOCK_CALENDAR,
    CLOCK_IN_SPEC,
    Segment,
    last_event_time,
)

# Estimator labels for gamma-hat, recorded on every result: they differ by whether the
# incomplete trailing gap contributes and, for the complete-gap form, by the divisor. At
# small N either choice moves the statistic.
#
# Source: Kvaloy and Lindqvist, Technometrics 62(1):101-115 (2020), Section 3.6, p. 104,
# define gamma_hat = sigma_hat/mu_hat from the "sample mean" and "sample standard
# deviation" of the complete gaps. Every number they print (Table 2 48.61 and 0.888, the
# Section 8.1 LR 0.681, the Section 8.2 Table 4 p-values) uses the sample 1/(N-1) divisor,
# which is what their R code computes (`sqrt(var(xvec))` in findCV). Their Appendix A.2,
# p. 114, writes the same estimator with 1/N(tau) inside a consistency argument, where the
# divisor is asymptotically irrelevant. Both forms therefore have a source. Each result's
# `notes` names the label that produced it; no ledger or bench row carries it.
GAMMA_COMPLETE = "complete_gaps"
GAMMA_COMPLETE_SAMPLE = "complete_gaps_sample"
GAMMA_TRUNCATED = "eq10_truncated"
GAMMA_ESTIMATORS = (GAMMA_COMPLETE, GAMMA_COMPLETE_SAMPLE, GAMMA_TRUNCATED)

# The estimator every check uses when the caller names none. One constant, so changing the
# default is one line and every signature follows it. The sample form is the default
# because every published number uses it, and because the population form inflates the
# asymptotic statistics: at m = 1, C1 by sqrt(N/(N-1)), C2 and CvM, which carry
# 1/gamma_hat^2, by N/(N-1); at m > 1 each segment carries its own factor. That is anti-conservative and largest exactly where segments are short. At
# m = 1 permutation p-values do not depend on this choice in exact arithmetic, since
# gamma_hat is invariant to reordering the gaps; in floats a tie between the observed and
# a permuted statistic can resolve either way under the two divisors.
GAMMA_DEFAULT = GAMMA_COMPLETE_SAMPLE

# Relative tolerance separating float cancellation from a genuinely negative eq (10)
# variance. Cancellation in E[x^2] - E[x]^2 lands ~1e-16 of scale; the real failure
# lands at 1e-1. Anything between is closer to noise than to signal, so it clamps.
_VAR_NEGATIVE_TOL = 1e-9


def gamma_hat(x: np.ndarray, tau: float, estimator: str = GAMMA_DEFAULT) -> float:
    """Estimated coefficient of variation of the gap distribution.

    `GAMMA_COMPLETE_SAMPLE`, the default, is the `gamma_hat` of Section 3.6 that eqs (4)
    to (7) scale by: the coefficient of variation of the complete gaps alone, with the
    sample 1/(N-1) divisor behind every number Kvaloy and Lindqvist print; see the note on
    the labels above. `GAMMA_COMPLETE` is the same estimator with the population 1/N
    divisor, matching eq (10)'s, so that it and `GAMMA_TRUNCATED` differ only by the
    residual term. The two complete-gap forms differ by the factor `sqrt(N/(N-1))`: at
    m = 1, C1, which scales by 1/gamma_hat, moves by its inverse, and C2 and CvM, which
    scale by 1/gamma_hat^2, by the inverse square. At m > 1 each segment has its own N.

    `GAMMA_TRUNCATED` is eq (10)'s `gamma_tilde` - `mu = tau/N`,
    `sigma^2 = (1/N)[sum(x^2) + (tau - T_N)^2] - mu^2` - the alternative that folds the
    incomplete trailing gap in as information rather than discarding it. It is consistent
    (renewal reward gives `(1/N)sum(x^2) -> E[X^2]` and `tau/N -> E[X]`) but it is a
    DIFFERENCE OF TWO LARGE TERMS and goes negative at small N. Worked example: a perfectly
    regular process with unit gaps truncated at tau = 5.5 has N = 5, mu = 1.1 and
    sigma^2 = 1.05 - 1.21 = -0.16. Measured consequence: at ~5 events per segment it
    failed 383 of 400 bench replicates, which is why the default is the complete-gap form
    and why segments are kept above `MIN_EXPECTED_EVENTS_PER_SEGMENT`.

    Their eq (11) `sigma*^2` is deliberately absent: the authors evaluate it and decline
    to use it ("less satisfactory significance level properties"), so shipping it would
    mean shipping a variant its own source rejects.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n < 2:
        raise ValueError(f"gamma_hat needs at least 2 gaps; got {n}")
    if estimator == GAMMA_COMPLETE:
        # Population (1/N) form, matching eq (10)'s divisor rather than the unbiased
        # (1/(N-1)) one, so the two estimators differ ONLY by the residual term.
        mean = float(np.mean(x))
        var = float(np.mean(x**2) - mean**2)
    elif estimator == GAMMA_COMPLETE_SAMPLE:
        mean = float(np.mean(x))
        var = float(np.mean(x**2) - mean**2) * n / (n - 1)
    elif estimator == GAMMA_TRUNCATED:
        residual = float(tau) - last_event_time(x)
        mean = float(tau) / n
        var = float((np.sum(x**2) + residual**2) / n - mean**2)
    else:
        raise ValueError(
            f"unknown gamma estimator {estimator!r}; known: {list(GAMMA_ESTIMATORS)}"
        )
    if mean <= 0.0:
        raise ValueError(f"gamma_hat got a non-positive mean gap ({mean:.6g})")
    # A MATERIALLY negative variance is a different fact from float cancellation, and
    # clamping both to zero reports the wrong cause: eq (10) is a difference of two
    # large terms and goes genuinely negative at small N (unit gaps truncated at tau = 5.5
    # give N = 5, mu = 1.1, sigma^2 = 1.05 - 1.21 = -0.16, i.e. 13% of scale), which is
    # nothing to do with "every gap is identical". Separate the two.
    if var < -_VAR_NEGATIVE_TOL * mean**2:
        raise ValueError(
            f"the eq (10) variance is negative ({var:.6g} against mu^2 = {mean**2:.6g}). "
            "It is a difference of two large terms and is not usable at this segment "
            f"size (N = {n}); use {GAMMA_DEFAULT!r}, which cannot go negative."
        )
    var = max(var, 0.0)
    if var == 0.0:
        raise ValueError(
            "gamma_hat is zero: every gap is identical, so the statistic it scales is "
            "undefined. A fully quantised duration vector (every window one read long) "
            "reaches this, and it is a degenerate cell rather than a result."
        )
    return float(np.sqrt(var) / mean)


def gamma_hat_batch(
    x_matrix: np.ndarray, tau: float, estimator: str = GAMMA_DEFAULT
) -> np.ndarray:
    """Row-wise `gamma_hat` over a `(B, N)` matrix of permuted gaps.

    Permuting gaps cannot change `sum(x)`, `sum(x^2)` or `N`, so every row shares the
    observed value - but it is computed rather than broadcast, because that invariance is
    the reason the permutation test is exact here and an assertion is cheaper than a
    comment.
    """
    x_matrix = np.asarray(x_matrix, dtype=float)
    n = x_matrix.shape[1]
    if n < 2:
        raise ValueError(f"gamma_hat needs at least 2 gaps; got {n}")
    if estimator == GAMMA_COMPLETE:
        mean = x_matrix.mean(axis=1)
        var = (x_matrix**2).mean(axis=1) - mean**2
    elif estimator == GAMMA_COMPLETE_SAMPLE:
        mean = x_matrix.mean(axis=1)
        var = ((x_matrix**2).mean(axis=1) - mean**2) * n / (n - 1)
    elif estimator == GAMMA_TRUNCATED:
        residual = float(tau) - np.cumsum(x_matrix, axis=1)[:, -1]
        mean = float(tau) / n
        var = ((x_matrix**2).sum(axis=1) + residual**2) / n - mean**2
        mean = np.full(x_matrix.shape[0], mean, dtype=float)
    else:
        raise ValueError(
            f"unknown gamma estimator {estimator!r}; known: {list(GAMMA_ESTIMATORS)}"
        )
    # Same split as the scalar `gamma_hat`: a materially negative eq (10) variance is a
    # different fact from float cancellation and must not be reported as "identical gaps".
    if np.any(var < -_VAR_NEGATIVE_TOL * np.asarray(mean) ** 2):
        raise ValueError(
            "a permuted row produced a negative eq (10) variance. It is a difference of "
            f"two large terms and is not usable at this segment size (N = {n}); use "
            f"{GAMMA_DEFAULT!r}, which cannot go negative."
        )
    var = np.maximum(var, 0.0)
    if np.any(mean <= 0.0) or np.any(var == 0.0):
        raise ValueError(
            "a permuted gap vector produced a degenerate gamma_hat (zero mean or zero "
            "variance). Permutation cannot change either, so the observed vector is "
            "already degenerate and the cell must be reported as such."
        )
    return np.sqrt(var) / mean


def _segment_starts(
    birth_types: list[str],
    t_birth_s: np.ndarray | None = None,
    gap_spans_s: list[tuple[float, float]] | None = None,
) -> list[int]:
    """Indices where a new segment begins.

    Birth types alone are NOT sufficient, which is why `gap_spans_s` exists.
    `analyzers/windows.carve` only emits `gap_resume` when the gap is adjacent to an
    in-spec read; if the record went out of spec before the gap, or stayed out across it,
    the carve emits `down_crossing` ... `up_crossing` and the gap leaves no trace in the
    window table at all. Two renewal processes separated by unobserved hours then became
    one segment, and an interior censored window got folded into `tau`.

    So the authoritative gap list from `WindowsResult.diagnostics["gap_spans_s"]` is used
    when available. It is whole-record and finite-filtered, and `carve` never lets a window
    span a gap, so every gap falls strictly between two windows and each one starts a
    segment at the first window born after it. Birth types remain a second source: they
    catch `scan_start` and agree with the gap list wherever the carve did record one.
    """
    starts = {
        i
        for i, birth in enumerate(birth_types)
        if birth in (BIRTH_SCAN_START, BIRTH_GAP_RESUME)
    }
    if gap_spans_s and t_birth_s is not None and len(t_birth_s):
        births = np.asarray(t_birth_s, dtype=float)
        for _t_before_s, t_after_s in gap_spans_s:
            # First window born at or after the gap ends. searchsorted keeps this O(log n)
            # and, being "left", lands on the window whose birth IS the resuming read.
            idx = int(np.searchsorted(births, float(t_after_s), side="left"))
            if 0 < idx < len(births):
                starts.add(idx)
    # A record whose first window was born by an up_crossing (the scan began out of
    # spec) still starts a segment at index 0.
    starts.add(0)
    return sorted(starts)


def segments_from_windows(
    windows: pd.DataFrame,
    *,
    clock: str,
    min_events: int = 2,
    observation_end_s: float | None = None,
    gap_spans_s: list[tuple[float, float]] | None = None,
) -> tuple[list[Segment], int]:
    """Split one threshold's in-spec window table into time-censored renewal segments.

    A `side` column with any value other than `SIDE_IN_SPEC` raises ValueError; a table
    without one (the bench carve) is taken as in-spec.

    Returns `(segments, n_segments_dropped)`. A segment with fewer than `min_events`
    complete gaps is dropped rather than padded - it carries no information about
    ordering, and silently keeping it would let a record of twenty one-window segments
    report `n_segments=20` while testing nothing.

    `CLOCK_IN_SPEC`: an event is a window dying by `down_crossing`, `x` is its duration,
    and `tau` is the segment's TOTAL in-spec time (complete windows plus the trailing
    censored one). Time advances only while in spec, so this asks whether successive
    in-spec lifetimes look like a renewal process.

    `CLOCK_CALENDAR`: an event is a window BIRTH, `x` is the wall-clock time between
    consecutive births, and `tau` is the segment's observed length. This asks whether
    failures arrive as a renewal process in real time, which is the question a
    maintenance schedule actually poses. It is honest only within a segment - across a
    gap, wall-clock time passed unobserved, which is why segments exist.

    `observation_end_s` is when watching STOPPED, and it applies to the final segment
    only (interior segments stop at the read gap that ends them, which is already in the
    table). It matters more than it looks: without it the final segment's `tau` is taken
    from the last window's death, which is an EVENT-DETERMINED boundary, and eq (7)
    requires a truncation time chosen independently of the events. Left at None, a record
    whose last window is a single read gets `tau == T_N` and C2 correctly refuses it -
    measured at 25% of replicates on an iid read series, so this is the common case, not
    a corner one.

    `gap_spans_s` is `WindowsResult.diagnostics["gap_spans_s"]` and does two jobs. It
    splits segments the birth taxonomy cannot see (see `_segment_starts`), and it supplies
    each interior segment's truncation time: observation of a segment stopped when the gap
    STARTED, not when its last window happened to die. Using the death is the same
    event-determined boundary `observation_end_s` exists to avoid for the final segment.
    Passing it is strongly recommended; omitted, both behaviours degrade to the old ones.
    """
    if clock not in (CLOCK_IN_SPEC, CLOCK_CALENDAR):
        raise ValueError(
            f"unknown clock {clock!r}; known: {CLOCK_IN_SPEC!r}, {CLOCK_CALENDAR!r}"
        )
    required = {"birth_type", "censored", "duration_s", "t_birth_s", "t_death_s"}
    missing = required - set(windows.columns)
    if missing:
        raise KeyError(f"segments_from_windows needs columns {sorted(missing)}")
    if "side" in windows.columns and (windows["side"] != SIDE_IN_SPEC).any():
        other = sorted({str(v) for v in windows["side"]} - {SIDE_IN_SPEC})
        raise ValueError(
            "segments_from_windows takes the in-spec window table; got rows with side "
            f"{other}. Both clocks take their events from the in-spec table (a death by "
            "down_crossing on the in-spec clock, a birth on the calendar clock), and an "
            "out-of-spec table would put every event at the opposite crossing."
        )
    if not len(windows):
        return [], 0

    frame = windows.reset_index(drop=True)
    births = [str(b) for b in frame["birth_type"]]
    t_birth_all = frame["t_birth_s"].to_numpy(dtype=float)
    if np.any(np.diff(t_birth_all) < 0.0):
        raise ValueError(
            "window births are not monotone, so this is not ONE threshold's window table. "
            "`WindowsResult.windows_in_spec` concatenates every threshold; select one with "
            "`windows[windows['threshold_label'] == label]` first. Passing the whole table "
            "silently produces segments split at the wrong places."
        )
    starts = _segment_starts(births, t_birth_all, gap_spans_s)
    bounds = [*starts, len(frame)]
    # Gap START times, for truncating interior segments at the moment watching stopped.
    gap_starts_s = sorted(float(before) for before, _after in (gap_spans_s or []))

    segments: list[Segment] = []
    n_dropped = 0
    n_blocks = len(bounds) - 1
    for position, (lo, hi) in enumerate(zip(bounds[:-1], bounds[1:])):
        block = frame.iloc[lo:hi]
        if not len(block):
            continue
        is_final_block = position == n_blocks - 1
        complete = ~block["censored"].to_numpy(dtype=bool)
        t_birth = block["t_birth_s"].to_numpy(dtype=float)
        t_death = block["t_death_s"].to_numpy(dtype=float)
        # When observation of THIS block stopped. For an interior block that is the start
        # of the gap that ended it - a boundary set by the record, not by an event.
        block_end_s: float | None = None
        after = [g for g in gap_starts_s if g >= float(t_death[-1])]
        if not is_final_block:
            if after:
                block_end_s = after[0]
        elif observation_end_s is not None:
            # The FINAL block can end at a trailing gap too - a record that stops
            # observing and never resumes. Taking `observation_end_s` unconditionally made
            # tau span the unobserved stretch: a gap 12->500 followed by out-of-spec reads
            # to 504 gave a residual of 495 s of which 488 s were never watched, handed to
            # eq (7) as the incomplete trailing gap. Whichever came first is the truth.
            block_end_s = (
                min(float(observation_end_s), after[0])
                if after
                else float(observation_end_s)
            )

        if clock == CLOCK_IN_SPEC:
            durations = block["duration_s"].to_numpy(dtype=float)
            # Only the LAST window of a block may be censored: blocks split at read gaps, so
            # a gap_start death lands at a block end. `x` drops a censored window's duration
            # while `tau` keeps it, so an interior one shifts every later T_i earlier - eq (4)
            # and eq (7) both run on the wrong ones, and the residual is reported as the
            # censored duration when it is really zero.
            #
            # Reachable: `check_ledger.make_inputs_from_windows` takes `gap_spans_s` from
            # `diagnostics.get(...)`, so a caller that omits it splits on birth types alone.
            # `Segment` documents `n_censored_dropped` as 0 or 1 and nothing enforces it.
            interior_censored = int((~complete[:-1]).sum())
            x = durations[complete]
            # In-spec time does not accrue during a gap, so an interior block's tau is its
            # own accumulated in-spec time regardless of when the gap started.
            tau = float(durations.sum())
            n_censored = int((~complete).sum())
        else:
            # Not applicable on this clock: the guard below is in-spec only, because what the
            # calendar clock's tau should be when a block ends at a gap is an open question.
            interior_censored = 0
            x = np.diff(t_birth)
            end_s = float(t_death[-1]) if block_end_s is None else block_end_s
            tau = float(end_s - t_birth[0])
            # On this clock the final window is always the residual, censored or not.
            n_censored = 1
        if len(x) < min_events:
            n_dropped += 1
            continue
        # After the drop, not before: `min_events` means "too small to say anything about",
        # so a block on its way to `n_dropped` must not abort the run.
        if clock == CLOCK_IN_SPEC and interior_censored:
            raise ValueError(
                f"block at position {position} has {interior_censored} censored window(s) "
                f"before its last. Only the final window of a block may be censored; without "
                f"gap_spans_s a gap flanked by out-of-spec reads leaves a censored death "
                f"mid-block, and every event time after it would be wrong."
            )
        segments.append(Segment(x=x, tau=tau, n_censored_dropped=n_censored))
    return segments, n_dropped
