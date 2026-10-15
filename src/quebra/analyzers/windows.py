"""Window carving on both sides of a threshold: gap policy, censoring, per-read table.

Ported from `monoliths/v13fig/v13_shape.py::spacing/carve`, with two deliberate
differences from that reference:

- **Death timestamp.** This module keeps the repo's convention - a window dies at the
  FIRST OUT-OF-SPEC read, so the interval spanning the crossing is inside the lifetime.
  `v13_shape` measured to the last in-spec read instead. Censored windows have no
  out-of-spec read to point at, so they die at their last in-spec read.
- **`min_reads` is not ported HERE.** `v13_shape` drops windows with fewer than 5 reads
  before any statistic; carving without that filter is what preserves parity with the
  pre-existing carve. The filter does exist downstream: `analyzers/distinguish_band.py`
  applies `shape_min_reads` (default 5) to decide which windows enter the excursion
  shape statistics, and nothing else. It never touches a boundary, duration or count.

Everything else - the gap threshold, the strict `>` test, the positive-only median, the
birth/death taxonomy, gap-wins-ties - is the monolith's, verbatim.

`margin` is the one definition of in spec for the whole package: in spec means
`margin >= 0`, so a read exactly at the threshold is in spec in both directions. Every
consumer classifies through `in_spec_mask` or reads the carve's `in_spec` column; none
compares a metric with a threshold itself.

Uncertainty is an ANNOTATION and never moves a window boundary. The carve is crisp
(`margin >= 0`), so `k` and `use_uncertainty` change the per-read `state` column and
nothing else. Deciding whether an uncertain read should break a window is a question
`jobs/bench/probe_unresolved.py` exists to answer; it must not be answered by assumption
here.

All times are SI seconds (`_s` suffixes). Values and thresholds share whatever unit the
caller supplies; only their comparison matters.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

DEFAULT_GAP_MULT = 10.0

# Crossings are named relative to spec: an up_crossing enters spec, a down_crossing leaves
# it. An in-spec window is born by an up_crossing and dies by a down_crossing; an
# out-of-spec window the other way round. A birth or death that is not a crossing was not
# observed. The window tables carry that as two booleans, `birth_observed` and `censored`,
# and consumers read the booleans rather than comparing names. `n_endurance_bags` is the
# legacy name for the count of windows whose birth was not observed.
BIRTH_UP_CROSSING = "up_crossing"
BIRTH_DOWN_CROSSING = "down_crossing"
BIRTH_SCAN_START = "scan_start"
BIRTH_GAP_RESUME = "gap_resume"

DEATH_DOWN_CROSSING = "down_crossing"
DEATH_UP_CROSSING = "up_crossing"
DEATH_GAP_START = "gap_start"
DEATH_SCAN_END = "scan_end"

# Per-read states. The two _uncertain states appear only when use_uncertainty is on.
STATE_IN_SPEC = "in_spec"
STATE_OUT_OF_SPEC = "out_of_spec"
STATE_IN_SPEC_UNCERTAIN = "in_spec_uncertain"
STATE_OUT_OF_SPEC_UNCERTAIN = "out_of_spec_uncertain"
# Not a state of the system - a state of the RECORD. A timeline that simply breaks at a
# gap is ambiguous (missing? in spec? a layout artefact?); an explicit grey band says
# "the instrument was not reporting here" and keeps the axis continuous so two timelines
# stay comparable read for read.
STATE_UNOBSERVED = "unobserved"

# The two sides of a threshold. The strings are the read states, so a side and the state of
# the reads in it cannot be spelled two ways.
SIDE_IN_SPEC = STATE_IN_SPEC
SIDE_OUT_OF_SPEC = STATE_OUT_OF_SPEC

# t_birth_s is the window's first read and t_last_s its last. An observed death is at the
# first read past the crossing (t_death_s); a censored window dies at t_last_s. The
# intervals are closed: an observed birth lies in [t_before_birth_s, t_birth_s] and an
# observed death in [t_last_s, t_death_s], so with both observed the length lies in
# [t_last_s - t_birth_s, t_death_s - t_before_birth_s]. A duplicate timestamp at a
# crossing makes its interval a single point, an exact value. t_before_birth_s is NaN
# when the birth was not observed.
# extreme_margin is the smallest margin over the window's reads: the closest approach to
# the threshold in spec, minus the depth out of spec.
WINDOW_COLUMNS = [
    "dataset_id",
    "threshold_label",
    "threshold_value",
    "big_values_good",
    "side",
    "window_index",
    "t_before_birth_s",
    "t_birth_s",
    "t_last_s",
    "t_death_s",
    "duration_s",
    "birth_type",
    "death_type",
    "birth_observed",
    "censored",
    "n_reads",
    "extreme_margin",
]

# window_index, window_age_s and forward_time_s refer to windows_in_spec only and are
# null on an out-of-spec read. Both window tables number from 0, so joining
# window_index onto windows_out_of_spec pairs reads with the wrong windows.
READ_COLUMNS = [
    "dataset_id",
    "threshold_label",
    "window_index",
    "t_read_s",
    "value",
    "margin",
    "window_age_s",
    "forward_time_s",
    "sigma_v",
    "sigma_known",
    "in_spec",
    "state",
]


@dataclass(slots=True)
class WindowsInputs:
    t_rel_s: np.ndarray
    values: np.ndarray
    thresholds: list[tuple[str, float, bool]]  # (label, value, big_values_good)
    sigma: np.ndarray | None = None  # per-read 1-sigma, same unit as values
    dataset_id: str = ""
    gap_mult: float = DEFAULT_GAP_MULT
    k: float = 1.0
    use_uncertainty: bool = False


@dataclass(slots=True)
class WindowsResult:
    windows_in_spec: pd.DataFrame
    windows_out_of_spec: pd.DataFrame
    reads: pd.DataFrame
    meta: dict[str, object]
    diagnostics: dict[str, object] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Gap policy
# ---------------------------------------------------------------------------


def spacing(t_rel_s: np.ndarray, gap_mult: float = DEFAULT_GAP_MULT):
    """Median read spacing and the derived gap threshold.

    The median is over POSITIVE differences only. A repeated timestamp gives a zero
    difference, and enough of them would drive the median to zero, making every positive
    step a gap and destroying every window. With fewer than two reads there is no spacing:
    the median is NaN and nothing is a gap.

    Returns (median_spacing_s, gap_threshold_s, n_gaps, n_nonpositive_steps). Raises
    ValueError unless gap_mult is finite and positive: an infinite one makes the gap
    threshold infinite, a zero one makes every step a gap.
    """
    if not (np.isfinite(gap_mult) and gap_mult > 0.0):
        raise ValueError(f"gap_mult must be finite and positive, got {gap_mult!r}")
    dt = np.diff(np.asarray(t_rel_s, dtype=float))
    if len(dt) == 0:
        return float("nan"), float("inf"), 0, 0
    n_nonpositive = int((dt <= 0).sum())
    positive = dt[dt > 0]
    if not len(positive):
        # Every step is a duplicate or a step backwards, so there is no observed
        # spacing to scale the gap threshold by. Inventing one (the reference
        # monolith falls back to 1.0) would silently put a fabricated physical
        # scale into every window boundary downstream.
        raise ValueError(
            f"cannot derive a read spacing: all {len(dt)} inter-read steps are "
            f"non-positive. A gap threshold needs at least one positive spacing."
        )
    median_s = float(np.median(positive))
    gap_threshold_s = gap_mult * median_s
    # Strict: a spacing exactly equal to the threshold is NOT a gap.
    n_gaps = int((dt > gap_threshold_s).sum())
    return median_s, gap_threshold_s, n_gaps, n_nonpositive


def gap_flags(t_rel_s: np.ndarray, gap_threshold_s: float) -> np.ndarray:
    """Read i opens a new segment. Index 0 is never a gap, so the first read of a
    record is always a scan_start and never a gap_resume."""
    dt = np.diff(np.asarray(t_rel_s, dtype=float))
    return np.r_[False, dt > gap_threshold_s]


@dataclass(slots=True)
class ObservedReads:
    """The reads a record observed, with the one gap rule every consumer shares.

    A *segment* is a maximal run of reads with no gap between neighbours. `is_gap[i]` is
    True when read i opens a segment after a gap; read 0 opens the first segment and is
    never a gap.
    """

    t_rel_s: np.ndarray
    values: np.ndarray
    sigma: np.ndarray | None
    is_gap: np.ndarray
    median_spacing_s: float
    gap_threshold_s: float
    n_gaps: int
    n_nonpositive_steps: int
    n_dropped_nonfinite: int

    def segments(self) -> list[tuple[int, int]]:
        """`[start, end)` read-index ranges of the segments, in order."""
        n = len(self.t_rel_s)
        if n == 0:
            return []
        starts = [0, *(int(i) for i in np.flatnonzero(self.is_gap))]
        return list(zip(starts, [*starts[1:], n]))

    def within_segment(self) -> np.ndarray:
        """`within[i]` is True when the spacing from read i to read i+1 is observed,
        that is, does not cross a gap. One entry per consecutive pair of reads."""
        return ~self.is_gap[1:]

    def effective_spacing_s(self) -> float:
        """sum(delta^2) / sum(delta) over the within-segment spacings, gaps excluded.

        The mean length of the interval a crossing lands in, when a stationary process
        crosses at a constant rate: the interval is chosen with probability proportional
        to its length. It equals the mean spacing times (1 + CV^2), with the CV taken with
        the 1/n variance. NaN when no within-segment spacing is positive.
        """
        delta = np.diff(self.t_rel_s)[self.within_segment()]
        total = float(np.sum(delta))
        return float(np.sum(delta * delta)) / total if total > 0.0 else float("nan")


def observed_reads(
    t_rel_s: np.ndarray,
    values: np.ndarray,
    *,
    sigma: np.ndarray | None = None,
    gap_mult: float = DEFAULT_GAP_MULT,
) -> ObservedReads:
    """Drop non-finite (time, value) pairs, then apply the gap rule.

    ORDER IS LOAD-BEARING: dropping non-finite pairs BEFORE measuring spacing lets a run of
    failed fits widen the interval into a real gap. The other way round either measures
    spacing on the wrong array or lets a NaN read stand inside a segment.

    Raises ValueError on mismatched lengths (values or sigma), on time that steps
    backwards, and, through `spacing`, when no step between reads is positive. Fewer than
    two finite reads leave no spacing: `median_spacing_s` is NaN and there is no gap.
    """
    t_all = np.asarray(t_rel_s, dtype=float)
    v_all = np.asarray(values, dtype=float)
    if len(t_all) != len(v_all):
        raise ValueError(
            f"t_rel_s and values must be the same length; got {len(t_all)} and "
            f"{len(v_all)}"
        )
    sigma_all = None
    if sigma is not None:
        sigma_all = np.asarray(sigma, dtype=float)
        if len(sigma_all) != len(t_all):
            raise ValueError(
                f"sigma must be the same length as t_rel_s; got {len(sigma_all)} and "
                f"{len(t_all)}"
            )
    finite = np.isfinite(t_all) & np.isfinite(v_all)
    t = t_all[finite]

    # Durations are t_death - t_birth, so time running backwards yields NEGATIVE
    # lifetimes that flow straight onto the survival curve - finite, plausible and
    # wrong. Loaders sort by timestamp; a caller that did not has a bug upstream.
    if len(t) > 1 and bool(np.any(np.diff(t) < 0)):
        raise ValueError(
            "t_rel_s must be non-decreasing; got a time array that steps backwards. "
            "Sort the reads by timestamp before carving."
        )

    median_s, gap_threshold_s, n_gaps, n_nonpositive = spacing(t, gap_mult)
    return ObservedReads(
        t_rel_s=t,
        values=v_all[finite],
        sigma=sigma_all[finite] if sigma_all is not None else None,
        is_gap=gap_flags(t, gap_threshold_s) if len(t) else np.zeros(0, dtype=bool),
        median_spacing_s=median_s,
        gap_threshold_s=gap_threshold_s,
        n_gaps=n_gaps,
        n_nonpositive_steps=n_nonpositive,
        n_dropped_nonfinite=int(len(t_all) - len(t)),
    )


def carve(
    t_rel_s: np.ndarray,
    in_spec: np.ndarray,
    is_gap: np.ndarray,
) -> list[dict[str, object]]:
    """Segment a boolean mask into runs of True, with birth and death types.

    Returns dicts of {s, e, birth_type, death_type} where `s` is the index of the first
    True read and `e` is one past the last - so for a down_crossing `e` indexes the first
    False read, which is the death timestamp. The crossing names describe the MASK turning
    on (up) and off (down); `_side_windows` renames them relative to spec for the
    out-of-spec side.

    The gap check runs BEFORE the mask check, so a read that is both post-gap and False
    produces gap_start, not down_crossing. Gap wins ties.
    """
    n = len(in_spec)
    windows: list[dict[str, object]] = []
    i = 0
    while i < n:
        if not in_spec[i]:
            i += 1
            continue
        s = i
        if s == 0:
            birth = BIRTH_SCAN_START
        elif is_gap[s]:
            birth = BIRTH_GAP_RESUME
        else:
            birth = BIRTH_UP_CROSSING
        j = s
        death: str | None = None
        while j + 1 < n:
            if is_gap[j + 1]:
                death = DEATH_GAP_START
                break
            if not in_spec[j + 1]:
                death = DEATH_DOWN_CROSSING
                break
            j += 1
        if death is None:
            death = DEATH_SCAN_END
        windows.append({"s": s, "e": j + 1, "birth_type": birth, "death_type": death})
        i = j + 1
    return windows


def mark_gaps_in_segments(
    segments: list[tuple[float, float, str]],
    gap_spans_h: list[tuple[float, float]],
) -> list[tuple[float, float, str]]:
    """Replace the part of any timeline segment inside a read gap with STATE_UNOBSERVED.

    A compliance bar drawn across a gap asserts a state through hours the instrument
    was not reporting - the same claim the carve refuses to make when it terminates a
    window with `gap_start`. But a bar that simply STOPS is just as bad: a blank reads
    as "in spec" or as a rendering artefact. The gap gets its own colour instead, so the
    timeline stays continuous and unobserved time is stated rather than implied.

    Preconditions, held by the only producer (`run` emits one span per consecutive read
    pair): `gap_spans_h` are disjoint, and `segments` tile one contiguous stretch.
    Overlapping spans would double-count the overlap.
    """
    if not gap_spans_h:
        return sorted(segments)
    out: list[tuple[float, float, str]] = []
    for start, end, state in segments:
        pieces = [(start, end)]
        for lo, hi in gap_spans_h:
            nxt: list[tuple[float, float]] = []
            for a, b in pieces:
                if b <= lo or a >= hi:
                    nxt.append((a, b))
                    continue
                if a < lo:
                    nxt.append((a, lo))
                if b > hi:
                    nxt.append((hi, b))
            pieces = nxt
        out.extend((a, b, state) for a, b in pieces if b > a)
    # the gap itself, once per span, clipped to the timeline's own extent
    if segments:
        lo_edge = min(a for a, _, _ in segments)
        hi_edge = max(b for _, b, _ in segments)
        for lo, hi in gap_spans_h:
            a, b = max(lo, lo_edge), min(hi, hi_edge)
            if b > a:
                out.append((a, b, STATE_UNOBSERVED))
    return sorted(out)


def margin(
    values: np.ndarray, threshold_value: float, big_values_good: bool
) -> np.ndarray:
    """Signed distance from the threshold, positive on the good side.

    `values - threshold_value` when big values are good (T2*), `threshold_value - values`
    when small values are good (infidelity). In the caller's units.

    Raises ValueError on a non-finite threshold or value: neither has a side, and a NaN
    compared with zero would read as out of spec. `run` drops non-finite reads first.
    """
    values = np.asarray(values, dtype=float)
    threshold_value = float(threshold_value)
    if not math.isfinite(threshold_value):
        raise ValueError(f"threshold_value must be finite; got {threshold_value!r}")
    if not np.all(np.isfinite(values)):
        raise ValueError(
            f"margin needs finite values; {int(np.sum(~np.isfinite(values)))} of "
            f"{values.size} are not. Drop non-finite reads before classifying them."
        )
    if big_values_good:
        return values - threshold_value
    return threshold_value - values


def check_k(k: float) -> float:
    """Return the confidence level k, or raise ValueError unless it is finite and positive:
    k = 0 would make every read confident. The one check: `confident_mask` and
    `windows.run` both ask it."""
    if not (np.isfinite(k) and k > 0.0):
        raise ValueError(f"k must be finite and positive, got {k!r}")
    return float(k)


def confident_mask(margin_v: np.ndarray, sigma: np.ndarray, k: float) -> np.ndarray:
    """A read is confident at level k when its sigma is finite and |margin| >= k sigma.

    A read whose sigma is unknown is never confident. With `use_uncertainty`, the carve's
    uncertain reads are exactly the reads that are not confident, so a read with no sigma
    is uncertain, its value picking the side. Raises ValueError on a bad k (`check_k`) and
    on a negative sigma, which is no error bar.
    """
    check_k(k)
    m = np.asarray(margin_v, dtype=float)
    s = np.asarray(sigma, dtype=float)
    known = np.isfinite(s)
    if np.any(s[known] < 0.0):
        raise ValueError("a reported sigma is negative")
    out = np.zeros(m.shape, dtype=bool)
    out[known] = np.abs(m[known]) >= k * s[known]
    return out


def in_spec_mask(
    values: np.ndarray, threshold_value: float, big_values_good: bool
) -> np.ndarray:
    """In spec means `margin >= 0`: a read exactly at the threshold is in spec, either way.

    Raises, through `margin`, on a non-finite threshold or value.
    """
    return margin(values, threshold_value, big_values_good) >= 0.0


# The crossing that bears and the one that kills a window on each side, in spec terms.
_SIDE_CROSSINGS = {
    SIDE_IN_SPEC: (BIRTH_UP_CROSSING, DEATH_DOWN_CROSSING),
    SIDE_OUT_OF_SPEC: (BIRTH_DOWN_CROSSING, DEATH_UP_CROSSING),
}


def _side_windows(
    t: np.ndarray,
    on_side: np.ndarray,
    is_gap: np.ndarray,
    margin_v: np.ndarray,
    side: str,
    base: dict[str, object],
) -> list[dict[str, object]]:
    """Window rows for one side. `_s` and `_e` carry the read indices for the checks."""
    born_by, dies_by = _SIDE_CROSSINGS[side]
    rows: list[dict[str, object]] = []
    for index, w in enumerate(carve(t, on_side, is_gap)):
        s, e = int(w["s"]), int(w["e"])
        birth_observed = w["birth_type"] == BIRTH_UP_CROSSING
        death_observed = w["death_type"] == DEATH_DOWN_CROSSING
        t_birth_s = float(t[s])
        t_last_s = float(t[e - 1])
        # An observed death is at the first read past the crossing (index e). A censored
        # window has no such read, so it dies at its last read (index e-1).
        t_death_s = float(t[e]) if death_observed else t_last_s
        rows.append(
            {
                **base,
                "side": side,
                "window_index": index,
                "t_before_birth_s": float(t[s - 1]) if birth_observed else np.nan,
                "t_birth_s": t_birth_s,
                "t_last_s": t_last_s,
                "t_death_s": t_death_s,
                "duration_s": t_death_s - t_birth_s,
                "birth_type": born_by if birth_observed else str(w["birth_type"]),
                "death_type": dies_by if death_observed else str(w["death_type"]),
                "birth_observed": bool(birth_observed),
                "censored": not death_observed,
                "n_reads": e - s,
                "extreme_margin": float(np.min(margin_v[s:e])),
                "_s": s,
                "_e": e,
            }
        )
    return rows


def _check_tiling_and_accounting(
    t: np.ndarray,
    is_gap: np.ndarray,
    ins: np.ndarray,
    rows: list[dict[str, object]],
    label: str,
) -> None:
    """Raise unless both sides' windows and the gaps tile the record, and each side's
    durations add up to the observed time its reads were on that side.

    Tiling is checked by float equality of shared endpoints, in read-index order, so it
    involves no summation. Accounting compares two computations that share nothing: the
    windows' durations, and the observed spacing charged to each read by its own state
    (a gap interval charged to nobody). On one side, n_windows times the mean duration
    must equal that side's fraction times the observed time.
    """

    def fail(what: str) -> RuntimeError:
        return RuntimeError(f"carve invariant violated at threshold {label!r}: {what}")

    next_index, next_time = 0, float(t[0])
    for row in sorted(rows, key=lambda r: int(r["_s"])):
        s, e = int(row["_s"]), int(row["_e"])
        if s != next_index or row["t_birth_s"] != next_time:
            raise fail(f"a window starts at read {s}, expected read {next_index}")
        if not row["censored"]:
            if e >= len(t) or is_gap[e] or row["t_death_s"] != t[e]:
                raise fail(f"an observed death at read {e} is not the next read")
            next_time = float(t[e])
        elif row["death_type"] == DEATH_GAP_START:
            if e >= len(t) or not is_gap[e] or row["t_death_s"] != t[e - 1]:
                raise fail(f"a gap death at read {e} does not meet a gap")
            next_time = float(t[e])
        elif e != len(t) or row["t_death_s"] != t[-1]:
            raise fail("a scan-end death is not at the last read")
        next_index = e
    if next_index != len(t):
        raise fail(f"the windows cover {next_index} of {len(t)} reads")

    observed = ~is_gap[1:]
    dt = np.diff(t)
    total = math.fsum(dt[observed])
    for side, on_side in ((SIDE_IN_SPEC, ins), (SIDE_OUT_OF_SPEC, ~ins)):
        from_reads = math.fsum(dt[observed & on_side[:-1]])
        from_windows = math.fsum(r["duration_s"] for r in rows if r["side"] == side)
        if not math.isclose(
            from_reads, from_windows, rel_tol=1e-12, abs_tol=1e-12 * total
        ):
            raise fail(
                f"{side} windows last {from_windows!r} s but its reads were observed for "
                f"{from_reads!r} s"
            )


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


def make_inputs_from_frame(
    frame: pd.DataFrame,
    *,
    time_col: str,
    value_col: str,
    thresholds: list[tuple[str, float, bool]],
    dataset_id: str = "",
    sigma_col: str | None = None,
    gap_mult: float = DEFAULT_GAP_MULT,
    k: float = 1.0,
    use_uncertainty: bool = False,
) -> WindowsInputs:
    for col in (time_col, value_col):
        if col not in frame.columns:
            raise KeyError(f"windows requires column {col!r} in the frame")
    sigma = None
    if sigma_col is not None:
        if sigma_col not in frame.columns:
            raise KeyError(
                f"windows was asked for sigma column {sigma_col!r}, which the frame "
                f"does not have. Columns: {list(frame.columns)}"
            )
        sigma = frame[sigma_col].to_numpy(dtype=float)
    return WindowsInputs(
        t_rel_s=frame[time_col].to_numpy(dtype=float),
        values=frame[value_col].to_numpy(dtype=float),
        thresholds=list(thresholds),
        sigma=sigma,
        dataset_id=dataset_id,
        gap_mult=gap_mult,
        k=k,
        use_uncertainty=use_uncertainty,
    )


def make_inputs_from_norm(
    norm: Mapping[str, object],
    *,
    value_key: str,
    thresholds: list[tuple[str, float, bool]],
    sigma_key: str | None = None,
    gap_mult: float = DEFAULT_GAP_MULT,
    k: float = 1.0,
    use_uncertainty: bool = False,
) -> WindowsInputs:
    if "t_rel_s" not in norm:
        raise KeyError("windows requires 't_rel_s' in the normalized mapping.")
    if value_key not in norm:
        raise KeyError(f"windows requires {value_key!r} in the normalized mapping.")
    meta = norm.get("meta", {}) if isinstance(norm.get("meta", {}), Mapping) else {}
    sigma = None
    if sigma_key is not None:
        if sigma_key not in norm:
            raise KeyError(
                f"windows was asked for sigma key {sigma_key!r}, which the norm does "
                f"not have."
            )
        sigma = np.asarray(norm[sigma_key], dtype=float)
    return WindowsInputs(
        t_rel_s=np.asarray(norm["t_rel_s"], dtype=float),
        values=np.asarray(norm[value_key], dtype=float),
        thresholds=list(thresholds),
        sigma=sigma,
        dataset_id=str(meta.get("dataset_id", "")),
        gap_mult=gap_mult,
        k=k,
        use_uncertainty=use_uncertainty,
    )


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------


def _window_frame(rows: list[dict[str, object]]) -> pd.DataFrame:
    frame = pd.DataFrame(rows, columns=WINDOW_COLUMNS)
    if len(frame):
        frame = frame.astype(
            {
                "window_index": "int64",
                "n_reads": "int64",
                "birth_observed": bool,
                "censored": bool,
            }
        )
    return frame


def run(inputs: WindowsInputs) -> WindowsResult:
    """Carve every threshold on the ladder and emit the window and read tables.

    Non-finite (time, value) pairs are dropped before carving and counted in
    `diagnostics`. Raises ValueError on mismatched lengths (values or sigma),
    `use_uncertainty` without sigma, time that steps backwards, no positive step between
    reads (`spacing`), a repeated threshold label or a non-finite threshold value; and
    RuntimeError when the windows fail the tiling or accounting check.
    """
    t_all = np.asarray(inputs.t_rel_s, dtype=float)
    # k is recorded on the windows and read downstream, so it is checked whether or not
    # this carve draws uncertainty.
    check_k(inputs.k)
    if inputs.use_uncertainty and inputs.sigma is None:
        raise ValueError(
            "use_uncertainty=True requires a sigma array. A dataset with no error "
            "column is legal, but then use_uncertainty must be False - the uncertain "
            "states cannot be silently inferred."
        )
    # Every table and `per_threshold` key a threshold by its label. Finiteness is also
    # checked in `margin`, which a record of fewer than two reads never reaches.
    labels = [label for label, _, _ in inputs.thresholds]
    repeated = sorted({label for label in labels if labels.count(label) > 1})
    if repeated:
        raise ValueError(
            f"threshold labels must be unique; {repeated} appear more than once, and "
            f"their windows would merge under one label."
        )
    for label, thr_value, _ in inputs.thresholds:
        if not math.isfinite(float(thr_value)):
            raise ValueError(f"threshold {label!r} has non-finite value {thr_value!r}")

    # A NaN left in place would evaluate NaN >= threshold as False and fabricate a
    # down-crossing at every failed fit; `observed_reads` drops it first.
    obs = observed_reads(
        t_all, inputs.values, sigma=inputs.sigma, gap_mult=inputs.gap_mult
    )
    t, v, sigma, is_gap = obs.t_rel_s, obs.values, obs.sigma, obs.is_gap
    median_s, gap_threshold_s = obs.median_spacing_s, obs.gap_threshold_s
    n_gaps, n_nonpositive = obs.n_gaps, obs.n_nonpositive_steps

    sigma_known = (
        np.isfinite(sigma) if sigma is not None else np.zeros(len(t), dtype=bool)
    )

    window_rows: list[dict[str, object]] = []
    out_window_rows: list[dict[str, object]] = []
    read_rows: list[dict[str, object]] = []
    per_threshold: dict[str, dict[str, object]] = {}

    for label, thr_value, big_values_good in inputs.thresholds:
        if len(t) < 2:
            per_threshold[label] = {
                "n_windows": 0,
                "n_censored": 0,
                "n_endurance_bags": 0,
                "n_windows_out_of_spec": 0,
                "n_censored_out_of_spec": 0,
            }
            continue

        ins = in_spec_mask(v, thr_value, big_values_good)
        margin_v = margin(v, thr_value, big_values_good)
        base = {
            "dataset_id": inputs.dataset_id,
            "threshold_label": label,
            "threshold_value": float(thr_value),
            "big_values_good": bool(big_values_good),
        }
        in_rows = _side_windows(t, ins, is_gap, margin_v, SIDE_IN_SPEC, base)
        out_rows = _side_windows(t, ~ins, is_gap, margin_v, SIDE_OUT_OF_SPEC, base)
        _check_tiling_and_accounting(t, is_gap, ins, in_rows + out_rows, label)
        window_rows.extend(in_rows)
        out_window_rows.extend(out_rows)

        window_of_read = np.full(len(t), -1, dtype=int)
        birth_of_read = np.full(len(t), np.nan, dtype=float)
        death_of_read = np.full(len(t), np.nan, dtype=float)
        for row in in_rows:
            s, e = int(row["_s"]), int(row["_e"])
            window_of_read[s:e] = int(row["window_index"])
            birth_of_read[s:e] = row["t_birth_s"]
            death_of_read[s:e] = row["t_death_s"]

        state = np.where(ins, STATE_IN_SPEC, STATE_OUT_OF_SPEC).astype(object)
        if inputs.use_uncertainty:
            # Not confident, including every read whose sigma is unknown.
            uncertain = ~confident_mask(margin_v, sigma, inputs.k)
            state[uncertain & ins] = STATE_IN_SPEC_UNCERTAIN
            state[uncertain & ~ins] = STATE_OUT_OF_SPEC_UNCERTAIN

        for i in range(len(t)):
            has_window = window_of_read[i] >= 0
            read_rows.append(
                {
                    "dataset_id": inputs.dataset_id,
                    "threshold_label": label,
                    "window_index": int(window_of_read[i]) if has_window else None,
                    "t_read_s": float(t[i]),
                    "value": float(v[i]),
                    "margin": float(margin_v[i]),
                    "window_age_s": float(t[i] - birth_of_read[i])
                    if has_window
                    else None,
                    "forward_time_s": float(death_of_read[i] - t[i])
                    if has_window
                    else None,
                    "sigma_v": float(sigma[i]) if sigma is not None else None,
                    "sigma_known": bool(sigma_known[i]),
                    "in_spec": bool(ins[i]),
                    "state": str(state[i]),
                }
            )

        per_threshold[label] = {
            # n_censored and n_endurance_bags OVERLAP (a scan_start + scan_end window is
            # both). Never sum them. The unsuffixed counts are the in-spec side.
            "n_windows": len(in_rows),
            "n_censored": sum(1 for r in in_rows if r["censored"]),
            "n_endurance_bags": sum(1 for r in in_rows if not r["birth_observed"]),
            "n_windows_out_of_spec": len(out_rows),
            "n_censored_out_of_spec": sum(1 for r in out_rows if r["censored"]),
        }

    windows_df = _window_frame(window_rows)
    out_windows_df = _window_frame(out_window_rows)
    reads_df = pd.DataFrame(read_rows, columns=READ_COLUMNS)
    if len(reads_df):
        reads_df["window_index"] = reads_df["window_index"].astype("Int64")

    # (t_before, t_after) for every gap, so consumers can break a line across one
    # instead of drawing a segment through unobserved time.
    gap_spans_s: list[tuple[float, float]] = []
    if len(t) > 1:
        for idx in np.flatnonzero(is_gap[1:]) + 1:
            gap_spans_s.append((float(t[idx - 1]), float(t[idx])))

    n_sigma_unknown = int((~sigma_known).sum()) if sigma is not None else len(t)
    diagnostics: dict[str, object] = {
        "n_reads_raw": int(len(t_all)),
        "n_reads_finite": int(len(t)),
        "n_reads_dropped_nonfinite": int(len(t_all) - len(t)),
        "n_reads_sigma_unknown": n_sigma_unknown,
        "median_spacing_s": median_s,
        "gap_threshold_s": gap_threshold_s,
        "n_gaps": n_gaps,
        "gap_spans_s": gap_spans_s,
        "n_nonpositive_steps": n_nonpositive,
        "per_threshold": per_threshold,
    }
    print(
        f"[windows] reads={len(t)}/{len(t_all)} thresholds={len(inputs.thresholds)} "
        f"windows={len(windows_df)} gaps={n_gaps} "
        f"gap_thr={gap_threshold_s:.1f}s median_dt={median_s:.1f}s",
        flush=True,
    )
    return WindowsResult(
        windows_in_spec=windows_df,
        windows_out_of_spec=out_windows_df,
        reads=reads_df,
        meta={
            "dataset_id": inputs.dataset_id,
            "gap_mult": inputs.gap_mult,
            "k": inputs.k,
            "use_uncertainty": inputs.use_uncertainty,
        },
        diagnostics=diagnostics,
    )
