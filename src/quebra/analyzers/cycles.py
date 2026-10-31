"""Calibration cycles of a record of spells, the share of time per status, and the
segments the check battery reads.

A record of spells is a list of rows, each a start, an end and a status, between a declared
record start and end, and optionally a list of unobserved spans: stretches no sample saw.
Nothing here reads a file. A status history of point samples (`schemas/qi_status_history`)
becomes a record of spells through `spells_from_samples`, whose rules are in its own
docstring.
Two statuses are named: `calibrating` and `offline`. Any other status a file declares (idle,
executing, benchmarking) is carried through to the status shares and plays no part in the
cycles: a calibration is a run of calibrating rows only, and a cycle starts when calibrating
ends, whatever status follows.

A **cycle** runs from the end of one calibration to the start of the next. Offline rows and
unobserved spans cut the record into **segments**: in time order, each ends the current
stretch at its start and opens the next at its end. Only stretches of positive length are
segments; two cuts that touch, or one at the record's start or end, leave none between them.
An offline row of zero length is not an outage: it cuts nothing, and it is dropped and
counted (`n_zero_outages_dropped`).
Inside a segment:

- a cycle between two calibrations is complete (its death is observed);
- the last calibration's end to the segment's end (an outage start, or the record end) is
  a censored cycle;
- the segment's start (an outage end, or the record start) to its first calibration has an
  unknown age, so it is not a cycle: it is dropped and counted
  (`n_unobserved_birth_dropped`), as a window born at a gap resume is;
- a segment with no calibration is dropped and counted the same way.

Every cycle has positive length. Calibrating rows that touch (one ends where the next
starts) are one calibration, merged before cycles are built and counted
(`n_calibrations_merged`), so the Kaplan-Meier curve, the MTBC and the battery see one
cycle list. A stretch of zero length before the first or after the last calibration of a
segment is neither a cycle nor counted.

The cycles become an event table (`event_table.CYCLES`), so the estimators are the ones the
windows use. `battery_segments` gives the check battery one `checks.result.Segment` per
segment that holds a cycle: the complete cycles in time order, and tau their sum plus the
censored tail. Calibration time is not in x or tau, so the battery's clock stops while a
calibration runs; the ledger rows still carry the calendar clock's label, because the
cycles are wall-clock durations.

The status shares partition the record because the rows are refused unless they are in
time order, disjoint and inside the record: the time per status and the time between rows
then sum to the record length by construction.

`cycle_checks` runs the battery on those segments and scores each row with the ledger's own
rule (`check_ledger.score_segments`), so a row on cycles means what a row on windows means.
The rows' `threshold_label` column holds `CYCLES_LABEL`, since cycles have no threshold.
`build_cycles_panel_data` builds `CyclesPanelData`, the typed artifact a cycles figure
draws: the curve, the MTBC with the check outcome attached, the verdict rows, the status
shares in hours and every count; `panels/calibration_cycles.py` draws it.

Validity assumptions

- Assumption: a stretch from a segment's start to its first calibration is a residual, not
  a cycle.
  Diagnostic: `n_unobserved_birth_dropped` travels with the table to the panel.
  Consequence of violation: entered as cycles, residuals would pull the MTBC down.
  Reference: the same rule for windows, `analyzers/kaplan_meier.py` module docstring.
- Assumption: calibrating rows that touch are one calibration.
  Diagnostic: `n_calibrations_merged` travels to the panel.
  Consequence of violation: each merge hides a cycle of length 0, so the MTBC is slightly
  high and the cycle count slightly low.
  Reference: no source located; the rule is this project's.
- Assumption: every battery segment ends inside a cycle, so tau > T_N.
  Diagnostic: a segment whose last calibration ends exactly at the segment's end has
  tau = T_N; C2 refuses it, and the ledger's retry drops C1 and C2 for the whole record,
  whose rows then read `not computed` with the note "C1/C2 undefined here".
  Consequence of violation: C1 and C2 give no answer on that record, until the checks read
  a segment that ends at an event in its failure-truncated form; the cycles reach the
  battery only through the shared segment path, so they take that form unchanged.
  Reference: `checks.result.validate_segment`, `check_ledger.score_segments`.
"""

from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

import quebra.analyzers.checks.c3_serial_copula as c3
from quebra.analyzers.calibration_summary import bench_acceptance_at_n
from quebra.analyzers.check_attachment import A2_RENEWAL_CYCLES_ID, CheckAttachment
from quebra.analyzers.check_ledger import LEDGER_COLUMNS, score_segments
from quebra.analyzers.check_outcome import (
    PERMUTATION_KEYS,
    CheckOutcome,
    RowKey,
    battery_flags,
    build_check_outcome,
    check_attachment_of,
    resolve,
)
from quebra.analyzers.checks.battery import drop_short_segments
from quebra.analyzers.checks.result import CLOCK_CALENDAR, Segment
from quebra.analyzers.event_table import CYCLES, EventTable
from quebra.analyzers.kaplan_meier import (
    CONF_LEVEL,
    KaplanMeierCurve,
    curve_from_event_table,
)
from quebra.analyzers.mtbf import Mtbc, mtbc
from quebra.core._artifact_guard import StaleArtifactGuard
from quebra.core.types import Norm

STATUS_CALIBRATING = "calibrating"
STATUS_OFFLINE = "offline"
# Time no row covers, in the status shares. A row may not carry it as its status.
NOT_LISTED = "not listed"

# The `threshold_label` of the ledger rows on cycles: the row schema is the windows', and
# cycles have no threshold, so the column names what was scored.
CYCLES_LABEL = "cycles"

# Carried on `CyclesPanelData.bench_note`: the rows match the bench's calendar cells, which
# were measured on gaps between window births, not on cycles.
BENCH_TRANSFER_NOTE = (
    "bench cells measured on window-birth gaps; their transfer to cycles is unresolved"
)

# A step between samples longer than this many spacings is time no sample saw.
GAP_FACTOR = 2.0

SPELL_KEYS = ("spell_start_unix_s", "spell_end_unix_s", "spell_status")
UNOBSERVED_KEYS = ("unobserved_start_unix_s", "unobserved_end_unix_s")
# Counts the samples step leaves on a record of spells; a record built some other way has
# none to report and reads 0.
SAMPLE_COUNT_KEYS = (
    "n_duplicates_identical_dropped",
    "n_duplicates_conflicting_dropped",
)
RECORD_KEYS = ("record_start_unix_s", "record_end_unix_s", "dataset_id")


def spells_from_samples(samples: Norm, *, spacing_s: float) -> Norm:
    """A status history of point samples as a record of spells.

    Rules, in the order they apply:

    - Duplicate timestamps: rows that agree are one sample; rows that disagree keep the last
      in file order. The rows dropped in each case are counted
      (`n_duplicates_identical_dropped`, `n_duplicates_conflicting_dropped`).
    - A step between two samples longer than `GAP_FACTOR` x `spacing_s` is a polling gap:
      time no sample saw. It becomes an unobserved span from the last sample before it to
      the first after it, which the status shares count as not listed and which cuts the
      record like an outage.
    - A spell is a run of samples of one status with no gap inside: it starts at its first
      sample and ends at the first sample of the next run, or at its own last sample when a
      gap or the record's end follows. A run that a gap splits is two spells.

    The record runs from the first sample to the last.

    Validity assumptions

    - Assumption: a status holds from one sample until the next.
      Diagnostic: none; every spell boundary is known to within one sample spacing.
      Consequence of violation: a status shorter than the spacing can be missed between
      two samples. Each end of a cycle is placed late by less than one spacing, so its
      length is off by less than one spacing.
      Reference: no source located.
    - Assumption: a step longer than `GAP_FACTOR` spacings is time the poller did not see,
      not a status.
      Diagnostic: the unobserved spans travel with the spells and their count to the panel.
      Consequence of violation: a status held across a gap is censored at the gap's start
      and its continuation dropped as a residual, so a long cycle reads as two shorter
      pieces.
      Reference: no source located; the rule is this project's.
    - Assumption: of conflicting rows at one timestamp, the last in file order holds.
      Diagnostic: `n_duplicates_conflicting_dropped` travels to the panel.
      Consequence of violation: one sample's status is misread.
      Reference: no source located; the rule is this project's.
    """
    if not (math.isfinite(spacing_s) and spacing_s > 0.0):
        raise ValueError(f"spacing_s must be finite and positive; got {spacing_s!r}")
    meta = samples["meta"] if "meta" in samples else {}
    missing = [k for k in ("sample_time_unix_s", "sample_status") if k not in samples]
    missing += [k for k in ("dataset_id",) if k not in meta]
    if missing:
        raise KeyError(f"a status history needs {missing}")
    raw_status = list(samples["sample_status"])
    # Checked before any conversion, as `spells_from_norm` does: `dtype=str` would turn a
    # missing value into a status "nan" or "None".
    if not all(isinstance(value, str) and value for value in raw_status):
        raise ValueError("every sample status must be a non-empty string")
    time = np.asarray(samples["sample_time_unix_s"], dtype=float)
    status = np.asarray(raw_status, dtype=str)
    if len(time) == 0:
        raise ValueError("a status history with no samples has no record")
    if len(time) != len(status):
        raise ValueError(
            f"sample arrays differ in length: {len(time)} times, {len(status)} statuses"
        )
    if not np.all(np.isfinite(time)):
        raise ValueError("sample times must be finite")
    if np.any(np.diff(time) < 0.0):
        raise ValueError("samples must be in time order")

    # Duplicate timestamps: one sample per time, the last row's status.
    keep_time: list[float] = []
    keep_status: list[str] = []
    n_identical = 0
    n_conflicting = 0
    first = 0
    while first < len(time):
        last = first
        while last + 1 < len(time) and time[last + 1] == time[first]:
            last += 1
        group = status[first : last + 1]
        if last > first:
            if np.all(group == group[0]):
                n_identical += last - first
            else:
                n_conflicting += last - first
        keep_time.append(float(time[first]))
        keep_status.append(str(group[-1]))
        first = last + 1

    starts: list[float] = []
    ends: list[float] = []
    statuses: list[str] = []
    gap_starts: list[float] = []
    gap_ends: list[float] = []
    run_first = 0
    for k in range(1, len(keep_time)):
        gap = keep_time[k] - keep_time[k - 1] > GAP_FACTOR * spacing_s
        if gap or keep_status[k] != keep_status[run_first]:
            starts.append(keep_time[run_first])
            ends.append(keep_time[k - 1] if gap else keep_time[k])
            statuses.append(keep_status[run_first])
            if gap:
                gap_starts.append(keep_time[k - 1])
                gap_ends.append(keep_time[k])
            run_first = k
    starts.append(keep_time[run_first])
    ends.append(keep_time[-1])
    statuses.append(keep_status[run_first])

    return Norm(
        {
            "spell_start_unix_s": np.asarray(starts, dtype=float),
            "spell_end_unix_s": np.asarray(ends, dtype=float),
            "spell_status": np.asarray(statuses, dtype=str),
            "unobserved_start_unix_s": np.asarray(gap_starts, dtype=float),
            "unobserved_end_unix_s": np.asarray(gap_ends, dtype=float),
            "meta": {
                **dict(meta),
                "record_start_unix_s": keep_time[0],
                "record_end_unix_s": keep_time[-1],
                "spacing_s": float(spacing_s),
                "n_duplicates_identical_dropped": n_identical,
                "n_duplicates_conflicting_dropped": n_conflicting,
            },
        }
    )


@dataclass
class Spells:
    """The rows of a record of spells, validated: in time order, disjoint, inside the
    record, every status a non-empty string; the unobserved spans in time order, inside
    the record, of positive length and overlapping no spell."""

    dataset_id: str
    start_unix_s: np.ndarray
    end_unix_s: np.ndarray
    status: np.ndarray
    record_start_unix_s: float
    record_end_unix_s: float
    unobserved_start_unix_s: np.ndarray
    unobserved_end_unix_s: np.ndarray
    n_duplicates_identical_dropped: int
    n_duplicates_conflicting_dropped: int


def spells_from_norm(norm: Norm) -> Spells:
    """Read and validate the spell contract; raise on any breach, guess nothing."""
    meta = norm["meta"] if "meta" in norm else {}
    missing = [key for key in SPELL_KEYS if key not in norm]
    missing += [key for key in RECORD_KEYS if key not in meta]
    if missing:
        raise KeyError(f"a record of spells needs {missing}")
    raw_status = list(norm["spell_status"])
    # Checked before any conversion: `np.asarray(..., dtype=str)` would turn a missing value
    # into the status "nan", and an offline row so mislabelled would stop cutting the record.
    if not all(isinstance(value, str) and value for value in raw_status):
        raise ValueError("every spell status must be a non-empty string")
    if NOT_LISTED in raw_status:
        raise ValueError(
            f"{NOT_LISTED!r} names the time no row covers; a row cannot carry it as a status"
        )
    start = np.asarray(norm["spell_start_unix_s"], dtype=float)
    end = np.asarray(norm["spell_end_unix_s"], dtype=float)
    status = np.asarray(raw_status, dtype=str)
    record_start = float(meta["record_start_unix_s"])
    record_end = float(meta["record_end_unix_s"])
    if not (len(start) == len(end) == len(status)):
        raise ValueError(
            f"spell arrays differ in length: {len(start)}, {len(end)}, {len(status)}"
        )
    if not (np.all(np.isfinite(start)) and np.all(np.isfinite(end))):
        raise ValueError("spell starts and ends must be finite")
    if not (math.isfinite(record_start) and math.isfinite(record_end)):
        raise ValueError("the record start and end must be finite")
    if record_end <= record_start:
        raise ValueError("the record must end after it starts")
    if np.any(end < start):
        raise ValueError(f"{int(np.sum(end < start))} spell(s) end before they start")
    if np.any(start[1:] < end[:-1]):
        raise ValueError(
            "spells must be in time order and must not overlap: "
            f"{int(np.sum(start[1:] < end[:-1]))} row(s) start before the previous ends"
        )
    # In time order and disjoint, so the first row starts earliest and the last ends last.
    if len(start) and (start[0] < record_start or end[-1] > record_end):
        raise ValueError("a spell lies outside the record's start and end")
    gap_start, gap_end = _unobserved(norm, start, end, record_start, record_end)
    return Spells(
        dataset_id=str(meta["dataset_id"]),
        start_unix_s=start,
        end_unix_s=end,
        status=status,
        record_start_unix_s=record_start,
        record_end_unix_s=record_end,
        unobserved_start_unix_s=gap_start,
        unobserved_end_unix_s=gap_end,
        n_duplicates_identical_dropped=int(meta.get(SAMPLE_COUNT_KEYS[0], 0)),
        n_duplicates_conflicting_dropped=int(meta.get(SAMPLE_COUNT_KEYS[1], 0)),
    )


def _unobserved(
    norm: Norm,
    start: np.ndarray,
    end: np.ndarray,
    record_start: float,
    record_end: float,
) -> tuple[np.ndarray, np.ndarray]:
    """The unobserved spans, validated; empty when the record declares none."""
    present = [key in norm for key in UNOBSERVED_KEYS]
    if not any(present):
        return np.array([], dtype=float), np.array([], dtype=float)
    if not all(present):
        raise KeyError(f"unobserved spans need both {list(UNOBSERVED_KEYS)}")
    lo = np.asarray(norm[UNOBSERVED_KEYS[0]], dtype=float)
    hi = np.asarray(norm[UNOBSERVED_KEYS[1]], dtype=float)
    if len(lo) != len(hi):
        raise ValueError("unobserved span starts and ends differ in length")
    if not (np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))):
        raise ValueError("unobserved spans must be finite")
    if np.any(hi <= lo):
        raise ValueError("an unobserved span must have positive length")
    if np.any(lo[1:] < hi[:-1]):
        raise ValueError("unobserved spans must be in time order and must not overlap")
    if len(lo) and (lo[0] < record_start or hi[-1] > record_end):
        raise ValueError("an unobserved span lies outside the record's start and end")
    overlap = (start[:, None] < hi[None, :]) & (end[:, None] > lo[None, :])
    if np.any(overlap):
        raise ValueError(
            f"{int(np.any(overlap, axis=1).sum())} spell(s) overlap an unobserved span"
        )
    return lo, hi


@dataclass
class CycleTable:
    """The cycles of one record, in time order, with what was merged and dropped.

    `event` is True for a cycle ended by a calibration and False for one censored at an
    outage, at an unobserved span or at the record end. `segment` numbers the segment a cycle
    lies in, counting every segment of the record (stretches of positive length between
    cuts), including those that hold no cycle; `n_segments` is their number. The last four
    counts come with the record: offline rows of zero length dropped, unobserved spans, and
    the duplicate sample rows the samples step dropped.
    """

    dataset_id: str
    start_unix_s: np.ndarray
    end_unix_s: np.ndarray
    length_s: np.ndarray
    event: np.ndarray
    segment: np.ndarray
    n_segments: int
    n_unobserved_birth_dropped: int
    n_calibrations_merged: int
    n_zero_outages_dropped: int
    n_unobserved_spans: int
    n_duplicates_identical_dropped: int
    n_duplicates_conflicting_dropped: int

    @property
    def n_cycles(self) -> int:
        return int(len(self.length_s))


_CUT = "unobserved"


def _segments(
    spells: Spells,
) -> tuple[list[tuple[float, float, list[tuple[float, float]]]], int]:
    """The segments in time order, as (start, end, calibrating rows), only those of
    positive length, and the number of zero-length offline rows dropped. Touching
    calibrations are merged later, by `_merge`."""
    rows = [
        (float(s), 0, float(e), str(status))
        for s, e, status in zip(spells.start_unix_s, spells.end_unix_s, spells.status)
    ]
    # A span sorts after a row that starts with it: such a row is a zero-length spell at the
    # last sample before the span, and belongs to the stretch the span closes.
    rows += [
        (float(s), 1, float(e), _CUT)
        for s, e in zip(spells.unobserved_start_unix_s, spells.unobserved_end_unix_s)
    ]
    rows.sort(key=lambda row: (row[0], row[1]))
    segments: list[tuple[float, float, list[tuple[float, float]]]] = []
    n_zero_outages = 0
    lo = spells.record_start_unix_s
    current: list[tuple[float, float]] = []
    for start, _order, end, status in rows:
        if status == STATUS_OFFLINE and end == start:
            n_zero_outages += 1
        elif status in (STATUS_OFFLINE, _CUT):
            if start > lo:
                segments.append((lo, start, current))
            lo, current = end, []
        elif status == STATUS_CALIBRATING:
            current.append((start, end))
    if spells.record_end_unix_s > lo:
        segments.append((lo, spells.record_end_unix_s, current))
    return segments, n_zero_outages


def _merge(
    calibrations: list[tuple[float, float]],
) -> tuple[list[tuple[float, float]], int]:
    """Touching calibrations as one, and how many rows were merged into a predecessor."""
    merged: list[tuple[float, float]] = []
    n_merged = 0
    for start, end in calibrations:
        if merged and start == merged[-1][1]:
            merged[-1] = (merged[-1][0], end)
            n_merged += 1
        else:
            merged.append((start, end))
    return merged, n_merged


def cycles(norm: Norm) -> CycleTable:
    """The cycles of a record of spells, by the rule in the module docstring."""
    spells = spells_from_norm(norm)
    starts: list[float] = []
    ends: list[float] = []
    events: list[bool] = []
    segment_of: list[int] = []
    n_dropped = 0
    n_merged = 0
    segments, n_zero_outages = _segments(spells)
    for index, (lo, hi, rows) in enumerate(segments):
        calibrations, merges = _merge(rows)
        n_merged += merges
        if not calibrations:
            n_dropped += 1
            continue
        n_dropped += int(calibrations[0][0] > lo)
        for (_s0, e0), (s1, _e1) in zip(calibrations, calibrations[1:]):
            starts.append(e0)
            ends.append(s1)
            events.append(True)
            segment_of.append(index)
        last_end = calibrations[-1][1]
        if hi > last_end:
            starts.append(last_end)
            ends.append(hi)
            events.append(False)
            segment_of.append(index)
    start_arr = np.asarray(starts, dtype=float)
    end_arr = np.asarray(ends, dtype=float)
    return CycleTable(
        dataset_id=spells.dataset_id,
        start_unix_s=start_arr,
        end_unix_s=end_arr,
        length_s=end_arr - start_arr,
        event=np.asarray(events, dtype=bool),
        segment=np.asarray(segment_of, dtype=int),
        n_segments=len(segments),
        n_unobserved_birth_dropped=n_dropped,
        n_calibrations_merged=n_merged,
        n_zero_outages_dropped=n_zero_outages,
        n_unobserved_spans=int(len(spells.unobserved_start_unix_s)),
        n_duplicates_identical_dropped=spells.n_duplicates_identical_dropped,
        n_duplicates_conflicting_dropped=spells.n_duplicates_conflicting_dropped,
    )


def cycle_event_table(table: CycleTable) -> EventTable:
    """The cycles as an event table, the one input every survival estimator reads."""
    return EventTable(
        dataset_id=table.dataset_id,
        threshold_label="",
        side=CYCLES,
        age_s=table.length_s,
        event=table.event,
        n_unobserved_birth_dropped=table.n_unobserved_birth_dropped,
    )


def battery_segments(table: CycleTable, min_events: int) -> tuple[list[Segment], int]:
    """One battery segment per segment that holds a cycle, short ones dropped and counted.

    x is the segment's complete cycles in time order; tau is their sum plus the censored
    cycle that ends the segment, if there is one. A segment with fewer than `min_events`
    complete cycles is dropped by `checks.battery.drop_short_segments` and counted in the
    returned number; a segment holding no cycle at all is not a battery segment, so that
    number and `CycleTable.n_segments` count different things.
    """
    segments: list[Segment] = []
    for index in np.unique(table.segment):
        mask = table.segment == index
        lengths = table.length_s[mask]
        complete = table.event[mask]
        x = lengths[complete]
        tail = float(np.sum(lengths[~complete]))
        tau = (float(np.cumsum(x)[-1]) if len(x) else 0.0) + tail
        segments.append(
            Segment(x=x, tau=tau, n_censored_dropped=int(np.count_nonzero(~complete)))
        )
    return drop_short_segments(segments, min_events)


@dataclass
class StatusShares:
    """Seconds and fraction of the record per status, `NOT_LISTED` for time no row covers.

    Statuses appear in the order the record first declares them, `NOT_LISTED` last. The
    seconds sum to `record_length_s` by construction (module docstring).
    """

    dataset_id: str
    statuses: tuple[str, ...]
    seconds: tuple[float, ...]
    fractions: tuple[float, ...]
    record_length_s: float


def status_shares(norm: Norm) -> StatusShares:
    """Time per status, and the time no row covers, summed gap by gap."""
    spells = spells_from_norm(norm)
    order: list[str] = []
    for status in spells.status:
        if status not in order:
            order.append(str(status))
    durations = spells.end_unix_s - spells.start_unix_s
    seconds = [
        math.fsum(durations[spells.status == status].tolist()) for status in order
    ]
    if len(spells.start_unix_s):
        gaps = [
            float(spells.start_unix_s[0]) - spells.record_start_unix_s,
            *(spells.start_unix_s[1:] - spells.end_unix_s[:-1]).tolist(),
            spells.record_end_unix_s - float(spells.end_unix_s[-1]),
        ]
    else:
        gaps = [spells.record_end_unix_s - spells.record_start_unix_s]
    order.append(NOT_LISTED)
    seconds.append(math.fsum(gaps))
    length = spells.record_end_unix_s - spells.record_start_unix_s
    return StatusShares(
        dataset_id=spells.dataset_id,
        statuses=tuple(order),
        seconds=tuple(seconds),
        fractions=tuple(s / length for s in seconds),
        record_length_s=length,
    )


# ---------------------------------------------------------------------------
# The check battery on cycles
# ---------------------------------------------------------------------------


@dataclass
class CycleChecks:
    """The battery's rows on one record's cycles, and the grid they reshape into.

    `rows` has the ledger's columns (`check_ledger.LEDGER_COLUMNS`), `threshold_label` set to
    `CYCLES_LABEL` and `clock` to the calendar clock. `n_segments_dropped` counts the
    segments holding a cycle but too few complete ones for the battery. `battery_note` is
    the reason the scoring gave when it wrote no row per check (one `(all)` row instead),
    and empty otherwise. The scoring knobs are kept, so a row can be read without the job.
    """

    dataset_id: str
    rows: pd.DataFrame
    outcome: CheckOutcome
    n_segments_dropped: int
    battery_note: str
    alpha: float
    min_events_pass: int
    tie_cutoff_distinct: int
    lag_max: int
    n_permutations: int
    seed: int
    min_events_per_segment: int


def cycle_checks(
    table: CycleTable,
    bench_size_table: pd.DataFrame,
    *,
    run_set: tuple[RowKey, ...],
    alpha: float,
    min_events_pass: int,
    tie_cutoff_distinct: int,
    lag_max: int,
    n_permutations: int,
    seed: int,
    min_events_per_segment: int,
) -> CycleChecks:
    """The check battery on the cycles, scored by the ledger's own rule.

    The run-set must be permutation rows (`check_outcome.PERMUTATION_KEYS`), the rows the
    a2 record names as its diagnostic. It sets the battery's switches
    (`check_outcome.battery_flags`), which are coarse, so `rows` can hold rows outside it;
    the grid shows the run-set only. Never raises on degenerate cycles: with no segment
    left, the one row is `not computed` with its reason (`battery_note`), and every asked
    row reads as unanswered.
    """
    outside = [key for key in run_set if tuple(key) not in PERMUTATION_KEYS]
    if outside:
        raise ValueError(
            f"cycle checks run permutation rows only; got {[list(k) for k in outside]}"
        )
    run, _display = resolve(run_set, run_set)
    flags = battery_flags(run)
    segments, n_dropped = battery_segments(table, min_events_per_segment)
    acceptance = bench_acceptance_at_n(bench_size_table, alpha=alpha)
    grid = sorted({int(n) for n in acceptance["n_target"]})
    rows = score_segments(
        segments,
        n_dropped=n_dropped,
        n_used=table.n_cycles,
        label=CYCLES_LABEL,
        clock=CLOCK_CALENDAR,
        dataset_id=table.dataset_id,
        acceptance=acceptance,
        grid=grid,
        alpha=alpha,
        min_events_pass=min_events_pass,
        tie_cutoff_distinct=tie_cutoff_distinct,
        lag_max=lag_max,
        n_permutations=n_permutations,
        seed=seed,
        min_events_per_segment=min_events_per_segment,
        include_tau_checks=flags["include_tau_checks"],
        include_c2_asymptotic=flags["include_c2_asymptotic"],
        include_c3=flags["include_c3"],
        c3_n_null_sim=c3.N_NULL_SIM,
    )
    frame = pd.DataFrame(rows).reindex(columns=LEDGER_COLUMNS)
    blank = frame[frame["check_id"] == "(all)"]
    outcome = build_check_outcome(
        frame,
        threshold_label=CYCLES_LABEL,
        run_set=run,
        display_set=run,
        datasets=(table.dataset_id,),
        clocks=(CLOCK_CALENDAR,),
        alpha=alpha,
    )
    return CycleChecks(
        dataset_id=table.dataset_id,
        rows=frame,
        outcome=outcome,
        n_segments_dropped=n_dropped,
        battery_note="; ".join(str(note) for note in blank["notes"]),
        alpha=alpha,
        min_events_pass=min_events_pass,
        tie_cutoff_distinct=tie_cutoff_distinct,
        lag_max=lag_max,
        n_permutations=n_permutations,
        seed=seed,
        min_events_per_segment=min_events_per_segment,
    )


def attachment(checks: CycleChecks) -> CheckAttachment:
    """The check outcome as the fields the cycles' MTBC and band carry, naming a2."""
    return check_attachment_of(checks.outcome, A2_RENEWAL_CYCLES_ID)


# ---------------------------------------------------------------------------
# The cycles figure's typed artifact
# ---------------------------------------------------------------------------


@dataclass
class CyclesPanelData(StaleArtifactGuard):
    """Everything a cycles figure draws, every number computed here.

    The curve and the MTBC are in minutes on the figure's axis (`mtbc_min` is point, lower,
    upper). `checks` is the check outcome the band and the MTBC interval rest on, and
    `check_summary` its one line. `verdict_rows` holds (check, p-value, verdict, note) for
    each check the run-set asked for, on the calendar clock; `battery_note` is the reason
    when the scoring wrote no row per check, so a column of absent cells says why. The
    counts: cycles, censored cycles, residuals dropped, calibrations merged, segments,
    segments the battery dropped, zero-length outages dropped, unobserved spans, and
    duplicate sample rows dropped (identical, conflicting).
    """

    dataset_id: str
    curve: KaplanMeierCurve
    mtbc: Mtbc
    mtbc_min: tuple[float, float, float]
    tau_min: float
    checks: CheckAttachment
    check_summary: str
    bench_note: str
    battery_note: str
    verdict_rows: tuple[tuple[str, float | None, str, str], ...]
    shares: StatusShares
    share_h: tuple[float, ...]
    n_cycles: int
    n_censored: int
    n_unobserved_birth_dropped: int
    n_calibrations_merged: int
    n_segments: int
    n_segments_dropped: int
    n_zero_outages_dropped: int
    n_unobserved_spans: int
    n_duplicates_identical_dropped: int
    n_duplicates_conflicting_dropped: int


def build_cycles_panel_data(
    table: CycleTable,
    checks: CycleChecks,
    shares: StatusShares,
    *,
    conf_level: float = CONF_LEVEL,
) -> CyclesPanelData:
    """The complete `CyclesPanelData` of one record, from its cycles, checks and shares.
    Raises if the three describe different records."""
    ids = {table.dataset_id, checks.dataset_id, shares.dataset_id}
    if len(ids) != 1:
        raise ValueError(f"cycles, checks and shares describe different records: {ids}")
    event_table = cycle_event_table(table)
    attached = attachment(checks)
    summary = dataclasses.replace(
        mtbc(event_table, conf_level=conf_level), checks=attached
    )
    curve = curve_from_event_table(
        event_table, label=table.dataset_id, conf_level=conf_level
    )
    verdict_rows = tuple(
        (grid.label, grid.p_values[0][0], grid.verdicts[0][0], grid.notes[0][0])
        for grid in checks.outcome.grids
    )
    return CyclesPanelData(
        dataset_id=table.dataset_id,
        curve=curve,
        mtbc=summary,
        mtbc_min=(
            summary.rmst_s / 60.0,
            summary.lower_s / 60.0,
            summary.upper_s / 60.0,
        ),
        tau_min=summary.tau_s / 60.0,
        checks=attached,
        check_summary=attached.summary(),
        bench_note=BENCH_TRANSFER_NOTE,
        battery_note=checks.battery_note,
        verdict_rows=verdict_rows,
        shares=shares,
        share_h=tuple(seconds / 3600.0 for seconds in shares.seconds),
        n_cycles=table.n_cycles,
        n_censored=int(np.count_nonzero(~table.event)),
        n_unobserved_birth_dropped=table.n_unobserved_birth_dropped,
        n_calibrations_merged=table.n_calibrations_merged,
        n_segments=table.n_segments,
        n_segments_dropped=checks.n_segments_dropped,
        n_zero_outages_dropped=table.n_zero_outages_dropped,
        n_unobserved_spans=table.n_unobserved_spans,
        n_duplicates_identical_dropped=table.n_duplicates_identical_dropped,
        n_duplicates_conflicting_dropped=table.n_duplicates_conflicting_dropped,
    )
