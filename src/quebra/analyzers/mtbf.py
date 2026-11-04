"""Inter-event intervals of a calibration log, and the mean time between calibrations.

`run` takes a list of CalibrationEvent objects (from a CalibrationLogSchema Norm) and
computes inter-event intervals. No filtering is applied - all events in the log are used
regardless of chi_squared or other quality metrics.

`mtbc` is the mean time between calibrations (MTBC): the restricted mean of the
product-limit curve on an event table of cycles (`event_table.CYCLES`), with its standard
error and a normal interval. The point and the SE are `survival.curve_summaries`, the code
the window estimators use; only the interval is added here. `intervals_event_table` turns
the intervals of `run` into such a table, every interval a completed cycle.

Validity assumptions

- Assumption: the cycles of one record are draws from one law, independent of one another.
  Diagnostic: the battery's permutation rows on the cycles, where a job runs them and
  attaches the outcome; none for a log's inter-event intervals.
  Consequence of violation: the SE and the interval do not hold their level; the point is
  still the area under this record's curve.
  Reference: R `survival` 3.8.6 `summary(fit, rmean = tau)` for the restricted mean and its
  SE (pinned in `tests/test_survival_r_reference.py`); no source located for the interval.
- Assumption: the interval rmst +/- z SE is close to its level.
  Diagnostic: none here; a restricted mean is bounded below by 0, so a lower end at or
  below 0 says the normal approximation is poor for this record.
  Consequence of violation: coverage away from `conf_level`.
  Reference: no source located.
- Assumption: every cycle at risk at age t shares one hazard, and censoring is independent
  of the cycle's future (the multiplicative intensity model every estimator on an event
  table rests on).
  Diagnostic: the attached check outcome; inter-event intervals of a log are never censored.
  Consequence of violation: the curve estimates a mixture, and with censored cycles the
  point moves too.
  Reference: the `survival` module docstring, which states it for every estimator there.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from quebra.analyzers.check_attachment import (
    A2_RENEWAL_CYCLES_ID,
    DESCRIPTIVE,
    RESTS_ON_EXCHANGEABILITY,
    CheckAttachment,
)
from quebra.analyzers.event_table import CYCLES, EventTable
from quebra.analyzers.kaplan_meier import CONF_LEVEL, z_two_sided
from quebra.analyzers.survival import CurveSummaries, curve_summaries
from quebra.core.types import CalibrationEvent, Norm


@dataclass
class MtbfInputs:
    events: list[CalibrationEvent]
    meta: dict = field(default_factory=dict)


@dataclass
class MtbfResult:
    intervals_s: np.ndarray  # N-1 inter-event intervals (seconds)
    event_times_unix_s: (
        np.ndarray
    )  # unix times at which each interval ENDS (N-1 values)
    stats: dict  # mean, std, count, min, max (all in seconds)
    diagnostics: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)


def run(inputs: MtbfInputs) -> MtbfResult:
    if len(inputs.events) < 2:
        raise ValueError(
            f"inter-event intervals require at least 2 events; got {len(inputs.events)}"
        )
    times = np.array([e.t_event_unix_s for e in inputs.events], dtype=float)
    times.sort()
    intervals = np.diff(times)
    stats = {
        "mean_s": float(np.mean(intervals)),
        "std_s": float(np.std(intervals)),
        "count": int(len(intervals)),
        "min_s": float(np.min(intervals)),
        "max_s": float(np.max(intervals)),
    }
    print(
        f"[mtbf] events={len(inputs.events)} intervals={len(intervals)} "
        f"mean={stats['mean_s'] / 3600:.2f}h min={stats['min_s']:.1f}s max={stats['max_s'] / 3600:.1f}h"
    )
    return MtbfResult(
        intervals_s=intervals,
        event_times_unix_s=times[1:],
        stats=stats,
        meta=inputs.meta,
    )


def make_inputs_from_norm(norm: Norm) -> MtbfInputs:
    return MtbfInputs(
        events=norm.events["calibration"],
        meta=dict(norm["meta"]),
    )


# ---------------------------------------------------------------------------
# Mean time between calibrations
# ---------------------------------------------------------------------------


@dataclass
class Mtbc:
    """The restricted mean of the product-limit curve on cycles, with its SE and interval.

    `rmst_s` is the area under the curve on [0, `tau_s`], with `tau_s` the longest cycle in
    the table, censored or not (`survival.curve_summaries`' default, R's for one curve).
    `lower_s` and `upper_s` are rmst_s -/+ z rmst_se_s at `conf_level`; a lower end at or
    below 0 is reported as computed. The point is descriptive, with `n_censored` beside it;
    the SE and the interval rest on exchangeable cycles (`a2_renewal_cycles`), and
    `checks` carries the check outcome; empty means NOT ASSESSED, which is what a log's
    inter-event intervals always read, since no check runs on them.

    The SE is 0 when every death sits at tau, which is the case for a single completed
    cycle: R reports the same 0, and the interval then has no width and says nothing about
    the spread.
    """

    dataset_id: str
    rmst_s: float
    rmst_se_s: float
    lower_s: float
    upper_s: float
    tau_s: float
    conf_level: float
    n_cycles: int
    n_deaths: int
    n_censored: int
    n_unobserved_birth_dropped: int
    point_standing: str = DESCRIPTIVE
    interval_standing: str = RESTS_ON_EXCHANGEABILITY
    checks: CheckAttachment = field(
        default_factory=lambda: CheckAttachment(assumption_id=A2_RENEWAL_CYCLES_ID)
    )


def rmst_interval(summary: CurveSummaries, conf_level: float) -> tuple[float, float]:
    """The Wald interval of a restricted mean: rmst -/+ z SE, z two-sided at `conf_level`.

    The one definition: the MTBC (cycles) and the record figures' curve summaries (windows
    of one bag at one rung) both call it. A lower end at or below 0 is returned as
    computed. NaN in, NaN out.

    Validity assumptions

    - Assumption: the restricted mean is close to normal with the standard error
      `survival.curve_summaries` gives, which rests on exchangeable durations (cycles:
      `a2_renewal_cycles`; windows: `a1_renewal_durations`).
      Diagnostic: the check outcome the summary carries; an SE of 0, which a curve with no
      death before tau gives, makes an interval of no width that says nothing about the
      spread.
      Consequence of violation: coverage away from `conf_level`, worst with few deaths.
      Reference: no source located.
    """
    half_s = z_two_sided(conf_level) * summary.rmst_se_s
    return summary.rmst_s - half_s, summary.rmst_s + half_s


def mtbc(table: EventTable, *, conf_level: float = CONF_LEVEL) -> Mtbc:
    """MTBC from an event table of cycles: `curve_summaries`, plus rmst +/- z SE.

    Raises unless the table holds cycles: a restricted mean of windows on a side of a
    threshold is a mean window length, not a time between calibrations.
    """
    if table.side != CYCLES:
        raise ValueError(
            f"mtbc needs an event table of cycles (side={CYCLES!r}); got {table.side!r}"
        )
    summary = curve_summaries(table, conf_level=conf_level)
    lower_s, upper_s = rmst_interval(summary, conf_level)
    return Mtbc(
        dataset_id=table.dataset_id,
        rmst_s=summary.rmst_s,
        rmst_se_s=summary.rmst_se_s,
        lower_s=lower_s,
        upper_s=upper_s,
        tau_s=summary.rmst_tau_s,
        conf_level=conf_level,
        n_cycles=summary.n_windows,
        n_deaths=summary.n_deaths,
        n_censored=summary.n_censored,
        n_unobserved_birth_dropped=summary.n_unobserved_birth_dropped,
    )


def intervals_event_table(
    intervals_s: np.ndarray, *, dataset_id: str = ""
) -> EventTable:
    """Inter-event intervals as an event table of cycles, each one completed.

    A log of event times declares no record end, so no interval is censored. A 0 interval
    (two events at one instant) stays in, as a cycle of age 0.
    """
    age_s = np.asarray(intervals_s, dtype=float)
    return EventTable(
        dataset_id=dataset_id,
        threshold_label="",
        side=CYCLES,
        age_s=age_s,
        event=np.ones(len(age_s), dtype=bool),
    )


# ---------------------------------------------------------------------------
# Interval histogram
#
# One definition of the binning, used by BOTH the across-calibration panel's histogram
# subplot (via panels/_across_calibration_compute._interval_histogram) and the standalone
# mean-time-between-calibrations figure. Two copies would be two chances for the
# same bar to mean different things in two figures of the same log.
# ---------------------------------------------------------------------------


def log_interval_histogram(
    intervals_s: np.ndarray, n_bins: int = 50
) -> tuple[np.ndarray, np.ndarray]:
    """Log-spaced histogram of POSITIVE intervals; returns (counts, edges_s).

    Non-positive intervals cannot be placed on a log axis and are dropped here. Callers
    that render the result must state how many - `IntervalHistogramResult` carries the
    count for exactly that reason. Empty (counts, edges) when nothing is positive.
    """
    ivs = np.asarray(intervals_s, dtype=float)
    valid = ivs[ivs > 0]
    if len(valid) == 0:
        return np.array([]), np.array([])
    log_bins = np.logspace(
        np.log10(float(np.min(valid))), np.log10(float(np.max(valid))), n_bins
    )
    counts, edges = np.histogram(valid, bins=log_bins)
    return counts, edges


@dataclass
class IntervalHistogramResult:
    """Complete typed artifact for the single-axes interval histogram.

    Hours throughout (`_h`), because a calibration log spans seconds to weeks and hours
    is the only unit that keeps both ends of that range readable as plain numbers. The
    conversion happens HERE, in the step, not at draw time.

    `mean_h` is the arithmetic mean of every positive interval - the "mean time between
    calibrations" the figure is named for. It is reported beside `median_h` because on a
    log axis the two land in visibly different places, and a reader given only the mean
    cannot tell that.

    It leaves 0 s intervals out, as the bars do. `mtbc`, the MTBC the across-calibration
    panel shows, counts them as cycles of age 0, so the two agree only on a log with no 0 s
    interval.
    """

    counts: np.ndarray
    edges_h: np.ndarray
    mean_h: float
    median_h: float
    n_intervals: int
    n_nonpositive_dropped: int
    label: str = ""
    dataset_id: str = ""
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if len(self.counts) and len(self.edges_h) != len(self.counts) + 1:
            raise ValueError(
                "incomplete IntervalHistogramResult: edges_h must have one more entry "
                f"than counts ({len(self.edges_h)} vs {len(self.counts)})"
            )


def make_interval_histogram(
    result: MtbfResult, *, label: str = "", n_bins: int = 50
) -> IntervalHistogramResult:
    """Adapter: MtbfResult -> the typed artifact the histogram figure draws.

    The mean and median are taken over the SAME positive intervals the bars count, so the
    marked mean cannot sit at a value no bar contributed to.
    """
    intervals_s = np.asarray(result.intervals_s, dtype=float)
    positive = intervals_s[intervals_s > 0]
    if len(positive) == 0:
        raise ValueError(
            "interval histogram needs at least one positive inter-event interval; "
            f"got {len(intervals_s)} intervals, none positive"
        )
    counts, edges_s = log_interval_histogram(intervals_s, n_bins=n_bins)
    meta = dict(result.meta)
    return IntervalHistogramResult(
        counts=counts,
        edges_h=edges_s / 3600.0,
        mean_h=float(np.mean(positive)) / 3600.0,
        median_h=float(np.median(positive)) / 3600.0,
        n_intervals=int(len(positive)),
        n_nonpositive_dropped=int(len(intervals_s) - len(positive)),
        label=label,
        dataset_id=str(meta.get("dataset_id", "")),
        meta=meta,
    )
