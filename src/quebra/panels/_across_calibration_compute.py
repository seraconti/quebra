"""Output-builder for AcrossCalibrationPanel: computes the COMPLETE typed panel data.

Owns every data-derived quantity the panel plots - elapsed days, intervals in
hours, 14-day binned interval stats, and the log-spaced histogram (counts + edges)
- so the materialized AcrossCalibrationPanelData artifact is complete and the renderer
(panels/across_calibration.py) is a pure function of it. Elapsed days, hours and the
histogram are the pre-split draw-time arithmetic. The binned statistics are not: they keep
every bin, an empty one as NaN, and close the final bin on the right.
"""

from __future__ import annotations

import numpy as np

from quebra.analyzers.mtbf import log_interval_histogram
from quebra.panels.across_calibration import AcrossCalibrationPanelData


def _binned_interval_stats(
    elapsed_days: np.ndarray,
    intervals_h: np.ndarray,
    bin_days: float = 14.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Bin intervals by elapsed days; return (centers, median, q1, q3, p90).

    Every bin between the first and last event is returned, and a bin with no events carries
    NaN in all four statistics. A line joining two points asserts something about the interval
    between them, so an empty bin must reach the renderer as a hole, where matplotlib breaks
    the median, the IQR fill and the p90 line, rather than be skipped and bridged.

    The final bin is closed on the right, as in `np.histogram`. `np.digitize` puts a value on
    or past the last edge one past the final bin, which would drop the latest event whenever
    the span is an exact multiple of `bin_days`, or the float `arange` lands the last edge a
    rounding error below it.
    """
    if len(elapsed_days) == 0:
        return (np.array([]),) * 5
    t0, t1 = float(np.min(elapsed_days)), float(np.max(elapsed_days))
    edges = np.arange(t0, t1 + bin_days, bin_days)
    if len(edges) < 2:
        edges = np.array([t0, t1 + bin_days])
    idx = np.minimum(np.digitize(elapsed_days, edges) - 1, len(edges) - 2)
    centers, medians, q1s, q3s, p90s = [], [], [], [], []
    for i in range(len(edges) - 1):
        vals = intervals_h[idx == i]
        centers.append(0.5 * (edges[i] + edges[i + 1]))
        if len(vals) == 0:
            for column in (medians, q1s, q3s, p90s):
                column.append(float("nan"))
            continue
        medians.append(float(np.median(vals)))
        q1s.append(float(np.percentile(vals, 25)))
        q3s.append(float(np.percentile(vals, 75)))
        p90s.append(float(np.percentile(vals, 90)))
    return (
        np.asarray(centers),
        np.asarray(medians),
        np.asarray(q1s),
        np.asarray(q3s),
        np.asarray(p90s),
    )


def _interval_histogram(
    intervals_s: np.ndarray, n_bins: int = 50
) -> tuple[np.ndarray, np.ndarray]:
    """Log-spaced histogram of positive intervals; returns (counts, edges).

    Matches the pre-split ax.hist(valid, bins=logspace(...)) exactly (ax.hist counts
    identically to np.histogram). Empty (counts, edges) when no positive intervals.

    Delegates to `analyzers.mtbf.log_interval_histogram`, which is the single definition:
    the standalone mean-time-between-calibrations figure bins the same intervals, and two
    copies of a binning rule are two chances for the panel and the figure to disagree
    about what a bar means.
    """
    return log_interval_histogram(intervals_s, n_bins=n_bins)


def build_across_calibration_panel_data(
    *,
    intervals_s: np.ndarray,
    event_times_unix_s: np.ndarray,
    stats: dict[str, object],
    meta: dict[str, object],
    bin_days: float = 14.0,
    hist_bins: int = 50,
) -> AcrossCalibrationPanelData:
    """Compute every data-derived quantity and return a COMPLETE AcrossCalibrationPanelData.

    Sole constructor path: the renderer assumes the derived fields are populated.
    """
    intervals_arr = np.asarray(intervals_s, dtype=float)
    times_arr = np.asarray(event_times_unix_s, dtype=float)

    if len(times_arr) > 0:
        t0 = float(times_arr[0])
        elapsed_days = (times_arr - t0) / 86400.0
    else:
        elapsed_days = np.array([])
    intervals_h = intervals_arr / 3600.0

    counts, edges = _interval_histogram(intervals_arr, n_bins=hist_bins)

    return AcrossCalibrationPanelData(
        intervals_s=intervals_arr,
        event_times_unix_s=times_arr,
        stats=stats,
        meta=meta,
        elapsed_days=elapsed_days,
        intervals_h=intervals_h,
        binned_interval_stats=_binned_interval_stats(
            elapsed_days, intervals_h, bin_days=bin_days
        ),
        histogram_counts=counts,
        histogram_edges=edges,
    )
