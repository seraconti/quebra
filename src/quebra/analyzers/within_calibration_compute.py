"""Output-builder for WithinCalibrationPanel: computes the COMPLETE typed panel data.

This is the compute half of the panel. It owns every data-derived quantity -
cumulative time out of spec, cumulative damage, TTF, per-threshold window stats,
in-spec window survival, 30-min binned stats, CV, and in-spec fractions - so that
the materialized WithinCalibrationPanelData artifact is complete and the renderer
(panels/within_calibration.py) is a pure function of it (no data arithmetic at draw
time; only axis/theme concerns stay there).

Arithmetic here is ported verbatim from the pre-split draw-time methods; parity is
asserted array-equal against golden baselines captured on the known-good commit.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd

from quebra.analyzers import distinguish_band, reliability_band, signal_band
from quebra.analyzers.within_calibration_data import WithinCalibrationPanelData

# ---------------------------------------------------------------------------
# Direction-agnostic threshold primitives
# ---------------------------------------------------------------------------


def _out_of_spec_mask(
    series: np.ndarray, threshold_value: float, big_values_good: bool
) -> np.ndarray:
    if not big_values_good:
        return series > threshold_value
    return series < threshold_value


def _excess(
    series: np.ndarray, threshold_value: float, big_values_good: bool
) -> np.ndarray:
    if not big_values_good:
        return np.maximum(series - threshold_value, 0.0)
    return np.maximum(threshold_value - series, 0.0)


# ---------------------------------------------------------------------------
# Scalar / distribution helpers (no domain assumptions)
# ---------------------------------------------------------------------------


def _compute_cv(series: np.ndarray) -> float:
    s = np.asarray(series, dtype=float)
    s = s[np.isfinite(s)]
    if len(s) == 0:
        return np.nan
    mean = float(np.mean(s))
    if mean == 0.0:
        return np.nan
    return float(np.std(s) / abs(mean))


def _analyze_threshold_windows(
    t_h: np.ndarray, series: np.ndarray, threshold_value: float
) -> dict[str, object]:
    """Factual above/below window statistics; direction-agnostic."""
    t = np.asarray(t_h, dtype=float)
    s = np.asarray(series, dtype=float)
    mask = np.isfinite(t) & np.isfinite(s)
    t, s = t[mask], s[mask]
    if len(t) < 2:
        return {"above": {}, "below": {}}

    dt = np.diff(t) * 60.0  # minutes
    above = s >= threshold_value
    windows_above: list[float] = []
    windows_below: list[float] = []
    current = 0.0
    in_above = bool(above[0]) if len(above) > 0 else False

    for i in range(len(above) - 1):
        current += float(dt[i])
        if above[i] != above[i + 1]:
            (windows_above if in_above else windows_below).append(current)
            current = 0.0
            in_above = above[i + 1]
    # NO tail accumulation here. The loop above already consumes ALL of dt
    # (len(dt) == len(above) - 1), so adding dt[-1] again inflated every trace's
    # final window by one read interval - an all-in-spec 2-hour span reported 180
    # minutes instead of 120.
    (windows_above if in_above else windows_below).append(current)

    def _stats(windows: list[float]) -> dict[str, object]:
        if not windows:
            return {
                "longest": np.nan,
                "mean": np.nan,
                "median": np.nan,
                "p90": np.nan,
                "count": 0,
            }
        w = np.asarray(windows, dtype=float)
        return {
            "longest": float(np.max(w)),
            "mean": float(np.mean(w)),
            "median": float(np.median(w)),
            "p90": float(np.percentile(w, 90)),
            "count": len(windows),
        }

    return {"above": _stats(windows_above), "below": _stats(windows_below)}


def _finite_stat(values: np.ndarray | None, fn: Callable[[np.ndarray], float]) -> float:
    """`fn` over the finite entries; NaN when there are none (never silently zero)."""
    if values is None:
        return float("nan")
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    return float(fn(v)) if len(v) else float("nan")


def _hist_view_limit(
    counts: np.ndarray, edges: np.ndarray, quantile: float = 0.99
) -> tuple[float, int]:
    """Upper x-limit holding `quantile` of the mass, and the count left outside it.

    A handful of failed fits carry errors orders of magnitude above the bulk and would
    flatten the histogram against the left edge. The count is data - it is rendered as
    text and a reader treats it as a fact - so it is computed here, not at draw time.
    """
    if len(counts) == 0:
        return float("nan"), 0
    cumulative = np.cumsum(counts)
    total = int(cumulative[-1])
    if total <= 0:
        return float("nan"), 0
    inside = int(np.searchsorted(cumulative, quantile * total, side="left")) + 1
    x_max = float(edges[min(inside, len(edges) - 1)])
    if x_max <= float(edges[0]):
        return float("nan"), 0
    n_outside = total - int(cumulative[min(inside - 1, len(cumulative) - 1)])
    return x_max, int(n_outside)


def _value_histogram(
    values: np.ndarray | None, n_bins: int = 40
) -> tuple[np.ndarray, np.ndarray]:
    """Histogram of the finite entries of `values`; returns (counts, edges).

    Empty (counts, edges) when there is nothing finite to bin - the renderer treats
    that as "no distribution to draw" rather than crashing on an empty axis.
    """
    if values is None:
        return np.array([], dtype=int), np.array([], dtype=float)
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) == 0 or float(np.max(v)) == float(np.min(v)):
        return np.array([], dtype=int), np.array([], dtype=float)
    counts, edges = np.histogram(v, bins=n_bins)
    return counts, edges


def _window_survival(windows_min: list[float]) -> list[tuple[float, float]]:
    w = np.asarray(windows_min, dtype=float)
    w = w[np.isfinite(w)]
    if len(w) == 0:
        return []
    unique_w = np.unique(np.sort(w))
    return [
        (
            float(length),
            float(np.clip(np.round(np.sum(w >= length) / len(w), 10), 0.0, 1.0)),
        )
        for length in unique_w
    ]


# ---------------------------------------------------------------------------
# Window-table consumers
#
# Carving lives in analyzers/windows.py and reaches the builder as two DataFrames.
# The panel never carves: it reads what the step produced, so the gap policy,
# censoring and window identity are the same facts everywhere they are used.
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Per-threshold derived series/scalars
# ---------------------------------------------------------------------------


def _observed_dt_h(
    t_f: np.ndarray, gap_spans_h: list[tuple[float, float]]
) -> np.ndarray:
    """Inter-read intervals with GAP intervals zeroed out.

    An interval that spans a read gap was not observed: the instrument was not
    reporting. Counting it credits unobserved hours to whichever state happened to
    hold at the left edge, which for a 30-minute gap opening on an in-spec read means
    30 in-spec minutes the record never contains. The figure already refuses to draw
    across a gap; this makes the numbers agree with it.
    """
    dt_h = np.diff(t_f)
    if not gap_spans_h or len(dt_h) == 0:
        return dt_h
    observed = dt_h.copy()
    for lo_h, hi_h in gap_spans_h:
        # Zero by OVERLAP, not by index identity. Matching the gap's left edge against
        # a read timestamp failed two ways: the boundary read may have been dropped by
        # the finite mask (a failed fit is exactly what tends to precede an instrument
        # gap), and np.isclose's relative tolerance is ~12 minutes at this project's
        # time base (t ~ 21000 h), which could zero the wrong interval entirely.
        observed[(t_f[:-1] < hi_h) & (t_f[1:] > lo_h)] = 0.0
    return observed


def _to_full_length(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Put a finite-subset result back on the full row index, NaN where rows were dropped.

    Callers hold the unmasked row index, so a subset-length return is silently misaligned
    against it. NaN rather than 0.0 because these are cumulative quantities: a zero would
    drop the curve back to the origin at each missing read.
    """
    full = np.full(mask.shape, np.nan, dtype=float)
    full[mask] = values
    return full


def _cumulative_time_out_of_spec(
    t_h: np.ndarray,
    primary_series: np.ndarray,
    thresholds: list[tuple[str, float, bool]],
    gap_spans_h: list[tuple[float, float]],
) -> dict[str, np.ndarray]:
    """Left-Riemann cumulative time out of spec per threshold (OBSERVED hours)."""
    t = np.asarray(t_h, dtype=float)
    s = np.asarray(primary_series, dtype=float)
    mask = np.isfinite(t) & np.isfinite(s)
    t_f, s_f = t[mask], s[mask]

    result: dict[str, np.ndarray] = {}
    for label, thr_val, big_values_good in thresholds:
        if len(t_f) < 2:
            result[label] = _to_full_length(np.zeros(len(t_f)), mask)
            continue
        oos = _out_of_spec_mask(s_f, thr_val, big_values_good)
        dt_h = _observed_dt_h(t_f, gap_spans_h)
        increments = oos[:-1].astype(float) * dt_h
        cum = np.zeros(len(t_f))
        cum[1:] = np.cumsum(increments)
        result[label] = _to_full_length(cum, mask)
    return result


def _cumulative_damage(
    t_h: np.ndarray,
    primary_series: np.ndarray,
    thresholds: list[tuple[str, float, bool]],
    damage_fn: Callable[[np.ndarray], np.ndarray] | None,
    gap_spans_h: list[tuple[float, float]],
) -> dict[str, np.ndarray]:
    """Trapezoidal cumulative damage per threshold (primary_unit · h).

    With the default (damage_fn=None → identity) this is the cumulative integral of
    the excess-over-threshold - a real derived curve, NOT a no-op, distinct from the
    cumulative time out of spec (which integrates a 0/1 mask).

    # EXTENSION: future DamageModel. `damage_fn` is the seam for a user-defined,
    # job-declared damage meaning. It is a plain builder parameter ONLY - never routed
    # through node kwargs / identity / provenance labels, because a callable cannot be
    # hashed deterministically or labeled stably. A future modular DamageModel would
    # carry a stable versioned id (e.g. "linear_v1") that folds into identity instead.
    """
    t = np.asarray(t_h, dtype=float)
    s = np.asarray(primary_series, dtype=float)
    mask = np.isfinite(t) & np.isfinite(s)
    t_f, s_f = t[mask], s[mask]

    apply_damage: Callable[[np.ndarray], np.ndarray] = (
        damage_fn if damage_fn is not None else (lambda x: x)
    )

    result: dict[str, np.ndarray] = {}
    for label, thr_val, big_values_good in thresholds:
        if len(t_f) < 2:
            result[label] = _to_full_length(np.zeros(len(t_f)), mask)
            continue
        excess = _excess(s_f, thr_val, big_values_good)
        damage_rate = apply_damage(excess)
        # Gap intervals are zeroed here too: trapezoidal integration would otherwise
        # accrue damage across hours the instrument was not reporting.
        dt_h = _observed_dt_h(t_f, gap_spans_h)
        trap_steps = 0.5 * (damage_rate[:-1] + damage_rate[1:]) * dt_h
        cum = np.zeros(len(t_f))
        cum[1:] = np.cumsum(trap_steps)
        result[label] = _to_full_length(cum, mask)
    return result


def _ttf(
    t_h: np.ndarray,
    primary_series: np.ndarray,
    thresholds: list[tuple[str, float, bool]],
) -> dict[str, float | None]:
    """Time to FIRST threshold crossing (elapsed hours from t[0]), per threshold.

    TTF, not MTTF: this is one crossing of one trace, not a mean over a population.
    Calling it a mean overclaimed a statistic that was never computed.
    """
    t = np.asarray(t_h, dtype=float)
    s = np.asarray(primary_series, dtype=float)
    mask = np.isfinite(t) & np.isfinite(s)
    t_f, s_f = t[mask], s[mask]

    result: dict[str, float | None] = {}
    for label, thr_val, big_values_good in thresholds:
        oos = _out_of_spec_mask(s_f, thr_val, big_values_good)
        indices = np.where(oos)[0]
        if len(indices) == 0 or len(t_f) == 0:
            result[label] = None
        else:
            result[label] = float(t_f[indices[0]]) - float(t_f[0])
    return result


def _threshold_in_spec_frac(
    t_h: np.ndarray,
    primary_series: np.ndarray,
    thresholds: list[tuple[str, float, bool]],
    gap_spans_h: list[tuple[float, float]],
) -> dict[str, float]:
    """Occupancy: fraction of OBSERVED time in spec, per threshold.

    Denominator is observed time, not wall clock - see _observed_dt_h. This is the
    single occupancy definition in the repo; the reliability band exposes it as
    `occupancy` and the renderer's >=5% timeline cull reads the same number.
    """
    t = np.asarray(t_h, dtype=float)
    s = np.asarray(primary_series, dtype=float)
    mask = np.isfinite(t) & np.isfinite(s)
    t_f, s_f = t[mask], s[mask]

    result: dict[str, float] = {}
    for label, thr_val, big_values_good in thresholds:
        if len(t_f) < 2:
            result[label] = 0.0
            continue
        dt = _observed_dt_h(t_f, gap_spans_h)
        observed_h = float(np.sum(dt))
        if observed_h == 0.0:
            result[label] = 0.0
            continue
        oos = _out_of_spec_mask(s_f, thr_val, big_values_good)
        oos_h = float(np.sum(dt[oos[:-1]]))
        result[label] = 1.0 - oos_h / observed_h
    return result


def _threshold_summary(
    t_h: np.ndarray,
    primary_series: np.ndarray,
    thresholds: list[tuple[str, float, bool]],
    gap_spans_h: list[tuple[float, float]],
) -> dict[str, dict[str, float] | None]:
    """Per-threshold out-of-spec summary (time_oos_h, frac_oos_pct).

    Value is None where the summary skips the threshold (fewer than 2 finite points),
    matching the pre-split summary's `continue`.
    """
    t = np.asarray(t_h, dtype=float)
    s = np.asarray(primary_series, dtype=float)
    mask = np.isfinite(t) & np.isfinite(s)
    t_f, s_f = t[mask], s[mask]

    result: dict[str, dict[str, float] | None] = {}
    for label, thr_val, big_values_good in thresholds:
        if len(t_f) < 2:
            result[label] = None
            continue
        oos = _out_of_spec_mask(s_f, thr_val, big_values_good)
        # Same observed-time denominator as _threshold_in_spec_frac: these two land
        # on the same band and a reader compares them, so they must not disagree.
        dt = _observed_dt_h(t_f, gap_spans_h)
        total_h = float(np.sum(dt))
        time_oos_h = float(np.sum(dt[oos[:-1]])) if len(dt) > 0 else 0.0
        frac_oos = 100.0 * time_oos_h / total_h if total_h > 0 else 0.0
        result[label] = {"time_oos_h": time_oos_h, "frac_oos_pct": frac_oos}
    return result


# ---------------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------------


def build_within_calibration_panel_data(
    *,
    t_h: np.ndarray,
    primary_series: np.ndarray,
    primary_label: str,
    thresholds: list[tuple[str, float, bool]],
    meta: dict[str, object],
    windows: pd.DataFrame,
    reads: pd.DataFrame,
    gap_spans_s: list[tuple[float, float]],
    primary_sigma: np.ndarray | None = None,
    traces: list[tuple[str, np.ndarray]] | None = None,
    use_log_scale: bool = False,
    color: object = None,
    damage_fn: Callable[[np.ndarray], np.ndarray] | None = None,
    include_cumulative_time: bool = True,
    include_cumulative_damage: bool = True,
    include_ttf: bool = True,
    shape_min_reads: int = 5,
    xi_seed: int = 0,
    k: float = 1.0,
    use_uncertainty: bool = False,
) -> WithinCalibrationPanelData:
    """Compute every data-derived quantity and return a COMPLETE WithinCalibrationPanelData.

    This is the sole constructor path: the renderer assumes the derived fields are
    populated. All derived quantities are computed regardless of the include_* flags
    (completeness); those flags only gate what the renderer draws.

    `windows` and `reads` are the tables from `analyzers/windows.py::run`, carved on the
    SAME series in the same units as `thresholds`. They are required: the panel does not
    carve, so there is exactly one carve in the repo and the gap policy, censoring and
    per-read state cannot diverge between the artifact and the figure.

    `gap_spans_s` is the carve's `WindowsResult.diagnostics["gap_spans_s"]`, and it is
    required too: the tables do not hold every gap, and a record read as gap-free draws
    through its holes and counts unobserved hours as observed. Pass `[]` for a record with
    no gaps.

    `primary_sigma` is the per-read 1-sigma on `primary_series` (same units), used for
    error bars. None when the dataset carries no uncertainty.

    `damage_fn` is the descoped damage seam - see _cumulative_damage.
    """
    if gap_spans_s is None:
        raise ValueError(
            "gap_spans_s is required: pass WindowsResult.diagnostics['gap_spans_s'], "
            "or [] for a record with no gaps"
        )
    t_arr = np.asarray(t_h, dtype=float)
    s_arr = np.asarray(primary_series, dtype=float)
    sigma_arr = (
        np.asarray(primary_sigma, dtype=float) if primary_sigma is not None else None
    )
    if sigma_arr is not None and len(sigma_arr) != len(s_arr):
        raise ValueError(
            f"primary_sigma must match primary_series in length; got "
            f"{len(sigma_arr)} and {len(s_arr)}"
        )
    # The tables are joined to the panel by threshold_label. A carve run on a
    # different ladder than the panel renders would otherwise pass silently, showing
    # n_windows=0, an empty survival curve and a blank timeline instead of failing.
    if len(reads):
        carved_labels = set(reads["threshold_label"].unique())
        missing = [label for label, _, _ in thresholds if label not in carved_labels]
        if missing:
            raise ValueError(
                f"the window tables were carved on a different ladder than this panel "
                f"renders: no carved reads for threshold(s) {missing}. Carve and render "
                f"must use the same threshold labels."
            )

    signal = signal_band.run(
        signal_band.make_inputs_from_windows(
            t_h=t_arr,
            values=s_arr,
            sigma=sigma_arr,
            reads=reads,
            gap_spans_s=gap_spans_s,
        )
    )
    distinguish = distinguish_band.run(
        distinguish_band.make_inputs_from_windows(
            reads=reads,
            windows=windows,
            thresholds=list(thresholds),
            median_read_spacing_s=signal.median_read_spacing_s,
            gap_spans_h=signal.gap_spans_h,
            sigma_display=sigma_arr,
            shape_min_reads=shape_min_reads,
            xi_seed=xi_seed,
            k=k,
            use_uncertainty=use_uncertainty,
        )
    )
    reliability = reliability_band.run(
        reliability_band.make_inputs_from_windows(
            t_h=t_arr,
            values=s_arr,
            reads=reads,
            windows=windows,
            thresholds=list(thresholds),
            gap_spans_h=signal.gap_spans_h,
            damage_fn=damage_fn,
        )
    )
    return WithinCalibrationPanelData(
        signal=signal,
        distinguish=distinguish,
        reliability=reliability,
        meta=meta,
        thresholds=list(thresholds),
        primary_label=primary_label,
        traces=traces,
        use_log_scale=use_log_scale,
        color=color,
        include_cumulative_time=include_cumulative_time,
        include_cumulative_damage=include_cumulative_damage,
        include_ttf=include_ttf,
    )
