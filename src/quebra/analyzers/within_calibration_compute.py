"""Output-builder for WithinCalibrationPanel: computes the COMPLETE typed panel data.

The compute half of the panel. `build_within_calibration_panel_data` takes the carve's
window and read tables and the job's one Kaplan-Meier node, and returns a
WithinCalibrationPanelData with every data-derived field filled, so the renderer
(panels/within_calibration.py) draws from it and does no data arithmetic.

Computed here: the in-spec fraction of observed time, cumulative time out of spec, time
to first crossing, cumulative damage, the value, error and relative-error histograms and
the coefficient of variation; the distinguish band adds the 4-state timeline and the
shape statistics. Not computed here:

- classification: every in/out flag is the carve's `in_spec` column
  (`windows.in_spec_mask`);
- survival: the reliability band carries the Kaplan-Meier set it is given.

The carve and the panel must describe one record and one ladder. The panel's series may
be the carve's values times one positive display scale (µs for a T2* carve in s); the
builder raises when read times, per-label read counts, threshold values or directions
disagree beyond that scale, and when an in-spec window is not the run of reads the read
table gives it.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd

from quebra.analyzers import distinguish_band, reliability_band, signal_band
from quebra.analyzers import windows as _windows
from quebra.analyzers.kaplan_meier import KaplanMeierSet
from quebra.analyzers.within_calibration_data import WithinCalibrationPanelData

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


def _carved_rows(reads: pd.DataFrame, label: str, t_f_h: np.ndarray) -> pd.DataFrame:
    """The carve's read rows for `label`, checked to line up with the finite display reads.

    The carve drops the same non-finite pairs, so its rows for one label line up with the
    finite display reads one for one; a table that does not raises rather than
    misaligning.
    """
    sub = reads[reads["threshold_label"] == label]
    if len(sub) != len(t_f_h):
        raise ValueError(
            f"threshold {label!r}: the carve has {len(sub)} reads but the display series "
            f"has {len(t_f_h)} finite reads; they must come from the same record"
        )
    t_read_h = sub["t_read_s"].to_numpy(dtype=float) / 3600.0
    if not np.allclose(t_read_h, t_f_h, rtol=1e-12, atol=1e-12):
        raise ValueError(
            f"threshold {label!r}: the carve's read times do not match the display "
            f"series; they must come from the same record in the same order"
        )
    return sub


def _carved_in_spec(reads: pd.DataFrame, label: str, t_f_h: np.ndarray) -> np.ndarray:
    """The carve's per-read `in_spec` for `label`, aligned to the finite display reads.

    Classification happens once, in the carve's own units (`windows.in_spec_mask`). The
    panel's series may be rescaled for display, and reclassifying a rescaled value can
    move a read that sits on the threshold, so this reads the flags instead.
    """
    return _carved_rows(reads, label, t_f_h)["in_spec"].to_numpy(dtype=bool)


# Relative tolerance of the carve-to-panel comparison: well above the rounding of a
# display rescale (double precision, ~1e-16 relative), well below a mistyped threshold.
_DISPLAY_SCALE_RTOL = 1e-9


def _require_one_ladder(
    t_f_h: np.ndarray,
    s_f: np.ndarray,
    reads: pd.DataFrame,
    thresholds: list[tuple[str, float, bool]],
) -> None:
    """Raise unless the carve and the panel use one ladder, up to one display scale.

    The in/out flags come from the carve, while threshold lines and damage use the
    panel's ladder, so a label carved at one value and drawn at another would report one
    threshold's numbers against another's line. The panel may rescale for display, so
    the contract is one positive scale `c`: the finite display series is `c` times the
    carve's `value`, and for every label the panel's `windows.margin` is `c` times the
    carve's `margin`. Margin is affine in the value, so this pins each threshold's value
    and direction, on every label with reads, whichever side its windows fall on. Each
    entry is compared relative to its own operands, so one outlying read cannot widen
    the tolerance for the rest.
    """
    if len(s_f) < 2 or not thresholds:
        return  # below two finite reads the carve emits no read rows
    first = _carved_rows(reads, thresholds[0][0], t_f_h)
    parts = [("the series", first["value"].to_numpy(dtype=float), s_f, np.abs(s_f))]
    for label, thr_val, big_values_good in thresholds:
        rows = _carved_rows(reads, label, t_f_h)
        parts.append(
            (
                f"threshold {label!r}",
                rows["margin"].to_numpy(dtype=float),
                _windows.margin(s_f, thr_val, big_values_good),
                np.abs(s_f) + abs(float(thr_val)),
            )
        )
    # Least squares for `c` on the series alone, so a wrong rung is named alone. A record
    # whose carved values are all zero fixes nothing there and the margins stand in; with
    # every carved entry zero any `c` fits, and 1.0 stands in.
    basis = parts[:1] if np.any(parts[0][1]) else parts
    carve = np.concatenate([p[1] for p in basis])
    panel = np.concatenate([p[2] for p in basis])
    norm = float(np.dot(carve, carve))
    scale = float(np.dot(carve, panel)) / norm if norm > 0.0 else 1.0
    if not (np.isfinite(scale) and scale > 0.0):
        raise ValueError(
            f"the panel's series and ladder are no positive multiple of the carve's "
            f"(best scale {scale!r}); carve and render must use the same record and ladder"
        )
    bad = [
        name
        for name, carved, drawn, magnitude in parts
        if np.any(np.abs(drawn - scale * carved) > _DISPLAY_SCALE_RTOL * magnitude)
    ]
    if bad:
        raise ValueError(
            f"the carve and the panel disagree on {bad} at display scale {scale:.6g}: "
            f"carve and render must use the same series and the same ladder (values and "
            f"directions), up to one positive display scale"
        )


def _require_windows_match_reads(windows: pd.DataFrame, reads: pd.DataFrame) -> None:
    """Raise unless every in-spec window is the run of reads the read table gives it.

    The bands count windows from `windows` and occupancy and the timeline from `reads`,
    so the two must be one carve. Per (label, window_index), the two key sets must be
    equal, and the window's `t_birth_s`, `t_last_s` and `n_reads` must equal the first
    and last `t_read_s` and the read count of the reads carrying that index: one carve
    gives these from the same floats, so that comparison is exact. Its `t_death_s` must
    equal every such read's `t_read_s + forward_time_s` up to the rounding of that sum
    (4 ulp of the larger time), which pins the death and so the censoring and gap policy:
    two carves of the same reads under different `gap_mult` tile them alike but end
    windows differently.
    """
    keys = ["threshold_label", "window_index"]
    in_window = reads[reads["window_index"].notna()].astype({"window_index": "int64"})
    from_reads = in_window.groupby(keys)["t_read_s"].agg(
        t_birth_s="min", t_last_s="max", n_reads="size"
    )
    from_windows = windows.set_index(keys)[
        ["t_birth_s", "t_last_s", "n_reads", "t_death_s"]
    ].sort_index()
    tiled = ["t_birth_s", "t_last_s", "n_reads"]
    same = from_windows.index.equals(from_reads.index) and np.array_equal(
        from_windows[tiled].to_numpy(dtype=float), from_reads.to_numpy(dtype=float)
    )
    if same and len(in_window):
        per_read = in_window.join(from_windows["t_death_s"], on=keys)
        t_read = per_read["t_read_s"].to_numpy(dtype=float)
        t_death = per_read["t_death_s"].to_numpy(dtype=float)
        rebuilt = t_read + per_read["forward_time_s"].to_numpy(dtype=float)
        tol = 4.0 * np.finfo(float).eps * np.maximum(np.abs(t_read), np.abs(t_death))
        same = bool(np.all(np.abs(rebuilt - t_death) <= tol))
    if not same:
        raise ValueError(
            "the in-spec window table and the read table are not one carve: some "
            "window's first read, last read, read count or death differs from the reads "
            "that carry its window_index. Pass both tables from the same windows.run"
        )


def _cumulative_time_out_of_spec(
    t_h: np.ndarray,
    primary_series: np.ndarray,
    reads: pd.DataFrame,
    thresholds: list[tuple[str, float, bool]],
    gap_spans_h: list[tuple[float, float]],
) -> dict[str, np.ndarray]:
    """Left-Riemann cumulative time out of spec per threshold (OBSERVED hours)."""
    t = np.asarray(t_h, dtype=float)
    s = np.asarray(primary_series, dtype=float)
    mask = np.isfinite(t) & np.isfinite(s)
    t_f = t[mask]

    result: dict[str, np.ndarray] = {}
    for label, _, _ in thresholds:
        if len(t_f) < 2:
            result[label] = _to_full_length(np.zeros(len(t_f)), mask)
            continue
        oos = ~_carved_in_spec(reads, label, t_f)
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
        excess = np.maximum(-_windows.margin(s_f, thr_val, big_values_good), 0.0)
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
    reads: pd.DataFrame,
    thresholds: list[tuple[str, float, bool]],
) -> dict[str, float | None]:
    """Time to FIRST threshold crossing (elapsed hours from t[0]), per threshold.

    TTF, not MTTF: this is one crossing of one trace, not a mean over a population.
    Calling it a mean overclaimed a statistic that was never computed.
    """
    t = np.asarray(t_h, dtype=float)
    s = np.asarray(primary_series, dtype=float)
    mask = np.isfinite(t) & np.isfinite(s)
    t_f = t[mask]

    result: dict[str, float | None] = {}
    for label, _, _ in thresholds:
        if len(t_f) < 2:
            result[label] = None
            continue
        oos = ~_carved_in_spec(reads, label, t_f)
        indices = np.where(oos)[0]
        if len(indices) == 0:
            result[label] = None
        else:
            result[label] = float(t_f[indices[0]]) - float(t_f[0])
    return result


def _threshold_in_spec_frac(
    t_h: np.ndarray,
    primary_series: np.ndarray,
    reads: pd.DataFrame,
    thresholds: list[tuple[str, float, bool]],
    gap_spans_h: list[tuple[float, float]],
) -> dict[str, float]:
    """Occupancy: fraction of OBSERVED time in spec, per threshold.

    Denominator is observed time, not wall clock - see _observed_dt_h. This is the
    single occupancy definition in the repo; the reliability band exposes it as
    `occupancy` and the renderer's >=5% timeline cull reads the same number.

    NaN when there is no observed time to divide by (fewer than two finite reads, or
    every interval a gap or of zero length): not measured, which 0.0 would misstate as
    never in spec.
    """
    t = np.asarray(t_h, dtype=float)
    s = np.asarray(primary_series, dtype=float)
    mask = np.isfinite(t) & np.isfinite(s)
    t_f = t[mask]

    result: dict[str, float] = {}
    for label, _, _ in thresholds:
        if len(t_f) < 2:
            result[label] = float("nan")
            continue
        dt = _observed_dt_h(t_f, gap_spans_h)
        observed_h = float(np.sum(dt))
        if observed_h == 0.0:
            result[label] = float("nan")
            continue
        oos = ~_carved_in_spec(reads, label, t_f)
        oos_h = float(np.sum(dt[oos[:-1]]))
        result[label] = 1.0 - oos_h / observed_h
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
    kaplan_meier: KaplanMeierSet,
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
    (completeness); those flags only decide what the renderer draws.

    `windows` and `reads` are the in-spec window table and the read table from
    `analyzers/windows.py::run`, carved on the same reads under the same threshold labels
    as `primary_series` and `thresholds`. Values may differ by one positive display scale
    shared by the series and the whole ladder (the T2* carve runs in s, the panel in µs).
    The builder raises on a different series or ladder (`_require_one_ladder`) and on a
    window table that is not the read table's carve (`_require_windows_match_reads`).
    They are required: the panel does not carve, so there is exactly one carve in the
    repo and the gap policy, censoring and per-read state cannot diverge between the
    artifact and the figure.

    `gap_spans_s` is the carve's `WindowsResult.diagnostics["gap_spans_s"]`, and it is
    required too: the tables do not hold every gap, and a record read as gap-free draws
    through its holes and counts unobserved hours as observed. Pass `[]` for a record with
    no gaps.

    `kaplan_meier` is the job's one Kaplan-Meier node, estimated from the same carve; the
    builder never computes survival itself. The reliability band raises unless the set
    names the tables' record, holds every label and counts the same in-spec windows.

    `meta["dataset"]` is required and must be a non-empty string: the caption names the
    record it draws.

    `primary_sigma` is the per-read 1-sigma on `primary_series` (same units), used for
    error bars. None when the dataset carries no uncertainty.

    `damage_fn` is the descoped damage seam - see _cumulative_damage.
    """
    if gap_spans_s is None:
        raise ValueError(
            "gap_spans_s is required: pass WindowsResult.diagnostics['gap_spans_s'], "
            "or [] for a record with no gaps"
        )
    dataset = meta.get("dataset")
    if not isinstance(dataset, str) or not dataset.strip():
        raise ValueError(
            f"meta['dataset'] must name the record, and the caption shows it; got "
            f"{dataset!r}"
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
    finite = np.isfinite(t_arr) & np.isfinite(s_arr)
    _require_one_ladder(t_arr[finite], s_arr[finite], reads, list(thresholds))
    _require_windows_match_reads(windows, reads)

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
            kaplan_meier=kaplan_meier,
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
