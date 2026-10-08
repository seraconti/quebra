"""Band 1 of the composed within-calibration panel: the signal and its distributions.

What was measured, before any threshold is applied. The raw series against the scan
clock, the marginal distribution of the metric, the distribution of the per-read fit
error, and the distribution of the relative error. Nothing here knows what "in spec"
means.

Units: `values` and the histogram edges derived from it are in the panel's DISPLAY
units - µs for T2*, matching `primary_label` - because this band is what the series axis
draws. Times carry SI suffixes (`t_h` in hours, `median_read_spacing_s`, `max_gap_s`),
so the two are never confused.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from quebra.core._artifact_guard import StaleArtifactGuard


@dataclass
class SignalBand(StaleArtifactGuard):
    """Complete contract for the signal band. No `slots=True` - see StaleArtifactGuard."""

    # the series the panel draws, in display units
    t_h: np.ndarray
    values: np.ndarray
    # (t_before_h, t_after_h) per read gap; the trace is broken across these. Required:
    # a band read as gap-free draws through its holes, so there is no empty default.
    gap_spans_h: list[tuple[float, float]]
    # per-read 1-sigma in the same units; None when the dataset carries no error column
    sigma: np.ndarray | None = None

    # marginal distribution of the metric; mean and population std over finite reads, in
    # display units, for the figure's caption
    value_hist_counts: np.ndarray = field(default_factory=lambda: np.array([]))
    value_hist_edges: np.ndarray = field(default_factory=lambda: np.array([]))
    value_mean: float = field(default_factory=lambda: float("nan"))
    value_std: float = field(default_factory=lambda: float("nan"))
    # distribution of the per-read fit error
    sigma_hist_counts: np.ndarray = field(default_factory=lambda: np.array([]))
    sigma_hist_edges: np.ndarray = field(default_factory=lambda: np.array([]))
    sigma_hist_view_x_max: float = field(default_factory=lambda: float("nan"))
    sigma_hist_n_above_view: int = 0
    sigma_mean: float = field(default_factory=lambda: float("nan"))
    sigma_median: float = field(default_factory=lambda: float("nan"))
    # distribution of sigma / value, dimensionless
    relerr_hist_counts: np.ndarray = field(default_factory=lambda: np.array([]))
    relerr_hist_edges: np.ndarray = field(default_factory=lambda: np.array([]))
    relerr_mean: float = field(default_factory=lambda: float("nan"))
    relerr_median: float = field(default_factory=lambda: float("nan"))
    # reads excluded from the relative-error histogram: value <= 0 or sigma not finite
    relerr_n_excluded: int = 0

    n_reads: int = 0
    n_unknown_error: int = 0
    median_read_spacing_s: float = field(default_factory=lambda: float("nan"))
    max_gap_s: float = field(default_factory=lambda: float("nan"))
    n_gaps: int = 0

    cv: float = field(default_factory=lambda: float("nan"))

    def __post_init__(self) -> None:
        if len(self.t_h) != len(self.values):
            raise ValueError(
                f"signal band: t_h and values must be the same length; got "
                f"{len(self.t_h)} and {len(self.values)}"
            )
        if self.sigma is not None and len(self.sigma) != len(self.values):
            raise ValueError(
                f"signal band: sigma must match values in length; got "
                f"{len(self.sigma)} and {len(self.values)}"
            )


@dataclass(slots=True)
class SignalBandInputs:
    t_h: np.ndarray
    values: np.ndarray
    sigma: np.ndarray | None
    reads: pd.DataFrame
    gap_spans_s: list[tuple[float, float]]
    median_read_spacing_s: float
    n_bins: int = 40


def median_read_spacing_s(reads: pd.DataFrame) -> float:
    """Median positive spacing between distinct read timestamps.

    Derived from the read table rather than taken from the carve diagnostics so the
    band needs nothing but its two tables. Positive-only for the same reason the carve
    is: duplicate timestamps would otherwise drive it to zero.
    """
    if not len(reads):
        return float("nan")
    t = np.unique(reads["t_read_s"].to_numpy(dtype=float))
    if len(t) < 2:
        return float("nan")
    dt = np.diff(t)
    positive = dt[dt > 0]
    return float(np.median(positive)) if len(positive) else float("nan")


def make_inputs_from_windows(
    *,
    t_h: np.ndarray,
    values: np.ndarray,
    sigma: np.ndarray | None,
    reads: pd.DataFrame,
    gap_spans_s: list[tuple[float, float]],
    n_bins: int = 40,
) -> SignalBandInputs:
    """Build inputs from the display series plus `analyzers.windows.run` output."""
    return SignalBandInputs(
        t_h=np.asarray(t_h, dtype=float),
        values=np.asarray(values, dtype=float),
        sigma=None if sigma is None else np.asarray(sigma, dtype=float),
        reads=reads,
        gap_spans_s=list(gap_spans_s),
        median_read_spacing_s=median_read_spacing_s(reads),
        n_bins=n_bins,
    )


def run(inputs: SignalBandInputs) -> SignalBand:
    from quebra.analyzers.within_calibration_compute import (
        _compute_cv,
        _finite_stat,
        _hist_view_limit,
        _value_histogram,
    )

    values = inputs.values
    sigma = inputs.sigma

    value_hist = _value_histogram(values, inputs.n_bins)
    sigma_hist = _value_histogram(sigma, inputs.n_bins)
    sigma_view = _hist_view_limit(*sigma_hist)

    # Relative error is only defined where the metric is positive and the error is
    # known. Excluded reads are counted, never silently turned into inf.
    if sigma is None:
        relerr = None
        n_excluded = int(len(values))
    else:
        usable = np.isfinite(sigma) & np.isfinite(values) & (values > 0.0)
        relerr = sigma[usable] / values[usable]
        n_excluded = int((~usable).sum())
    relerr_hist = _value_histogram(relerr, inputs.n_bins)

    max_gap_s = (
        max((hi - lo) for lo, hi in inputs.gap_spans_s) if inputs.gap_spans_s else 0.0
    )
    n_unknown = (
        int((~inputs.reads["sigma_known"]).sum())
        if len(inputs.reads)
        else int(len(values))
    )
    # The read table repeats every read once per threshold; the unknown-error count is
    # a property of the read, not of the ladder.
    n_thresholds = (
        int(inputs.reads["threshold_label"].nunique()) if len(inputs.reads) else 1
    )
    n_unknown = n_unknown // max(n_thresholds, 1)

    print(
        f"[signal_band] reads={len(values)} unknown_error={n_unknown} "
        f"gaps={len(inputs.gap_spans_s)} max_gap={max_gap_s:.1f}s",
        flush=True,
    )
    return SignalBand(
        t_h=np.asarray(inputs.t_h, dtype=float),
        values=values,
        sigma=sigma,
        gap_spans_h=[(lo / 3600.0, hi / 3600.0) for lo, hi in inputs.gap_spans_s],
        value_hist_counts=value_hist[0],
        value_hist_edges=value_hist[1],
        value_mean=_finite_stat(values, np.mean),
        value_std=_finite_stat(values, np.std),
        sigma_hist_counts=sigma_hist[0],
        sigma_hist_edges=sigma_hist[1],
        sigma_hist_view_x_max=sigma_view[0],
        sigma_hist_n_above_view=sigma_view[1],
        sigma_mean=_finite_stat(sigma, np.mean),
        sigma_median=_finite_stat(sigma, np.median),
        relerr_hist_counts=relerr_hist[0],
        relerr_hist_edges=relerr_hist[1],
        relerr_mean=_finite_stat(relerr, np.mean),
        relerr_median=_finite_stat(relerr, np.median),
        relerr_n_excluded=n_excluded,
        n_reads=int(len(values)),
        n_unknown_error=n_unknown,
        median_read_spacing_s=inputs.median_read_spacing_s,
        max_gap_s=float(max_gap_s),
        n_gaps=len(inputs.gap_spans_s),
        cv=_compute_cv(values),
    )
