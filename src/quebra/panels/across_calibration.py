"""Panel for the across-calibration intervals of a calibration log (repairable-system tier).

VOCABULARY. This project calls these two tiers `within-calibration` and
`across-calibration`, not `non-repairable` and `repairable`. The rename is ours and it is
deliberate: our tiers are separated by a CALIBRATION BOUNDARY, which is what the data
actually records, while the reliability literature separates them by whether a system is
restored after failure. The two distinctions coincide here - a calibration restores the
qubit, so an across-calibration record is a repairable-system record - but they are not the
same idea, and naming ours after the boundary keeps the code honest about what it measured.

The literature's terms are NOT renamed where they name the literature. `repairable system`
is standard usage from Ascher and Feingold and from Rigdon and Basu; rewriting it inside a
sentence about that theory would make the sentence false. So the term survives in prose
that cites the field, and nowhere else.


Two halves:
  - AcrossCalibrationPanelData is the COMPLETE typed artifact - raw inputs plus every
    data-derived quantity (elapsed days, intervals in hours, 14-day binned stats,
    log-spaced histogram, the MTBC with its interval). Built solely by
    panels._across_calibration_compute.build_across_calibration_panel_data.
  - AcrossCalibrationPanel is a PURE renderer: it reads fields and draws. No data
    arithmetic at draw time.

Renders a 4-subplot figure:
  - Top:           inter-event interval vs. elapsed days (scatter, log y), with the MTBC
                   and its interval
  - Middle:        14-day binned moving average (median, IQR, p90)
  - Bottom-left:   histogram on log x-axis (multimodal distribution); with
                   `show_dropped_count`, the count of 0-length intervals it cannot place

DEVIATION from docs/FIGURE_STANDARD.md, "Dropped data": with `show_dropped_count` False (the
default) the 0 s intervals the MTBC counts appear neither on the log-y scatter nor in the
histogram, and their count is only on the artifact (`n_nonpositive_dropped`).
  - Bottom-right:  summary text in human-readable time units
"""

from __future__ import annotations

from dataclasses import dataclass, field

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go

from quebra.analyzers.mtbf import Mtbc, MtbfResult
from quebra.core._artifact_guard import StaleArtifactGuard
from quebra.plots.base import BasePlot
from quebra.plots.theme import (
    ANNOTATION,
    apply_common_style,
    qubit_color,
    style_context,
)


def _empty() -> np.ndarray:
    return np.array([])


@dataclass
class AcrossCalibrationPanelData(StaleArtifactGuard):
    """Complete typed contract for AcrossCalibrationPanel.

    Built by build_across_calibration_panel_data (panels/_across_calibration_compute.py), the sole
    constructor path. `__post_init__` enforces completeness: constructing with
    derived arrays whose lengths are inconsistent with the raw inputs raises
    ValueError, so an incomplete artifact never materializes - go through the
    builder. Unpickling a stale pre-split artifact (one missing any field) raises
    ValueError via StaleArtifactGuard (that path bypasses __post_init__).
    """

    # raw inputs
    intervals_s: np.ndarray
    event_times_unix_s: np.ndarray
    stats: dict[str, object]
    meta: dict[str, object]

    # derived (populated by build_across_calibration_panel_data)
    elapsed_days: np.ndarray = field(default_factory=_empty)
    intervals_h: np.ndarray = field(default_factory=_empty)
    binned_interval_stats: tuple[np.ndarray, ...] = field(
        default_factory=lambda: (_empty(),) * 5
    )
    histogram_counts: np.ndarray = field(default_factory=_empty)
    histogram_edges: np.ndarray = field(default_factory=_empty)
    mtbc: Mtbc | None = None
    # (point, lower, upper) of the MTBC in hours, the scatter's axis unit.
    mtbc_h: tuple[float, float, float] = (float("nan"),) * 3
    # Intervals of 0 s (two events at one instant): counted in the MTBC, but a log axis
    # has no place for them, so the histogram leaves them out.
    n_nonpositive_dropped: int = 0
    # Draw that count on the histogram. A render choice, set as a step kwarg so it reaches
    # the provenance label.
    show_dropped_count: bool = False
    # The check outcome the MTBC interval rests on, as the line the panel prints.
    mtbc_check_summary: str = ""

    def __post_init__(self) -> None:
        # Completeness contract: the builder's derived arrays are length-consistent
        # with the raw inputs; direct construction (empty defaults) alongside
        # non-empty raw inputs is incomplete and must not materialize.
        # (__setstate__ handles the unpickle path and bypasses __post_init__.)
        if len(self.elapsed_days) != len(self.event_times_unix_s):
            raise ValueError(
                "incomplete AcrossCalibrationPanelData: elapsed_days length "
                f"{len(self.elapsed_days)} != event_times_unix_s length "
                f"{len(self.event_times_unix_s)} - construct via "
                "build_across_calibration_panel_data()"
            )
        if len(self.intervals_h) != len(self.intervals_s):
            raise ValueError(
                "incomplete AcrossCalibrationPanelData: intervals_h length "
                f"{len(self.intervals_h)} != intervals_s length "
                f"{len(self.intervals_s)} - construct via build_across_calibration_panel_data()"
            )
        if len(self.binned_interval_stats) != 5:
            raise ValueError(
                "incomplete AcrossCalibrationPanelData: binned_interval_stats must be a "
                "5-tuple (centers, median, q1, q3, p90)"
            )
        if len(self.histogram_counts) and (
            len(self.histogram_edges) != len(self.histogram_counts) + 1
        ):
            raise ValueError(
                "incomplete AcrossCalibrationPanelData: histogram_edges must have one more "
                "entry than histogram_counts"
            )
        if len(self.intervals_s) and self.mtbc is None:
            raise ValueError(
                "incomplete AcrossCalibrationPanelData: mtbc is missing - construct via "
                "build_across_calibration_panel_data()"
            )


def make_mtbf_panel_data(
    result: MtbfResult, *, show_dropped_count: bool = False
) -> AcrossCalibrationPanelData:
    """Adapter: MtbfResult → complete AcrossCalibrationPanelData (via the output-builder)."""
    # Local import breaks the module cycle (compute imports AcrossCalibrationPanelData here).
    from quebra.panels._across_calibration_compute import (
        build_across_calibration_panel_data,
    )

    return build_across_calibration_panel_data(
        intervals_s=result.intervals_s,
        event_times_unix_s=result.event_times_unix_s,
        stats=result.stats,
        meta=result.meta,
        show_dropped_count=show_dropped_count,
    )


def _human_time(seconds: float) -> str:
    # The unit follows the magnitude, so a negative interval end (a normal interval on a
    # short record) reads in the same unit as its positive twin.
    size = abs(seconds)
    if size < 120:
        return f"{seconds:.1f} s"
    if size < 7200:
        return f"{seconds / 60:.1f} min"
    if size < 172800:
        return f"{seconds / 3600:.2f} h"
    return f"{seconds / 86400:.2f} days"


class AcrossCalibrationPanel(BasePlot):
    """Interval and MTBC panel for calibration-event logs (pure renderer)."""

    def build_matplotlib(
        self, result: AcrossCalibrationPanelData, style: str = "default"
    ) -> plt.Figure:
        if not isinstance(result, AcrossCalibrationPanelData):
            raise TypeError("AcrossCalibrationPanel expects AcrossCalibrationPanelData")
        pd_ = result
        color = qubit_color(meta=pd_.meta)

        with style_context(style):
            fig = plt.figure(
                figsize=(14, 13), constrained_layout=True, facecolor="white"
            )
            gs = fig.add_gridspec(3, 2, height_ratios=[1.6, 1.2, 1.0])
            ax_scatter = fig.add_subplot(gs[0, :])
            ax_roll = fig.add_subplot(gs[1, :])
            ax_hist = fig.add_subplot(gs[2, 0])
            ax_text = fig.add_subplot(gs[2, 1])

            for ax in (ax_scatter, ax_roll, ax_hist):
                apply_common_style(ax)

            self._draw_scatter(ax_scatter, pd_, color)
            self._draw_rolling(ax_roll, pd_, color)
            self._draw_histogram(ax_hist, pd_, color)
            self._draw_summary(ax_text, pd_)

        return fig

    def build_plotly(self, result: object) -> go.Figure:
        raise NotImplementedError(f"{self.__class__.__name__} has no plotly backend")

    @staticmethod
    def _draw_scatter(
        ax: plt.Axes, pd_: AcrossCalibrationPanelData, color: str
    ) -> None:
        ax.scatter(
            pd_.elapsed_days,
            pd_.intervals_h,
            s=4,
            alpha=0.8,
            color=color,
            linewidths=0,
        )
        ax.set_yscale("log")
        AcrossCalibrationPanel._draw_mtbc(ax, pd_)
        ax.set_xlabel("Elapsed time (days)")
        ax.set_ylabel("Inter-event interval (h)")
        ax.set_title(
            f"Calibration intervals over time "
            f"(qubit {pd_.meta.get('qubit', '?')}, {pd_.meta.get('device', '')})"
        )
        ax.legend(frameon=False, fontsize=8)
        ax.grid(True, which="both", color="lightgray", alpha=0.4)

    @staticmethod
    def _draw_mtbc(ax: plt.Axes, pd_: AcrossCalibrationPanelData) -> None:
        """The MTBC line and its interval span, on the scatter's log-hours axis."""
        mtbc = pd_.mtbc
        if mtbc is None:
            return
        point_h, lower_h, upper_h = pd_.mtbc_h
        ax.axhline(
            point_h,
            color="gray",
            linestyle="--",
            linewidth=1.2,
            label=f"MTBC = {_human_time(mtbc.rmst_s)}",
        )
        label = (
            f"{mtbc.conf_level:.0%} interval {_human_time(mtbc.lower_s)} - "
            f"{_human_time(mtbc.upper_s)}"
        )
        if pd_.mtbc_check_summary:
            label += f"; {pd_.mtbc_check_summary}"
        if lower_h <= 0.0:
            # A log axis has no place for a lower end at or below 0: the span starts at the
            # axis floor, and the label says so. The limits are frozen first, so the floor
            # the span starts from is the floor drawn. The number itself is in the artifact.
            floor_h, top_h = ax.get_ylim()
            ax.set_ylim(floor_h, max(top_h, upper_h))
            label += " (lower end <= 0, drawn from the axis floor)"
            lower_h = floor_h
        ax.axhspan(lower_h, upper_h, color="gray", alpha=0.15, linewidth=0, label=label)

    @staticmethod
    def _draw_rolling(
        ax: plt.Axes, pd_: AcrossCalibrationPanelData, color: str
    ) -> None:
        centers, medians, q1s, q3s, p90s = pd_.binned_interval_stats
        if len(centers) == 0:
            ax.axis("off")
            return
        ax.plot(
            centers, medians, "-", linewidth=1.4, color=color, label="Median", zorder=2
        )
        # A line needs two consecutive finite points, so a populated bin with no populated
        # neighbour draws nothing on the median or p90 line (its IQR fill is a hairline).
        # Mark both there, so no populated bin vanishes from either line.
        finite = np.isfinite(medians)
        lone = finite & ~np.r_[False, finite[:-1]] & ~np.r_[finite[1:], False]
        if lone.any():
            ax.plot(
                centers[lone], medians[lone], "o", markersize=3, color=color, zorder=2
            )
            ax.plot(
                centers[lone],
                p90s[lone],
                "o",
                markersize=2.5,
                color=color,
                alpha=0.55,
                zorder=2,
            )
        ax.fill_between(
            centers, q1s, q3s, color=color, alpha=0.2, zorder=1, label="IQR"
        )
        ax.plot(
            centers,
            p90s,
            "--",
            linewidth=0.8,
            color=color,
            alpha=0.55,
            label="p90",
            zorder=1,
        )
        ax.set_xlabel("Elapsed time (days)")
        ax.set_ylabel("Interval (h)")
        ax.set_title("14-day binned statistics (median, IQR, p90)")
        ax.legend(frameon=False, fontsize=8)
        ax.grid(True, alpha=0.25)

    @staticmethod
    def _draw_histogram(
        ax: plt.Axes, pd_: AcrossCalibrationPanelData, color: str
    ) -> None:
        counts, edges = pd_.histogram_counts, pd_.histogram_edges
        if len(counts) == 0:
            ax.axis("off")
            return
        ax.stairs(counts, edges, fill=True, color=color, alpha=0.75)
        if pd_.show_dropped_count:
            ax.text(
                0.98,
                0.95,
                f"excluded: {pd_.n_nonpositive_dropped} of {len(pd_.intervals_s)} "
                "intervals (0 s, no place on a log axis)",
                transform=ax.transAxes,
                ha="right",
                va="top",
                **ANNOTATION,
            )
        ax.set_xscale("log")
        ax.set_xlabel("Interval (s)")
        ax.set_ylabel("Count")
        ax.set_title("Distribution of inter-event intervals")
        ax.grid(True, which="both", color="lightgray", alpha=0.8)

    @staticmethod
    def _draw_summary(ax: plt.Axes, pd_: AcrossCalibrationPanelData) -> None:
        s = pd_.stats
        mtbc = pd_.mtbc
        mtbc_lines = (
            [
                f"MTBC:        {_human_time(mtbc.rmst_s)}",
                f"{mtbc.conf_level:.0%} interval: {_human_time(mtbc.lower_s)} - "
                f"{_human_time(mtbc.upper_s)}",
                f"tau:         {_human_time(mtbc.tau_s)}",
            ]
            if mtbc is not None
            else []
        )
        lines = [
            f"dataset:     {pd_.meta.get('dataset_id', '?')}",
            f"qubit:       {pd_.meta.get('qubit', '?')}",
            f"device:      {pd_.meta.get('device', '?')}",
            f"n events:    {pd_.meta.get('n_events', s['count'] + 1)}",
            f"n intervals: {s['count']}",
            "",
            *mtbc_lines,
            f"std:         {_human_time(s['std_s'])}",
            f"min:         {_human_time(s['min_s'])}",
            f"max:         {_human_time(s['max_s'])}",
        ]
        ax.axis("off")
        ax.text(
            0.05,
            0.95,
            "\n".join(lines),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.5,
            family="monospace",
            bbox={
                "facecolor": "lightyellow",
                "edgecolor": "gray",
                "boxstyle": "round,pad=0.5",
            },
        )
