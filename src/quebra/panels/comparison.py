"""Generic N-series comparison panel.

Overlays an arbitrary number of labeled ``(x_h, y)`` series on one axis - for
comparing a metric series (e.g. T2*) across datasets/qubits. Deliberately
minimal (overlay + legend only); richer comparison panels can come later. The
``series`` list is arbitrary length, so all-qubit comparisons need no API change.
"""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go

from quebra.plots.base import BasePlot
from quebra.plots.theme import apply_common_style, qubit_color, style_context


@dataclass
class CompareSeriesData:
    """Input contract for CompareSeriesPanel.

    series  : list of (label, x_h, y) - one per dataset/qubit, arbitrary length.
              x_h is elapsed time in hours; y is the metric (same unit across series).
    x_label : x-axis label (default "Elapsed time (h)").
    y_label : y-axis label (metric + unit).
    title   : figure title.
    """

    series: list[tuple[str, np.ndarray, np.ndarray]]
    x_label: str = "Elapsed time (h)"
    y_label: str = ""
    title: str = ""


class CompareSeriesPanel(BasePlot):
    """Overlay comparison of a metric series across N datasets/qubits."""

    def build_matplotlib(
        self, result: CompareSeriesData, style: str = "default"
    ) -> plt.Figure:
        if not isinstance(result, CompareSeriesData):
            raise TypeError("CompareSeriesPanel expects CompareSeriesData")
        with style_context(style):
            fig = plt.figure(
                figsize=(14, 6), constrained_layout=True, facecolor="white"
            )
            ax = fig.add_subplot(1, 1, 1)
            apply_common_style(ax)
            for label, x_h, y in result.series:
                ax.plot(
                    np.asarray(x_h, dtype=float),
                    np.asarray(y, dtype=float),
                    linewidth=1.0,
                    alpha=0.85,
                    label=label,
                    color=qubit_color(dataset_id=label),
                )
            ax.set_xlabel(result.x_label)
            ax.set_ylabel(result.y_label)
            ax.set_title(result.title)
            ax.grid(True, alpha=0.3)
            if result.series:
                # Outside-right legend (captured by bbox_inches="tight" at save).
                ax.legend(
                    frameon=False,
                    fontsize=8,
                    loc="upper left",
                    bbox_to_anchor=(1.01, 1.0),
                    borderaxespad=0.0,
                )
        return fig

    def build_plotly(self, result: object) -> go.Figure:
        raise NotImplementedError(f"{self.__class__.__name__} has no plotly backend")
