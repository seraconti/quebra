"""The blocks of a record figure, one `BasePlot` per block, each a pure renderer.

Each block draws one typed artifact of `analyzers/record_figure.py` through a function of a
figure or subfigure (`draw_<block>(fig, data, style)`), so the block's own PDF and any page
that places it beside others draw it the same way. What is left here is a function of axes
and theme: limits, ticks, colours and where text sits.

**Survival** (`SurvivalBlockPlot`). Small multiples, one column per rung: for each bag, the
Kaplan-Meier curves with their bands on one log age axis shared by every column, a line of
counts under each (windows, deaths, censored, and the windows excluded because their birth
is unobserved, named in the caption), then the check outcome of exactly that band, checks down and clocks across,
each cell in its verdict's colour. A key column names the verdicts. The caption says what
is drawn and, for a record the paper's dataset tables shade, carries its flag.
"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure, SubFigure
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import NullFormatter

from quebra.analyzers.record_figure import SurvivalBlockData
from quebra.plots import theme
from quebra.plots.base import BasePlot
from quebra.plots.km_survival_plot import band_vertices, step_arrays

# Heights in inches of one bag's two rows: the curves with their counts, and the check
# grid, whose rows add a fixed pitch per check.
_CURVE_H = 2.4
_GRID_BASE_H = 0.3
_GRID_ROW_H = 0.19
# Width in inches of one rung's column, and of the key column.
_COLUMN_W = 1.9
_KEY_W = 1.7
# The id each check cell carries, so a test can find the cells of one band.
CELL_GID = "check_cell"


def _short_check(label: str) -> str:
    """`C5 rank autocorrelation (studentized)` -> `C5`: the check's code, enough beside
    the key and the caption that name each one."""
    return label.split()[0]


def _age_limits(data: SurvivalBlockData) -> tuple[float, float]:
    """Decade limits that hold every drawn curve's first positive age and its longest."""
    firsts: list[float] = []
    lasts: list[float] = []
    for curve in data.curves.values():
        if not curve.n_windows:
            continue
        positive = curve.time_min[curve.time_min > 0]
        if len(positive):
            firsts.append(float(np.min(positive)))
        lasts.append(float(curve.max_observed_min))
    lo = 10.0 ** math.floor(math.log10(min(firsts))) if firsts else 0.1
    hi = 10.0 ** math.ceil(math.log10(max(lasts))) if lasts and max(lasts) > 0 else 10.0
    return lo, max(hi, lo * 10.0)


def _decade_ticks(lo: float, hi: float) -> list[float]:
    """At most three decade ticks, so a narrow column's labels do not collide."""
    first, last = round(math.log10(lo)), round(math.log10(hi))
    step = max(1, math.ceil((last - first) / 2))
    return [10.0**d for d in range(first, last + 1, step)]


def _counts_lines(data: SurvivalBlockData, label: str, side: str) -> list[str]:
    curve = data.curves[(label, side)]
    return [
        f"n {curve.n_windows}, d {curve.n_deaths}, c {curve.n_censored}",
        f"excluded {curve.n_unobserved_birth_dropped}",
    ]


def _draw_grid(
    ax: plt.Axes,
    data: SurvivalBlockData,
    label: str,
    side: str,
    names: bool,
    style: str,
) -> None:
    """The check outcome of the band above: one row per check, one column per clock."""
    outcome = data.outcomes[(label, side)]
    n = len(outcome.grids)
    for i, grid in enumerate(outcome.grids):
        y = n - 1 - i
        for j, verdict in enumerate(grid.verdicts[0]):
            ax.add_patch(
                Rectangle(
                    (j, y),
                    1.0,
                    1.0,
                    facecolor=theme.verdict_color(verdict),
                    edgecolor=theme.CELL_EDGE,
                    linewidth=0.5,
                    gid=CELL_GID,
                )
            )
    ax.set_xlim(0.0, len(data.clocks))
    ax.set_ylim(0.0, max(n, 1))
    ax.set_xticks([])
    ax.set_yticks([n - 1 - i + 0.5 for i in range(n)])
    ax.set_yticklabels(
        [_short_check(grid.label) for grid in outcome.grids] if names else [],
        **theme.scaled_text(theme.ANNOTATION, style),
    )
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)


def _draw_key(ax: plt.Axes, data: SurvivalBlockData, style: str) -> None:
    ax.set_axis_off()
    ax.legend(
        handles=[
            Patch(facecolor=theme.verdict_color(v), edgecolor=theme.CELL_EDGE, label=v)
            for v in data.verdict_key
        ],
        loc="center left",
        frameon=False,
        title="check outcome\ncolumns: " + ", ".join(f"{c} clock" for c in data.clocks),
        **theme.scaled_text(theme.LEGEND_TEXT, style),
    )


def _caption(data: SurvivalBlockData) -> str:
    parts = [
        f"Kaplan-Meier survival of the windows of each bag at each rung of {data.dataset_id}, "
        f"with its {data.conf_level:.0%} log-log band; ticks mark censored windows.",
        "Under each band, the check outcome of that rung and bag: checks down, clocks "
        "across.",
        "Under each curve: n windows, d deaths, c censored, and the windows excluded "
        "because their birth is unobserved.",
    ]
    if data.flag:
        parts.append(f"This record is {data.flag}.")
    return " ".join(parts)


def draw_survival(fig: Figure | SubFigure, data: SurvivalBlockData, style: str) -> None:
    """The survival block of one record on `fig` (module docstring)."""
    caption = theme.scaled_text(theme.CAPTION, style)
    if data.failure or not data.ladder:
        ax = fig.add_subplot(1, 1, 1)
        ax.set_axis_off()
        reason = data.failure or "no rung on the ladder"
        ax.text(
            0.5,
            0.5,
            f"not computed: {reason}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            wrap=True,
            **theme.scaled_text(theme.ANNOTATION, style),
        )
        fig.suptitle(f"Window survival by rung, {data.dataset_id}")
        return
    n = len(data.ladder)
    bag_rows = [_CURVE_H, _GRID_BASE_H + _GRID_ROW_H * len(data.run_set)]
    grid = fig.add_gridspec(
        2 * len(data.sides),
        n + 1,
        width_ratios=[_COLUMN_W] * n + [_KEY_W],
        height_ratios=bag_rows * len(data.sides),
    )
    x_lo, x_hi = _age_limits(data)
    ticks = _decade_ticks(x_lo, x_hi)
    annotation = theme.scaled_text(theme.ANNOTATION, style)
    first: plt.Axes | None = None
    for s, side in enumerate(data.sides):
        for i, label in enumerate(data.ladder):
            ax = fig.add_subplot(grid[2 * s, i], sharex=first, sharey=first)
            first = first or ax
            curve = data.curves[(label, side)]
            color = theme.threshold_color(i, n)
            if curve.n_windows:
                x, surv, lo, hi = step_arrays(curve, x_lo)
                ax.fill_between(
                    *band_vertices(x, lo, hi),
                    color=color,
                    alpha=theme.BAND_STYLE["fill_alpha"],
                    linewidth=0.0,
                )
                ax.step(x, surv, where="post", color=color)
                if len(curve.censor_time_min):
                    ax.plot(
                        np.maximum(curve.censor_time_min, x_lo),
                        curve.censor_survival,
                        linestyle="none",
                        color=color,
                        **theme.CENSOR_TICK,
                    )
            else:
                ax.text(
                    0.5,
                    0.5,
                    "no window",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    **annotation,
                )
            ax.set_xscale("log")
            ax.set_xlim(x_lo, x_hi)
            ax.set_xticks(ticks)
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.set_ylim(-0.02, 1.05)
            ax.set_yticks([0.0, 0.5, 1.0])
            # The counts sit under the age label, so the layout sizes them with it.
            ax.set_xlabel(
                "\n".join(["Window age (min)", *_counts_lines(data, label, side)])
            )
            if i:
                ax.tick_params(axis="y", labelleft=False)
            else:
                ax.set_ylabel(f"Surviving fraction,\n{theme.BAG_LABELS[side]}")
            ax.set_title(label)
            theme.apply_common_style(ax)
            _draw_grid(
                fig.add_subplot(grid[2 * s + 1, i]), data, label, side, i == 0, style
            )
    _draw_key(fig.add_subplot(grid[:, n]), data, style)
    fig.suptitle(f"Window survival by rung, {data.dataset_id}")
    fig.supxlabel(_caption(data), wrap=True, **caption)


class SurvivalBlockPlot(BasePlot):
    """The survival block of one record: both bags' curves at every rung, each band with
    its check outcome under it."""

    def build_matplotlib(self, result: object, style: str = "default") -> plt.Figure:
        if not isinstance(result, SurvivalBlockData):
            raise TypeError(
                f"{type(self).__name__} draws a SurvivalBlockData; got "
                f"{type(result).__name__}"
            )
        data = result
        n = max(len(data.ladder), 1)
        bag_h = _CURVE_H + _GRID_BASE_H + _GRID_ROW_H * len(data.run_set)
        with theme.style_context(style):
            fig = plt.figure(
                figsize=(_COLUMN_W * n + _KEY_W, bag_h * len(data.sides) + 1.2),
                layout="constrained",
                facecolor="white",
            )
            draw_survival(fig, data, style)
        return fig
