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

**Hazard** (`HazardBlockPlot`). One axes per bag: every rung's Nelson-Aalen cumulative
hazard with its band, on log axes shared by both bags, the counts per rung in the legend.
Under each bag, the check outcome of each rung's band: checks down, and per rung its two
clocks side by side. The key and the caption read as in the survival block.
"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure, SubFigure
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import NullFormatter

from quebra.analyzers.record_figure import HazardBlockData, SurvivalBlockData
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


def _draw_key(
    ax: plt.Axes, data: SurvivalBlockData | HazardBlockData, style: str
) -> None:
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


def _draw_failure(
    fig: Figure | SubFigure, data: SurvivalBlockData | HazardBlockData, style: str
) -> bool:
    """Draw why a block has nothing to draw, and say whether it did."""
    if not (data.failure or not data.ladder):
        return False
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
    return True


def draw_survival(fig: Figure | SubFigure, data: SurvivalBlockData, style: str) -> None:
    """The survival block of one record on `fig` (module docstring)."""
    caption = theme.scaled_text(theme.CAPTION, style)
    if _draw_failure(fig, data, style):
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


# The hazard block: inches per bag's hazard axes and per check row of its grid, and per
# clock column of a rung in the grid.
_HAZARD_H = 3.0
# Points between the x axis and its tick labels on the log-log hazard axes, so the first x
# label clears the lowest y label at the corner.
_HAZARD_X_TICK_PAD = 6.0
# The legend row under each bag's hazard axes: per rung, two legend columns.
_HAZARD_LEGEND_ROW_H = 0.17
_HAZARD_GRID_BASE_H = 0.55
_HAZARD_CELL_W = 0.34
_HAZARD_MIN_W = 4.5


def _hazard_limits(
    data: HazardBlockData,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Decade limits holding every drawn hazard and its band, both bags, ages from the
    first death (H is 0 before it, which a log axis cannot draw)."""
    ages: list[float] = []
    levels: list[float] = []
    for key, item in data.hazards.items():
        if not item.n_deaths:
            continue
        ages += [float(data.time_min[key][1]), float(data.time_min[key][-1])]
        drawn = np.concatenate(
            [item.cumulative_hazard[1:], item.band_lower[1:], item.band_upper[1:]]
        )
        drawn = drawn[np.isfinite(drawn) & (drawn > 0)]
        levels += [float(drawn.min()), float(drawn.max())]

    def decades(
        values: list[float], default: tuple[float, float]
    ) -> tuple[float, float]:
        if not values:
            return default
        lo = 10.0 ** math.floor(math.log10(min(values)))
        hi = 10.0 ** math.ceil(math.log10(max(values)))
        return lo, max(hi, lo * 10.0)

    return decades(ages, (0.1, 10.0)), decades(levels, (0.01, 10.0))


def _draw_hazard_grid(
    ax: plt.Axes, data: HazardBlockData, side: str, names: bool, style: str
) -> None:
    """Each rung's band outcome: checks down, the rung's clocks side by side."""
    n_clocks = len(data.clocks)
    n_checks = len(data.run_set)
    for i, label in enumerate(data.ladder):
        outcome = data.outcomes[(label, side)]
        for k, grid in enumerate(outcome.grids):
            for j, verdict in enumerate(grid.verdicts[0]):
                ax.add_patch(
                    Rectangle(
                        (n_clocks * i + j, n_checks - 1 - k),
                        1.0,
                        1.0,
                        facecolor=theme.verdict_color(verdict),
                        edgecolor=theme.CELL_EDGE,
                        linewidth=0.5,
                        gid=CELL_GID,
                    )
                )
    ax.set_xlim(0.0, n_clocks * len(data.ladder))
    ax.set_ylim(0.0, max(n_checks, 1))
    ax.set_xticks([n_clocks * (i + 0.5) for i in range(len(data.ladder))])
    ax.set_xticklabels(list(data.ladder), **theme.scaled_text(theme.ANNOTATION, style))
    first = data.outcomes[(data.ladder[0], side)]
    ax.set_yticks([n_checks - 1 - k + 0.5 for k in range(n_checks)])
    ax.set_yticklabels(
        [_short_check(grid.label) for grid in first.grids] if names else [],
        **theme.scaled_text(theme.ANNOTATION, style),
    )
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)


def _hazard_caption(data: HazardBlockData) -> str:
    parts = [
        f"Nelson-Aalen cumulative hazard of the windows of each bag at each rung of "
        f"{data.dataset_id}, with its {data.conf_level:.0%} band symmetric on log H; a "
        "rung with no death draws no line.",
        "Under each bag, the check outcome of each rung's band: checks down, and for each "
        "rung its side clock, then its calendar clock.",
        "In the legend: n windows, d deaths, and the windows excluded because their birth "
        "is unobserved.",
    ]
    if data.flag:
        parts.append(f"This record is {data.flag}.")
    return " ".join(parts)


def draw_hazard(fig: Figure | SubFigure, data: HazardBlockData, style: str) -> None:
    """The hazard block of one record on `fig` (module docstring)."""
    title = f"Nelson-Aalen cumulative hazard by rung, {data.dataset_id}"
    if _draw_failure(fig, data, style):
        fig.suptitle(title)
        return
    n = len(data.ladder)
    grid = fig.add_gridspec(
        3,
        len(data.sides) + 1,
        width_ratios=[1.0] * len(data.sides) + [0.35],
        height_ratios=[
            _HAZARD_H,
            _HAZARD_LEGEND_ROW_H * math.ceil(n / 2),
            _HAZARD_GRID_BASE_H + _GRID_ROW_H * len(data.run_set),
        ],
    )
    (x_lo, x_hi), (y_lo, y_hi) = _hazard_limits(data)
    legend = theme.scaled_text(theme.LEGEND_TEXT, style)
    first: plt.Axes | None = None
    for s, side in enumerate(data.sides):
        ax = fig.add_subplot(grid[0, s], sharex=first, sharey=first)
        first = first or ax
        for i, label in enumerate(data.ladder):
            item = data.hazards[(label, side)]
            color = theme.threshold_color(i, n)
            ax.plot(
                [],
                [],
                color=color,
                label=(
                    f"{label}: n {item.n_windows}, d {item.n_deaths}, "
                    f"excluded {item.n_unobserved_birth_dropped}"
                ),
            )
            if not item.n_deaths:
                continue
            x = data.time_min[(label, side)][1:]
            ax.fill_between(
                *band_vertices(x, item.band_lower[1:], item.band_upper[1:]),
                color=color,
                alpha=theme.BAND_STYLE["fill_alpha"],
                linewidth=0.0,
            )
            ax.step(x, item.cumulative_hazard[1:], where="post", color=color)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.tick_params(axis="x", pad=_HAZARD_X_TICK_PAD)
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        ax.set_xlabel("Window age (min)")
        if s:
            ax.tick_params(axis="y", labelleft=False)
        else:
            ax.set_ylabel("Cumulative hazard")
        ax.set_title(theme.BAG_LABELS[side], loc="left")
        theme.apply_common_style(ax)
        # The legend in a row of its own, so it never covers a curve.
        key = fig.add_subplot(grid[1, s])
        key.set_axis_off()
        handles, labels = ax.get_legend_handles_labels()
        key.legend(handles, labels, loc="upper left", ncol=2, frameon=False, **legend)
        _draw_hazard_grid(fig.add_subplot(grid[2, s]), data, side, s == 0, style)
    _draw_key(fig.add_subplot(grid[:, len(data.sides)]), data, style)
    fig.suptitle(title)
    fig.supxlabel(
        _hazard_caption(data), wrap=True, **theme.scaled_text(theme.CAPTION, style)
    )


class HazardBlockPlot(BasePlot):
    """The hazard block of one record: both bags' cumulative hazards at every rung, each
    band with its check outcome under its bag."""

    def build_matplotlib(self, result: object, style: str = "default") -> plt.Figure:
        if not isinstance(result, HazardBlockData):
            raise TypeError(
                f"{type(self).__name__} draws a HazardBlockData; got "
                f"{type(result).__name__}"
            )
        data = result
        n = max(len(data.ladder), 1)
        bag_w = max(_HAZARD_MIN_W, _HAZARD_CELL_W * len(data.clocks) * n)
        height = (
            _HAZARD_H
            + _HAZARD_LEGEND_ROW_H * math.ceil(n / 2)
            + _HAZARD_GRID_BASE_H
            + _GRID_ROW_H * len(data.run_set)
            + 0.9
        )
        with theme.style_context(style):
            fig = plt.figure(
                figsize=(bag_w * len(data.sides) + _KEY_W, height),
                layout="constrained",
                facecolor="white",
            )
            draw_hazard(fig, data, style)
        return fig
