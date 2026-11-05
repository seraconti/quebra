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

**Fraction** (`FractionBlockPlot`). The time-weighted in-spec fraction at each rung as a
bar in the in-spec colour, its bounds as a whisker, drawn as the compliance survey draws
them: a dotted whisker from 0 to 1 where the fraction is computed and its bounds are not,
a grey cell where the fraction itself is not computed. The gaps left out of the observed
time are counted in the axes.

**Timelines** (`TimelineBlockPlot`). Read state along the scan clock, one row per drawn
rung: the four states with uncertainty above, the two the carve used below, a key of the
states over both and each rung's first read out of spec beside the lower one. The culled
rungs are named on the upper title when they fit on one line, otherwise counted there and
listed under the lower axes.

**Record figure** (`RecordFigurePlot`). One page per record: one constrained `SubFigure`
per block, in the page's order, each drawn by the block's own `draw_*` function at the
height its own figure has (`<block>_size`), the width the widest block's. One subfigure per
row, rather than one gridspec shared across rows, which collapsed constrained layout on a
measured prototype of this page.

**Curve summaries** (`SummariesTablePlot`). No drawing: `build_table` hands the step-built
frame to the `table` target unchanged.

**Resolvability** (`ResolvabilityBlockPlot`). Three panels. The spectrum on log-log axes
with the band edge (dotted), the edge level over the top tenth of the band (heavy), the noise
floor (dashed) and the top decade G10 reads (shaded), and no fitted line. The count ratio of
each rung across the decimation ladder, with its two limits as labelled guides. The crossing
placement counts of each rung at k = 1 and 2, stacked by kind. A group the tellability did
not compute shows its status instead: a panel in its place; a mark of the spectrum in the
caption, which names only the marks drawn and gives the status of each one not drawn. A
count ratio of exactly 0 has no place on the log axis: it is drawn as a down triangle on the
axis floor and named in the caption. `ResolvabilityTablePlot` hands the job's two stacked
tables to the `table` target unchanged.
"""

from __future__ import annotations

import math
from collections.abc import Callable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure, SubFigure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import NullFormatter

from quebra.analyzers.record_figure import (
    LIMIT_INDEPENDENT,
    LIMIT_SMOOTH,
    OCCUPANCY_DRAWN,
    FractionBlockData,
    RecordFigureData,
    HazardBlockData,
    ResolvabilityBlockData,
    SurvivalBlockData,
    TimelineBlockData,
)
from quebra.analyzers.windows import STATE_UNOBSERVED
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
    fig: Figure | SubFigure,
    data: SurvivalBlockData | HazardBlockData | FractionBlockData | TimelineBlockData,
    style: str,
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


def survival_size(data: SurvivalBlockData) -> tuple[float, float]:
    """The survival block's size in inches, its own figure's and its page row's."""
    n = max(len(data.ladder), 1)
    bag_h = _CURVE_H + _GRID_BASE_H + _GRID_ROW_H * len(data.run_set)
    return _COLUMN_W * n + _KEY_W, bag_h * len(data.sides) + 1.2


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
        with theme.style_context(style):
            fig = plt.figure(
                figsize=survival_size(data), layout="constrained", facecolor="white"
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


def hazard_size(data: HazardBlockData) -> tuple[float, float]:
    """The hazard block's size in inches, its own figure's and its page row's."""
    n = max(len(data.ladder), 1)
    bag_w = max(_HAZARD_MIN_W, _HAZARD_CELL_W * len(data.clocks) * n)
    height = (
        _HAZARD_H
        + _HAZARD_LEGEND_ROW_H * math.ceil(n / 2)
        + _HAZARD_GRID_BASE_H
        + _GRID_ROW_H * len(data.run_set)
        + 0.9
    )
    return bag_w * len(data.sides) + _KEY_W, height


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
        with theme.style_context(style):
            fig = plt.figure(
                figsize=hazard_size(data), layout="constrained", facecolor="white"
            )
            draw_hazard(fig, data, style)
        return fig


# The fraction block: inches per rung, the axes' height, and the bar's width in rung units.
_FRACTION_RUNG_W = 0.75
_FRACTION_H = 3.2
_FRACTION_BAR_W = 0.62
# The ids a fraction bar, a bounds whisker, a whisker of bounds not computed, and a
# not-computed cell carry, so a test can find each.
FRACTION_GID = "fraction"
BOUNDS_GID = "bounds"
NO_BOUNDS_GID = "bounds_not_computed"
NOT_COMPUTED_GID = "fraction_not_computed"
# Bounds not computed: the compliance survey's dotted whisker.
_NO_BOUNDS = {**theme.BOUNDS_WHISKER, "linestyle": ":"}


def _fraction_caption(data: FractionBlockData) -> str:
    k = f"{data.unsure_k:g}" if data.unsure_k is not None else "the carve's"
    parts = [
        f"Time-weighted in-spec fraction of the observed time of {data.dataset_id} at "
        "each rung, the gaps left out.",
        f"The whisker bounds it, counting every stretch the reads cannot settle at {k} "
        "sigma either way; dotted from 0 to 1 where the bounds are not computed, grey "
        "where the fraction is not.",
    ]
    if data.flag:
        parts.append(f"This record is {data.flag}.")
    return " ".join(parts)


def draw_fraction(fig: Figure | SubFigure, data: FractionBlockData, style: str) -> None:
    """The fraction block of one record on `fig` (module docstring)."""
    title = f"In-spec fraction by rung, {data.dataset_id}"
    if _draw_failure(fig, data, style):
        fig.suptitle(title)
        return
    ax = fig.add_subplot(1, 1, 1)
    grey = theme.verdict_color("not computed")
    for j, (fraction, lower, upper) in enumerate(
        zip(data.fraction, data.lower, data.upper, strict=True)
    ):
        if fraction is None:
            ax.add_patch(
                Rectangle(
                    (j - _FRACTION_BAR_W / 2, 0.0),
                    _FRACTION_BAR_W,
                    1.0,
                    facecolor=grey,
                    edgecolor="none",
                    gid=NOT_COMPUTED_GID,
                )
            )
        else:
            (bar,) = ax.bar(
                j,
                fraction,
                width=_FRACTION_BAR_W,
                color=theme.STATE_COLORS["in_spec"],
                linewidth=0,
            )
            bar.set_gid(FRACTION_GID)
        if lower is not None and upper is not None:
            ax.vlines(j, lower, upper, zorder=3, gid=BOUNDS_GID, **theme.BOUNDS_WHISKER)
        elif fraction is not None:
            ax.vlines(j, 0.0, 1.0, zorder=3, gid=NO_BOUNDS_GID, **_NO_BOUNDS)
    ax.set_xlim(-0.5, len(data.ladder) - 0.5)
    ax.set_ylim(0.0, 1.0)
    ax.set_xticks(range(len(data.ladder)))
    ax.set_xticklabels(list(data.ladder))
    ax.set_xlabel("Rung")
    ax.set_ylabel("In-spec fraction of observed time")
    gap = (
        f"gaps left out: {data.n_gaps} ({data.gap_h:.3g} h)"
        if data.n_gaps is not None and data.gap_h is not None
        else "gaps left out: not counted"
    )
    ax.set_title(gap, loc="right", **theme.scaled_text(theme.ANNOTATION, style))
    ax.legend(
        handles=[
            Patch(facecolor=theme.STATE_COLORS["in_spec"], label="in-spec fraction"),
            Line2D([], [], **theme.BOUNDS_WHISKER, label="bounds"),
            Line2D([], [], **_NO_BOUNDS, label="bounds not computed"),
            Patch(facecolor=grey, label="not computed"),
        ],
        loc="lower left",
        bbox_to_anchor=(0.0, 1.0),
        ncol=4,
        frameon=False,
        borderaxespad=0.2,
        **theme.scaled_text(theme.LEGEND_TEXT, style),
    )
    theme.apply_common_style(ax)
    fig.suptitle(title)
    fig.supxlabel(
        _fraction_caption(data), wrap=True, **theme.scaled_text(theme.CAPTION, style)
    )


def fraction_size(data: FractionBlockData) -> tuple[float, float]:
    """The fraction block's size in inches, its own figure's and its page row's."""
    width = max(_FRACTION_RUNG_W * max(len(data.ladder), 1) + 1.5, 6.0)
    return width, _FRACTION_H + 1.2


class FractionBlockPlot(BasePlot):
    """The fraction block of one record: the in-spec fraction at every rung with its
    bounds."""

    def build_matplotlib(self, result: object, style: str = "default") -> plt.Figure:
        if not isinstance(result, FractionBlockData):
            raise TypeError(
                f"{type(self).__name__} draws a FractionBlockData; got "
                f"{type(result).__name__}"
            )
        data = result
        with theme.style_context(style):
            fig = plt.figure(
                figsize=fraction_size(data), layout="constrained", facecolor="white"
            )
            draw_fraction(fig, data, style)
        return fig


# The timeline block: inches per drawn rung's row, per axes' frame and labels, the key
# row, the first-out column, and the figure's width. Culled names longer than
# `_CULL_TITLE_CHARS` move off the title into a list under the lower axes.
_TIMELINE_ROW_H = 0.28
_TIMELINE_FRAME_H = 0.8
_TIMELINE_KEY_H = 0.45
_TIMELINE_FIRST_OUT_W = 1.1
_TIMELINE_W = 12.0
_CULL_TITLE_CHARS = 70
_BAR_H = 0.75
# The states a timeline draws, in the key's order, with the key's words.
_TIMELINE_STATES = (
    ("in_spec", "in spec"),
    ("in_spec_uncertain", "in spec, error overlaps threshold, or no error"),
    ("out_of_spec", "out of spec"),
    ("out_of_spec_uncertain", "out of spec, error overlaps threshold, or no error"),
)


def _cull_lines(data: TimelineBlockData) -> tuple[str, str]:
    """The title's cull clause and, when the names do not fit there, the list under the
    lower axes."""
    parts = []
    if data.culled_below:
        parts.append(
            f"below {OCCUPANCY_DRAWN:.0%} in spec: {', '.join(data.culled_below)}"
        )
    if data.culled_unmeasured:
        parts.append(f"not measured: {', '.join(data.culled_unmeasured)}")
    if not parts:
        return "", ""
    named = f"culled, {'; '.join(parts)}"
    if len(named) <= _CULL_TITLE_CHARS:
        return f" ({named})", ""
    n = len(data.culled_below) + len(data.culled_unmeasured)
    return f" ({n} culled, listed below)", f"Culled, {'; '.join(parts)}"


def _draw_states(
    ax: plt.Axes,
    data: TimelineBlockData,
    states: dict[str, tuple[tuple[float, float, str], ...]],
) -> None:
    n = len(data.drawn)
    for y, label in zip(range(n - 1, -1, -1), data.drawn, strict=True):
        for t0, t1, state in states[label]:
            unobserved = state == STATE_UNOBSERVED
            ax.barh(
                y,
                t1 - t0,
                left=t0,
                height=_BAR_H,
                color=theme.state_color(state),
                edgecolor=theme.UNOBSERVED_EDGE if unobserved else "none",
                hatch=theme.UNOBSERVED_HATCH if unobserved else None,
                linewidth=0.4 if unobserved else 0.0,
                gid=f"state:{label}",
            )
    ax.set_ylim(-0.5, n - 0.5)
    ax.set_yticks(range(n))
    ax.set_yticklabels(list(data.drawn)[::-1])
    ax.set_xlim(data.t_min_h, data.t_max_h)
    ax.grid(True, axis="x", alpha=0.3)


def draw_timelines(
    fig: Figure | SubFigure, data: TimelineBlockData, style: str
) -> None:
    """The timeline block of one record on `fig` (module docstring)."""
    title = f"Read state along the scan clock by rung, {data.dataset_id}"
    if _draw_failure(fig, data, style):
        fig.suptitle(title)
        return
    annotation = theme.scaled_text(theme.ANNOTATION, style)
    rows = _TIMELINE_ROW_H * max(len(data.drawn), 1) + _TIMELINE_FRAME_H
    grid = fig.add_gridspec(
        3,
        2,
        width_ratios=[_TIMELINE_W - _TIMELINE_FIRST_OUT_W, _TIMELINE_FIRST_OUT_W],
        height_ratios=[_TIMELINE_KEY_H, rows, rows],
    )
    key = fig.add_subplot(grid[0, :])
    key.set_axis_off()
    handles = [
        Patch(facecolor=theme.state_color(state), label=name)
        for state, name in _TIMELINE_STATES
    ]
    handles.append(
        Patch(
            facecolor=theme.state_color(STATE_UNOBSERVED),
            edgecolor=theme.UNOBSERVED_EDGE,
            hatch=theme.UNOBSERVED_HATCH,
            label="unobserved (read gap)",
        )
    )
    key.legend(
        handles=handles,
        loc="center left",
        ncol=3,
        frameon=False,
        **theme.scaled_text(theme.LEGEND_TEXT, style),
    )
    upper = fig.add_subplot(grid[1, 0])
    lower = fig.add_subplot(grid[2, 0], sharex=upper)
    first_out = fig.add_subplot(grid[2, 1], sharey=lower)
    fig.add_subplot(grid[1, 1]).set_axis_off()
    clause, listed = _cull_lines(data)
    upper.set_title(f"Four states, with uncertainty{clause}", loc="left")
    lower.set_title("Two states, as carved", loc="left")
    first_out.set_axis_off()
    if not data.drawn:
        for ax in (upper, lower):
            ax.set_yticks([])
        upper.text(
            0.5,
            0.5,
            f"no rung reaches {OCCUPANCY_DRAWN:.0%} in spec",
            transform=upper.transAxes,
            ha="center",
            va="center",
            **annotation,
        )
    else:
        _draw_states(upper, data, data.states_4)
        _draw_states(lower, data, data.states_2)
        first_out.set_title("first out", loc="left", **annotation)
        for y, label in zip(
            range(len(data.drawn) - 1, -1, -1), data.drawn, strict=True
        ):
            hour = data.first_out_h[label]
            first_out.text(
                0.0,
                y,
                f"{hour:.1f} h" if hour is not None else "never",
                ha="left",
                va="center",
                **annotation,
            )
    upper.tick_params(axis="x", labelbottom=False)
    lower.set_xlabel("\n".join(["Scan clock (h)", *([listed] if listed else [])]))
    for ax in (upper, lower):
        theme.apply_common_style(ax)
    fig.suptitle(title)
    caption = (
        f"Each read's state at each rung of {data.dataset_id}, the four states above and "
        f"the two the carve used below; a rung is drawn when it is in spec at least "
        f"{OCCUPANCY_DRAWN:.0%} of the observed time, the rest are named as culled. "
        "Hatched spans are read gaps."
    )
    if data.flag:
        caption += f" This record is {data.flag}."
    fig.supxlabel(caption, wrap=True, **theme.scaled_text(theme.CAPTION, style))


def timeline_size(data: TimelineBlockData) -> tuple[float, float]:
    """The timeline block's size in inches, its own figure's and its page row's."""
    rows = _TIMELINE_ROW_H * max(len(data.drawn), 1) + _TIMELINE_FRAME_H
    return _TIMELINE_W, _TIMELINE_KEY_H + 2 * rows + 1.2


class TimelineBlockPlot(BasePlot):
    """The timeline block of one record: read state by rung, four states and two."""

    def build_matplotlib(self, result: object, style: str = "default") -> plt.Figure:
        if not isinstance(result, TimelineBlockData):
            raise TypeError(
                f"{type(self).__name__} draws a TimelineBlockData; got "
                f"{type(result).__name__}"
            )
        data = result
        with theme.style_context(style):
            fig = plt.figure(
                figsize=timeline_size(data), layout="constrained", facecolor="white"
            )
            draw_timelines(fig, data, style)
        return fig


class SummariesTablePlot(BasePlot):
    """The curve summaries of a record-figure job, for the `table` target only."""

    def build_table(self, result: object) -> pd.DataFrame:
        if not isinstance(result, pd.DataFrame):
            raise TypeError(
                f"{type(self).__name__} writes a DataFrame; got {type(result).__name__}"
            )
        return result


# Width and height in inches of one panel of the resolvability block.
_PANEL_W = 3.4
_PANEL_H = 2.8
# The kinds of down-crossing, stacked in this order, and the k each pair of bars reads.
CROSSING_KINDS = ("clean", "merged", "unbracketed")
# The id each crossing bar carries, so a test can read the bars back.
CROSSING_GID = "crossing_bar"
# The dash of each count-ratio limit, so the two guides differ in more than their label.
LIMIT_DASHES = {LIMIT_INDEPENDENT: "--", LIMIT_SMOOTH: "-."}
# The id the marker of each count ratio of exactly 0 carries, so a test can find it.
ZERO_RATIO_GID = "zero_ratio"


def _status_text(ax: plt.Axes, text: str, style: str) -> None:
    ax.text(
        0.5,
        0.5,
        text,
        transform=ax.transAxes,
        ha="center",
        va="center",
        wrap=True,
        **theme.scaled_text(theme.ANNOTATION, style),
    )


def _spectrum_marks(data: ResolvabilityBlockData) -> tuple[list[str], list[str]]:
    """The marks of the spectrum panel the block holds, by name, and a note with the status
    of each one it does not."""
    drawn: list[str] = []
    notes: list[str] = []
    if data.band_edge_hz is not None:
        drawn.append("band edge")
    else:
        notes.append(f"band edge ({data.band_edge_status})")
    if data.edge_stretch_hz is not None and data.edge_level_per_hz is not None:
        drawn.append("edge level")
    else:
        notes.append("edge level (no band edge to place it on)")
    if data.noise_floor_per_hz is not None:
        drawn.append("noise floor")
    else:
        notes.append(f"noise floor ({data.noise_floor_status})")
    if data.top_decade_hz is not None:
        drawn.append("the top decade G10 reads")
    else:
        notes.append(f"the top decade G10 reads ({data.g10_status})")
    return drawn, notes


def _draw_spectrum(ax: plt.Axes, data: ResolvabilityBlockData, style: str) -> None:
    if not data.frequency_hz.size:
        ax.set_axis_off()
        _status_text(ax, f"spectrum {data.spectrum_status}", style)
        return
    f, s = data.frequency_hz, data.spectrum_per_hz
    if data.top_decade_hz is not None:
        ax.axvspan(
            *data.top_decade_hz, label="top decade (G10)", **theme.TOP_DECADE_SHADE
        )
    ax.plot(f, s, label=f"S of {data.value_key}", **theme.SPECTRUM_LINE)
    if data.band_edge_hz is not None:
        ax.axvline(data.band_edge_hz, label="band edge", **theme.BAND_EDGE_LINE)
    if data.edge_stretch_hz is not None and data.edge_level_per_hz is not None:
        ax.plot(
            data.edge_stretch_hz,
            (data.edge_level_per_hz, data.edge_level_per_hz),
            label="edge level",
            **theme.EDGE_LEVEL_LINE,
        )
    if data.noise_floor_per_hz is not None:
        ax.axhline(
            data.noise_floor_per_hz, label="noise floor", **theme.NOISE_FLOOR_LINE
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel(f"S of {data.value_key} (per Hz)")
    ax.set_title("Aliasing: the spectrum and its band")
    ax.legend(frameon=False, **theme.scaled_text(theme.LEGEND_TEXT, style))
    theme.apply_common_style(ax)


def _zero_ratios(data: ResolvabilityBlockData) -> list[tuple[str, list[int]]]:
    """Each rung with a count ratio of exactly 0, and the decimations it is 0 at."""
    zeros = []
    for label in data.ladder:
        ms = [
            m
            for m, r in zip(
                data.decimation_ladder, data.count_ratio[label], strict=False
            )
            if r == 0.0
        ]
        if ms:
            zeros.append((label, ms))
    return zeros


def _draw_counts(ax: plt.Axes, data: ResolvabilityBlockData, style: str) -> None:
    m = np.asarray(data.decimation_ladder, dtype=float)
    drawn = [label for label in data.ladder if data.count_ratio[label]]
    if not drawn:
        ax.set_axis_off()
        _status_text(ax, "count ratio not computed at any rung", style)
        return
    for name, values in data.count_limits.items():
        ax.plot(
            m,
            values,
            label=f"limit: {name}",
            **{**theme.REFERENCE_LINE, "linestyle": LIMIT_DASHES[name]},
        )
    n = len(data.ladder)
    for i, label in enumerate(data.ladder):
        ratio = data.count_ratio[label]
        if not ratio:
            continue
        kept = [
            (mm, r) for mm, r in zip(m, ratio, strict=True) if r is not None and r > 0
        ]
        ax.plot(
            [mm for mm, _ in kept],
            [r for _, r in kept],
            marker="o",
            color=theme.threshold_color(i, n),
            label=label,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    zeros = _zero_ratios(data)
    if zeros:
        # Held where the positive values put it, so the markers sit on the floor.
        floor, top = ax.get_ylim()
        ax.set_ylim(floor, top)
        for label, ms in zeros:
            ax.plot(
                ms,
                [floor] * len(ms),
                marker="v",
                linestyle="none",
                color=theme.threshold_color(data.ladder.index(label), n),
                clip_on=False,
                gid=ZERO_RATIO_GID,
            )
    ax.set_xlabel("Decimation m")
    ax.set_ylabel("Count ratio per kept pair")
    ax.set_title("Detection: crossings as reads thin")
    ax.legend(frameon=False, **theme.scaled_text(theme.LEGEND_TEXT, style))
    theme.apply_common_style(ax)


def _draw_crossings(ax: plt.Axes, data: ResolvabilityBlockData, style: str) -> None:
    if not any(data.crossing_counts[label] for label in data.ladder):
        ax.set_axis_off()
        _status_text(ax, "crossing counts not computed at any rung", style)
        return
    width = 0.38
    names: list[str] = []
    for i, label in enumerate(data.ladder):
        counts = data.crossing_counts[label]
        if not counts:
            ax.text(
                i,
                0.0,
                data.crossing_status[label].replace("not computed: ", ""),
                rotation=90,
                ha="center",
                va="bottom",
                **theme.scaled_text(theme.ANNOTATION, style),
            )
        for j, count in enumerate(counts):
            x = i + (j - (len(counts) - 1) / 2.0) * width
            bottom = 0
            for kind in CROSSING_KINDS:
                value = int(getattr(count, kind))
                ax.bar(
                    x,
                    value,
                    width=width * 0.9,
                    bottom=bottom,
                    color=theme.CROSSING_COLORS[kind],
                    label=kind if i == 0 and j == 0 else None,
                    gid=CROSSING_GID,
                )
                bottom += value
        names.append(label)
    ax.set_xticks(range(len(data.ladder)))
    ax.set_xticklabels(names, rotation=45, ha="right")
    ax.set_ylabel("Down-crossings")
    levels = ", ".join(f"{k:g}" for k in data.confidence_levels)
    ax.set_title(f"Placement: each rung at k = {levels}")
    ax.legend(frameon=False, **theme.scaled_text(theme.LEGEND_TEXT, style))
    theme.apply_common_style(ax)


def _listed(items: list[str]) -> str:
    return items[0] if len(items) == 1 else f"{', '.join(items[:-1])} and {items[-1]}"


def _resolvability_caption(data: ResolvabilityBlockData) -> str:
    drawn, notes = _spectrum_marks(data) if data.frequency_hz.size else ([], [])
    marks = f" with its {_listed(drawn)}" if drawn else ""
    parts = [
        f"Event resolvability of {data.dataset_id} ({data.metric_key}): the spectrum of "
        f"{data.value_key}{marks}; per rung of {data.threshold_value_key}, the count ratio "
        "across the decimation ladder with its limits, and the down-crossings clean, "
        "merged or unbracketed at each confidence level k, left to right in increasing k.",
        "G10 and its halves are in the record's table.",
    ]
    if notes:
        parts.append(f"Not drawn on the spectrum: {_listed(notes)}.")
    zeros = _zero_ratios(data)
    if zeros:
        where = "; ".join(
            f"{label} at m = {', '.join(str(m) for m in ms)}" for label, ms in zeros
        )
        parts.append(
            "A count ratio of 0, no down-crossing between kept reads, sits on the axis "
            f"floor as a down triangle: {where}."
        )
    if data.flag:
        parts.append(f"This record is {data.flag}.")
    return " ".join(parts)


def draw_resolvability(
    fig: Figure | SubFigure, data: ResolvabilityBlockData, style: str
) -> None:
    """The resolvability block of one record on `fig` (module docstring)."""
    title = f"Event resolvability, {data.dataset_id}"
    if data.failure:
        ax = fig.add_subplot(1, 1, 1)
        ax.set_axis_off()
        _status_text(ax, f"not computed: {data.failure}", style)
        fig.suptitle(title)
        return
    axes = fig.subplots(1, 3)
    _draw_spectrum(axes[0], data, style)
    _draw_counts(axes[1], data, style)
    _draw_crossings(axes[2], data, style)
    fig.suptitle(f"{title} ({data.metric_key})")
    fig.supxlabel(
        _resolvability_caption(data),
        wrap=True,
        **theme.scaled_text(theme.CAPTION, style),
    )


class ResolvabilityBlockPlot(BasePlot):
    """The resolvability block of one record: aliasing, detection and placement."""

    def build_matplotlib(self, result: object, style: str = "default") -> plt.Figure:
        if not isinstance(result, ResolvabilityBlockData):
            raise TypeError(
                f"{type(self).__name__} draws a ResolvabilityBlockData; got "
                f"{type(result).__name__}"
            )
        with theme.style_context(style):
            fig = plt.figure(
                figsize=(3 * _PANEL_W, _PANEL_H + 1.2),
                layout="constrained",
                facecolor="white",
            )
            draw_resolvability(fig, result, style)
        return fig


class ResolvabilityTablePlot(BasePlot):
    """One of the job's two stacked resolvability tables, for the `table` target only."""

    def build_table(self, result: object) -> pd.DataFrame:
        if not isinstance(result, pd.DataFrame):
            raise TypeError(
                f"{type(self).__name__} writes a DataFrame; got {type(result).__name__}"
            )
        return result


def resolvability_size(data: ResolvabilityBlockData) -> tuple[float, float]:
    """The resolvability block's size in inches, as `ResolvabilityBlockPlot` draws it."""
    return 3 * _PANEL_W, _PANEL_H + 1.2


# The page's rows: each block's artifact type, the function that draws it, and its size.
_PAGE_ROWS: dict[
    type, tuple[Callable[..., None], Callable[..., tuple[float, float]]]
] = {
    SurvivalBlockData: (draw_survival, survival_size),
    HazardBlockData: (draw_hazard, hazard_size),
    FractionBlockData: (draw_fraction, fraction_size),
    ResolvabilityBlockData: (draw_resolvability, resolvability_size),
    TimelineBlockData: (draw_timelines, timeline_size),
}
# Inches the page adds above its rows for its own title.
_PAGE_TITLE_H = 0.5


def _page_row(
    block: object,
) -> tuple[Callable[..., None], Callable[..., tuple[float, float]]]:
    if type(block) not in _PAGE_ROWS:
        raise TypeError(f"no page row draws a {type(block).__name__}")
    return _PAGE_ROWS[type(block)]


def page_size(data: RecordFigureData) -> tuple[float, float]:
    """The page's size in inches: the widest block's width, every block's height summed."""
    sizes = [_page_row(block)[1](block) for block in data.blocks]
    return max(w for w, _h in sizes), sum(h for _w, h in sizes) + _PAGE_TITLE_H


def draw_record_figure(fig: Figure, data: RecordFigureData, style: str) -> None:
    """One record's page on `fig` (module docstring)."""
    heights = [_page_row(block)[1](block)[1] for block in data.blocks]
    rows = fig.subfigures(len(data.blocks), 1, height_ratios=heights)
    for row, block in zip(np.atleast_1d(rows), data.blocks, strict=True):
        _page_row(block)[0](row, block, style)
    fig.suptitle(f"Record figure, {data.dataset_id}")


class RecordFigurePlot(BasePlot):
    """One record's page: every block, one row each, in the page's order."""

    def build_matplotlib(self, result: object, style: str = "default") -> plt.Figure:
        if not isinstance(result, RecordFigureData):
            raise TypeError(
                f"{type(self).__name__} draws a RecordFigureData; got "
                f"{type(result).__name__}"
            )
        with theme.style_context(style):
            fig = plt.figure(
                figsize=page_size(result), layout="constrained", facecolor="white"
            )
            draw_record_figure(fig, result, style)
        return fig
