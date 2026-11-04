"""One figure per instrument: is the independence assumption safe across every dataset?

One instrument, one figure, because a reader asks about one instrument at a time and a
single composite of all of them would make each grid unreadably small. `SURVEY_PLOTS` holds
one class per entry of `independence_survey.SURVEY_KEYS`, in that order, so the job can
wire them in a loop, and `GROUP_SURVEY_PLOTS` one per group of `check_ledger.GROUP_KEYS`,
whose caption says the group reads no bench cell.

Each figure is one instrument: its grids are the clocks down and the bags across, each grid
datasets down and threshold ladder across, each headed on a strip of its bag's colour
(`theme.bag_header`).

PURE RENDERER. Every value drawn is a field of `IndependenceSurveyData`. Colour comes from
`plots.theme.verdict_color` - the SAME palette the per-dataset ledger panel uses, so a
reader moving between the two figures does not have to relearn what green means.

Follows `docs/FIGURE_STANDARD.md`: the title names the instrument and does not argue, the
caption states what the instrument rejects and what is excluded, spaced hyphens throughout.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

from quebra.analyzers.independence_survey import (
    VERDICT_ORDER,
    IndependenceSurveyData,
    _fmt_p,
)
from quebra.analyzers.check_ledger import GROUP_CALIBRATION_NOTE, GROUP_KEYS
from quebra.analyzers.independence_survey import SURVEY_KEYS
from quebra.plots import theme
from quebra.plots.base import BasePlot


def _check(result: object) -> IndependenceSurveyData:
    if not isinstance(result, IndependenceSurveyData):
        raise TypeError(
            f"{type(result).__name__} is not IndependenceSurveyData; this figure draws "
            "only the typed artifact its builder returns."
        )
    return result


# A group grid's members are those that answer in each cell, which can differ by cell.
GROUP_MEMBERS_NOTE = (
    "A cell calibrates the members that answer in it, so they can differ from cell to "
    "cell; each ledger row names them."
)


class _InstrumentSurveyPlot(BasePlot):
    """Base: subclasses set `KEY` to one entry of `battery.ROW_KEYS`, or `GROUP` and `KEY`
    to one entry of `check_ledger.GROUP_KEYS`."""

    KEY: tuple[str, str, str] = ()
    # A group row reads no bench cell, so its caption says what its size rests on instead.
    GROUP: bool = False

    def build_matplotlib(self, result: object, style: str = "default") -> plt.Figure:
        pd_ = _check(result)
        grids = [g for g in pd_.grids if g.key == self.KEY]
        if not grids:
            raise ValueError(
                f"no grid for {self.KEY} - the survey artifact does not carry this "
                "instrument, so there is nothing to draw and a blank figure would lie"
            )

        sides = [side for side in pd_.sides if any(g.side == side for g in grids)]
        clocks = [clock for clock in pd_.clocks if any(g.clock == clock for g in grids)]
        n_rows = len(pd_.datasets)
        n_cols = len(pd_.thresholds)
        with theme.style_context(style):
            fig, axes = plt.subplots(
                len(clocks),
                len(sides),
                figsize=(
                    max(7.0, 0.62 * n_cols + 4.2) * len(sides),
                    max(3.0, 0.30 * n_rows) * len(clocks) + 2.2,
                ),
                squeeze=False,
                facecolor="white",
            )
            for row, clock in enumerate(clocks):
                for col, side in enumerate(sides):
                    grid = pd_.grid(self.KEY, side, clock)
                    if grid is None:
                        axes[row][col].axis("off")
                        continue
                    self._draw_grid(axes[row][col], grid, pd_)

            first = grids[0]
            fig.suptitle(first.label)
            legend = [
                Patch(facecolor=theme.verdict_color(v), label=v) for v in VERDICT_ORDER
            ]
            fig.legend(
                handles=legend,
                loc="lower center",
                ncol=len(VERDICT_ORDER),
                frameon=False,
                **theme.LEGEND_TEXT,
            )
            decided = [
                f"{theme.BAG_LABELS[g.side]}, {g.clock} clock: {g.counts['fail']} of "
                f"{g.counts['pass'] + g.counts['fail']} decided cells reject, "
                f"{g.counts['not computed']} of {g.verdicts.size} not computed"
                for g in grids
            ]
            pale = (
                "too few events"
                if self.GROUP
                else "too few events, a size the bench found wrong, or no bench cell"
            )
            calibration = (
                GROUP_CALIBRATION_NOTE[0].upper()
                + GROUP_CALIBRATION_NOTE[1:]
                + " "
                + GROUP_MEMBERS_NOTE
                if self.GROUP
                else "Every instrument here is scored against a measured bench size."
            )
            fig.text(
                0.5,
                0.035,
                f"This instrument {first.null_statement}. A red cell is a statement about "
                f"the DATA, not a failed test. A grid title's colour names its bag, not a "
                f"verdict. Pale cells carry no information: an "
                f"inconclusive non-rejection ({pale}) is not evidence of independence. "
                f"{'; '.join(decided)}. {calibration}",
                ha="center",
                wrap=True,
                **theme.CAPTION,
            )
            fig.tight_layout(rect=(0, 0.10, 1, 0.96))
        return fig

    def _draw_grid(self, ax, grid, pd_) -> None:
        verdicts = grid.verdicts.to_numpy(dtype=object)
        p_values = grid.p_values.to_numpy(dtype=object)
        n_rows, n_cols = verdicts.shape

        for i in range(n_rows):
            for j in range(n_cols):
                verdict = verdicts[i, j]
                colour = (
                    theme.verdict_color(str(verdict))
                    if isinstance(verdict, str) and verdict in VERDICT_ORDER
                    # A cell with no row at all is not a verdict. It reuses the
                    # `not computed` tone because that is what it means: nothing ran.
                    else theme.verdict_color("not computed")
                )
                ax.add_patch(
                    plt.Rectangle(
                        (j, i), 1, 1, facecolor=colour, edgecolor="white", linewidth=0.8
                    )
                )
                text = _fmt_p(p_values[i, j])
                if text:
                    ax.text(
                        j + 0.5,
                        i + 0.5,
                        text,
                        ha="center",
                        va="center",
                        color="white"
                        if verdict in ("pass", "fail")
                        else theme.ON_FILL_TEXT["color"],
                        **{"fontsize": theme.ON_FILL_TEXT["fontsize"] - 1},
                    )

        ax.set_xlim(0, n_cols)
        ax.set_ylim(0, n_rows)
        ax.set_xticks(np.arange(n_cols) + 0.5)
        ax.set_xticklabels(
            list(grid.verdicts.columns), rotation=45, ha="right", **theme.LABEL_TEXT
        )
        ax.set_yticks(np.arange(n_rows) + 0.5)
        ax.set_yticklabels(list(grid.verdicts.index), **theme.LABEL_TEXT)
        ax.invert_yaxis()
        ax.set_xlabel("Threshold")
        ax.set_ylabel("Dataset")
        ax.set_title(
            f"{theme.BAG_LABELS[grid.side]} - {grid.clock} clock",
            loc="left",
            **theme.LABEL_TEXT,
            **theme.bag_header(grid.side),
        )
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(length=0)


def _make(
    key: tuple[str, str, str], name: str, group: bool = False
) -> type[_InstrumentSurveyPlot]:
    return type(name, (_InstrumentSurveyPlot,), {"KEY": key, "GROUP": group})


def _class_name(key: tuple[str, str, str]) -> str:
    """`("c5_rank_autocorr", "permutation", "studentized")` -> `C5RankAutocorrPermutationStudentizedSurveyPlot`."""
    parts = [p for part in key for p in part.split("_") if p]
    return "".join(p.capitalize() for p in parts) + "SurveyPlot"


# GENERATED from SURVEY_KEYS, one class per entry, rather than hand-listed by index. A
# hand-written list goes stale by omission: an entry appended to the tuple draws nothing and
# an entry inserted shifts every index after it. Deriving the list means adding a check adds
# its figure here, and `test_there_is_exactly_one_figure_per_surveyed_instrument` pins it.
SURVEY_PLOTS: list[type[_InstrumentSurveyPlot]] = [
    _make(key, _class_name(key)) for key in SURVEY_KEYS
]

# One per group of `check_ledger.GROUP_KEYS`, generated the same way.
GROUP_SURVEY_PLOTS: list[type[_InstrumentSurveyPlot]] = [
    _make(key, _class_name(key), group=True) for key in GROUP_KEYS
]

# Bound at module level too, so an importer can name one directly.
for _cls in (*SURVEY_PLOTS, *GROUP_SURVEY_PLOTS):
    globals()[_cls.__name__] = _cls


class IndependenceSurveyOverviewPlot(BasePlot):
    """EVERY single check in one image: one row per check, columns grouped by bag, one
    column per clock within each bag, each column headed on a strip of its bag's colour. The
    two group grids are not drawn here; each has its own figure (`GROUP_SURVEY_PLOTS`).

    The per-instrument figures above are for reading one instrument closely. This one is
    for the question they cannot answer between them - does the picture change when you
    swap the instrument, or only when you swap the BAG or the CLOCK? Laid out so that
    comparison is a left-right glance rather than a memory exercise across as many files.

    TRANSPOSED relative to the per-instrument figures: thresholds on y (the ladder's rungs) and
    datasets on x (28 on the device-wide survey). Keeping datasets on y would make every
    panel that many rows tall and the image unusable; the ladder is the shorter axis and belongs on the
    short side. p-values are dropped here for the same reason - at this density they would
    be unreadable, and the per-instrument figure is where a number is read.
    """

    def build_matplotlib(self, result: object, style: str = "default") -> plt.Figure:
        pd_ = _check(result)
        columns = [(side, clock) for side in pd_.sides for clock in pd_.clocks]
        keys = [k for k in SURVEY_KEYS if any(g.key == k for g in pd_.grids)]
        if not keys:
            raise ValueError("the survey artifact carries no instrument grids")

        n_datasets = len(pd_.datasets)
        with theme.style_context(style):
            fig, axes = plt.subplots(
                len(keys),
                len(columns),
                figsize=(
                    0.30 * n_datasets * len(columns) + 3.0,
                    1.55 * len(keys) + 2.6,
                ),
                squeeze=False,
                facecolor="white",
            )
            for row, key in enumerate(keys):
                for col, (side, clock) in enumerate(columns):
                    ax = axes[row][col]
                    grid = pd_.grid(key, side, clock)
                    if grid is None:
                        ax.axis("off")
                        continue
                    self._draw_transposed(ax, grid)
                    if row == 0:
                        ax.set_title(
                            f"{theme.BAG_LABELS[side]} - {clock} clock",
                            **theme.LABEL_TEXT,
                            **theme.bag_header(side),
                        )
                    if col == 0:
                        ax.set_ylabel(
                            grid.label.replace(" (", "\n("),
                            rotation=0,
                            ha="right",
                            va="center",
                            labelpad=8,
                            **theme.LABEL_TEXT,
                        )
                    # Dataset names only on the bottom row: one copy of a
                    # long dataset axis per row is noise, not information.
                    if row == len(keys) - 1:
                        ax.set_xticks(np.arange(n_datasets) + 0.5)
                        ax.set_xticklabels(
                            list(grid.verdicts.index), rotation=90, **theme.LEGEND_TEXT
                        )
                    else:
                        ax.set_xticks([])

            legend = [
                Patch(facecolor=theme.verdict_color(v), label=v) for v in VERDICT_ORDER
            ]
            fig.legend(
                handles=legend,
                loc="lower center",
                ncol=len(VERDICT_ORDER),
                frameon=False,
                **theme.LEGEND_TEXT,
            )
            fig.suptitle(
                "Window independence across every dataset - every single check"
            )
            fig.text(
                0.5,
                0.028,
                "Rows are the single checks (the two group grids have their own "
                "figures); columns are grouped by bag, the windows on each "
                "side of the threshold, each group headed in its colour and holding the "
                "two clocks; a header's colour names its bag, not a verdict. Within each "
                "grid the y axis is the threshold ladder and the x "
                "axis is the dataset. Red is a rejection - a statement about the DATA, not "
                "a failed test. Grey carries no information. The two clocks ask different "
                "questions of the same bag: the side clock advances only while the record "
                "is on the bag's side and its events are window DEATHS, the crossings that "
                "leave the side; calendar time is wall clock and its events are window "
                "BIRTHS - so they can disagree without either being wrong. Read a p-value "
                "off the per-instrument figures, not this one.",
                ha="center",
                wrap=True,
                **theme.CAPTION,
            )
            fig.tight_layout(rect=(0, 0.075, 1, 0.965))
        return fig

    def _draw_transposed(self, ax, grid) -> None:
        """Thresholds down, datasets across - the transpose of the per-instrument view."""
        verdicts = grid.verdicts.to_numpy(dtype=object).T
        n_thresholds, n_datasets = verdicts.shape
        for i in range(n_thresholds):
            for j in range(n_datasets):
                verdict = verdicts[i, j]
                colour = (
                    theme.verdict_color(str(verdict))
                    if isinstance(verdict, str) and verdict in VERDICT_ORDER
                    else theme.verdict_color("not computed")
                )
                ax.add_patch(
                    plt.Rectangle(
                        (j, i), 1, 1, facecolor=colour, edgecolor="white", linewidth=0.4
                    )
                )
        ax.set_xlim(0, n_datasets)
        ax.set_ylim(0, n_thresholds)
        ax.set_yticks(np.arange(n_thresholds) + 0.5)
        ax.set_yticklabels(list(grid.verdicts.columns), **theme.LEGEND_TEXT)
        ax.invert_yaxis()
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(length=0)
