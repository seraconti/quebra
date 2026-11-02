"""One verdict grid per displayed check, as a separate figure from the survival curve.

Separate for a structural reason, not an aesthetic one: `job.figure` binds one plot class to
one node, and `plots/km_survival_plot.py` carries no in-panel text on the author's
instruction. Nothing here goes into the survival axes.

PURE RENDERER. Every value drawn is a field of `CheckOutcome`. Colour comes from
`theme.verdict_color`, the same palette the ledger panel and the independence survey use, so
a reader moving between those figures does not relearn what green means.
"""

from __future__ import annotations

import matplotlib.pyplot as plt

from matplotlib.patches import Patch

from quebra.analyzers.check_ledger import VERDICT_FAIL, VERDICT_PASS
from quebra.analyzers.check_outcome import CheckOutcome, OutcomeGrid
from quebra.analyzers.independence_survey import VERDICT_ORDER, _fmt_p
from quebra.plots import theme
from quebra.plots.base import BasePlot

# Wide enough for a dataset label, tall enough that a single-row grid is still readable.
_CELL_W = 1.5
_CELL_H = 0.45


class CheckOutcomePlot(BasePlot):
    """One row of grids, one per displayed check, sharing a dataset axis."""

    def build_matplotlib(self, result: object, style: str = "default") -> plt.Figure:
        if not isinstance(result, CheckOutcome):
            raise TypeError(
                f"{type(self).__name__} draws a CheckOutcome; got "
                f"{type(result).__name__}"
            )

        grids = result.grids
        if not grids:
            # An empty display-set is a legitimate configuration, so this renders rather
            # than raising: the checks ran, the job asked for none of them to be shown.
            fig, ax = plt.subplots(figsize=(6.0, 1.4))
            ax.text(
                0.5,
                0.5,
                "no check selected for display\n"
                f"({len(result.run_set)} ran and are in the artifact)",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.axis("off")
            return fig

        n_rows = max(len(g.datasets) for g in grids)
        n_cols = max(len(g.clocks) for g in grids)
        width = max(4.0, _CELL_W * n_cols * len(grids) + 2.2)
        height = max(2.0, _CELL_H * n_rows + 1.8)

        fig, axes = plt.subplots(1, len(grids), figsize=(width, height), squeeze=False)
        for ax, grid in zip(axes[0], grids, strict=True):
            self._draw_grid(ax, grid, show_dataset_labels=ax is axes[0][0])

        # Colour is the only encoding of a verdict, so the key is not optional.
        shown = [v for v in VERDICT_ORDER if any(v in g.counts for g in grids)]
        if shown:
            fig.legend(
                handles=[
                    Patch(facecolor=theme.verdict_color(v), label=v) for v in shown
                ],
                loc="lower center",
                ncol=len(shown),
                frameon=False,
                **theme.LEGEND_TEXT,
            )
        # The bag is named in the title: a column header reads only "side".
        bag = f", {theme.BAG_LABELS[result.side]}" if result.side else ""
        fig.suptitle(f"Independence checks at {result.threshold_label}{bag}")
        fig.text(
            0.5,
            0.02,
            self._caption(result),
            ha="center",
            wrap=True,
            **theme.CAPTION,
        )
        fig.tight_layout(rect=(0, 0.10, 1, 0.94))
        return fig

    @staticmethod
    def _caption(result: CheckOutcome) -> str:
        """States what was asked for and what did not answer. Description, not argument."""
        parts = [
            f"{len(result.display_set)} of {len(result.run_set)} checks shown; "
            f"all of them are in the artifact.",
            "A cell is a verdict, not a p-value: 'did not reject' is not 'passed', and the "
            "desaturated tones carry no information.",
            "The band this accompanies is drawn in every case - these checks annotate it "
            "and never suppress it.",
        ]
        if result.asked_but_unanswered:
            names = ", ".join(
                " ".join(p for p in key if p) for key in result.asked_but_unanswered
            )
            parts.append(f"Asked for but produced no answer anywhere: {names}.")
        return " ".join(parts)

    @staticmethod
    def _draw_grid(ax, grid: OutcomeGrid, *, show_dataset_labels: bool) -> None:
        """The whole drawing, as a function of axes, theme and one grid."""
        n_rows, n_cols = len(grid.datasets), len(grid.clocks)
        for i in range(n_rows):
            for j in range(n_cols):
                verdict = grid.verdicts[i][j]
                ax.add_patch(
                    plt.Rectangle(
                        (j, n_rows - 1 - i),
                        1,
                        1,
                        facecolor=theme.verdict_color(verdict),
                        edgecolor="white",
                        linewidth=0.8,
                    )
                )
                p_value = grid.p_values[i][j]
                if p_value is not None:
                    # `_fmt_p`, not `:.2f`: a permutation p-value floors at 1/(B+1) = 0.001
                    # at B = 999, and `:.2f` prints that as "0.00", a value it cannot take.
                    # White on the two dark fills, as the sibling grid figure does.
                    ax.text(
                        j + 0.5,
                        n_rows - 1 - i + 0.5,
                        _fmt_p(p_value),
                        ha="center",
                        va="center",
                        color="white"
                        if verdict in (VERDICT_PASS, VERDICT_FAIL)
                        else theme.ON_FILL_TEXT["color"],
                        **{"fontsize": theme.ON_FILL_TEXT["fontsize"] - 1},
                    )

        ax.set_xlim(0, n_cols)
        ax.set_ylim(0, n_rows)
        ax.set_xticks([j + 0.5 for j in range(n_cols)])
        ax.set_xticklabels(list(grid.clocks), **theme.LEGEND_TEXT)
        if show_dataset_labels:
            ax.set_yticks([n_rows - 1 - i + 0.5 for i in range(n_rows)])
            ax.set_yticklabels(list(grid.dataset_labels), **theme.LABEL_TEXT)
        else:
            ax.set_yticks([])
        ax.set_title(grid.label, loc="left", **theme.LABEL_TEXT)
        if grid.null_statement:
            ax.set_xlabel(grid.null_statement, wrap=True, **theme.CAPTION)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(length=0)
