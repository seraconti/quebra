"""CheckLedgerPanel: the ledger as a table, over the statistic across the ladder.

A pure renderer. Every verdict arrives already decided by `analyzers/check_ledger.py`; this
file chooses colours and positions and nothing else.

Two parts, one figure, because they answer two halves of one question. The table says what
each check concluded and whether that conclusion is worth anything. The panels below say
how the evidence changes down the threshold ladder - and because the event count is drawn
on a second axis, a reader can see power collapsing at the loose thresholds rather than
having to infer it.

Per `docs/FIGURE_STANDARD.md`: the title names the object and does not argue, the locked
vocabulary is used (window, read, check, band), units are in the axis labels, hyphens are
spaced, and the count of thresholds excluded from the table is stated IN the panel.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

from quebra.analyzers.check_ledger import CheckLedger, VERDICTS
from quebra.analyzers.checks.result import CLOCK_CALENDAR, CLOCK_IN_SPEC
from quebra.panels._check_ledger_render import (
    CHECK_SHORT,
    format_p,
    ledger_grid,
    statistic_series,
)
from quebra.plots.base import BasePlot
from quebra.plots.theme import (
    REFERENCE_LINE,
    UNOBSERVED_HATCH,
    ordered_color,
    style_context,
    verdict_color,
)

# Reuses the record-state hatch: both mean "this cell is not what it appears to be".
MISCALIBRATED_HATCH = UNOBSERVED_HATCH

__all__ = ["CheckLedgerPanel"]


class CheckLedgerPanel(BasePlot):
    def build_matplotlib(
        self, result: CheckLedger, style: str = "default"
    ) -> plt.Figure:
        if not isinstance(result, CheckLedger):
            raise TypeError("CheckLedgerPanel expects a CheckLedger artifact")
        with style_context(style):
            clocks = [
                clock
                for clock in (CLOCK_IN_SPEC, CLOCK_CALENDAR)
                if len(result.rows[result.rows["clock"] == clock])
            ]
            checks = [
                check for check in CHECK_SHORT if check in set(result.rows["check_id"])
            ]
            fig = plt.figure(figsize=(16, 4.2 + 3.4 * len(clocks)))
            grid = fig.add_gridspec(
                len(clocks) + 1,
                1,
                height_ratios=[3.4] * len(clocks) + [4.2],
                hspace=0.45,
            )

            for index, clock in enumerate(clocks):
                self._draw_table(fig.add_subplot(grid[index, 0]), result, clock)
            self._draw_ladder(fig.add_subplot(grid[len(clocks), 0]), result, checks)

            fig.suptitle(
                f"Check ledger - {result.dataset_id}, alpha {result.alpha:g}, "
                f"B = {result.n_permutations}, seed {result.seed}"
            )
            fig.tight_layout(rect=(0, 0, 1, 0.97))
            return fig

    def build_plotly(self, result: CheckLedger):
        raise NotImplementedError("CheckLedgerPanel has no plotly backend")

    # --- private drawing methods (read precomputed fields; no data arithmetic) ---

    def _draw_table(self, ax, result: CheckLedger, clock: str) -> None:
        thresholds, columns, grid = ledger_grid(result.rows, clock)
        ax.set_xlim(0, max(len(columns), 1))
        ax.set_ylim(0, max(len(thresholds), 1))
        ax.invert_yaxis()
        ax.set_xticks(np.arange(len(columns)) + 0.5)
        ax.set_xticklabels(columns, rotation=30, ha="left")
        ax.xaxis.set_ticks_position("top")
        ax.set_yticks(np.arange(len(thresholds)) + 0.5)
        ax.set_yticklabels(thresholds)
        ax.set_ylabel("Threshold")
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(length=0)
        ax.grid(False)

        for row_index, row in enumerate(grid):
            for column_index, cell in enumerate(row):
                if cell is None:
                    continue
                verdict, p_value, bench_ok = cell
                ax.add_patch(
                    plt.Rectangle(
                        (column_index, row_index),
                        1,
                        1,
                        facecolor=verdict_color(verdict),
                        edgecolor="white",
                        linewidth=1.0,
                        # A rejection from a check the bench found miscalibrated at this
                        # event count is still a rejection, but it is not the same
                        # evidence as one from a calibrated check. The hatch says so at a
                        # glance; burying it in a notes column would not.
                        #
                        # BOTH non-calibrated states hatch, not just one. `None` means the
                        # bench has no cell for this check at all, which is a weaker
                        # footing than `False`, not a stronger one. Hatching only `False`
                        # drew an unbenched cell identically to a calibrated one.
                        hatch=None if bench_ok is True else MISCALIBRATED_HATCH,
                    )
                )
                ax.text(
                    column_index + 0.5,
                    row_index + 0.5,
                    format_p(p_value),
                    ha="center",
                    va="center",
                    color="white" if verdict in ("pass", "fail") else "black",
                )
        ax.set_title(f"{clock} clock", loc="left")
        if clock == CLOCK_IN_SPEC:
            ax.legend(
                handles=[
                    Patch(facecolor=verdict_color(v), edgecolor="white", label=v)
                    for v in VERDICTS
                ]
                + [
                    Patch(
                        facecolor="white",
                        edgecolor="black",
                        hatch=MISCALIBRATED_HATCH,
                        label="bench: miscalibrated at this n, or no bench cell",
                    )
                ],
                loc="upper left",
                bbox_to_anchor=(1.01, 1.0),
                frameon=False,
            )

    def _draw_ladder(self, ax, result: CheckLedger, checks: list[str]) -> None:
        labels = [label for label, _v, _b in result.thresholds]
        positions = {label: index for index, label in enumerate(labels)}
        for index, check in enumerate(checks):
            series = statistic_series(result.rows, CLOCK_IN_SPEC, check, labels)
            if all(np.isnan(series)):
                continue
            ax.plot(
                range(len(labels)),
                series,
                marker="o",
                markersize=4.0,
                color=ordered_color(index, max(len(checks), 1)),
                label=CHECK_SHORT.get(check, check),
            )
        ax.axhline(result.alpha, **REFERENCE_LINE)
        ax.set_yscale("log")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels)
        ax.set_xlabel("Threshold")
        ax.set_ylabel("p-value")
        ax.set_title("Evidence across the threshold ladder - in-spec clock", loc="left")
        ax.legend(loc="lower left", frameon=False, ncol=2)

        # Event count on a second axis: the reader sees power collapsing down the ladder
        # instead of having to infer it from the p-values.
        counts = (
            result.rows[result.rows["clock"] == CLOCK_IN_SPEC]
            .groupby("threshold_label")["n_events"]
            .max()
        )
        twin = ax.twinx()
        twin.step(
            [positions[t] for t in counts.index if t in positions],
            [counts[t] for t in counts.index if t in positions],
            where="mid",
            color=REFERENCE_LINE["color"],
            linewidth=REFERENCE_LINE["linewidth"],
            alpha=REFERENCE_LINE["alpha"],
        )
        twin.set_ylabel("Events (windows)")
        twin.set_yscale("log")
        twin.grid(False)

        # FIGURE_STANDARD: say what was dropped, in the panel.
        empty = [t for t in labels if t not in set(counts[counts > 0].index)]
        if empty:
            ax.text(
                0.99,
                0.03,
                f"no usable windows: {', '.join(empty)}",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
            )
