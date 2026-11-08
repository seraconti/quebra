"""InfidelityCheckLedgerPanel: the check ledger of infidelity windows, captioned with its bench.

A pure renderer. It draws `CheckLedgerPanel`'s figure for one bag and adds one caption line
the artifact carries: the record's name in the paper's dataset tables, and the statement that
the bench's size cells and ties cutoff are not measured on an |x| carve
(`fidelity_ledger_recipe.InfidelityCheckLedger`).
"""

from __future__ import annotations

import textwrap
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt

from quebra.panels.check_ledger import CheckLedgerPanel
from quebra.plots.theme import CAPTION, style_context

if TYPE_CHECKING:
    from quebra.fidelity_ledger_recipe import InfidelityCheckLedger

__all__ = ["InfidelityCheckLedgerPanel"]

# Figure-height fraction one caption line takes, and the base figure's own bottom margin.
LINE_HEIGHT = 0.016
BASE_BOTTOM = 0.02


class InfidelityCheckLedgerPanel(CheckLedgerPanel):
    def build_matplotlib(
        self, result: InfidelityCheckLedger, style: str = "default"
    ) -> plt.Figure:
        paper_name = getattr(result, "paper_name", "")
        bench_transfer = getattr(result, "bench_transfer", "")
        if not paper_name or not bench_transfer:
            raise TypeError(
                "InfidelityCheckLedgerPanel expects an InfidelityCheckLedger carrying its "
                "paper name and bench statement"
            )
        fig = super().build_matplotlib(result, style)
        # Wrapped to the figure's width at the caption's size, as the compliance survey's
        # caption is: matplotlib's own wrap lets a line run past the edge.
        per_line = int(fig.get_figwidth() * 72 / (0.55 * CAPTION["fontsize"]) * 0.92)
        caption = textwrap.fill(
            f"{paper_name} in the paper's dataset tables. {bench_transfer}", per_line
        )
        n_lines = caption.count("\n") + 1
        with style_context(style):
            fig.text(
                0.01, BASE_BOTTOM + 0.004, caption, ha="left", va="bottom", **CAPTION
            )
            fig.tight_layout(
                rect=(0, BASE_BOTTOM + LINE_HEIGHT * n_lines + 0.006, 1, 0.97)
            )
        return fig
