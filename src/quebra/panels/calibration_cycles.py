"""Panel for the calibration cycles of one status history: survival, check outcome, shares.

A pure renderer of `analyzers/cycles.CyclesPanelData`, which `build_cycles_panel_data`
builds with every number drawn here. What is left to this module is a function of axes and
theme: limits, colours, where text sits, and clipping the MTBC interval at the axis floor.

Three axes:

- Left: the Kaplan-Meier survival of the cycles on a linear minutes axis, its band, the
  censored cycles as ticks, the MTBC as a line with its interval shaded behind it, and tau,
  the horizon the restricted mean runs to. The counts of what was dropped, merged and cut
  sit in this axes, as `docs/FIGURE_STANDARD.md` asks.
- Top right: the check outcome on the cycles, one row per check the run-set asked for, the
  verdict as the cell's colour (`theme.verdict_color`), the bench-transfer note under it,
  and the reason when the battery scored nothing.
- Bottom right: the hours per status, time no row covers drawn white and hatched as the
  other panels draw unobserved time.

The caption carries the check outcome's one line (`CheckAttachment.summary`), so the band and
the MTBC interval are never shown without it.
"""

from __future__ import annotations

import math
import re

import matplotlib.pyplot as plt
import numpy as np

from quebra.analyzers.check_ledger import (
    VERDICT_FAIL,
    VERDICT_INCONCLUSIVE,
    VERDICT_PASS,
)
from quebra.analyzers.cycles import NOT_LISTED, CyclesPanelData
from quebra.analyzers.independence_survey import _fmt_p
from quebra.plots import theme
from quebra.plots.base import BasePlot
from quebra.plots.km_survival_plot import band_vertices, step_arrays

# The x axis runs this far past the longest cycle or tau, so the last step is not drawn
# on the frame.
_X_MARGIN = 1.03
# The id the MTBC interval's shaded span carries, so a test can find the one drawn patch.
MTBC_SPAN_GID = "mtbc_interval"
# Verdicts that answer; any other cell's note says why it did not.
ANSWERED = (VERDICT_PASS, VERDICT_FAIL, VERDICT_INCONCLUSIVE)
# The check table's columns, in its own units (the axes run 0 to 4): the p-value's centre,
# then the verdict cell's left edge and width. The label column starts at 0.
_P_X = 2.4
_CELL_X = 2.75
_CELL_W = 1.25


class CalibrationCyclesPanel(BasePlot):
    """Survival of calibration cycles with the MTBC, the check outcome and status shares."""

    def build_matplotlib(self, result: object, style: str = "default") -> plt.Figure:
        if not isinstance(result, CyclesPanelData):
            raise TypeError(
                f"{type(self).__name__} draws a CyclesPanelData; got "
                f"{type(result).__name__}"
            )
        data = result
        with theme.style_context(style):
            fig = plt.figure(figsize=(13.0, 7.5), facecolor="white")
            grid = fig.add_gridspec(2, 2, width_ratios=[1.6, 1.0], height_ratios=[1, 1])
            ax_km = fig.add_subplot(grid[:, 0])
            ax_checks = fig.add_subplot(grid[0, 1])
            ax_shares = fig.add_subplot(grid[1, 1])
            self._draw_survival(ax_km, data, style)
            self._draw_checks(ax_checks, data, style)
            self._draw_shares(ax_shares, data, style)
            fig.suptitle(f"Calibration cycle survival, {data.dataset_id}")
            fig.text(
                0.5,
                0.015,
                _caption(data),
                ha="center",
                wrap=True,
                **theme.scaled_text(theme.CAPTION, style),
            )
            fig.tight_layout(rect=(0, 0.07, 1, 0.96))
        return fig

    @staticmethod
    def _draw_survival(ax: plt.Axes, data: CyclesPanelData, style: str) -> None:
        curve = data.curve
        ink = theme.POSTER_INK_SURVIVAL
        x, s, lo, hi = step_arrays(curve, 0.0)
        x_hi = max(float(x[-1]), data.tau_min) * _X_MARGIN
        ax.set_xlim(0.0, x_hi)
        ax.set_ylim(-0.02, 1.04)
        ax.fill_between(
            *band_vertices(x, lo, hi),
            color=ink,
            alpha=theme.BAND_STYLE["fill_alpha"],
            linewidth=0.0,
            zorder=2,
            label=f"{curve.conf_level:.0%} band",
        )
        ax.step(
            x,
            s,
            where="post",
            color=ink,
            zorder=3,
            label=f"{data.n_cycles} cycles, {data.n_censored} censored",
        )
        if len(curve.censor_time_min):
            ax.plot(
                curve.censor_time_min,
                curve.censor_survival,
                linestyle="none",
                color=ink,
                zorder=4,
                **theme.CENSOR_TICK,
            )
        point, lower, upper = data.mtbc_min
        label = f"MTBC {point:.1f} min"
        if math.isfinite(lower) and math.isfinite(upper):
            label += f" [{lower:.1f}, {upper:.1f}]"
            # The axis starts at 0, so an interval reaching below it is drawn from 0; the
            # number on the label is kept.
            if lower <= 0.0:
                label += " (lower end <= 0, drawn from 0)"
            ax.axvspan(
                max(lower, 0.0),
                upper,
                gid=MTBC_SPAN_GID,
                zorder=1,
                linewidth=0.0,
                **theme.HIGHLIGHT_SHADE,
            )
        else:
            label += " (no interval)"
        ax.axvline(
            point,
            color=theme.GUIDE_LABEL["color"],
            linewidth=1.2,
            zorder=5,
            label=label,
        )
        ax.axvline(
            data.tau_min,
            zorder=1,
            label=f"tau {data.tau_min:.1f} min",
            **theme.REFERENCE_LINE,
        )
        ax.text(
            0.98,
            0.98,
            _counts(data),
            transform=ax.transAxes,
            ha="right",
            va="top",
            **theme.scaled_text(theme.ANNOTATION, style),
        )
        ax.set_xlabel("Cycle length (min)")
        ax.set_ylabel("Surviving fraction of cycles")
        ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
        ax.legend(
            loc="center right",
            frameon=False,
            **theme.scaled_text(theme.LEGEND_TEXT, style),
        )
        theme.apply_common_style(ax)

    @staticmethod
    def _draw_checks(ax: plt.Axes, data: CyclesPanelData, style: str) -> None:
        rows = data.verdict_rows
        n = len(rows)
        on_fill = theme.scaled_text(theme.ON_FILL_TEXT, style)
        label_text = theme.scaled_text(theme.LABEL_TEXT, style)
        for i, (label, p_value, verdict, _note) in enumerate(rows):
            y = n - 1 - i
            ax.text(0.0, y + 0.5, label, ha="left", va="center", **label_text)
            ax.text(
                _P_X, y + 0.5, _fmt_p(p_value), ha="center", va="center", **label_text
            )
            ax.add_patch(
                plt.Rectangle(
                    (_CELL_X, y),
                    _CELL_W,
                    1,
                    facecolor=theme.verdict_color(verdict),
                    edgecolor="white",
                    linewidth=0.8,
                )
            )
            text_color = (
                "white" if verdict in (VERDICT_PASS, VERDICT_FAIL) else on_fill["color"]
            )
            ax.text(
                _CELL_X + _CELL_W / 2,
                y + 0.5,
                verdict,
                ha="center",
                va="center",
                **{**on_fill, "color": text_color},
            )
        if not rows:
            ax.text(0.0, 0.5, "no check asked", ha="left", va="center", **label_text)
        ax.text(_P_X, max(n, 1) + 0.25, "p", ha="center", va="bottom", **label_text)
        ax.text(
            _CELL_X + _CELL_W / 2,
            max(n, 1) + 0.25,
            "verdict",
            ha="center",
            va="bottom",
            **label_text,
        )
        ax.set_xlabel(
            "\n".join(_check_notes(data)),
            wrap=True,
            **theme.scaled_text(theme.CAPTION, style),
        )
        ax.set_xlim(0.0, 4.0)
        ax.set_ylim(0.0, max(n, 1) + 1.0)
        ax.set_title(
            f"Check outcome on cycles ({data.checks.assumption_id}), calendar clock",
            loc="left",
        )
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    @staticmethod
    def _draw_shares(ax: plt.Axes, data: CyclesPanelData, style: str) -> None:
        statuses = data.shares.statuses
        listed = [s for s in statuses if s != NOT_LISTED]
        y = np.arange(len(statuses))[::-1]
        annotation = theme.scaled_text(theme.ANNOTATION, style)
        for position, status, hours, fraction in zip(
            y, statuses, data.share_h, data.shares.fractions, strict=True
        ):
            if status == NOT_LISTED:
                ax.barh(
                    position,
                    hours,
                    color=theme.STATE_COLORS["unobserved"],
                    edgecolor=theme.UNOBSERVED_EDGE,
                    hatch=theme.UNOBSERVED_HATCH,
                )
            else:
                ax.barh(
                    position,
                    hours,
                    color=theme.ordered_color(listed.index(status), len(listed)),
                )
            ax.text(
                hours,
                position,
                f"  {hours:.1f} h ({fraction:.1%})",
                ha="left",
                va="center",
                **annotation,
            )
        ax.set_yticks(y)
        ax.set_yticklabels(list(statuses), **theme.scaled_text(theme.LABEL_TEXT, style))
        # Room on the right for the longest bar's label.
        ax.set_xlim(0.0, max(data.share_h, default=0.0) * 1.35 or 1.0)
        ax.set_xlabel("Time in status (h)")
        ax.set_title("Status shares of the record", loc="left")
        theme.apply_common_style(ax)


def _check_notes(data: CyclesPanelData) -> list[str]:
    """Under the table: the bench-transfer note, then why the checks without an answer gave
    none. When the battery scored nothing, its one reason covers every check; otherwise
    each such check's note, to its first clause (the whole note is on the artifact)."""
    notes = [data.bench_note]
    if data.battery_note:
        return [*notes, data.battery_note]
    for label, _p_value, verdict, note in data.verdict_rows:
        if verdict not in ANSWERED and note:
            notes.append(f"{label}: {re.split(r'[:;]', note, maxsplit=1)[0]}")
    return notes


def _counts(data: CyclesPanelData) -> str:
    """What the cycles keep, merge and drop, each count from the artifact."""
    return "\n".join(
        [
            f"cycles: {data.n_cycles} ({data.n_censored} censored)",
            f"excluded, start unknown: {data.n_unobserved_birth_dropped}",
            f"calibrations merged: {data.n_calibrations_merged}",
            f"segments: {data.n_segments}, too few cycles for the checks: "
            f"{data.n_segments_dropped}, all cycles equal: {data.n_segments_degenerate}",
            f"zero-length outages dropped: {data.n_zero_outages_dropped}",
            f"gaps: {data.n_unobserved_spans}",
            "duplicate timestamps dropped: "
            f"{data.n_duplicates_identical_dropped} identical, "
            f"{data.n_duplicates_conflicting_dropped} conflicting",
        ]
    )


def _caption(data: CyclesPanelData) -> str:
    """What is drawn and what is excluded, then the check outcome's one line."""
    return " ".join(
        [
            f"Kaplan-Meier survival of the calibration cycles of {data.dataset_id}, with its "
            f"{data.mtbc.conf_level:.0%} log-log band.",
            "A cycle runs from the end of one calibration to the start of the next; an "
            "outage, a gap or the record end censors it.",
            "MTBC is the restricted mean up to tau, the longest cycle, with its Wald "
            "interval.",
            f"Excluded: {data.n_unobserved_birth_dropped} stretches whose start is unknown "
            "(after the record start, an outage or a gap).",
            f"Check outcome: {data.check_summary}",
        ]
    )
