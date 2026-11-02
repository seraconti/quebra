"""The compliance survey: the in-spec fraction of every record at every rung of the ladder.

One panel per device, side by side on a shared ladder axis; one row per record of the
paper's dataset table. Each row is a small bar chart: a bar per rung, its height the
time-weighted in-spec fraction (the row spans 0 to 1), its whisker the lower and upper
bounds. A record the paper shades sits on a shaded row, and the legend carries the flag
text. A rung or record that was not computed is a grey cell, never an empty slot; a rung
whose fraction was computed but whose bounds were not carries a dotted whisker across
its whole row, never no whisker, which would read as bounds equal to the fraction. The
panel prints both counts from the artifact. A carved record outside the paper's table is
listed under its panel, not drawn.

PURE RENDERER. Every value drawn is a field of `ComplianceSurveyData`. Each bar and whisker
carries a `gid` naming its record and rung, so a test can read the figure back.

Follows `docs/FIGURE_STANDARD.md`: the title names the object, units sit in the axis label,
the caption states what is drawn and what is not, spaced hyphens throughout.
"""

from __future__ import annotations

import textwrap

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

from quebra.analyzers.compliance_survey import ComplianceSurveyData, SurveyRecord
from quebra.plots import theme
from quebra.plots.base import BasePlot

# The height of a fraction of 1 within a record's row, in row units; the rest is space.
ROW_HEIGHT = 0.78
BAR_WIDTH = 0.62
# The rules a row is read against, its 0 and its 1: the neutral leader, thinner.
_GUIDE = {**theme.LEADER_LINE, "linewidth": 0.6, "zorder": 0}
# Bounds not computed: the whisker dotted, spanning the row from 0 to 1. The not-computed
# grey is too faint on white for a line this thin.
_NO_BOUNDS = {**theme.BOUNDS_WHISKER, "linestyle": ":"}


def _check(result: object) -> ComplianceSurveyData:
    if not isinstance(result, ComplianceSurveyData):
        raise TypeError(
            f"{type(result).__name__} is not ComplianceSurveyData; this figure draws "
            "only the typed artifact its builder returns."
        )
    return result


def bar_gid(record: str, label: str) -> str:
    return f"bar|{record}|{label}"


def bounds_gid(record: str, label: str) -> str:
    return f"bounds|{record}|{label}"


def no_bounds_gid(record: str, label: str) -> str:
    return f"no_bounds|{record}|{label}"


class ComplianceSurveyPlot(BasePlot):
    def build_matplotlib(self, result: object, style: str = "default") -> plt.Figure:
        data = _check(result)
        panels = [(d, data.drawn(d)) for d in data.devices]
        if not any(rows for _, rows in panels):
            raise ValueError(
                "no record of the paper's table is on the artifact; a blank figure would lie"
            )
        n_rows = max(len(rows) for _, rows in panels)
        width_in = 4.4 * len(panels) + 1.2
        with theme.style_context(style):
            fig, axes = plt.subplots(
                1,
                len(panels),
                sharex=True,
                squeeze=False,
                figsize=(width_in, 0.27 * n_rows + 3.0),
                facecolor="white",
            )
            for ax, (device, rows) in zip(axes[0], panels):
                self._draw_panel(ax, data, device, rows, n_rows)

            k = data.unsure_k
            legend = [
                Patch(
                    facecolor=theme.STATE_COLORS["in_spec"],
                    label="in-spec fraction of observed time",
                ),
                Line2D(
                    [],
                    [],
                    **theme.BOUNDS_WHISKER,
                    label=f"lower and upper bounds at k = {k:g}"
                    if k is not None
                    else "lower and upper bounds",
                ),
            ]
            if any(row.flag for _, rows in panels for row in rows):
                legend.append(
                    Patch(
                        facecolor=theme.HIGHLIGHT_SHADE["color"],
                        alpha=theme.HIGHLIGHT_SHADE["alpha"],
                        label=data.flag_text,
                    )
                )
            legend.append(
                Patch(
                    facecolor=theme.verdict_color("not computed"), label="not computed"
                )
            )
            if any(data.counts[d].bounds_not_computed for d, _ in panels):
                legend.append(Line2D([], [], **_NO_BOUNDS, label="bounds not computed"))
            fig.legend(
                handles=legend,
                loc="lower center",
                ncol=2,
                frameon=False,
                bbox_to_anchor=(0.5, 0.05),
                **theme.LEGEND_TEXT,
            )
            gap = (
                f" (gap_mult = {data.gap_mult:g})" if data.gap_mult is not None else ""
            )
            caption = (
                "Each row spans an in-spec fraction of 0 to 1. The fraction is of observed "
                f"time: an interval longer than the gap rule is not observed{gap}. The "
                "lower bound counts every interval not settled at k out of spec, the "
                "upper bound counts it in spec. Records and shading: "
                f"{data.table_source}."
            )
            # Wrapped to the figure's width at the caption's size: matplotlib's own wrap
            # lets a centred line run past the edge.
            per_line = int(width_in * 72 / (0.55 * theme.CAPTION["fontsize"]) * 0.92)
            fig.text(
                0.5,
                0.008,
                textwrap.fill(caption, per_line),
                ha="center",
                va="bottom",
                **theme.CAPTION,
            )
            fig.suptitle("In-spec fraction on the T2* ladder")
            fig.tight_layout(rect=(0, 0.13, 1, 0.97))
        return fig

    def _draw_panel(
        self,
        ax: plt.Axes,
        data: ComplianceSurveyData,
        device: str,
        rows: tuple[SurveyRecord, ...],
        n_rows: int,
    ) -> None:
        n_rungs = len(data.thresholds)
        grey = theme.verdict_color("not computed")
        for i, row in enumerate(rows):
            base = n_rows - 1 - i
            if row.flag:
                ax.axhspan(base - 0.11, base + 0.89, zorder=0, **theme.HIGHLIGHT_SHADE)
            ax.hlines([base, base + ROW_HEIGHT], -0.5, n_rungs - 0.5, **_GUIDE)
            if row.compliance.record_failure:
                # One band across the row, so the failure's note reads on it.
                ax.add_patch(
                    Rectangle(
                        (-0.5 + (1 - BAR_WIDTH) / 2, base),
                        n_rungs - (1 - BAR_WIDTH),
                        ROW_HEIGHT,
                        facecolor=grey,
                        edgecolor="none",
                    )
                )
                ax.text(
                    -0.4,
                    base + ROW_HEIGHT / 2,
                    f"not computed: {row.compliance.record_failure}",
                    va="center",
                    clip_on=True,
                    **theme.ANNOTATION,
                )
                continue
            for j, rung in enumerate(row.compliance.rungs):
                if rung.in_spec_fraction is None:
                    ax.add_patch(
                        Rectangle(
                            (j - BAR_WIDTH / 2, base),
                            BAR_WIDTH,
                            ROW_HEIGHT,
                            facecolor=grey,
                            edgecolor="none",
                        )
                    )
                else:
                    (bar,) = ax.bar(
                        j,
                        ROW_HEIGHT * rung.in_spec_fraction,
                        bottom=base,
                        width=BAR_WIDTH,
                        color=theme.STATE_COLORS["in_spec"],
                        linewidth=0,
                    )
                    bar.set_gid(bar_gid(row.record, rung.threshold_label))
                if rung.in_spec_fraction_lower is not None:
                    whisker = ax.vlines(
                        j,
                        base + ROW_HEIGHT * rung.in_spec_fraction_lower,
                        base + ROW_HEIGHT * rung.in_spec_fraction_upper,
                        zorder=3,
                        **theme.BOUNDS_WHISKER,
                    )
                    whisker.set_gid(bounds_gid(row.record, rung.threshold_label))
                elif rung.in_spec_fraction is not None:
                    whisker = ax.vlines(
                        j, base, base + ROW_HEIGHT, zorder=3, **_NO_BOUNDS
                    )
                    whisker.set_gid(no_bounds_gid(row.record, rung.threshold_label))

        ax.set_xlim(-0.5, n_rungs - 0.5)
        ax.set_ylim(n_rows - len(rows) - 0.3, n_rows)
        ax.set_xticks(range(n_rungs))
        ax.set_xticklabels(
            [f"{value * 1e6:g}" for _, value, _ in data.thresholds], **theme.LABEL_TEXT
        )
        ax.set_xlabel("T2* threshold (µs)")
        ax.set_yticks([n_rows - 1 - i + ROW_HEIGHT / 2 for i in range(len(rows))])
        ax.set_yticklabels([row.label for row in rows], **theme.LABEL_TEXT)
        ax.set_title(f"{device}, {len(rows)} records", loc="left", **theme.LABEL_TEXT)
        ax.grid(False)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        ax.tick_params(axis="y", length=0)

        counts = data.counts[device]
        notes = [
            f"not computed: {counts.fraction_not_computed} of {counts.n_rungs} rungs",
            f"bounds not computed beside a drawn fraction: {counts.bounds_not_computed} "
            "rungs",
        ]
        extra = data.not_drawn(device)
        if extra:
            notes.append(
                f"not drawn, not in the paper's dataset table: {len(extra)} records ("
                + ", ".join(r.label for r in extra)
                + ")"
            )
        ax.text(
            0.0,
            -0.12 if n_rows < 12 else -0.06,
            "\n".join(notes),
            transform=ax.transAxes,
            va="top",
            wrap=True,
            **theme.ANNOTATION,
        )
