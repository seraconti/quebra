"""The comparison figures: a pair across the ladder, a forest at one rung, the pooled ladder.

Renderers only. Every number drawn is on the artifact a builder step made in
`analyzers/record_comparison.py` (`PairLadderData`, `ForestData`, `PooledLadderData`);
what stays here is a function of axes and theme. Ages, restricted means and tau are stored
in seconds and drawn in minutes, the unit of the Kaplan-Meier curves.

What each figure states in the image, per docs/FIGURE_STANDARD.md:

- the comparison tau and the record that set it, per panel, per group or per rung;
- what was dropped: windows whose birth was not observed (pair panels), and records with
  no windows or with no death before tau, by name (forest) or by count (ladder);
- the check outcome of what is drawn, which reads NOT ASSESSED while no step fills it,
  and, in each pair panel, each record's check and group-row tallies when filled;
- the evidence line of the pair or of each group, and how a pair was chosen when declared;
- in the caption, that a pair's p-values are unadjusted for the other thresholds and
  sides and assume the windows within each record are independent, and that groups drawn
  together are described, not tested.

A group with fewer than `record_comparison.MIN_POOLED` records with a variance has no
pooled numbers, so it is drawn as a table of its records, never as a forest. A flagged
record (`flags` on the artifact) carries a dagger beside its name, and the line of text
the flag holds is printed with the names it marks: in the group's notes for a forest, in
the caption for a pair.
"""

from __future__ import annotations

import textwrap

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from quebra.analyzers.event_table import SIDES
from quebra.analyzers.record_comparison import (
    TAU_ACROSS_GROUPS,
    ComparisonTau,
    ForestData,
    PairLadderData,
    PairPanel,
    PooledLadderData,
    RandomEffectsSummary,
)
from quebra.analyzers.windows import SIDE_IN_SPEC, SIDE_OUT_OF_SPEC
from quebra.plots import theme
from quebra.plots.base import BasePlot
from quebra.plots.km_survival_plot import band_vertices, step_arrays

_S_PER_MIN = 60.0
_SIDE_WORDS = {SIDE_IN_SPEC: "in spec", SIDE_OUT_OF_SPEC: "out of spec"}
# Record a solid, record b dash-dot: the pair separates by line style as well as by colour,
# and neither is dashed, the comparison tau's line (`theme.REFERENCE_LINE`).
_LINESTYLES = ("-", "-.")
# The dark and middle tones of the ordered ramp: its light end is too faint on white.
_PAIR_INKS = (theme.ordered_color(0, 3), theme.ordered_color(1, 3))
# One ink for forests and ladders: each group already has a panel of its own.
_GROUP_INK = theme.ordered_color(0, 3)
_WRAP = 58
# Rungs to a row of the pair figure; each rung is a column of two panels, in spec over
# out of spec, and the next rungs wrap onto a new pair of rows.
_COLUMNS = 4
DAGGER = "†"

UNADJUSTED = (
    "Each panel's log-rank p is unadjusted for the other thresholds and sides, and "
    "assumes the windows within each record are independent."
)
DESCRIBED_NOT_TESTED = "The groups are described, not tested: no test compares them."


def _known(value: float | None) -> float:
    """A number the drawn result guarantees; None means a malformed artifact, so raise."""
    if value is None:
        raise ValueError("a number the drawn result guarantees is None")
    return value


def _min(seconds: float | None) -> float:
    return _known(seconds) / _S_PER_MIN


def _minutes(seconds: float | None) -> str:
    return f"{_min(seconds):.3g}"


def _fill(text: str, width: int) -> str:
    # Breaks at spaces only: "log-rank" and a record name stay whole on one line.
    return textwrap.fill(text, width, break_long_words=False, break_on_hyphens=False)


def _wrapped(lines: list[str]) -> str:
    return "\n".join(_fill(line, _WRAP) for line in lines)


def _tau_line(tau: ComparisonTau, scope: str | None = None) -> str:
    where = "" if scope is None else f" ({_scope_words(scope)})"
    if tau.rmst_tau_s is None:
        return f"τ: no record with windows{where}"
    return f"τ = {_minutes(tau.rmst_tau_s)} min, set by {tau.tau_set_by}{where}"


def _marked(record: str, flags: dict[str, str]) -> str:
    return f"{record} {DAGGER}" if record in flags else record


def _flag_lines(records: list[str], flags: dict[str, str]) -> list[str]:
    """One line per flag text among `records`: the text, then the records it marks."""
    by_text: dict[str, list[str]] = {}
    for record in records:
        if record in flags:
            by_text.setdefault(flags[record], []).append(record)
    return [f"{DAGGER} {text}: {', '.join(names)}" for text, names in by_text.items()]


def _scope_words(scope: str) -> str:
    return "one τ for every group" if scope == TAU_ACROSS_GROUPS else "τ per group"


# ---------------------------------------------------------------------------
# A pair across the ladder
# ---------------------------------------------------------------------------


def _pair_lines(panel: PairPanel) -> list[str]:
    pair = panel.pair
    lines = [_tau_line(pair.tau)]
    if pair.not_computed_reason is not None:
        lines.append(f"not computed: {pair.not_computed_reason}")
    else:
        # A descriptive pair has no SE; its log-rank line below says why.
        difference = f"a - b = {_minutes(pair.rmst_difference_s)} min"
        if pair.rmst_difference_se_s is not None:
            difference += f", SE {_minutes(pair.rmst_difference_se_s)}"
        lines.append(difference)
        lines.append(
            f"log-rank p = {pair.logrank_p:.2g}"
            if pair.logrank_p is not None
            else f"log-rank not computed: {pair.logrank_not_computed_reason}"
        )
    lines.append(
        f"windows (deaths): a {pair.n_windows_a} ({pair.n_deaths_a}), "
        f"b {pair.n_windows_b} ({pair.n_deaths_b})"
    )
    lines.append(
        f"birth unseen, excluded: a {panel.curve_a.n_unobserved_birth_dropped}, "
        f"b {panel.curve_b.n_unobserved_birth_dropped}"
    )
    return lines + list(panel.check_lines)


def _pair_panel(ax: plt.Axes, panel: PairPanel, style: str) -> None:
    x_hi = 0.0
    for index, curve in enumerate((panel.curve_a, panel.curve_b)):
        if curve.n_windows == 0:
            continue
        ink = _PAIR_INKS[index]
        x, s, lo, hi = step_arrays(curve, 0.0)
        ax.fill_between(
            *band_vertices(x, lo, hi),
            color=ink,
            alpha=theme.BAND_STYLE["fill_alpha"],
            linewidth=0.0,
        )
        ax.step(x, s, where="post", color=ink, linestyle=_LINESTYLES[index])
        if len(curve.censor_time_min):
            ax.plot(
                curve.censor_time_min,
                curve.censor_survival,
                linestyle="none",
                color=ink,
                **theme.CENSOR_TICK,
            )
        x_hi = max(x_hi, float(x[-1]))
    tau = panel.pair.tau.rmst_tau_s
    if tau is not None:
        ax.axvline(_min(tau), **theme.REFERENCE_LINE)
        x_hi = max(x_hi, _min(tau))
    ax.set_xlim(0.0, 1.05 * x_hi if x_hi > 0.0 else 1.0)
    ax.set_ylim(-0.02, 1.04)
    ax.text(
        0.03,
        0.03,
        _wrapped(_pair_lines(panel)),
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85},
        **theme.scaled_text(theme.ANNOTATION, style),
    )


class PairLadderPlot(BasePlot):
    """Two records' window survival at each drawn threshold, in spec over out of spec.

    The thresholds run along a row, `_COLUMNS` to a row, then wrap; which thresholds are
    drawn is the artifact's `ladder` (`record_comparison.select_rungs`).
    """

    def build_matplotlib(
        self, data: PairLadderData, style: str = "default"
    ) -> plt.Figure:
        if not isinstance(data, PairLadderData):
            raise TypeError(
                f"PairLadderPlot draws a PairLadderData, not {type(data).__name__}"
            )
        n = len(data.ladder)
        cols = min(n, _COLUMNS)
        blocks = -(-n // _COLUMNS)
        rows = len(SIDES) * blocks
        with theme.style_context(style):
            fig, axes = plt.subplots(
                rows, cols, figsize=(3.4 * cols, 3.7 * rows), squeeze=False, sharey=True
            )
            for index, label in enumerate(data.ladder):
                block, col = divmod(index, _COLUMNS)
                for row, side in enumerate(SIDES):
                    ax = axes[block * len(SIDES) + row][col]
                    _pair_panel(ax, data.panels[(label, side)], style)
                    ax.set_title(f"{_SIDE_WORDS[side]}, {label}")
                    if row == len(SIDES) - 1:
                        ax.set_xlabel("Window age (min)")
            # The last row of thresholds may be short: its empty slots are not drawn.
            for index in range(n, blocks * cols):
                block, col = divmod(index, _COLUMNS)
                for row in range(len(SIDES)):
                    axes[block * len(SIDES) + row][col].remove()
            for row in range(rows):
                axes[row][0].set_ylabel("Surviving fraction of windows")
            axes[0][0].legend(
                handles=[
                    Line2D(
                        [],
                        [],
                        color=_PAIR_INKS[index],
                        linestyle=_LINESTYLES[index],
                        label=f"{name}: {_marked(record, data.flags)}",
                    )
                    for index, (name, record) in enumerate(
                        (("a", data.record_a), ("b", data.record_b))
                    )
                ]
                + [Line2D([], [], **theme.REFERENCE_LINE, label="comparison tau")],
                loc="upper right",
                frameon=False,
                **theme.scaled_text(theme.LEGEND_TEXT, style),
            )
            fig.suptitle(f"Window survival, {data.record_a} vs {data.record_b}")
            panels = list(data.panels.values())
            checks = (
                "each panel's tallies, every cell in the table"
                if any(p.check_lines for p in panels)
                else " | ".join(dict.fromkeys(p.pair.check_summary() for p in panels))
            )
            rule = panels[0].pair.tau.tau_rule
            fig.text(
                0.01,
                0.01,
                _fill(
                    f"{UNADJUSTED} Dashed line: the comparison tau ({rule}); a - b is "
                    f"the difference of the restricted means up to it. Bands: "
                    f"{panels[0].pair.conf_level:.0%} pointwise. Evidence: "
                    f"{data.evidence}. Check outcome: {checks}. "
                    + (f"{data.declared}. " if data.declared else "")
                    + " ".join(
                        f"{line}."
                        for line in _flag_lines(
                            [data.record_a, data.record_b], data.flags
                        )
                    ),
                    40 * cols,
                ),
                **theme.scaled_text(theme.CAPTION, style),
            )
            fig.tight_layout(rect=(0.0, 0.05, 1.0, 0.96))
        return fig


# ---------------------------------------------------------------------------
# A forest at one threshold and side
# ---------------------------------------------------------------------------


def _not_pooled(summary: RandomEffectsSummary) -> list[str]:
    return [row.record for row in summary.rows if not row.pooled]


def _group_lines(
    group: str, data: ForestData, summary: RandomEffectsSummary
) -> list[str]:
    conf = f"{data.conf_level:.0%}"
    lines = [
        group,
        f"evidence: {data.evidence[group]}",
        _tau_line(summary.tau, data.tau_scope),
    ]
    if summary.pooled_s is None:
        lines.append(f"listed, not pooled: {summary.not_pooled_reason}")
    else:
        min2 = _S_PER_MIN**2
        lines += [
            f"pooled, k = {summary.k_pooled}: {_minutes(summary.pooled_s)} min, "
            f"{conf} CI [{_minutes(summary.pooled_lower_s)}, "
            f"{_minutes(summary.pooled_upper_s)}]",
            f"prediction interval [{_minutes(summary.prediction_lower_s)}, "
            f"{_minutes(summary.prediction_upper_s)}] min",
            f"between-record variance {_known(summary.between_var_s2) / min2:.3g} "
            f"min², {conf} Q-profile "
            f"[{_known(summary.between_var_lower_s2) / min2:.3g}, "
            f"{_known(summary.between_var_upper_s2) / min2:.3g}]; I² = "
            f"{_known(summary.i2_fraction):.0%}; Q = {_known(summary.q_statistic):.3g} on "
            f"{summary.q_df} df",
        ]
    no_windows = [record for record, _ in summary.tau.excluded]
    no_death = _not_pooled(summary)
    lines += [
        f"no windows ({len(no_windows)}): {', '.join(no_windows) or 'none'}",
        f"no death before τ, not pooled ({len(no_death)}): "
        f"{', '.join(no_death) or 'none'}",
        "birth unseen, excluded: {} of {} windows of the records with windows".format(
            *data.birth_unseen[group]
        ),
        f"check outcome: {summary.check_line or summary.check_summary()}",
    ]
    drawn = [row.record for row in summary.rows] + no_windows
    return lines + _flag_lines(drawn, data.flags)


def _forest(
    ax: plt.Axes, summary: RandomEffectsSummary, ink: str, flags: dict[str, str]
) -> None:
    rows = summary.rows
    y = np.arange(len(rows), 0, -1, dtype=float)
    for yi, row in zip(y, rows, strict=True):
        mean = _min(row.summary.rmst_s)
        if row.pooled:
            ax.plot(
                [_min(row.lower_s), _min(row.upper_s)],
                [yi, yi],
                color=ink,
            )
            ax.plot(mean, yi, marker="s", color=ink, linestyle="none")
        else:
            ax.plot(
                mean,
                yi,
                marker="o",
                markerfacecolor="white",
                color=ink,
                linestyle="none",
            )
    mu = _min(summary.pooled_s)
    lo, hi = _min(summary.pooled_lower_s), _min(summary.pooled_upper_s)
    ax.fill([lo, mu, hi, mu], [0.0, 0.3, 0.0, -0.3], color=ink)
    ax.plot(
        [
            _min(summary.prediction_lower_s),
            _min(summary.prediction_upper_s),
        ],
        [0.0, 0.0],
        color=ink,
        linestyle=":",
    )
    ax.axvline(_min(summary.tau.rmst_tau_s), **theme.REFERENCE_LINE)
    ax.set_yticks(
        [*y, 0.0],
        [
            *(
                _marked(row.record, flags)
                if row.pooled
                else f"{_marked(row.record, flags)} (no death before τ)"
                for row in rows
            ),
            f"pooled, k = {summary.k_pooled}",
        ],
    )
    ax.set_ylim(-1.0, len(rows) + 1.0)
    ax.set_xlabel("Restricted mean window age (min)")


def _record_table(
    ax: plt.Axes, summary: RandomEffectsSummary, flags: dict[str, str], style: str
) -> None:
    """The records of a group too small to pool, as text, since no forest is drawn."""
    ax.axis("off")
    lines = []
    for row in summary.rows:
        mean = _minutes(row.summary.rmst_s)
        name = _marked(row.record, flags)
        lines.append(
            f"{name}: {mean} min, SE {_minutes(row.summary.rmst_se_s)}"
            if row.pooled
            else f"{name}: {mean} min = τ, no death before τ"
        )
    lines += [
        f"{_marked(record, flags)}: no windows" for record, _ in summary.tau.excluded
    ]
    ax.text(
        0.0,
        1.0,
        _wrapped(["records, not pooled:", *lines]),
        transform=ax.transAxes,
        ha="left",
        va="top",
        **theme.scaled_text(theme.ANNOTATION, style),
    )


def _forest_height(summary: RandomEffectsSummary) -> float:
    # Rows of the forest or of the table, and room for the group's notes beside it.
    return float(max(len(summary.rows) + len(summary.tau.excluded) + 3, 6))


class ForestPlot(BasePlot):
    """Each group's records and pooled restricted mean at one threshold and side."""

    def build_matplotlib(self, data: ForestData, style: str = "default") -> plt.Figure:
        if not isinstance(data, ForestData):
            raise TypeError(f"ForestPlot draws a ForestData, not {type(data).__name__}")
        groups = list(data.summaries)
        heights = [_forest_height(data.summaries[group]) for group in groups]
        with theme.style_context(style):
            fig = plt.figure(figsize=(13.0, 1.4 + 0.3 * sum(heights)))
            grid = fig.add_gridspec(
                len(groups), 2, width_ratios=(3.0, 2.4), height_ratios=heights
            )
            shared: plt.Axes | None = None
            for index, group in enumerate(groups):
                summary = data.summaries[group]
                pooled = summary.pooled_s is not None
                ax = fig.add_subplot(
                    grid[index, 0],
                    sharex=shared if data.tau_scope == TAU_ACROSS_GROUPS else None,
                )
                if pooled:
                    _forest(ax, summary, _GROUP_INK, data.flags)
                    shared = shared or ax
                else:
                    _record_table(ax, summary, data.flags, style)
                ax.set_title(group, loc="left")
                note = fig.add_subplot(grid[index, 1])
                note.axis("off")
                note.text(
                    0.0,
                    1.0,
                    _wrapped(_group_lines(group, data, summary)),
                    transform=note.transAxes,
                    ha="left",
                    va="top",
                    **theme.scaled_text(theme.ANNOTATION, style),
                )
            fig.suptitle(
                f"Restricted mean window age, {_SIDE_WORDS[data.side]} at "
                f"{data.threshold_label}"
            )
            caption = [
                f"Squares: each record's restricted mean with its {data.conf_level:.0%} "
                f"interval; diamond: the pooled mean and its interval; dotted: the "
                f"prediction interval for a new record's true restricted mean; dashed: "
                f"the comparison tau ({_scope_words(data.tau_scope)}). Hollow circles: "
                f"records with no death before tau, left out of the pool."
            ]
            if len(groups) > 1:
                caption.append(DESCRIBED_NOT_TESTED)
            fig.text(
                0.01,
                0.005,
                _fill(" ".join(caption), 160),
                **theme.scaled_text(theme.CAPTION, style),
            )
            fig.tight_layout(rect=(0.0, 0.03, 1.0, 0.97))
        return fig


# ---------------------------------------------------------------------------
# The pooled ladder: every threshold and side
# ---------------------------------------------------------------------------


def _ladder_lines(data: PooledLadderData, group: str, side: str) -> list[str]:
    lines = [f"evidence: {data.evidence[group]}"]
    checks: dict[str, None] = {}
    for label in data.ladder:
        summary = data.summaries[(group, label, side)]
        tau = summary.tau
        if tau.rmst_tau_s is None:
            lines.append(f"{label}: no record with windows")
        else:
            pooled = (
                f"pooled {summary.k_pooled} of {len(summary.rows)}"
                if summary.pooled_s is not None
                else f"not pooled, {summary.k_pooled} with a variance"
            )
            dropped, carved = data.birth_unseen[(group, label, side)]
            lines.append(
                f"{label}: τ {_minutes(tau.rmst_tau_s)} min, set by {tau.tau_set_by}; "
                f"{pooled}, {len(tau.excluded)} without windows; birth unseen, "
                f"excluded {dropped} of {carved} windows"
            )
        line = f"{label}: {summary.check_line}" if summary.check_line else None
        checks[line or summary.check_summary()] = None
    lines.append(f"check outcome: {' | '.join(checks)}")
    return lines


class PooledLadderPlot(BasePlot):
    """Each group's pooled restricted mean against threshold, one row per group."""

    def build_matplotlib(
        self, data: PooledLadderData, style: str = "default"
    ) -> plt.Figure:
        if not isinstance(data, PooledLadderData):
            raise TypeError(
                f"PooledLadderPlot draws a PooledLadderData, not {type(data).__name__}"
            )
        x = np.arange(len(data.ladder), dtype=float)
        with theme.style_context(style):
            fig = plt.figure(figsize=(18.0, 1.4 + 3.4 * len(data.groups)))
            grid = fig.add_gridspec(len(data.groups), 4, width_ratios=(3, 2.2, 3, 2.2))
            for index, group in enumerate(data.groups):
                for column, side in enumerate(SIDES):
                    ax = fig.add_subplot(grid[index, 2 * column])
                    summaries = [
                        data.summaries[(group, label, side)] for label in data.ladder
                    ]
                    taus = [s.tau.rmst_tau_s for s in summaries]
                    with_tau = [i for i, tau in enumerate(taus) if tau is not None]
                    twin = ax.twinx()
                    twin.plot(
                        x[with_tau],
                        [_min(taus[i]) for i in with_tau],
                        **theme.REFERENCE_LINE,
                    )
                    twin.set_ylabel("Comparison tau (min)")
                    for i, summary in enumerate(summaries):
                        if summary.pooled_s is None:
                            continue
                        mu = _min(summary.pooled_s)
                        for lower, upper, width in (
                            (
                                summary.prediction_lower_s,
                                summary.prediction_upper_s,
                                0.8,
                            ),
                            (summary.pooled_lower_s, summary.pooled_upper_s, 2.4),
                        ):
                            ax.plot(
                                [x[i], x[i]],
                                [_min(lower), _min(upper)],
                                color=_GROUP_INK,
                                linewidth=width,
                            )
                        ax.plot(
                            x[i], mu, marker="D", color=_GROUP_INK, linestyle="none"
                        )
                    ax.set_xticks(x, data.ladder)
                    ax.set_title(f"{group}, {_SIDE_WORDS[side]}", loc="left")
                    ax.set_ylabel("Pooled restricted mean (min)")
                    ax.set_xlabel("Threshold")
                    note = fig.add_subplot(grid[index, 2 * column + 1])
                    note.axis("off")
                    note.text(
                        0.0,
                        1.0,
                        _wrapped(_ladder_lines(data, group, side)),
                        transform=note.transAxes,
                        ha="left",
                        va="top",
                        **theme.scaled_text(theme.ANNOTATION, style),
                    )
            fig.suptitle("Pooled restricted mean window age by threshold")
            caption = [
                f"Dashed, on the right axis: the comparison tau at each threshold "
                f"({_scope_words(data.tau_scope)}), so two thresholds' restricted means "
                f"are not on one restriction time. Diamonds: the pooled mean; thick bar: "
                f"its {data.conf_level:.0%} interval; thin bar: the prediction interval. "
                f"A threshold with no diamond was not pooled."
            ]
            if len(data.groups) > 1:
                caption.append(DESCRIBED_NOT_TESTED)
            fig.text(
                0.01,
                0.005,
                _fill(" ".join(caption), 200),
                **theme.scaled_text(theme.CAPTION, style),
            )
            fig.tight_layout(rect=(0.0, 0.03, 1.0, 0.97))
        return fig
