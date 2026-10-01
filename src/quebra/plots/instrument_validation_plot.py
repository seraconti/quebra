"""Figures for the instrument report. Four views of one artifact, drawn separately.

Separate figures rather than one panel, on purpose: these answer four different questions
and a reader normally wants one of them at a time. The tier matrix is the index; the other
three are the evidence behind three of its columns.

PURE RENDERER. Every number drawn here is a field of `InstrumentValidationData`, computed
by `analyzers.instrument_validation`. Only axis and colour concerns live in this file.

Follows `docs/FIGURE_STANDARD.md`: titles name the object and do not argue, units live in
the axis label and not the ticks, spaced hyphens throughout, and every view that drops or
omits data says so in the panel. Colours come from `plots.theme` only - the style ratchet
is 17 and a new hex would fail it.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from quebra.analyzers.instrument_validation import (
    TIER_ABSENT,
    TIER_FAIL,
    TIER_PARTIAL,
    TIER_PASS,
    InstrumentValidationData,
)
from quebra.plots import theme
from quebra.plots.base import BasePlot

# The tier vocabulary maps onto the ledger's verdict palette rather than inventing a second
# one. "absent" reuses `not computed` because that is exactly what it means here - no
# evidence was gathered - and the ledger already desaturates it so it cannot be misread as
# a result.
_TIER_VERDICT_COLOR = {
    TIER_PASS: "pass",
    TIER_FAIL: "fail",
    TIER_PARTIAL: "underpowered",
    TIER_ABSENT: "not computed",
}

_TIER_TITLES = {
    2: "Tier 2 - published values",
    3: "Tier 3 - calibration",
    4: "Tier 4 - cross-implementation",
}


def _check(result: object) -> InstrumentValidationData:
    if not isinstance(result, InstrumentValidationData):
        raise TypeError(
            f"{type(result).__name__} is not InstrumentValidationData; the figure draws "
            "only the typed artifact its builder returns."
        )
    return result


def _tier_color(verdict: str) -> str:
    return theme.verdict_color(_TIER_VERDICT_COLOR[verdict])


class TierMatrixPlot(BasePlot):
    """The index figure: which instrument has which evidence.

    Reads as a status table, and that is a deliberate exception to this repo's rule against
    status tables in DOCS - the rule exists because a table in a doc goes stale silently.
    This one is regenerated from the artifact on every run, so it cannot.
    """

    def build_matplotlib(self, result: object, style: str = "default") -> plt.Figure:
        pd_ = _check(result)
        instruments = pd_.instruments
        tiers = [2, 3, 4]

        with theme.style_context(style):
            fig, ax = plt.subplots(
                figsize=(9.5, 0.52 * len(instruments) + 2.4), facecolor="white"
            )
            for row, instrument in enumerate(instruments):
                for col, tier in enumerate(tiers):
                    entry = pd_.tier_verdict(instrument, tier)
                    verdict = entry.verdict if entry else TIER_ABSENT
                    ax.add_patch(
                        plt.Rectangle(
                            (col, row),
                            1,
                            1,
                            facecolor=_tier_color(verdict),
                            edgecolor="white",
                            linewidth=1.6,
                        )
                    )
                    ax.text(
                        col + 0.5,
                        row + 0.5,
                        verdict,
                        ha="center",
                        va="center",
                        # Unpack the whole style, then override only the colour on the
                        # two saturated fills where dark text would not read.
                        **{
                            **theme.ON_FILL_TEXT,
                            "color": "white"
                            if verdict in (TIER_PASS, TIER_FAIL)
                            else theme.ON_FILL_TEXT["color"],
                        },
                    )
            ax.set_xlim(0, len(tiers))
            ax.set_ylim(0, len(instruments))
            ax.set_xticks(np.arange(len(tiers)) + 0.5)
            ax.set_xticklabels([_TIER_TITLES[t] for t in tiers], **theme.LABEL_TEXT)
            ax.set_yticks(np.arange(len(instruments)) + 0.5)
            ax.set_yticklabels(instruments, **theme.LABEL_TEXT)
            ax.invert_yaxis()
            ax.set_title("Validation evidence by instrument and tier")
            for spine in ax.spines.values():
                spine.set_visible(False)
            ax.tick_params(length=0)
            # FIGURE_STANDARD: say what is not here, in the panel.
            fig.text(
                0.5,
                0.015,
                "Tier 1 (equation review) is not drawn - it has no number. "
                "'absent' means no evidence was gathered, not evidence of failure. "
                "None of this scores power: see jobs/bench/results/promotion_report.md.",
                ha="center",
                **theme.CAPTION,
                wrap=True,
            )
            fig.tight_layout(rect=(0, 0.06, 1, 1))
        return fig


class PublishedValuesPlot(BasePlot):
    """Tier 2: our numbers against the ones Kvaloy and Lindqvist print, on their data."""

    def build_matplotlib(self, result: object, style: str = "default") -> plt.Figure:
        pd_ = _check(result)
        rows = pd_.published
        labels = [r.quantity for r in rows]
        rel = [
            100.0 * r.relative_difference if np.isfinite(r.relative_difference) else 0.0
            for r in rows
        ]
        colors = [
            theme.verdict_color("pass" if r.agrees else "underpowered") for r in rows
        ]

        with theme.style_context(style):
            fig, (ax_val, ax_rel) = plt.subplots(
                1,
                2,
                figsize=(13.5, 0.52 * len(rows) + 2.6),
                width_ratios=[1.25, 1.0],
                facecolor="white",
            )
            y = np.arange(len(rows))

            # Left: the two numbers together, on a log axis because they span 0.6 to 55.
            # MARKERS, not bars: on a log axis a bar grows from an arbitrary baseline, so
            # its length encodes nothing. The connector between the pair is the quantity a
            # reader actually wants, and on four of seven rows it is invisibly short.
            for i, row in enumerate(rows):
                ax_val.plot(
                    [row.published, row.ours],
                    [i, i],
                    zorder=1,
                    **theme.PAIR_CONNECTOR,
                )
            ax_val.scatter(
                [r.published for r in rows],
                y,
                s=70,
                marker="o",
                facecolor="white",
                edgecolor=theme.REFERENCE_LINE["color"],
                linewidth=1.6,
                zorder=3,
                label="published",
            )
            ax_val.scatter(
                [r.ours for r in rows],
                y,
                s=48,
                marker="D",
                color=colors,
                zorder=4,
                label="this tool",
            )
            ax_val.set_xscale("log")
            ax_val.set_yticks(y)
            ax_val.set_yticklabels(labels, **theme.LABEL_TEXT)
            ax_val.invert_yaxis()
            ax_val.set_xlabel("Value (mixed units - see quantity)")
            ax_val.set_title("Published value against ours")
            ax_val.legend(loc="lower right", **theme.LEGEND_TEXT)
            ax_val.grid(axis="x", alpha=0.25)

            # Right: the relative difference, which is the only comparable axis across
            # quantities that do not share a unit.
            ax_rel.barh(y, rel, color=colors, height=0.55)
            ax_rel.axvline(0.0, **theme.REFERENCE_LINE)
            ax_rel.set_yticks(y)
            ax_rel.set_yticklabels([])
            ax_rel.invert_yaxis()
            ax_rel.set_xlabel("Difference from published (%)")
            ax_rel.set_title("Size of the disagreement")
            for i, row in enumerate(rows):
                ax_rel.text(
                    max(rel[i], 0.0) + 0.05,
                    i,
                    row.note,
                    va="center",
                    **theme.ANNOTATION,
                )
            ax_rel.set_xlim(left=-0.15)

            n_disagree = sum(1 for r in rows if not r.agrees)
            n_gaps = pd_.meta["load_haul_dump_complete_gaps"]
            tau_h = pd_.meta["load_haul_dump_tau_h"]
            n_censored = pd_.dropped["load_haul_dump_censored_gaps"]
            # A record ends either on a failure or at the censoring time; only the second
            # leaves a gap the complete-gap estimators drop.
            record = (
                f"time censored at {tau_h:.10g} h, {n_censored} censored gap dropped"
                if n_censored
                else f"observed to {tau_h:.10g} h, ending on a failure, so no gap is censored"
            )
            fig.suptitle(
                "Transcription against published values - load-haul-dump record, "
                "Kvaloy and Lindqvist (Technometrics 2020) Section 8.1"
            )
            fig.text(
                0.5,
                0.015,
                f"{n_gaps} complete gaps, {record}. "
                f"{len(rows) - n_disagree} of {len(rows)} quantities reproduce exactly"
                + (
                    f"; {n_disagree} differ, each row's note says why."
                    if n_disagree
                    else "."
                ),
                ha="center",
                **theme.CAPTION,
                wrap=True,
            )
            fig.tight_layout(rect=(0, 0.08, 1, 0.94))
        return fig


class CrossImplementationPlot(BasePlot):
    """Tier 4: agreement with the R packages, most of them the authors' own."""

    def build_matplotlib(self, result: object, style: str = "default") -> plt.Figure:
        pd_ = _check(result)
        rows = pd_.cross_implementation
        exact = [r for r in rows if not r.reference_is_random]
        random_rows = [r for r in rows if r.reference_is_random]

        with theme.style_context(style):
            fig, (ax_diff, ax_rand) = plt.subplots(
                1,
                2,
                figsize=(13.5, 0.5 * len(rows) + 3.0),
                width_ratios=[1.4, 1.0],
                facecolor="white",
            )

            # Left: |ours - R| on a log axis. Everything should sit far below the tolerance
            # line; the line is what makes "far below" legible rather than asserted.
            labels = [f"{r.statistic} - {r.case}" for r in exact]
            diffs = [max(r.abs_difference, 1e-18) for r in exact]
            y = np.arange(len(exact))
            # Markers for the same reason as the published-values view.
            ax_diff.hlines(y, 1e-18, diffs, **theme.LEADER_LINE)
            ax_diff.scatter(
                diffs, y, s=54, marker="D", color=theme.verdict_color("pass"), zorder=3
            )
            ax_diff.axvline(
                pd_.meta["agreement_tol"],
                label=f"agreement tolerance ({pd_.meta['agreement_tol']:.0e})",
                **theme.REFERENCE_LINE,
            )
            ax_diff.set_xscale("log")
            ax_diff.set_yticks(y)
            ax_diff.set_yticklabels(labels, **theme.LABEL_TEXT)
            ax_diff.invert_yaxis()
            ax_diff.set_xlabel("Absolute difference from the R value")
            ax_diff.set_title("Agreement where the reference is a number")
            ax_diff.legend(loc="lower right", **theme.LEGEND_TEXT)

            # Right: the case where the reference is NOT a number. Drawn as a distribution
            # with our deterministic value marked, because that is what the comparison is.
            if random_rows:
                row = random_rows[0]
                spread = row.reference_spread
                grid = np.linspace(
                    row.reference - 4 * spread, row.reference + 4 * spread, 200
                )
                density = np.exp(-0.5 * ((grid - row.reference) / spread) ** 2)
                ax_rand.fill_between(
                    grid,
                    density,
                    color=theme.calibration_color("permutation"),
                    alpha=0.30,
                )
                ax_rand.axvline(
                    row.reference,
                    label=f"XICOR mean over tie breaks ({row.reference:.4f})",
                    **theme.REFERENCE_LINE,
                )
                ax_rand.axvline(
                    row.ours,
                    color=theme.verdict_color("pass"),
                    linewidth=2.0,
                    label=f"ours, deterministic ({row.ours:.4f})",
                )
                ax_rand.set_yticks([])
                ax_rand.set_xlabel("Chatterjee xi")
                ax_rand.set_title("Where the reference is a DISTRIBUTION, not a value")
                ax_rand.legend(loc="upper right", **theme.LEGEND_TEXT)
                ax_rand.text(
                    0.02,
                    0.78,
                    f"Chatterjee (2021) breaks x-ties at random, so XICOR is a random variable here\n"
                    f"(sd {spread:.4f}). Equality is not the right test; membership is.",
                    transform=ax_rand.transAxes,
                    **theme.ANNOTATION,
                )

            fig.suptitle("Cross-implementation agreement - R reference values")
            fig.text(
                0.5,
                0.015,
                "XICOR is by Holmes and Chatterjee; energy is by Rizzo and Szekely. "
                "randtests::bartels.rank.test is NOT compared: it is the rank von Neumann "
                "ratio, a different functional from our studentized max-over-lags C5. "
                "The lag-1 rank autocorrelation the two share is pinned by test instead.",
                ha="center",
                **theme.CAPTION,
                wrap=True,
            )
            fig.tight_layout(rect=(0, 0.08, 1, 0.94))
        return fig


class TieExperimentPlot(BasePlot):
    """What response ties cost Chatterjee's xi - the three questions of `jobs/bench/xi_ties.py`."""

    def build_matplotlib(self, result: object, style: str = "default") -> plt.Figure:
        pd_ = _check(result)
        frame = pd_.tie_experiment
        n_values = sorted(frame["n"].unique())

        with theme.style_context(style):
            fig, axes = plt.subplots(1, 3, figsize=(16, 4.6), facecolor="white")
            ax_p, ax_size, ax_break = axes

            # X AXIS IS THE NUMBER OF DISTINCT RESPONSE LEVELS, not the tie fraction.
            # Tie fraction saturates: quantising to 20 levels already ties 83-99% of a
            # sample, so on that axis every cell piles up against 1.0 and the finding -
            # that the closed form is fine down to about 5 levels and fails at 2 to 3 -
            # is invisible. Levels is also the variable a reader controls: it is how many
            # distinct values their metric takes.
            levels = sorted(
                (int(v) for v in frame["y_levels"].unique() if v > 0), reverse=True
            )
            positions = {lv: i for i, lv in enumerate(levels)}
            positions[0] = len(levels)  # the tie-free control, at the right
            tick_labels = [str(lv) for lv in levels] + ["continuous"]

            for i, n in enumerate(n_values):
                sub = frame[frame["n"] == n].copy()
                sub["pos"] = sub["y_levels"].map(positions)
                sub = sub.sort_values("pos")
                color = theme.ordered_color(i, len(n_values))
                x = sub["pos"]
                ax_p.plot(
                    x,
                    sub["p_signed_diff_mean"],
                    marker="o",
                    color=color,
                    label=f"n = {int(n)}",
                )
                ax_p.fill_between(
                    x,
                    sub["p_signed_diff_mean"] - sub["p_signed_diff_se"],
                    sub["p_signed_diff_mean"] + sub["p_signed_diff_se"],
                    color=color,
                    alpha=0.20,
                )
                ax_size.plot(
                    x,
                    sub["type_i_closed"],
                    marker="o",
                    color=color,
                    label=f"n = {int(n)} closed form",
                )
                ax_size.plot(
                    x,
                    sub["type_i_perm"],
                    marker="s",
                    linestyle=":",
                    color=color,
                    alpha=0.65,
                )
                ax_break.plot(
                    x,
                    sub["xi_tie_break_sd"],
                    marker="o",
                    color=color,
                    label=f"n = {int(n)}",
                )

            threshold = pd_.meta["divergence_threshold"]
            for axis in (ax_p, ax_size, ax_break):
                axis.set_xticks(range(len(tick_labels)))
                axis.set_xticklabels(tick_labels, **theme.LABEL_TEXT)
                axis.set_xlabel("Distinct values in the response")

            for sign in (-1.0, 1.0):
                ax_p.axhline(sign * threshold, **theme.REFERENCE_LINE)
            ax_p.axhspan(-pd_.mc_floor, pd_.mc_floor, **theme.HIGHLIGHT_SHADE)
            ax_p.set_ylabel("Closed form minus permutation (p)")
            ax_p.set_title("Do the two p-values disagree")
            ax_p.legend(**theme.LEGEND_TEXT)
            ax_p.text(
                0.30,
                0.06,
                f"shaded: the permutation's own Monte Carlo error (+/-{pd_.mc_floor:.4f}) -\n"
                f"the reason this reads the SIGNED mean and not |dp|.\n"
                f"dashed: the {threshold} divergence threshold.",
                transform=ax_p.transAxes,
                **theme.ANNOTATION,
            )

            ax_size.axhline(0.05, **theme.REFERENCE_LINE)
            ax_size.set_ylabel("Rejection rate under independence")
            ax_size.set_title("Does the closed-form null survive ties")
            ax_size.legend(**theme.LEGEND_TEXT)
            ax_size.text(
                0.02,
                0.88,
                "solid: closed form.  dotted: permutation (the reference).\n"
                "dashed: nominal 0.05.",
                transform=ax_size.transAxes,
                **theme.ANNOTATION,
            )

            ax_break.set_ylabel("Spread of xi over random tie breaks (sd)")
            ax_break.set_title("Cost of a deterministic x-tie break")
            ax_break.legend(**theme.LEGEND_TEXT)
            ax_break.text(
                0.02,
                0.88,
                "x is quantised here TOO, unlike the other two panels -\n"
                "otherwise there are no x-ties to break.",
                transform=ax_break.transAxes,
                **theme.ANNOTATION,
            )

            fig.suptitle("Chatterjee xi under response ties")
            crossings = ", ".join(
                f"n={n}: " + ("never" if not np.isfinite(v) else f"{int(v)} levels")
                for n, v in sorted(pd_.divergence_levels.items())
            )
            fig.text(
                0.5,
                0.015,
                f"The two p-values agree until the response is nearly binary: systematic "
                f"divergence beyond {threshold} appears only at {crossings}. Read off the "
                f"SIGNED mean - the absolute difference is dominated by the permutation's "
                f"own resampling noise. The cost of ties shows up as LEVEL (centre panel), "
                f"not as p-value disagreement.",
                ha="center",
                **theme.CAPTION,
                wrap=True,
            )
            fig.tight_layout(rect=(0, 0.07, 1, 0.93))
        return fig
