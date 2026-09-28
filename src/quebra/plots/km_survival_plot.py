"""Single-axes Kaplan-Meier survival figure: two datasets, one threshold, nothing else.

A pure renderer. Every number it draws is already on the `KaplanMeierComparison` the step
produced - including WHICH two curves to draw, which is `comparison.pair`, ranked by
`log_time_separation` in the analyzer. Selecting a pair at draw time would make the figure
a function of something the provenance record does not contain.

The x-axis is logarithmic in window age, because the two records this figure exists to
contrast differ by more than a decade in median lifetime and a linear axis collapses the
shorter-lived one onto the origin. Consequences the code has to handle, not hide:

- S = 1 holds from age 0, and 0 has no place on a log axis. The flat head of each curve is
  drawn from the left limit, one full decade below the earliest step in either curve, so
  the eye reads the head as "before anything happened" rather than as a first event.
- A window of zero recorded age (birth and death share a timestamp, which the scan-
  structured record does produce) steps S down at t = 0 and cannot be placed. Those steps
  are clamped to the left limit.

DEVIATION FROM docs/FIGURE_STANDARD.md, on the author's instruction: this figure carries
no in-panel statement of what it excluded, and no caption. The standard requires both,
because a reader who cannot see an exclusion cannot judge it - and this estimate DOES
exclude windows whose birth was not observed. Those counts are on the materialized
`KaplanMeierCurve` and in the provenance record, and the legend's n is the number of
windows actually estimated from, so the exclusion is recoverable but is not visible
to someone holding only the image. Restoring it means putting `_exclusion_note` back.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter, LogLocator

from quebra.analyzers.kaplan_meier import KaplanMeierComparison, KaplanMeierCurve
from quebra.plots import theme
from quebra.plots.base import BasePlot

# The flat S = 1 head must be visibly a head rather than a plateau between two events, so
# the axis starts at the decade below the earliest step - unless that step sits almost on
# the decade boundary, which would leave no head at all, in which case it drops one more.
# Anchoring on decades rather than on the step itself keeps the axis limits stable when
# the data shifts slightly, and keeps the leftmost tick a labelled decade.
_MIN_HEAD_RATIO = 1.5

# One ink for the whole figure; the curves separate by line style alone. Solid reads as
# the reference and dotted as the comparison, and both survive greyscale and distance -
# which a two-hue pair does not.
_LINESTYLES = ("-", "--", ":", "-.")


def _plain_minutes(value: float, _pos: int) -> str:
    """Tick labels as plain numbers, never 10^n, and never carrying the unit.

    FIGURE_STANDARD puts the unit in the axis label alone. A log axis whose ticks read
    0.1 / 1 / 10 / 100 is also read faster on a poster than one reading 10^-1 / 10^0.
    """
    if value <= 0:
        return ""
    if value >= 1:
        return f"{value:.0f}"
    return f"{value:g}"


def _draw_arrays(
    curve: KaplanMeierCurve, x_lo: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Step arrays clamped into the drawable range and closed at the last observation.

    The final segment is extended to the curve's longest observed window age: the step
    function holds its last value there, and stopping the line at the last DEATH would
    leave the censored tail undrawn.
    """
    x = np.maximum(curve.time_min, x_lo)
    s, lo, hi = curve.survival, curve.band_lower, curve.band_upper
    x_end = max(curve.max_observed_min, float(x[-1]))
    if x_end > x[-1]:
        x = np.append(x, x_end)
        s = np.append(s, s[-1])
        lo = np.append(lo, lo[-1])
        hi = np.append(hi, hi[-1])
    return x, s, lo, hi


class KMSurvivalPlot(BasePlot):
    """Kaplan-Meier survival of in-spec windows for the comparison's selected pair."""

    def build_matplotlib(
        self, comparison: KaplanMeierComparison, style: str = "default"
    ) -> plt.Figure:
        if not isinstance(comparison, KaplanMeierComparison):
            raise TypeError(
                "KMSurvivalPlot draws a KaplanMeierComparison, not "
                f"{type(comparison).__name__}. The pair to draw is chosen in the step, "
                "not here."
            )
        curves = list(comparison.pair_curves())

        with theme.style_context(style):
            first_steps = [
                float(np.min(c.time_min[c.time_min > 0.0]))
                for c in curves
                if np.any(c.time_min > 0.0)
            ]
            if not first_steps:
                raise ValueError(
                    "every window in both curves has zero recorded age - there is no "
                    "log time axis to draw them on"
                )
            earliest = min(first_steps)
            decade = 10.0 ** np.floor(np.log10(earliest))
            x_lo = decade if earliest / decade >= _MIN_HEAD_RATIO else decade / 10.0
            x_hi = 10.0 ** np.ceil(np.log10(max(c.max_observed_min for c in curves)))

            fig, ax = plt.subplots(figsize=(12.0, 7.5))
            fig.patch.set_facecolor("white")
            ax.set_facecolor("white")

            ax.axhline(0.5, **theme.REFERENCE_LINE)

            color = theme.POSTER_INK_SURVIVAL
            for index, curve in enumerate(curves):
                x, s, lo, hi = _draw_arrays(curve, x_lo)
                linestyle = _LINESTYLES[index % len(_LINESTYLES)]
                ax.fill_between(
                    x,
                    lo,
                    hi,
                    step="post",
                    color=color,
                    alpha=theme.BAND_STYLE["fill_alpha"],
                    linewidth=0.0,
                    zorder=2,
                )
                ax.step(
                    x,
                    s,
                    where="post",
                    color=color,
                    linestyle=linestyle,
                    zorder=3,
                    label=_legend_label(curve),
                )
                if len(curve.censor_time_min):
                    ax.plot(
                        np.maximum(curve.censor_time_min, x_lo),
                        curve.censor_survival,
                        linestyle="none",
                        color=color,
                        zorder=4,
                        **theme.CENSOR_TICK,
                    )

            ax.set_xscale("log")
            ax.set_xlim(x_lo, x_hi)
            ax.set_ylim(-0.02, 1.04)
            ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
            ax.xaxis.set_major_locator(LogLocator(base=10.0))
            ax.xaxis.set_major_formatter(FuncFormatter(_plain_minutes))
            ax.set_xlabel("Window age (min)")
            ax.set_ylabel("Surviving fraction of windows")
            ax.set_title(_title(comparison))
            ax.grid(True, which="major")
            ax.grid(True, which="minor", alpha=theme.BAND_STYLE["fill_alpha"])
            # The legend is now the only text in the axes, so it takes the conventional
            # corner. Upper right is the one that is empty for a well-separated pair: the
            # flat heads run along the top LEFT, and the longer-lived curve has fallen
            # well below the top by the time it reaches the right-hand decades.
            ax.legend(loc="upper right", frameon=False)
            fig.tight_layout()
        return fig


def _legend_label(curve: KaplanMeierCurve) -> str:
    """Name, median, and the number of windows the estimate rests on.

    `n` is the only surviving trace in the image of how much data is behind each curve -
    see the deviation note in the module docstring - so it stays even though the legend is
    otherwise kept short.
    """
    median = curve.median_survival_min
    median_text = "not reached" if median is None else f"{median:.1f} min"
    return f"{curve.label}, n = {curve.n_windows}, median {median_text}"


def _title(comparison: KaplanMeierComparison) -> str:
    return f"In-spec window survival, T2* >= {comparison.threshold_label}"
