"""The four calibration figures: what the bench established, drawn.

Pure renderers. Every quantity comes from a typed artifact built by
`analyzers/calibration_summary.py`; nothing here derives a scientific value, and no colour
or font size is written down - both would fail `tests/test_style_baseline.py`, which pins
the count of hardcoded style literals across `panels/` and `plots/` and only lets it fall.

Each class draws ONE object and its title names that object without arguing about it, per
`docs/FIGURE_STANDARD.md`. Spaced hyphens throughout, units in axis labels in parentheses,
and where a figure rests on a subset of the grid it says so in the panel rather than in a
caption nobody will have.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from quebra.analyzers.calibration_summary import (
    PowerVsDependence,
    ReadDependence,
    SizeVsN,
    ValidationCurve,
)
from quebra.plots.base import BasePlot
from quebra.plots.theme import (
    REFERENCE_LINE,
    BAND_STYLE,
    calibration_color,
    ordered_color,
    style_context,
)


# Rows are labelled "check [calibration/variant]"; the calibration decides the colour.
def _calibration_of(row_label: str) -> str:
    inside = row_label.split("[", 1)[1].rstrip("]")
    return inside.split("/", 1)[0]


def _variant_of(row_label: str) -> str:
    inside = row_label.split("[", 1)[1].rstrip("]")
    return inside.split("/", 1)[1] if "/" in inside else ""


def _check_of(row_label: str) -> str:
    return row_label.split(" [", 1)[0]


def _short(check: str) -> str:
    return check.split("_", 1)[0].upper()


class CalibrationSizePlot(BasePlot):
    """Empirical size against event count, asymptotic against permutation."""

    def build_matplotlib(self, result: SizeVsN, style: str = "default") -> plt.Figure:
        if not isinstance(result, SizeVsN):
            raise TypeError("CalibrationSizePlot expects a SizeVsN artifact")
        with style_context(style):
            checks = sorted({_check_of(row) for row in result.rate_by_row})
            fig, axes = plt.subplots(
                1, len(checks), figsize=(4.1 * len(checks), 3.6), sharey=True
            )
            axes = np.atleast_1d(axes)
            n = np.asarray(result.n_values, dtype=float)

            for ax, check in zip(axes, checks):
                # The band a correctly calibrated check should sit inside: nominal, plus
                # or minus two Monte Carlo standard errors of the cell behind it.
                any_row = next(r for r in result.rate_by_row if _check_of(r) == check)
                se = np.asarray(result.se_by_row[any_row], dtype=float)
                ax.fill_between(
                    n,
                    result.alpha - 2 * se,
                    result.alpha + 2 * se,
                    color=REFERENCE_LINE["color"],
                    alpha=BAND_STYLE["fill_alpha"],
                    linewidth=0,
                )
                ax.axhline(result.alpha, **REFERENCE_LINE)

                for row in sorted(
                    r for r in result.rate_by_row if _check_of(r) == check
                ):
                    variant = _variant_of(row)
                    # The spread across the cells the mean pools, drawn so the reader can
                    # see how much averaging stands behind each point.
                    if row in result.rate_lo_by_row:
                        ax.fill_between(
                            n,
                            result.rate_lo_by_row[row],
                            result.rate_hi_by_row[row],
                            color=calibration_color(_calibration_of(row)),
                            alpha=BAND_STYLE["fill_alpha"],
                            linewidth=0,
                        )
                    ax.plot(
                        n,
                        result.rate_by_row[row],
                        marker="o",
                        markersize=3.5,
                        color=calibration_color(_calibration_of(row)),
                        linestyle="--" if variant == "unstudentized" else "-",
                        label=_calibration_of(row)
                        + (f" / {variant}" if variant else ""),
                    )
                ax.set_xscale("log")
                ax.set_xticks(result.n_values)
                ax.set_xticklabels([str(v) for v in result.n_values])
                ax.set_title(_short(check))
                ax.set_xlabel("Events in the record")
                ax.legend(loc="best", frameon=False)

            axes[0].set_ylabel("Empirical size")
            fig.suptitle(
                f"Empirical size against event count - {result.cell_description}"
            )
            # FIGURE_STANDARD: a panel that rests on a subset says so, in the panel.
            axes[-1].text(
                0.98,
                0.02,
                f"grey: nominal {result.alpha:g} +/- 2 MC SE\n"
                f"colour: spread over {result.n_cells_pooled} pooled cells",
                transform=axes[-1].transAxes,
                ha="right",
                va="bottom",
            )
            fig.tight_layout()
            return fig

    def build_plotly(self, result: SizeVsN):
        raise NotImplementedError("CalibrationSizePlot has no plotly backend")


class CalibrationPowerPlot(BasePlot):
    """Power against the duration dependence actually induced."""

    def build_matplotlib(
        self, result: PowerVsDependence, style: str = "default"
    ) -> plt.Figure:
        if not isinstance(result, PowerVsDependence):
            raise TypeError("CalibrationPowerPlot expects a PowerVsDependence artifact")
        with style_context(style):
            rows = list(result.rows_drawn) or sorted(result.power_by_row_and_n)
            fig, axes = plt.subplots(
                1, len(rows), figsize=(4.0 * len(rows), 3.9), sharey=True
            )
            axes = np.atleast_1d(axes)
            n_values = sorted(result.induced_lag1_by_n)

            for ax, row in zip(axes, rows):
                ax.axhline(result.alpha, **REFERENCE_LINE)
                for index, n in enumerate(n_values):
                    ax.plot(
                        result.induced_lag1_by_n[n],
                        result.power_by_row_and_n[row][n],
                        marker="o",
                        markersize=3.0,
                        color=ordered_color(index, len(n_values)),
                        label=f"n = {n}",
                    )
                ax.set_title(f"{_short(_check_of(row))} {_variant_of(row)}".strip())
                ax.set_xlabel("Induced duration lag-1 autocorrelation")
                ax.set_ylim(0.0, 1.0)
            axes[0].set_ylabel("Power")
            axes[-1].legend(loc="upper left", frameon=False)
            # FIGURE_STANDARD: a panel drawing a subset says which, in the panel.
            if result.selection_reason:
                fig.text(0.5, 0.005, f"shown: {result.selection_reason}", ha="center")
            fig.suptitle(
                "Power against induced duration dependence - Arm E, copula AR(1)"
            )
            fig.tight_layout(rect=(0, 0.045, 1, 1))
            return fig

    def build_plotly(self, result: PowerVsDependence):
        raise NotImplementedError("CalibrationPowerPlot has no plotly backend")


class CalibrationValidationPlot(BasePlot):
    """Eq (7) against its limiting null, with and without an estimated gamma.

    Eq (7) is Kvaloy and Lindqvist's Anderson-Darling statistic, Technometrics 62(1) 2020,
    p. 104.
    """

    def build_matplotlib(
        self, result: ValidationCurve, style: str = "default"
    ) -> plt.Figure:
        if not isinstance(result, ValidationCurve):
            raise TypeError(
                "CalibrationValidationPlot expects a ValidationCurve artifact"
            )
        with style_context(style):
            fig, axes = plt.subplots(
                1,
                len(result.n_values),
                figsize=(3.7 * len(result.n_values), 3.7),
                sharey=True,
            )
            axes = np.atleast_1d(axes)
            theoretical = np.asarray(result.theoretical, dtype=float)

            for ax, n in zip(axes, result.n_values):
                ax.plot(theoretical, theoretical, **REFERENCE_LINE)
                ax.plot(
                    theoretical,
                    result.empirical_gamma_one[n],
                    color=calibration_color("permutation"),
                    label="gamma = 1 (transcription)",
                )
                ax.plot(
                    theoretical,
                    result.empirical_gamma_hat[n],
                    color=calibration_color("asymptotic"),
                    label="estimated gamma (shipped)",
                )
                ax.set_title(f"n = {n}")
                ax.set_xlabel("Theoretical p-value")
                ax.set_xlim(0.0, 1.0)
                ax.set_ylim(0.0, 1.0)
                ax.set_aspect("equal", adjustable="box")
            axes[0].set_ylabel("Empirical p-value")
            axes[0].legend(loc="lower right", frameon=False)
            axes[-1].text(
                0.97,
                0.03,
                f"{result.n_replicates} replicates, seed {result.seed}",
                transform=axes[-1].transAxes,
                ha="right",
                va="bottom",
            )
            fig.suptitle(
                "Anderson-Darling p-values against the limiting null - "
                "Kvaloy and Lindqvist eq (7)"
            )
            fig.tight_layout()
            return fig

    def build_plotly(self, result: ValidationCurve):
        raise NotImplementedError("CalibrationValidationPlot has no plotly backend")


class CalibrationReadDependencePlot(BasePlot):
    """Duration dependence induced by read-level correlation."""

    def build_matplotlib(
        self, result: ReadDependence, style: str = "default"
    ) -> plt.Figure:
        if not isinstance(result, ReadDependence):
            raise TypeError(
                "CalibrationReadDependencePlot expects a ReadDependence artifact"
            )
        with style_context(style):
            fig, ax = plt.subplots(figsize=(6.4, 4.2))
            rho = np.asarray(result.read_rho, dtype=float)
            induced = np.asarray(result.induced_lag1, dtype=float)
            se = np.asarray(result.induced_lag1_se, dtype=float)

            ax.axhline(0.0, **REFERENCE_LINE)
            ax.plot(rho, rho, **REFERENCE_LINE)
            ax.errorbar(
                rho,
                induced,
                yerr=se,
                marker="o",
                markersize=4.0,
                capsize=2.5,
                color=calibration_color("permutation"),
                label="induced, measured",
            )
            ax.annotate(
                "if read correlation reached the durations,\nthe curve would follow this line",
                xy=(rho[len(rho) // 2], rho[len(rho) // 2]),
                xytext=(0.60, 0.30),
                textcoords="axes fraction",
                ha="center",
                va="top",
                arrowprops={
                    "arrowstyle": "->",
                    "color": REFERENCE_LINE["color"],
                    "linewidth": REFERENCE_LINE["linewidth"],
                },
            )
            ax.set_xlabel("Read-level AR(1) correlation")
            ax.set_ylabel("Induced duration lag-1 autocorrelation")
            ax.set_title(
                "Duration dependence induced by read correlation - Arm C, carved"
            )
            ax.legend(loc="upper left", frameon=False)
            # Headroom so the y = x reference and its note clear the top of the axes.
            ax.set_ylim(min(-0.06, float(induced.min()) - 0.03), max(rho) * 1.18)
            fig.tight_layout()
            return fig

    def build_plotly(self, result: ReadDependence):
        raise NotImplementedError("CalibrationReadDependencePlot has no plotly backend")


__all__ = [
    "CalibrationPowerPlot",
    "CalibrationReadDependencePlot",
    "CalibrationSizePlot",
    "CalibrationValidationPlot",
]
