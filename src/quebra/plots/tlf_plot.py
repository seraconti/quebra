from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go

from quebra.plots.base import BasePlot
from quebra.plots.theme import mix_with_white, qubit_color, style_context


class TLFPlot(BasePlot):
    def build_matplotlib(
        self, payload: dict[str, object], style: str = "default"
    ) -> plt.Figure:
        with style_context(style):
            result = payload["result"]
            values_hz = np.asarray(payload["values_hz"], dtype=float)
            meta = dict(payload.get("meta", {}))

            values_hz = values_hz[np.isfinite(values_hz)]
            if len(values_hz) == 0:
                raise ValueError("TLFPlot requires finite values_hz")

            # All x-units in this figure are kHz.
            values_khz = values_hz / 1e3

            color_main = qubit_color(
                dataset_id=str(meta.get("dataset_id", "")), meta=meta
            )
            color_alt = mix_with_white(color_main, amount=0.35)

            fig, (ax, ax_info) = plt.subplots(
                1,
                2,
                figsize=(12, 5),
                gridspec_kw={"width_ratios": [2, 1], "wspace": 0.18},
            )
            fig.patch.set_facecolor("white")
            ax.set_facecolor("white")

            n = len(values_khz)
            x_min = float(np.min(values_khz))
            x_max = float(np.max(values_khz))
            spread = x_max - x_min
            std_khz = float(np.std(values_khz, ddof=1)) if n > 1 else 0.0
            if n > 1 and std_khz > 0.0:
                bin_width = 3.5 * std_khz * (n ** (-1.0 / 3.0))
            else:
                bin_width = 0.0

            if spread <= 0.0:
                bins = 10
            elif bin_width > 0.0:
                bins = max(10, int(spread / bin_width))
            else:
                bins = 10

            density, bins, _ = ax.hist(
                values_khz,
                bins=bins,
                density=True,
                alpha=0.28,
                color=color_main,
                edgecolor="none",
                label="Data histogram",
            )

            if not np.isfinite(x_min) or not np.isfinite(x_max) or x_max <= x_min:
                x_min, x_max = -1.0, 1.0
            x_grid_khz = np.linspace(x_min, x_max, 800)

            means1 = np.asarray(result.gmm1.means_).ravel()
            covs1 = np.asarray(result.gmm1.covariances_).ravel()
            weights1 = np.asarray(result.gmm1.weights_).ravel()
            means1_axis = means1 / 1e3
            covs1_axis = covs1 / (1e3**2)
            y1 = np.zeros_like(x_grid_khz, dtype=float)
            for w, m, c in zip(weights1, means1_axis, covs1_axis):
                sigma = max(float(np.sqrt(c)), 1e-15)
                y1 += (
                    float(w)
                    * (1.0 / (sigma * np.sqrt(2.0 * np.pi)))
                    * np.exp(-0.5 * ((x_grid_khz - float(m)) / sigma) ** 2)
                )

            # No GMM(2) curve when the 2-component fit failed: a single-lobe curve under a
            # "GMM(2)" label reads as evidence against bimodality, not as an absent fit.
            fit_failed = bool(result.fit_failed) or result.gmm2 is None
            if fit_failed:
                means2_axis = np.asarray([], dtype=float)
                y2 = np.zeros_like(x_grid_khz, dtype=float)
            else:
                means2 = np.asarray(result.gmm2.means_).ravel()
                covs2 = np.asarray(result.gmm2.covariances_).ravel()
                weights2 = np.asarray(result.gmm2.weights_).ravel()
                means2_axis = means2 / 1e3
                covs2_axis = covs2 / (1e3**2)
                y2 = np.zeros_like(x_grid_khz, dtype=float)
                for w, m, c in zip(weights2, means2_axis, covs2_axis):
                    sigma = max(float(np.sqrt(c)), 1e-15)
                    y2 += (
                        float(w)
                        * (1.0 / (sigma * np.sqrt(2.0 * np.pi)))
                        * np.exp(-0.5 * ((x_grid_khz - float(m)) / sigma) ** 2)
                    )

            ax.plot(
                x_grid_khz, y1, "--", color=color_alt, linewidth=1.3, label="GMM(1)"
            )
            if not fit_failed:
                ax.plot(
                    x_grid_khz, y2, "-", color=color_main, linewidth=1.7, label="GMM(2)"
                )

            if len(means2_axis) >= 2:
                sorted_means = np.sort(means2_axis)
                y_marker = 0.92 * max(
                    float(np.max(y2)),
                    float(np.max(density)) if len(density) > 0 else 1.0,
                )
                left = float(sorted_means[0])
                right = float(sorted_means[-1])
                ax.annotate(
                    "",
                    xy=(right, y_marker),
                    xytext=(left, y_marker),
                    arrowprops={"arrowstyle": "<->", "color": "0.25", "lw": 1.0},
                )
                if result.lobe_separation_hz is not None:
                    ax.text(
                        0.5 * (left + right),
                        y_marker * 1.03,
                        f"lobe distance = {result.lobe_separation_hz / 1e3:.3f} kHz",
                        ha="center",
                        va="bottom",
                        fontsize=8,
                        color="0.25",
                    )

            ax_info.set_facecolor("white")
            ax_info.axis("off")
            lobe_khz = (
                float(result.lobe_separation_hz) / 1e3
                if result.lobe_separation_hz is not None
                else np.nan
            )
            sigma_khz = (
                float(result.within_lobe_spread_hz) / 1e3
                if result.within_lobe_spread_hz is not None
                else np.nan
            )
            if len(means2_axis) >= 2:
                order2 = np.argsort(means2_axis)
                means2_sorted = means2_axis[order2]
                stds2_sorted = np.sqrt(np.clip(covs2_axis[order2], 0.0, None))
                weights2_sorted = weights2[order2]
                mean1_khz = float(means2_sorted[0])
                std1_khz = float(stds2_sorted[0])
                w1 = float(weights2_sorted[0])
                mean2_khz = float(means2_sorted[-1])
                std2_khz = float(stds2_sorted[-1])
                w2 = float(weights2_sorted[-1])
                dominant_side = "left" if w1 > w2 else "right"
                dominant_pct = max(w1, w2) * 100.0
                lobe_lines = [
                    f"Lobe 1:  {mean1_khz:.1f} kHz   σ={std1_khz:.1f} kHz",
                    f"  occupation: {w1 * 100:.1f}%",
                    "",
                    f"Lobe 2:  {mean2_khz:.1f} kHz   σ={std2_khz:.1f} kHz",
                    f"  occupation: {w2 * 100:.1f}%",
                ]
            else:
                dominant_side = "n/a"
                dominant_pct = np.nan
                lobe_lines = [
                    "Lobe 1:  n/a",
                    "  occupation: n/a",
                    "",
                    "Lobe 2:  n/a",
                    "  occupation: n/a",
                ]
            # Build a more detailed textual summary including histogram metrics and dynamics (if available)
            bic_norm = (
                float(result.normalized_bic_delta)
                if getattr(result, "normalized_bic_delta", None) is not None
                else None
            )
            lobe_snr = (
                float(result.lobe_snr)
                if getattr(result, "lobe_snr", None) is not None
                else None
            )

            # Group 2: dynamics, within segments. Every count printed beside its statistic.
            n_trans = result.n_transitions
            sw_rate = result.switching_rate_per_hour
            dynamics_lines = [
                "",
                "Dynamics (within segments):",
                f"  n_transitions: {int(n_trans)}"
                if n_trans is not None
                else "  n_transitions: n/a",
                f"  switching_rate_per_hour: {sw_rate:.3f} /h"
                if sw_rate is not None
                else "  switching_rate_per_hour: n/a",
            ]
            for k in (0, 1):
                dwell = getattr(result, f"mean_dwell_s{k}")
                cv = getattr(result, f"dwell_cv_s{k}")
                exits = getattr(result, f"n_exits_s{k}")
                complete = getattr(result, f"n_complete_runs_s{k}")
                cut = getattr(result, f"n_cut_runs_s{k}")
                dynamics_lines.append(
                    f"  mean_dwell_s{k}: {dwell:.2f} s ({exits} exits)"
                    if dwell is not None
                    else f"  mean_dwell_s{k}: n/a"
                )
                dynamics_lines.append(
                    f"  dwell_cv_s{k}: {cv:.3f} ({complete} complete, {cut} cut runs)"
                    if cv is not None
                    else f"  dwell_cv_s{k}: n/a ({complete} complete, {cut} cut runs)"
                )

            n_used = result.diagnostics.get("n_reads")
            n_dropped = result.diagnostics.get("n_dropped_nonfinite")
            summary_lines = [
                "TLF Summary",
                f"reads: {n_used} used, {n_dropped} dropped (non-finite time or value)"
                if n_used is not None and n_dropped is not None
                else "reads: n/a",
                "",
                f"is_bimodal: {bool(result.is_bimodal)}"
                if not fit_failed
                else "is_bimodal: n/a (2-component fit failed)",
                "  True if BIC delta > 6",
                "",
                f"bic_delta: {float(result.bic_delta):.2f}"
                if result.bic_delta is not None
                else "bic_delta: n/a",
                "  BIC1 - BIC2",
                f"normalized_bic_delta: {bic_norm:.4f}"
                if bic_norm is not None
                else "normalized_bic_delta: n/a",
                "  bic_delta / N",
                "",
                f"lobe distance (kHz): {lobe_khz:.3f}"
                if np.isfinite(lobe_khz)
                else "lobe distance (kHz): n/a",
                "  |mu2 - mu1| from GMM(2)",
                "",
                *lobe_lines,
                "",
                f"within_lobe_spread: {sigma_khz:.3f} kHz"
                if np.isfinite(sigma_khz)
                else "within_lobe_spread: n/a",
                "  weighted Gaussian width",
                f"lobe_snr: {lobe_snr:.3f}"
                if lobe_snr is not None
                else "lobe_snr: n/a",
                f"Dominant lobe: {dominant_side}  ({dominant_pct:.1f}% occupied)"
                if np.isfinite(dominant_pct)
                else "Dominant lobe: n/a",
                *dynamics_lines,
            ]
            ax_info.text(
                0.02, 0.98, "\n".join(summary_lines), ha="left", va="top", fontsize=9
            )

            dataset_id = str(meta.get("dataset_id", ""))
            ax.set_title(f"TLF distribution and Gaussian fits - {dataset_id}")
            ax.set_xlabel("Delta frequency (kHz)")
            ax.set_ylabel("Normalized counts")
            ax.grid(True, alpha=0.25)
            ax.legend(frameon=False)
            fig.subplots_adjust(
                left=0.07, right=0.98, top=0.90, bottom=0.14, wspace=0.18
            )
            return fig

    def build_plotly(self, payload: dict[str, object]) -> go.Figure:
        raise NotImplementedError(f"{self.__class__.__name__} has no plotly backend")
