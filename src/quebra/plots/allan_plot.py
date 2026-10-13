from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go

from quebra.analyzers.allan import AllanResult
from quebra.plots.base import BasePlot
from quebra.plots.theme import CAPTION, mix_with_white, qubit_color, style_context


class AllanPlot(BasePlot):
    def build_matplotlib(
        self, result: AllanResult, style: str = "default"
    ) -> plt.Figure:
        with style_context(style):
            mode_names = list(result.modes.keys())
            if len(mode_names) == 0:
                raise ValueError("AllanResult contains no modes to plot")

            # Handle fractional_adev robustly: it may be None, a scalar, or an array-like.
            frac_raw = getattr(result, "fractional_adev", None)
            if frac_raw is None:
                frac = np.array([], dtype=float)
            else:
                frac = np.atleast_1d(np.asarray(frac_raw, dtype=float))
            has_fractional = frac.size > 0

            if has_fractional:
                fig, axes = plt.subplots(
                    3,
                    1,
                    figsize=(10, 8.0),
                    sharex=True,
                    gridspec_kw={"height_ratios": [2.3, 2.1, 1.5], "hspace": 0.28},
                )
                ax_main, ax_frac, ax_slope = axes
            else:
                fig, axes = plt.subplots(
                    2,
                    1,
                    figsize=(10, 6),
                    sharex=True,
                    gridspec_kw={"height_ratios": [3.0, 1.25], "hspace": 0.0},
                )
                ax_main, ax_slope = axes
                ax_frac = None

            # force white backgrounds regardless of active matplotlib style
            fig.patch.set_facecolor("white")
            ax_main.set_facecolor("white")
            ax_slope.set_facecolor("white")
            if ax_frac is not None:
                ax_frac.set_facecolor("white")

            meta = dict(result.meta)
            base_color = qubit_color(dataset_id=str(meta.get("dataset_id")), meta=meta)
            accent_color = mix_with_white(base_color, amount=0.25)

            first_mode = None
            for index, mode_name in enumerate(mode_names):
                frame = result.modes[mode_name]
                tau_s = frame["tau_s"].to_numpy(dtype=float)
                adev = frame["adev_hz"].to_numpy(dtype=float)
                color = base_color if index == 0 else accent_color
                ax_main.loglog(
                    tau_s,
                    adev,
                    "o-",
                    markersize=3.5,
                    linewidth=1.2,
                    color=color,
                    label=str(mode_name),
                )

                if first_mode is None and len(tau_s) > 0:
                    first_mode = (tau_s, adev)

            if first_mode is not None:
                tau_s, adev = first_mode
                tau_mask = (
                    np.isfinite(tau_s)
                    & np.isfinite(adev)
                    & (tau_s > 0.0)
                    & (adev > 0.0)
                )
                tau = tau_s[tau_mask]
                y = adev[tau_mask]
                if len(tau) > 0:
                    tau0 = float(tau[0])
                    y0 = float(y[0])
                    # Guide slopes of log ADEV on log tau, not spectral exponents.
                    for slope_guide in [0.5, 1.0, -0.5, -1.0]:
                        guide = y0 * (tau / tau0) ** slope_guide
                        ax_main.loglog(
                            tau,
                            guide,
                            "--",
                            color="gray",
                            linewidth=0.8,
                            alpha=0.3,
                            zorder=0,
                        )

                    max_tau = float(np.max(tau))
                    ax_main.axvline(max_tau, color="0.35", linestyle=":", linewidth=0.9)
                    ax_slope.axvline(
                        max_tau, color="0.35", linestyle=":", linewidth=0.9
                    )

                    y_min, y_max = ax_main.get_ylim()
                    label_y = (
                        y_min * (y_max / y_min) ** 0.96
                        if y_min > 0 and y_max > y_min
                        else y_max
                    )
                    ax_main.text(
                        max_tau,
                        label_y,
                        f"tau={max_tau:.3g}s",
                        rotation=90,
                        va="top",
                        ha="left",
                        fontsize="x-small",
                    )
                    ax_slope.text(
                        max_tau,
                        1.95,
                        f"tau={max_tau:.3g}s",
                        rotation=90,
                        va="top",
                        ha="left",
                        fontsize="x-small",
                    )

                    if len(tau) > 1:
                        tau_mid = np.sqrt(tau[:-1] * tau[1:])
                        slope = np.diff(np.log(y)) / np.diff(np.log(tau))
                        ax_slope.semilogx(
                            tau_mid, slope, "-", linewidth=1.3, color=base_color
                        )

                    if has_fractional and ax_frac is not None:
                        n = min(len(tau_s), len(frac))
                        tau_f = np.asarray(tau_s[:n], dtype=float)
                        frac_f = np.asarray(frac[:n], dtype=float)
                        frac_mask = (
                            np.isfinite(tau_f)
                            & np.isfinite(frac_f)
                            & (tau_f > 0.0)
                            & (frac_f > 0.0)
                        )
                        n_frac_dropped = int(np.sum(~frac_mask))
                        if np.any(frac_mask):
                            ax_frac.loglog(
                                tau_f[frac_mask],
                                frac_f[frac_mask],
                                "o--",
                                markersize=3.0,
                                linewidth=1.1,
                                color=base_color,
                                alpha=0.9,
                                label="Fractional ADEV",
                            )

                            nu0 = result.nu0_hz
                            if nu0 is not None and np.isfinite(float(nu0)):
                                ax_frac.annotate(
                                    f"normalization: ADEV / nu0 = {nu0 / 1e9:.5f} GHz"
                                    f"; dropped: {n_frac_dropped} non-positive or "
                                    f"non-finite points",
                                    xy=(0.98, 0.98),
                                    xycoords="axes fraction",
                                    ha="right",
                                    va="top",
                                    fontsize=7,
                                    color="gray",
                                )
                            ax_frac.set_title("Fractional Allan")
                            ax_frac.set_ylabel("Fractional ADEV")
                            ax_frac.grid(
                                True, which="major", color="lightgray", alpha=0.4
                            )
                            ax_frac.legend(frameon=False)

                tau0_s = meta.get("tau0_s")
                if tau0_s is not None and np.isfinite(float(tau0_s)):
                    title = f"Allan - {meta.get('dataset_id')} | tau0 = {float(tau0_s):.4g} s"
                else:
                    title = f"Allan - {meta.get('dataset_id')}"
                ax_main.set_title(title)
                diag = result.diagnostics
                ax_main.annotate(
                    f"{diag['n_segments']} segments, {diag['n_segments_too_short']} too "
                    f"short; spacing CV within segments "
                    f"{diag['spacing_cv_within_segments']:.2f} (tau is nominal when it is "
                    f"large); a tau is kept with >= {diag['min_terms_per_tau']} pooled terms"
                    f"\ndropped: {diag['n_dropped_nonfinite']} non-finite reads, "
                    f"{int(np.sum(~tau_mask))} non-positive or non-finite ADEV points",
                    xy=(0.02, 0.02),
                    xycoords="axes fraction",
                    ha="left",
                    va="bottom",
                    **CAPTION,
                )

            ax_main.set_ylabel("ADEV (Hz)")
            ax_main.grid(True, which="major", color="lightgray", alpha=0.4)
            ax_main.legend(frameon=False)

            # Keep the same x-scale across panels but redraw x labels/ticks on each panel for readability.
            ax_main.set_xlabel("Tau (s)")
            ax_main.tick_params(axis="x", which="both", labelbottom=True)
            if ax_frac is not None:
                ax_frac.set_xlabel("Tau (s)")
                ax_frac.tick_params(axis="x", which="both", labelbottom=True)

            ax_slope.set_xlabel("Tau (s)")
            ax_slope.set_ylabel("d log ADEV / d log tau")
            ax_slope.set_title(
                "Local slope of log ADEV on log tau (mu / 2, not the spectral alpha)"
            )
            ax_slope.axhline(0.0, color="gray", alpha=0.4, linewidth=0.8)
            ax_slope.grid(True, which="major", color="lightgray", alpha=0.4)
            if has_fractional:
                fig.subplots_adjust(
                    left=0.10, right=0.98, top=0.94, bottom=0.08, hspace=0.30
                )
            else:
                fig.tight_layout()
            return fig

    def build_plotly(self, result: AllanResult) -> go.Figure:
        raise NotImplementedError(f"{self.__class__.__name__} has no plotly backend")
