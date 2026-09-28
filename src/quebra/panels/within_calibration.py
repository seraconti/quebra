"""Generic within-calibration panel: degradation view with optional cumulative metrics.

WithinCalibrationPanelData (the contract) lives in _within_calibration_data and is re-exported
here, which stays the public import surface. Pure view arithmetic lives in
_within_calibration_render.

Artifacts materialized before the band split name
`panels.non_repairable.NonRepairablePanelData` in their pickle stream, and `output/` is
append-only.

`panels.non_repairable` does not exist, so those artifacts raise ModuleNotFoundError and no
re-export here can reach them: the module they name is gone, not the symbol. They are
recovered by re-running their jobs, not by editing this file.
`tests/test_artifact_guard.py` pins that outcome rather than hiding it.

Two halves:
  - WithinCalibrationPanelData is the COMPLETE typed artifact, composed of three bands
    (signal, distinguish, reliability) plus meta. Built solely by
    analyzers.within_calibration_compute.build_within_calibration_panel_data.
  - WithinCalibrationPanel is a PURE renderer: it reads fields and draws. It performs no
    data arithmetic; only axis/theme concerns (decade-guide ticks, bin geometry,
    colors) live here.

Accepts any monotonic or time-varying metric. No fidelity-specific logic lives here -
fidelity adaptation is in analyzers/fidelity.py::make_panel_data.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go

from quebra.panels import _within_calibration_render as render
from quebra.analyzers.reliability_band import estimator_name
from quebra.analyzers.windows import STATE_UNOBSERVED
from quebra.analyzers.within_calibration_data import WithinCalibrationPanelData
from quebra.plots import theme
from quebra.plots.base import BasePlot


# Per-threshold colours come from plots.theme.threshold_color: used on primary-axis
# dashed lines, cumulative subplots, and survival curves so all three can be visually
# correlated. Sampling a colormap across len(thresholds) replaces a fixed 8-entry list
# that wrapped modulo its length - a 10-entry ladder drew 1 us and 9 us in one blue.


class WithinCalibrationPanel(BasePlot):
    """Degradation analysis panel for the within-calibration tier.

    The reliability literature would call this the non-repairable case. See
    `panels/across_calibration.py` for why this project names its tiers after the
    calibration boundary instead, and where the literature's own terms are kept.

    Pure renderer over a complete WithinCalibrationPanelData. Three bands as rows:
      band 1  signal       series + value, fit-error and relative-error distributions
      band 2  distinguish  4-state read timeline, excursion shape, xi by threshold
      band 3  reliability  2-state compliance timeline, survival curve
    then cumulative time / damage [opt-in] and the summary text.

    The two timelines share an x-axis and sit vertically adjacent on purpose: if they
    look the same, resolvability costs nothing; if they differ, that is the finding.
    They are never merged.
    """

    def build_matplotlib(
        self, result: WithinCalibrationPanelData, style: str = "default"
    ) -> plt.Figure:
        if not isinstance(result, WithinCalibrationPanelData):
            raise TypeError("WithinCalibrationPanel expects WithinCalibrationPanelData")
        pd_ = result
        has_cum_time = pd_.include_cumulative_time
        has_cum_dmg = pd_.include_cumulative_damage
        has_extra_row = has_cum_time or has_cum_dmg

        with theme.style_context(style):
            n_thr = len(pd_.thresholds)
            # Both timelines get the SAME height so they read as a pair; a reader
            # compares them vertically, so anything that scales one must scale both.
            timeline_h = max(1.4, 0.30 * n_thr)
            rows = [
                ("signal", 3.6),
                ("distinguish_timeline", timeline_h),
                ("distinguish_detail", 4.4),
                ("reliability_timeline", timeline_h),
                ("survival", 3.2),
            ]
            if has_extra_row:
                rows.append(("cumulative", 2.4))
            rows.append(("summary", max(1.6, 0.28 * n_thr)))

            fig = plt.figure(
                figsize=(16, sum(h for _, h in rows)),
                constrained_layout=True,
                facecolor="white",
            )
            gs = fig.add_gridspec(len(rows), 2, height_ratios=[h for _, h in rows])
            row = {name: i for i, (name, _) in enumerate(rows)}

            # --- band 1: the signal and its distributions -----------------------
            # Column 2 is an empty spacer holding the threshold legend, so the signal
            # and its marginal sit flush (wspace 0.02) while the error histograms keep
            # their own breathing room.
            gs_signal = gs[row["signal"], :].subgridspec(
                1, 5, width_ratios=[6.4, 1.05, 0.5, 1.95, 1.95], wspace=0.03
            )
            ax_primary = fig.add_subplot(gs_signal[0, 0])
            ax_value_hist = fig.add_subplot(gs_signal[0, 1], sharey=ax_primary)
            ax_sigma_hist = fig.add_subplot(gs_signal[0, 3])
            ax_relerr_hist = fig.add_subplot(gs_signal[0, 4])

            # --- band 2: is the failure distinguishable -------------------------
            ax_state_timeline = fig.add_subplot(gs[row["distinguish_timeline"], :])
            # The two heatmaps STACK rather than sit side by side: ten threshold
            # columns each need the full half-width to stay legible, and stacking also
            # aligns them column-for-column so xi and rho are read together.
            gs_detail = gs[row["distinguish_detail"], :].subgridspec(
                1, 2, width_ratios=[1.0, 1.35], wspace=0.18
            )
            ax_shape = fig.add_subplot(gs_detail[0, 0])
            gs_stat = gs_detail[0, 1].subgridspec(
                2, 1, height_ratios=[1.0, 1.0], hspace=0.55
            )
            ax_xi = fig.add_subplot(gs_stat[0, 0])
            ax_rho = fig.add_subplot(gs_stat[1, 0])

            # --- band 3: reliability --------------------------------------------
            ax_compliance = fig.add_subplot(
                gs[row["reliability_timeline"], :], sharex=ax_state_timeline
            )
            ax_surv = fig.add_subplot(gs[row["survival"], :])

            axes = [
                ax_primary,
                ax_value_hist,
                ax_sigma_hist,
                ax_relerr_hist,
                ax_state_timeline,
                ax_shape,
                ax_xi,
                ax_rho,
                ax_compliance,
                ax_surv,
            ]
            for ax in axes:
                theme.apply_common_style(ax)

            if has_extra_row:
                if has_cum_time and has_cum_dmg:
                    ax_cum_time = fig.add_subplot(gs[row["cumulative"], 0])
                    ax_cum_dmg = fig.add_subplot(gs[row["cumulative"], 1])
                elif has_cum_time:
                    ax_cum_time = fig.add_subplot(gs[row["cumulative"], :])
                    ax_cum_dmg = None
                else:
                    ax_cum_time = None
                    ax_cum_dmg = fig.add_subplot(gs[row["cumulative"], :])
            else:
                ax_cum_time = None
                ax_cum_dmg = None
            ax_sum = fig.add_subplot(gs[row["summary"], :])

            color = pd_.color if pd_.color is not None else "C0"

            self._draw_primary(ax_primary, pd_, color, legend_ax=ax_value_hist)
            self._draw_primary_hist(ax_value_hist, pd_, color)
            self._draw_sigma_hist(ax_sigma_hist, pd_, color)
            self._draw_relerr_hist(ax_relerr_hist, pd_, color)

            self._draw_threshold_timeline(
                ax_state_timeline,
                pd_,
                pd_.distinguish.state_series_per_threshold,
                title="Read state by threshold (4-state, with uncertainty)",
                show_ttf=False,
                show_key=True,
            )
            self._draw_shape(ax_shape, pd_)
            # xi and Spearman side by side, never one alone: xi catches non-monotone
            # structure Spearman is blind to, Spearman is stronger against smooth
            # monotone dependence and carries a sign.
            self._draw_stat_heatmap(
                ax_xi,
                pd_,
                median_key="xi_m_to_f",
                pooled_lookup=lambda label: pd_.distinguish.pooled_xi_per_threshold.get(
                    label, (float("nan"),)
                )[0],
                title="Chatterjee xi, 0 = no dependence  (w = windows behind the cell)",
                cmap="YlOrBr",
                vmin=0.0,
                vmax=1.0,
                show_pooling_note=False,
                bar_label="xi",
            )
            self._draw_stat_heatmap(
                ax_rho,
                pd_,
                median_key="rho",
                pooled_lookup=lambda label: (
                    pd_.distinguish.pooled_rho_per_threshold.get(label, float("nan"))
                ),
                title="Spearman rho, sign = direction of the trend",
                cmap="PuOr_r",
                vmin=-1.0,
                vmax=1.0,
                show_pooling_note=True,
                bar_label="rho",
            )

            self._draw_threshold_timeline(
                ax_compliance,
                pd_,
                pd_.reliability.compliance_state_series,
                title="Compliance by threshold (2-state, as carved)",
            )
            self._draw_survival(ax_surv, pd_)

            if ax_cum_time is not None:
                theme.apply_common_style(ax_cum_time)
                self._draw_cumulative_time(ax_cum_time, pd_)
            if ax_cum_dmg is not None:
                theme.apply_common_style(ax_cum_dmg)
                self._draw_cumulative_damage(ax_cum_dmg, pd_)

            self._draw_summary(ax_sum, pd_)
        return fig

    def build_plotly(self, result: object) -> go.Figure:
        raise NotImplementedError(f"{self.__class__.__name__} has no plotly backend")

    # --- private drawing methods (read precomputed fields; no data arithmetic) ---

    def _draw_primary(
        self,
        ax: plt.Axes,
        pd_: WithinCalibrationPanelData,
        color: object,
        legend_ax: plt.Axes | None = None,
    ) -> None:
        plot_fn = ax.semilogy if pd_.use_log_scale else ax.plot
        # One call per observed stretch, split at the gaps the carve found, so no
        # segment is drawn across unobserved time.
        for lo, hi in render.observed_slices(pd_):
            plot_fn(
                pd_.signal.t_h[lo:hi],
                pd_.signal.values[lo:hi],
                color=color,
                linewidth=1.2,
                zorder=2,
            )
        # Extra traces overlay the same axis: the detail and binned subplots that used
        # to carry them are gone, and a computed series that is rendered nowhere is a
        # silent capability loss (fidelity supplies two).
        for offset, (trace_label, series) in enumerate(pd_.traces or []):
            if trace_label == pd_.primary_label:
                continue
            for lo, hi in render.observed_slices(pd_):
                plot_fn(
                    pd_.signal.t_h[lo:hi],
                    series[lo:hi],
                    linewidth=0.9,
                    alpha=0.75,
                    linestyle="--",
                    color=theme.mix_with_white(color, 0.35 + 0.15 * offset),
                    label=trace_label if lo == 0 else None,
                    zorder=2,
                )
        if pd_.signal.sigma is not None:
            # Error bars are drawn UNDER the trace and the threshold lines - the point
            # of showing them is to see which reads they overlap a threshold with.
            # The y-limits are snapshotted from the series and restored afterwards:
            # per-read sigma can be many times the signal range, and letting it drive
            # autoscale would squash the whole ladder into a few pixels.
            y_limits = ax.get_ylim()
            ax.errorbar(
                pd_.signal.t_h,
                pd_.signal.values,
                yerr=pd_.signal.sigma,
                fmt="none",
                ecolor=color,
                elinewidth=0.5,
                alpha=0.25,
                zorder=0,
            )
            ax.set_ylim(y_limits)
        if pd_.use_log_scale:
            render.draw_decade_guides(ax, pd_.signal.values)
        for i, (label, thr_val, _) in enumerate(pd_.thresholds):
            thr_color = theme.threshold_color(i, len(pd_.thresholds))
            ax.axhline(
                thr_val,
                color=thr_color,
                linestyle="--",
                linewidth=1.0,
                alpha=0.75,
                label=label,
                zorder=1,
            )
        if pd_.thresholds:
            # Anchored just right of the DISTRIBUTION, so it reads as belonging to the
            # signal block it labels rather than to the error histograms further right.
            (legend_ax or ax).legend(
                *ax.get_legend_handles_labels(),
                frameon=False,
                loc="upper left",
                bbox_to_anchor=(1.01, 1.0),
                borderaxespad=0.0,
            )
        ax.set_ylabel(pd_.primary_label)
        ax.set_xlabel("Scan clock (h)")
        ax.set_title(pd_.primary_label)
        ax.grid(True, which="both", color="lightgray", alpha=0.4)

    def _draw_primary_hist(
        self, ax: plt.Axes, pd_: WithinCalibrationPanelData, color: object
    ) -> None:
        """The series' own distribution, rotated to share the primary y axis."""
        if len(pd_.signal.value_hist_counts) == 0:
            ax.text(
                0.5,
                0.5,
                "no spread",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.axis("off")
            return
        centers, widths = render.bin_centers_and_widths(pd_.signal.value_hist_edges)
        # Drawn as a marginal of the series, not a chart in its own right: it shares the
        # y axis, carries no frame, title or ticks, and sits flush against the trace so
        # the two read as one image of the signal.
        ax.barh(
            centers,
            pd_.signal.value_hist_counts,
            height=widths,
            color=color,
            alpha=0.45,
            edgecolor="none",
        )
        ax.set_xlim(0, float(np.max(pd_.signal.value_hist_counts)) * 1.02)
        ax.tick_params(
            axis="both",
            which="both",
            labelleft=False,
            labelbottom=False,
            left=False,
            bottom=False,
        )
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_facecolor("none")
        ax.grid(False)

    def _draw_sigma_hist(
        self, ax: plt.Axes, pd_: WithinCalibrationPanelData, color: object
    ) -> None:
        """Spread of the per-read fit error itself, with its mean and median marked."""
        if len(pd_.signal.sigma_hist_counts) == 0:
            ax.text(
                0.5,
                0.5,
                "no per-read error",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.axis("off")
            return
        centers, widths = render.bin_centers_and_widths(pd_.signal.sigma_hist_edges)
        ax.bar(
            centers,
            pd_.signal.sigma_hist_counts,
            width=widths,
            color=color,
            alpha=0.6,
            edgecolor="none",
        )
        self._mark_centre(ax, pd_.signal.sigma_mean, pd_.signal.sigma_median)
        # A handful of failed fits carry errors orders of magnitude above the bulk and
        # would flatten the distribution against the left edge. Clip the VIEW, never
        # the data, and say how many reads fall outside it.
        if np.isfinite(pd_.signal.sigma_hist_view_x_max):
            ax.set_xlim(
                float(pd_.signal.sigma_hist_edges[0]), pd_.signal.sigma_hist_view_x_max
            )
        n_outside = pd_.signal.sigma_hist_n_above_view
        ax.set_xlabel(f"Fit error ({pd_.primary_label})")
        ax.set_ylabel("Count")
        # The clipped count goes in the title: this axis is too short to carry a second
        # annotation, and a panel that narrows its view has to say by how much.
        ax.set_title(
            "Per-read fit error"
            if not n_outside
            else f"Per-read fit error ({n_outside} above view)"
        )
        ax.grid(True, axis="y", color="lightgray", alpha=0.4)

    @staticmethod
    def _mark_centre(ax: plt.Axes, mean: float, median: float) -> None:
        """Mark mean and median, each labelled AT its own line.

        A shared corner block listing both left the reader matching two numbers to a
        solid and a dashed line by guesswork; on a skewed distribution, which is which
        is exactly the thing worth knowing.
        """
        # Staggered heights: on a skewed distribution mean and median sit close, and
        # two rotated labels at the same height overprint each other.
        # Opposite sides AND staggered heights: on a skewed distribution mean and
        # median sit close together, and any single placement rule overprints them.
        for value, style, name, y, side in (
            (mean, "-", "mean", 0.98, 1),
            (median, "--", "median", 0.72, -1),
        ):
            if not np.isfinite(value):
                continue
            ax.axvline(value, color="gray", linestyle=style, linewidth=1.1)
            ax.annotate(
                f"{name} {value:.3g}",
                xy=(value, y),
                xycoords=("data", "axes fraction"),
                xytext=(3 * side, 0),
                textcoords="offset points",
                ha="left" if side > 0 else "right",
                va="top",
                rotation=90,
                color="gray",
            )

    def _draw_relerr_hist(
        self, ax: plt.Axes, pd_: WithinCalibrationPanelData, color: object
    ) -> None:
        """Distribution of sigma / value: how large the error is RELATIVE to the metric."""
        if len(pd_.signal.relerr_hist_counts) == 0:
            ax.text(
                0.5,
                0.5,
                "no relative error",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.axis("off")
            return
        centers, widths = render.bin_centers_and_widths(pd_.signal.relerr_hist_edges)
        ax.bar(
            centers,
            pd_.signal.relerr_hist_counts,
            width=widths,
            color=color,
            alpha=0.6,
            edgecolor="none",
        )
        self._mark_centre(ax, pd_.signal.relerr_mean, pd_.signal.relerr_median)
        title = "Relative error (sigma/value)"
        if pd_.signal.relerr_n_excluded:
            title += f" ({pd_.signal.relerr_n_excluded} excluded)"
        ax.set_xlabel("sigma / value")
        ax.set_ylabel("Count")
        ax.set_title(title)
        ax.grid(True, axis="y", color="lightgray", alpha=0.4)

    def _draw_threshold_timeline(
        self,
        ax: plt.Axes,
        pd_: WithinCalibrationPanelData,
        segments_per_threshold: dict[str, list[tuple[float, float, str]]],
        title: str,
        show_ttf: bool = True,
        show_key: bool = False,
    ) -> None:
        """One Gantt row per threshold, coloured by per-read state.

        Takes the segment map rather than reading one, because the panel draws this
        TWICE: band 2 passes the 4-state distinguishability view, band 3 passes the
        2-state view the carve actually used. Comparing them is the panel's argument,
        so they must be drawn by the same code on the same x-scale.
        """
        if not pd_.thresholds:
            ax.text(
                0.5,
                0.5,
                "No thresholds defined",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.axis("off")
            return

        # Only plot thresholds with >=5% in-spec time; keep all in textual summary.
        # (Decision: keep the cull, preserving current figures.)
        plotted = [
            (label, thr_val, bvg)
            for label, thr_val, bvg in pd_.thresholds
            if pd_.reliability.occupancy.get(label, 0.0) >= 0.05
        ]

        n_culled = len(pd_.thresholds) - len(plotted)
        if n_culled:
            title = f"{title} ({n_culled} thresholds below 5% in spec)"
        if not plotted:
            ax.text(
                0.5,
                0.5,
                "No thresholds with >=5% in-spec time",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.axis("off")
            return

        ttf_map = (
            pd_.reliability.ttf_per_threshold if (pd_.include_ttf and show_ttf) else {}
        )
        n_plot = len(plotted)
        y_positions = np.arange(n_plot)[::-1]

        for y_pos, (label, _thr_val, _big_values_good) in zip(y_positions, plotted):
            # Segments and their states are computed by the builder from the read table
            # (direction and uncertainty already applied). Draw only.
            for t_start_h, t_end_h, state in segments_per_threshold.get(label, []):
                unobserved = state == STATE_UNOBSERVED
                ax.barh(
                    y_pos,
                    t_end_h - t_start_h,
                    left=t_start_h,
                    height=0.75,
                    color=theme.state_color(state),
                    alpha=1.0 if unobserved else 0.9,
                    edgecolor=theme.UNOBSERVED_EDGE if unobserved else "none",
                    hatch=theme.UNOBSERVED_HATCH if unobserved else None,
                    linewidth=0.4 if unobserved else 0.0,
                )

            if pd_.include_ttf and show_ttf:
                ttf = ttf_map.get(label)
                ttf_str = f"TTF={ttf:.1f}h" if ttf is not None else "-"
                # y fraction: assumes ylim = [-0.5, n_plot - 0.5] (set below)
                y_frac = (y_pos + 0.5) / n_plot
                ax.text(
                    1.01,
                    y_frac,
                    ttf_str,
                    transform=ax.transAxes,
                    ha="left",
                    va="center",
                    fontsize=7,
                    clip_on=False,
                )

        if show_key:
            # Five states reach these bars and one of them is white-on-white but for
            # its hatch; without a key a reader cannot tell "no data" from "not
            # applicable". Drawn once, on the upper timeline of the pair.
            from matplotlib.patches import Patch

            handles = [
                Patch(facecolor=theme.state_color(state), label=name)
                for state, name in (
                    ("in_spec", "in spec"),
                    ("in_spec_uncertain", "in spec, error overlaps threshold"),
                    ("out_of_spec", "out of spec"),
                    ("out_of_spec_uncertain", "out of spec, error overlaps threshold"),
                )
            ]
            handles.append(
                Patch(
                    facecolor=theme.state_color("unobserved"),
                    edgecolor=theme.UNOBSERVED_EDGE,
                    hatch=theme.UNOBSERVED_HATCH,
                    label="unobserved (read gap)",
                )
            )
            ax.legend(
                handles=handles,
                frameon=False,
                loc="upper left",
                bbox_to_anchor=(1.005, 1.0),
                borderaxespad=0.0,
            )

        ax.set_ylim(-0.5, n_plot - 0.5)
        ax.set_yticks(np.arange(n_plot))
        ax.set_yticklabels([label for label, _, _ in plotted][::-1], fontsize=8)
        ax.set_xlabel("Scan clock (h)")
        ax.set_ylabel("Threshold")
        ax.set_title(title)
        ax.set_xlim(float(np.min(pd_.signal.t_h)), float(np.max(pd_.signal.t_h)))
        ax.grid(True, axis="x", color="lightgray", alpha=0.4)

    def _draw_stat_heatmap(
        self,
        ax: plt.Axes,
        pd_: WithinCalibrationPanelData,
        *,
        median_key: str,
        pooled_lookup,
        title: str,
        cmap: str,
        vmin: float,
        vmax: float,
        show_pooling_note: bool,
        bar_label: str,
    ) -> None:
        """Chatterjee's xi per threshold, two ways: median over windows, and pooled.

        Replaces a zoom on one hand-picked threshold. That zoom showed the same reads as
        band 1 and privileged a threshold for no stated reason; this covers the whole
        ladder and makes the choice unnecessary.

        The two rows answer different questions and are shown together because either
        alone misleads: the median asks what a typical excursion looks like and is NaN
        wherever no window is long enough; the pooled value has coverage almost
        everywhere but can be driven by between-window structure rather than within.
        """
        labels = [label for label, _, _ in pd_.thresholds]
        if not labels:
            ax.text(
                0.5,
                0.5,
                "No thresholds defined",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.axis("off")
            return

        median_row, pooled_row, notes = [], [], []
        thin_pooling: list[str] = []
        for label in labels:
            median_row.append(
                pd_.distinguish.shape_medians_per_threshold.get(label, {}).get(
                    median_key, float("nan")
                )
            )
            # A 5-tuple: the trailing element is the calibration label that
            # produced the pooled p-value. The renderer uses only the counts.
            _xi, _p, n_reads, n_windows, _method = (
                pd_.distinguish.pooled_xi_per_threshold.get(
                    label, (float("nan"), float("nan"), 0, 0, "none")
                )
            )
            pooled_row.append(pooled_lookup(label))
            # The two rows count DIFFERENT things and must say so: the median rests on
            # windows that cleared shape_min_reads, the pooled value on every complete
            # window and its reads. Reads-per-window is shown for the pooled row because
            # a pooled xi over 475 windows of 1.6 reads each is not a within-excursion
            # statistic at all.
            notes.append(
                (
                    pd_.distinguish.shape_defined_counts_per_threshold.get(
                        label, {}
                    ).get("xi_m_to_f", 0),
                    n_windows,
                    n_reads,
                )
            )

        grid = np.array([median_row, pooled_row], dtype=float)
        # FIXED scale, not data-derived: on a dynamic scale 0.12 renders as the darkest
        # cell in the figure and reads as a strong effect. Against the statistic's own
        # full range it reads as what it is, and two runs become comparable by colour.
        image = ax.imshow(
            grid,
            aspect="auto",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
        )
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_yticks([0, 1])
        ax.set_yticklabels(
            [
                f"median over windows\n(>= {pd_.distinguish.shape_min_reads} reads)",
                "pooled over reads\n(all complete)",
            ]
        )
        # `label` MUST be re-bound from `labels` here. Reading the leaked loop variable
        # from the build loop above tested the LAST threshold on every column, so on the
        # shipped artifact - where 4-9 us are thin and 10 us is not - the list came out
        # empty and the caveat vanished from the figure entirely.
        for col, (label, (n_win_median, n_win_pooled, n_reads_pooled)) in enumerate(
            zip(labels, notes, strict=True)
        ):
            if pd_.distinguish.pooled_is_thin_per_threshold.get(label, False):
                thin_pooling.append(label)
            captions = (f"{n_win_median}w", f"{n_win_pooled}w")
            for row, caption in enumerate(captions):
                value = grid[row, col]
                shown = f"{value:.2f}" if np.isfinite(value) else "-"
                # Contrast against the fixed scale, not against the data. A DIVERGING
                # map is dark at both ends; a SEQUENTIAL one only at the top, so the
                # same midpoint rule would put white text on a pale cell.
                if vmin < 0.0:
                    dark = np.isfinite(value) and abs(
                        value - 0.5 * (vmin + vmax)
                    ) > 0.35 * (vmax - vmin)
                else:
                    dark = np.isfinite(value) and value > vmin + 0.6 * (vmax - vmin)
                ax.text(
                    col,
                    row,
                    f"{shown}\n{caption}",
                    ha="center",
                    va="center",
                    color="white" if dark else "black",
                )
        if thin_pooling and show_pooling_note:
            # The number that decides whether a pooled xi means anything: pooling over
            # windows of one or two reads correlates across the whole record, not
            # within excursions. Named once, rather than crowding ten cells.
            # Drawn into a shelf INSIDE the axes - imshow's extent is fixed, so text
            # below it is invisible to constrained_layout and collided with the next row.
            ax.set_ylim(2.6, -0.5)
            ax.text(
                -0.5,
                1.85,
                "pooled over windows averaging too few reads to be a "
                "within-excursion statistic:\n"
                f"{', '.join(thin_pooling)}",
                ha="left",
                va="center",
                color="gray",
            )
        ax.set_title(title)
        ax.grid(False)
        fig = ax.get_figure()
        bar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.01)
        bar.set_label(bar_label)

    def _draw_shape(self, ax: plt.Axes, pd_: WithinCalibrationPanelData) -> None:
        """Excursion shape: margin against normalised window age, supported thresholds only.

        A threshold with a handful of complete windows produces a median that looks like a
        measurement and is not one, so unsupported thresholds are named rather than drawn.
        """
        supported = [
            label
            for label, _, _ in pd_.thresholds
            if pd_.distinguish.shape_supported_per_threshold.get(label, False)
        ]
        skipped = [
            label
            for label, _, _ in pd_.thresholds
            if not pd_.distinguish.shape_supported_per_threshold.get(label, False)
        ]
        if not supported:
            ax.text(
                0.5,
                0.5,
                "no threshold has enough complete windows for a shape",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.axis("off")
            return

        annotations: list[str] = []
        for label in supported:
            grid, mean, stderr, n = pd_.distinguish.shape_curve_per_threshold.get(
                label, (None, None, None, 0)
            )
            if mean is None:
                continue
            color = theme.threshold_color(
                [lbl for lbl, _, _ in pd_.thresholds].index(label),
                len(pd_.thresholds),
            )
            ax.plot(grid, mean, color=color, linewidth=1.6, label=f"{label} (n={n})")
            if stderr is not None:
                ax.fill_between(
                    grid,
                    mean - stderr,
                    mean + stderr,
                    color=color,
                    alpha=theme.BAND_STYLE["fill_alpha"],
                    edgecolor="none",
                )
            annotations.append(
                f"{label}: "
                + render.shape_annotation(
                    pd_.distinguish.shape_medians_per_threshold.get(label, {}),
                    pd_.distinguish.shape_defined_counts_per_threshold.get(label, {}),
                    pd_.distinguish.peak_is_symmetric_per_threshold.get(label, False),
                )
            )
        if annotations:
            # Give the text a shelf INSIDE the axes by extending the y-limit downward,
            # then draw into it. Text placed outside the axes is invisible to
            # constrained_layout, so it collided with the next row's title.
            lo, hi = ax.get_ylim()
            shelf = 0.30 * (hi - lo)
            ax.set_ylim(lo - shelf, hi)
            ax.text(
                0.01,
                0.02,
                "\n".join(annotations),
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                color="gray",
                family="monospace",
            )
        ax.set_xlabel("Window age (normalised)")
        ax.set_ylabel("Margin / peak margin")
        title = "Excursion shape"
        if skipped:
            title += f" ({len(skipped)} thresholds without support)"
        ax.set_title(title)
        ax.legend(frameon=False, loc="upper right")
        ax.grid(True, color="lightgray", alpha=0.4)

    def _draw_survival(self, ax: plt.Axes, pd_: WithinCalibrationPanelData) -> None:
        if not pd_.thresholds:
            ax.text(
                0.5,
                0.5,
                "No thresholds defined",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.axis("off")
            return

        legend_kwargs = dict(
            frameon=False,
            fontsize=7,
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0.0,
        )
        plotted = 0
        excluded_only = 0
        for i, (label, _thr_val, _bvg) in enumerate(pd_.thresholds):
            color = theme.threshold_color(i, len(pd_.thresholds))
            # FIGURE_STANDARD: a panel that drops data says how much, in the panel. The
            # estimator drops censored windows, so each threshold's entry carries its count,
            # a threshold whose every window was censored included: it has no curve to label.
            dropped = pd_.reliability.n_censored_dropped[label]
            carved = pd_.reliability.n_windows[label]
            note = f"{label} - excluded: {dropped} of {carved} windows (censored)"
            survival = pd_.reliability.survival_curve_min.get(label, [])
            if not survival:
                if carved:
                    # No line sample: the entry reports a count, not a curve.
                    ax.plot(
                        [], [], color=color, linestyle="none", label=f"{note}, no curve"
                    )
                    excluded_only += 1
                continue
            surv_x, surv_y = zip(*survival)
            ax.semilogy(
                surv_x,
                surv_y,
                linewidth=1.2,
                markersize=4,
                markevery=max(1, len(surv_x) // 10),
                color=color,
                label=note,
                linestyle="-",
            )
            plotted += 1

        if plotted == 0:
            ax.text(
                0.5,
                0.5,
                "No complete in-spec windows for defined thresholds"
                if excluded_only
                else "No in-spec windows for defined thresholds",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            if excluded_only:
                ax.legend(**legend_kwargs)
            ax.axis("off")
            return

        ax.set_xlabel("Window length (minutes)")
        ax.set_ylabel("Fraction of windows lasting >= length")
        ax.set_title(
            f"In-spec window survival per threshold "
            f"({estimator_name(pd_.reliability.estimator)})"
        )
        ax.grid(True, which="both", color="lightgray", alpha=0.4)
        ax.legend(**legend_kwargs)

    def _draw_cumulative_time(
        self, ax: plt.Axes, pd_: WithinCalibrationPanelData
    ) -> None:
        """Cumulative time out of spec per threshold (hours)."""
        if not pd_.thresholds:
            ax.text(
                0.5,
                0.5,
                "No thresholds defined",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.axis("off")
            return

        cum_time = pd_.reliability.cumulative_time_per_threshold
        t = np.asarray(pd_.signal.t_h, dtype=float)
        slices = render.observed_slices(pd_)
        plotted = 0
        for i, (label, _, _) in enumerate(pd_.thresholds):
            arr = cum_time.get(label)
            if arr is None:
                continue
            if len(arr) != len(t):
                raise ValueError(
                    f"cumulative time for {label!r} has {len(arr)} points against {len(t)} "
                    "reads on the scan clock"
                )
            if len(arr) == 0:
                continue
            # One line per observed stretch: nothing is drawn across a gap. Only the first
            # stretch carries the legend label.
            for j, (lo, hi) in enumerate(slices):
                ax.plot(
                    t[lo:hi],
                    arr[lo:hi],
                    color=theme.threshold_color(i, len(pd_.thresholds)),
                    linewidth=1.2,
                    label=label if j == 0 else None,
                )
            plotted += 1

        if plotted == 0:
            ax.text(
                0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes
            )
            return

        ax.set_xlabel("Scan clock (h)")
        ax.set_ylabel("Cumulative time out of spec (h)")
        ax.set_title("Cumulative time out of spec")
        ax.grid(True, alpha=0.25)
        # Half-width axis: keep the legend inside, in the empty top-left corner
        # (cumulative curves rise toward top-right). Columns cap it at ≤5 rows.
        ax.legend(
            frameon=False,
            fontsize=7,
            loc="upper left",
            ncol=max(1, (len(pd_.thresholds) + 4) // 5),
        )

    def _draw_cumulative_damage(
        self, ax: plt.Axes, pd_: WithinCalibrationPanelData
    ) -> None:
        """Cumulative damage per threshold (primary_unit · h)."""
        if not pd_.thresholds:
            ax.text(
                0.5,
                0.5,
                "No thresholds defined",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.axis("off")
            return

        cum_dmg = pd_.reliability.cumulative_damage_per_threshold
        t = np.asarray(pd_.signal.t_h, dtype=float)
        slices = render.observed_slices(pd_)
        ylabel = f"Cumulative damage ({pd_.primary_label} · h)"
        plotted = 0
        for i, (label, _, _) in enumerate(pd_.thresholds):
            arr = cum_dmg.get(label)
            if arr is None:
                continue
            if len(arr) != len(t):
                raise ValueError(
                    f"cumulative damage for {label!r} has {len(arr)} points against {len(t)} "
                    "reads on the scan clock"
                )
            if len(arr) == 0:
                continue
            # One line per observed stretch: nothing is drawn across a gap. Only the first
            # stretch carries the legend label.
            for j, (lo, hi) in enumerate(slices):
                ax.plot(
                    t[lo:hi],
                    arr[lo:hi],
                    color=theme.threshold_color(i, len(pd_.thresholds)),
                    linewidth=1.2,
                    label=label if j == 0 else None,
                )
            plotted += 1

        if plotted == 0:
            ax.text(
                0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes
            )
            return

        ax.set_xlabel("Scan clock (h)")
        ax.set_ylabel(ylabel)
        ax.set_title("Cumulative damage")
        ax.grid(True, alpha=0.25)
        # Half-width axis: legend in the empty top-left corner, columnized.
        ax.legend(
            frameon=False,
            fontsize=7,
            loc="upper left",
            ncol=max(1, (len(pd_.thresholds) + 4) // 5),
        )

    def _draw_summary(self, ax: plt.Axes, pd_: WithinCalibrationPanelData) -> None:
        series = pd_.signal.values
        cv = pd_.signal.cv
        finite = series[np.isfinite(series)]

        lines: list[str] = [
            f"Metric: {pd_.primary_label}",
            f"Points: {len(finite)}",
            f"Initial value: {float(finite[0]):.6g}"
            if len(finite) > 0
            else "Initial value: N/A",
            f"Range: {float(np.min(finite)):.6g} to {float(np.max(finite)):.6g}"
            if len(finite) > 0
            else "Range: N/A",
            f"CV: {cv:.4f}" if np.isfinite(cv) else "CV: N/A",
        ]
        if pd_.meta:
            lines.append("")
            for k, v in list(pd_.meta.items())[:4]:
                lines.append(f"{k}: {v}")

        ttf_map = pd_.reliability.ttf_per_threshold if pd_.include_ttf else {}

        for label, thr_val, big_values_good in pd_.thresholds:
            summ = pd_.reliability.threshold_summary.get(label)
            if summ is None:
                continue
            time_oos_h = summ["time_oos_h"]
            frac_oos = summ["frac_oos_pct"]
            lines.append("")
            lines.append(f"{label} (thr={thr_val:.6g}):")
            lines.append(f"  Out of spec: {time_oos_h:.2f} h ({frac_oos:.1f}%)")
            w = pd_.reliability.threshold_window_stats.get(label, {}).get(
                "raw_series_runs", {"above": {}, "below": {}}
            )
            # big_values_good=False (infidelity): above threshold = oos, below = in-spec
            # big_values_good=True  (T2*):        below threshold = oos, above = in-spec
            oos_key = "above" if not big_values_good else "below"
            in_spec_key = "below" if not big_values_good else "above"
            ws_oos = w[oos_key]
            ws_in_spec = w[in_spec_key]
            if ws_oos.get("count", 0) > 0:
                lines.append(
                    f"  oos: count={ws_oos['count']}, mean={ws_oos['mean']:.1f} min, p90={ws_oos['p90']:.1f} min"
                )
            if ws_in_spec.get("count", 0) > 0:
                lines.append(
                    f"  in-spec: count={ws_in_spec['count']}, mean={ws_in_spec['mean']:.1f} min, p90={ws_in_spec['p90']:.1f} min"
                )
            if pd_.include_ttf and label in ttf_map:
                first_cross = ttf_map[label]
                if first_cross is None:
                    lines.append("  First crossing: none in dataset")
                else:
                    lines.append(f"  First crossing: {first_cross:.3f} h")

        ax.axis("off")
        ax.text(
            0.02,
            0.98,
            "\n".join(lines).strip(),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.5,
            family="monospace",
            bbox={
                "facecolor": "lightyellow",
                "edgecolor": "gray",
                "boxstyle": "round,pad=0.5",
            },
        )


__all__ = ["WithinCalibrationPanel", "WithinCalibrationPanelData"]
