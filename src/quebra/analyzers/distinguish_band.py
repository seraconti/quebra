"""Band 2 of the composed within-calibration panel: is the failure distinguishable.

Everything here answers one question - given the per-read uncertainty, can a reader tell
an in-spec read from an out-of-spec one at this threshold, and do the excursions have a
shape worth modelling. It carries the 4-state view of the reads; band 3 carries the
2-state view the carve actually used. The two are deliberately separate so the panel can
put them side by side.

Uncertainty is an annotation: `k` and `use_uncertainty` never move a window boundary.
Window carving is `analyzers/windows.py` and is untouched by anything in this module.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from quebra.analyzers.windows import (
    mark_gaps_in_segments,
    SIDE_IN_SPEC,
    STATE_IN_SPEC,
    STATE_IN_SPEC_UNCERTAIN,
    STATE_OUT_OF_SPEC,
    STATE_OUT_OF_SPEC_UNCERTAIN,
)
from quebra.core._artifact_guard import StaleArtifactGuard

# A shape statistic computed over a handful of windows is not a measurement. Thresholds with
# fewer surviving complete windows than this are reported as unsupported rather than
# drawn. On a chattering trace most thresholds yield single-read windows and clear
# nothing, which is the point: a median over a handful of windows is not a measurement.
SHAPE_SUPPORT_FLOOR = 20

# Below this, a pooled statistic has stopped describing excursions and is correlating
# across the record, because the "windows" are single reads.
POOLING_MIN_READS_PER_WINDOW = 2.0

# A peak this close to mid-window makes the excursion symmetric enough that Spearman is
# suppressed by construction, whatever the underlying predictability.
SYMMETRIC_PEAK_TOLERANCE = 0.08


@dataclass
class DistinguishBand(StaleArtifactGuard):
    """Complete contract for the distinguish band. Keyed by threshold label."""

    # (t_start_h, t_end_h, state) runs, 4-state
    state_series_per_threshold: dict[str, list[tuple[float, float, str]]] = field(
        default_factory=dict
    )
    # n_in, n_out, n_unresolved, n_unknown (+ the four-way split underneath)
    state_counts_per_threshold: dict[str, dict[str, int]] = field(default_factory=dict)
    # fraction of a window's reads that are resolved, one entry per window
    window_purity_per_threshold: dict[str, list[float]] = field(default_factory=dict)
    reads_per_window_per_threshold: dict[str, list[int]] = field(default_factory=dict)
    median_duration_over_read_spacing: dict[str, float] = field(default_factory=dict)
    # Half-width of the band around a threshold that the per-read errors cannot
    # separate: k * median sigma over the reads actually classified. Computed here so
    # the shaded region and the unresolved COUNT come from the same definition -
    # taking a median of the data at draw time made them disagree.
    unresolved_band_half_per_threshold: dict[str, float] = field(default_factory=dict)
    frac_windows_ge_3_reads: dict[str, float] = field(default_factory=dict)

    # per-window excursion shape statistics, complete windows only
    shape_stats_per_threshold: dict[str, pd.DataFrame] = field(default_factory=dict)
    # median of each statistic, and how many windows that median rests on
    shape_medians_per_threshold: dict[str, dict[str, float]] = field(
        default_factory=dict
    )
    shape_defined_counts_per_threshold: dict[str, dict[str, int]] = field(
        default_factory=dict
    )
    # (age_grid, mean_profile, stderr, n_windows) from the shape curve
    shape_curve_per_threshold: dict[str, tuple] = field(default_factory=dict)
    shape_supported_per_threshold: dict[str, bool] = field(default_factory=dict)
    # xi over all complete windows' reads pooled, per threshold: (xi, p, n_reads,
    # n_windows). Available at thresholds where no single window is long enough, which is
    # most of them - see the module docstring for what it does and does not say.
    # (xi, p_value, n_reads_pooled, n_windows_pooled, calibration_method)
    pooled_xi_per_threshold: dict[str, tuple[float, float, int, int, str]] = field(
        default_factory=dict
    )
    # Verdicts, not view choices: both get NAMED on the figure and a reader takes them
    # as results, so they are computed here and are visible to provenance and to any
    # other consumer of the artifact.
    pooled_reads_per_window_per_threshold: dict[str, float] = field(
        default_factory=dict
    )
    pooled_is_thin_per_threshold: dict[str, bool] = field(default_factory=dict)
    # Spearman over the same pooled sample; sign is meaningful, so it is not squared.
    pooled_rho_per_threshold: dict[str, float] = field(default_factory=dict)
    peak_is_symmetric_per_threshold: dict[str, bool] = field(default_factory=dict)

    # why windows were excluded from the shape statistics
    n_excluded_endurance: dict[str, int] = field(default_factory=dict)
    n_excluded_censored: dict[str, int] = field(default_factory=dict)
    n_excluded_short: dict[str, int] = field(default_factory=dict)
    shape_min_reads: int = 5
    k: float = 1.0
    # False means uncertainty was never evaluated, which is NOT the same as
    # "every read resolved" - without it a reader cannot tell those apart.
    use_uncertainty: bool = False

    def check_thresholds(self, labels: list[str]) -> None:
        """One entry per threshold in every per-threshold map.

        Called by the outer contract, which owns `thresholds`. Nesting moved these maps
        off the class that knows the ladder, so the sweep has to be handed the labels
        rather than reading them - if this is not called, the completeness contract
        silently drops to nothing while the unpickle guard stays green.
        """
        # The shape_* maps are deliberately NOT swept. run() writes an entry for every
        # threshold, empty where there is no support, so sweeping them would only
        # re-check what the loop already guarantees; the support FLAG is swept instead
        # and that is what tells a reader whether an entry means anything.
        maps = {
            "state_series_per_threshold": self.state_series_per_threshold,
            "state_counts_per_threshold": self.state_counts_per_threshold,
            "window_purity_per_threshold": self.window_purity_per_threshold,
            "reads_per_window_per_threshold": self.reads_per_window_per_threshold,
            "median_duration_over_read_spacing": self.median_duration_over_read_spacing,
            "frac_windows_ge_3_reads": self.frac_windows_ge_3_reads,
            "unresolved_band_half_per_threshold": self.unresolved_band_half_per_threshold,
            "n_excluded_endurance": self.n_excluded_endurance,
            "n_excluded_censored": self.n_excluded_censored,
            "n_excluded_short": self.n_excluded_short,
            "shape_supported_per_threshold": self.shape_supported_per_threshold,
            "pooled_xi_per_threshold": self.pooled_xi_per_threshold,
            "pooled_reads_per_window_per_threshold": (
                self.pooled_reads_per_window_per_threshold
            ),
            "pooled_is_thin_per_threshold": self.pooled_is_thin_per_threshold,
            "pooled_rho_per_threshold": self.pooled_rho_per_threshold,
            "peak_is_symmetric_per_threshold": self.peak_is_symmetric_per_threshold,
        }
        for name, mapping in maps.items():
            missing = [label for label in labels if label not in mapping]
            if missing:
                raise ValueError(
                    f"incomplete DistinguishBand: {name} is missing threshold(s) "
                    f"{missing} - construct via analyzers.distinguish_band.run()"
                )


@dataclass(slots=True)
class DistinguishBandInputs:
    reads: pd.DataFrame
    windows: pd.DataFrame
    thresholds: list[tuple[str, float, bool]]
    median_read_spacing_s: float
    gap_spans_h: list[tuple[float, float]]
    # per-read sigma in the panel's DISPLAY units (the read table's is SI)
    sigma_display: np.ndarray | None
    shape_min_reads: int = 5
    k: float = 1.0
    use_uncertainty: bool = False
    # Seed for xi's permutation calibration. Only consulted where the data is TIED - the
    # tie-free path is the closed form and uses no randomness at all - but it is a
    # constructor field rather than a default so that a job which ever hits tied data
    # declares it, and it reaches the provenance label like every other step kwarg.
    xi_seed: int = 0


def make_inputs_from_windows(
    *,
    reads: pd.DataFrame,
    windows: pd.DataFrame,
    thresholds: list[tuple[str, float, bool]],
    median_read_spacing_s: float,
    gap_spans_h: list[tuple[float, float]],
    sigma_display: np.ndarray | None,
    shape_min_reads: int = 5,
    k: float = 1.0,
    use_uncertainty: bool = False,
    xi_seed: int = 0,
) -> DistinguishBandInputs:
    return DistinguishBandInputs(
        reads=reads,
        windows=windows,
        thresholds=list(thresholds),
        median_read_spacing_s=median_read_spacing_s,
        gap_spans_h=list(gap_spans_h),
        sigma_display=sigma_display,
        shape_min_reads=shape_min_reads,
        k=k,
        use_uncertainty=use_uncertainty,
        xi_seed=xi_seed,
    )


def _timeline_segments(
    reads: pd.DataFrame, label: str
) -> list[tuple[float, float, str]]:
    """Contiguous (t_start_h, t_end_h, state) runs for one threshold.

    A run ends at the timestamp of the first read of the NEXT run, so the bars are
    contiguous. State is read from the table, never recomputed.
    """
    sub = reads[reads["threshold_label"] == label]
    if len(sub) == 0:
        return []
    t_h = sub["t_read_s"].to_numpy(dtype=float) / 3600.0
    states = sub["state"].to_numpy(dtype=object)
    segments: list[tuple[float, float, str]] = []
    start, state = float(t_h[0]), str(states[0])
    for idx in range(1, len(t_h)):
        if str(states[idx]) != state:
            end = float(t_h[idx])
            segments.append((start, end, state))
            start, state = end, str(states[idx])
    segments.append((start, float(t_h[-1]), state))
    return segments


def run(inputs: DistinguishBandInputs) -> DistinguishBand:
    from quebra.analyzers import shape_stats

    reads, windows = inputs.reads, inputs.windows
    # The shape statistics and window counts describe the in-spec side; the out-of-spec
    # table has the same schema and would be read without complaint.
    sides = sorted(set(windows["side"].unique()) - {SIDE_IN_SPEC})
    if sides:
        raise ValueError(
            f"expected the in-spec window table; got windows on side(s) {sides}. "
            f"Pass WindowsResult.windows_in_spec."
        )
    sigma_display = inputs.sigma_display
    finite_sigma = (
        sigma_display[np.isfinite(sigma_display)]
        if sigma_display is not None
        else np.array([])
    )
    half_display = (
        float(inputs.k * np.median(finite_sigma)) if len(finite_sigma) else float("nan")
    )
    band = DistinguishBand(
        shape_min_reads=inputs.shape_min_reads,
        k=inputs.k,
        use_uncertainty=inputs.use_uncertainty,
    )
    spacing_s = inputs.median_read_spacing_s

    for label, thr_value, big_values_good in inputs.thresholds:
        r = reads[reads["threshold_label"] == label]
        w = windows[windows["threshold_label"] == label]

        band.state_series_per_threshold[label] = mark_gaps_in_segments(
            _timeline_segments(reads, label), inputs.gap_spans_h
        )

        state = r["state"] if len(r) else pd.Series(dtype=object)
        n_unresolved = int(
            state.isin([STATE_IN_SPEC_UNCERTAIN, STATE_OUT_OF_SPEC_UNCERTAIN]).sum()
        )
        band.state_counts_per_threshold[label] = {
            "n_in": int((state == STATE_IN_SPEC).sum()),
            "n_out": int((state == STATE_OUT_OF_SPEC).sum()),
            "n_unresolved": n_unresolved,
            "n_unknown": int((~r["sigma_known"]).sum()) if len(r) else 0,
            "n_in_uncertain": int((state == STATE_IN_SPEC_UNCERTAIN).sum()),
            "n_out_uncertain": int((state == STATE_OUT_OF_SPEC_UNCERTAIN).sum()),
        }

        # DISPLAY units, not the read table's. The carve runs in SI (sigma_v is
        # seconds for T2*) while the panel plots µs, so taking the median off the read
        # table would put a band 1e6 too small on the figure - the unit-mismatch bug
        # class this repo is most exposed to.
        band.unresolved_band_half_per_threshold[label] = half_display

        in_spec_reads = r[r["in_spec"]] if len(r) else r
        if len(in_spec_reads):
            # A read with no error bar is NOT resolved - it is unclassifiable. Folding
            # it in with the crisp reads inflates purity with reads that were never
            # tested, which is the wrong-but-plausible failure mode.
            resolved = in_spec_reads["sigma_known"] & ~in_spec_reads["state"].isin(
                [STATE_IN_SPEC_UNCERTAIN, STATE_OUT_OF_SPEC_UNCERTAIN]
            )
            purity = (
                in_spec_reads.assign(_resolved=resolved)
                .groupby("window_index")["_resolved"]
                .mean()
            )
            band.window_purity_per_threshold[label] = [float(x) for x in purity]
        else:
            band.window_purity_per_threshold[label] = []

        n_reads_per_window = (
            w["n_reads"].to_numpy(dtype=int) if len(w) else np.array([])
        )
        band.reads_per_window_per_threshold[label] = [
            int(x) for x in n_reads_per_window
        ]
        band.frac_windows_ge_3_reads[label] = (
            float((n_reads_per_window >= 3).mean())
            if len(n_reads_per_window)
            else float("nan")  # no windows is not "no window has 3 reads"
        )
        durations = w["duration_s"].to_numpy(dtype=float) if len(w) else np.array([])
        band.median_duration_over_read_spacing[label] = (
            float(np.median(durations) / spacing_s)
            if len(durations) and np.isfinite(spacing_s) and spacing_s > 0
            else float("nan")
        )

        # Shape statistics: complete windows only, then shape_min_reads.
        complete = w[w["birth_observed"] & ~w["censored"]] if len(w) else w
        band.n_excluded_endurance[label] = (
            int((~w["birth_observed"]).sum()) if len(w) else 0
        )
        band.n_excluded_censored[label] = int(w["censored"].sum()) if len(w) else 0
        eligible = (
            complete[complete["n_reads"] >= inputs.shape_min_reads]
            if len(complete)
            else complete
        )
        band.n_excluded_short[label] = int(len(complete) - len(eligible))
        band.shape_supported_per_threshold[label] = len(eligible) >= SHAPE_SUPPORT_FLOOR

        # Pooled over COMPLETE windows, not the shape_min_reads subset: pooling exists
        # to give coverage where per-window statistics cannot.
        band.pooled_xi_per_threshold[label] = shape_stats.pooled_xi(
            r, complete, xi_seed=inputs.xi_seed
        )
        band.pooled_rho_per_threshold[label] = shape_stats.pooled_rho(r, complete)[0]
        _xi, _p, n_reads_pooled, n_win_pooled, _m = band.pooled_xi_per_threshold[label]
        per_window = n_reads_pooled / n_win_pooled if n_win_pooled else float("nan")
        band.pooled_reads_per_window_per_threshold[label] = per_window
        band.pooled_is_thin_per_threshold[label] = bool(
            n_win_pooled and per_window < POOLING_MIN_READS_PER_WINDOW
        )

        stats, medians, defined, curve = shape_stats.for_windows(
            reads=r,
            windows=eligible,
            threshold_value=float(thr_value),
            big_values_good=bool(big_values_good),
            xi_seed=inputs.xi_seed,
        )
        band.shape_stats_per_threshold[label] = stats
        band.shape_medians_per_threshold[label] = medians
        band.shape_defined_counts_per_threshold[label] = defined
        band.shape_curve_per_threshold[label] = curve
        peak = medians.get("peak_frac", float("nan"))
        band.peak_is_symmetric_per_threshold[label] = bool(
            np.isfinite(peak) and abs(peak - 0.5) <= SYMMETRIC_PEAK_TOLERANCE
        )

    supported = [k for k, v in band.shape_supported_per_threshold.items() if v]
    print(
        f"[distinguish_band] thresholds={len(inputs.thresholds)} "
        f"shape-supported={supported or 'none'} (floor {SHAPE_SUPPORT_FLOOR} windows, "
        f"shape_min_reads={inputs.shape_min_reads})",
        flush=True,
    )
    return band
