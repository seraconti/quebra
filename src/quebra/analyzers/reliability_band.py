"""Band 3 of the composed within-calibration panel: reliability.

The 2-state view - the one the carve actually used - and what follows from it: how long
windows last, how much of the observed record was in spec, and how many windows the
estimator had to discard.

`estimator` is a field, not a comment. The artifact declares which estimator produced its
survival curve, and the survival axis label is built from it, so a reader never has to
guess and the label cannot go stale. When Kaplan-Meier lands, only this module changes:
`estimator` flips to "kaplan_meier", `cumulative_hazard` and the band bounds populate,
and the panel is untouched. That is what the band split is for.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from quebra.analyzers.windows import (
    BIRTH_UP_CROSSING,
    DEATH_GAP_START,
    STATE_IN_SPEC,
    STATE_OUT_OF_SPEC,
    mark_gaps_in_segments,
)
from quebra.core._artifact_guard import StaleArtifactGuard

ESTIMATOR_CRUDE = "crude_empirical"

ESTIMATOR_LABELS = {
    ESTIMATOR_CRUDE: "empirical survival, censored windows dropped",
    "kaplan_meier": "Kaplan-Meier",
}


def estimator_name(estimator: str) -> str:
    """Human name of the estimator, derived from the field so it cannot go stale.

    Used in the survival TITLE, not the y-label: naming the estimator is essential but
    the string is too long to set sideways on an axis without colliding with its
    neighbours.
    """
    if estimator not in ESTIMATOR_LABELS:
        raise ValueError(
            f"unknown survival estimator: {estimator!r}. Known: {sorted(ESTIMATOR_LABELS)}"
        )
    return ESTIMATOR_LABELS[estimator]


@dataclass
class ReliabilityBand(StaleArtifactGuard):
    """Complete contract for the reliability band. Keyed by threshold label."""

    # (t_start_h, t_end_h, state) runs, 2-state: exactly what the carve saw
    compliance_state_series: dict[str, list[tuple[float, float, str]]] = field(
        default_factory=dict
    )
    occupancy: dict[str, float] = field(default_factory=dict)
    # (window length in MINUTES, surviving fraction); suffixed because every other
    # time on this artifact is hours or seconds.
    survival_curve_min: dict[str, list[tuple[float, float]]] = field(
        default_factory=dict
    )

    estimator: str = ESTIMATOR_CRUDE
    # Reserved for the Kaplan-Meier pass; None is the honest value until then.
    cumulative_hazard: dict[str, list[tuple[float, float]]] | None = None
    band_lower: dict[str, list[tuple[float, float]]] | None = None
    band_upper: dict[str, list[tuple[float, float]]] | None = None

    n_windows: dict[str, int] = field(default_factory=dict)
    n_complete: dict[str, int] = field(default_factory=dict)
    n_censored: dict[str, int] = field(default_factory=dict)
    n_endurance_bags: dict[str, int] = field(default_factory=dict)
    n_censored_dropped: dict[str, int] = field(default_factory=dict)
    n_gap_terminated: dict[str, int] = field(default_factory=dict)

    # retained from the pre-band panel; all per-threshold, all reliability-side
    cumulative_time_per_threshold: dict[str, np.ndarray] = field(default_factory=dict)
    cumulative_damage_per_threshold: dict[str, np.ndarray] = field(default_factory=dict)
    ttf_per_threshold: dict[str, float | None] = field(default_factory=dict)
    threshold_summary: dict[str, dict[str, float] | None] = field(default_factory=dict)
    threshold_window_stats: dict[str, dict[str, object]] = field(default_factory=dict)

    def check_thresholds(self, labels: list[str]) -> None:
        """One entry per threshold in every per-threshold map. See DistinguishBand."""
        maps = {
            "compliance_state_series": self.compliance_state_series,
            "occupancy": self.occupancy,
            "survival_curve_min": self.survival_curve_min,
            "n_windows": self.n_windows,
            "n_complete": self.n_complete,
            "n_censored": self.n_censored,
            "n_endurance_bags": self.n_endurance_bags,
            "n_censored_dropped": self.n_censored_dropped,
            "n_gap_terminated": self.n_gap_terminated,
            "cumulative_time_per_threshold": self.cumulative_time_per_threshold,
            "cumulative_damage_per_threshold": self.cumulative_damage_per_threshold,
            "ttf_per_threshold": self.ttf_per_threshold,
            "threshold_summary": self.threshold_summary,
            "threshold_window_stats": self.threshold_window_stats,
        }
        for name, mapping in maps.items():
            missing = [label for label in labels if label not in mapping]
            if missing:
                raise ValueError(
                    f"incomplete ReliabilityBand: {name} is missing threshold(s) "
                    f"{missing} - construct via analyzers.reliability_band.run()"
                )


@dataclass(slots=True)
class ReliabilityBandInputs:
    t_h: np.ndarray
    values: np.ndarray
    reads: pd.DataFrame
    windows: pd.DataFrame
    thresholds: list[tuple[str, float, bool]]
    gap_spans_h: list[tuple[float, float]]
    damage_fn: Callable[[np.ndarray], np.ndarray] | None = None


def make_inputs_from_windows(
    *,
    t_h: np.ndarray,
    values: np.ndarray,
    reads: pd.DataFrame,
    windows: pd.DataFrame,
    thresholds: list[tuple[str, float, bool]],
    gap_spans_h: list[tuple[float, float]],
    damage_fn: Callable[[np.ndarray], np.ndarray] | None = None,
) -> ReliabilityBandInputs:
    return ReliabilityBandInputs(
        t_h=np.asarray(t_h, dtype=float),
        values=np.asarray(values, dtype=float),
        reads=reads,
        windows=windows,
        thresholds=list(thresholds),
        gap_spans_h=list(gap_spans_h),
        damage_fn=damage_fn,
    )


def _compliance_segments(
    reads: pd.DataFrame, label: str
) -> list[tuple[float, float, str]]:
    """2-state (t_start_h, t_end_h, state) runs from `in_spec` alone.

    Deliberately NOT derived from `state`: this is the view the carve used, so it must
    ignore the uncertainty annotation entirely. Comparing it against the distinguish
    band's 4-state timeline is the argument the panel makes.
    """
    sub = reads[reads["threshold_label"] == label]
    if len(sub) == 0:
        return []
    t_h = sub["t_read_s"].to_numpy(dtype=float) / 3600.0
    in_spec = sub["in_spec"].to_numpy(dtype=bool)
    segments: list[tuple[float, float, str]] = []
    start, state = float(t_h[0]), bool(in_spec[0])
    for idx in range(1, len(t_h)):
        if bool(in_spec[idx]) != state:
            end = float(t_h[idx])
            segments.append((start, end, STATE_IN_SPEC if state else STATE_OUT_OF_SPEC))
            start, state = end, bool(in_spec[idx])
    segments.append(
        (start, float(t_h[-1]), STATE_IN_SPEC if state else STATE_OUT_OF_SPEC)
    )
    return segments


def run(inputs: ReliabilityBandInputs) -> ReliabilityBand:
    from quebra.analyzers.within_calibration_compute import (
        _analyze_threshold_windows,
        _cumulative_damage,
        _cumulative_time_out_of_spec,
        _ttf,
        _threshold_in_spec_frac,
        _threshold_summary,
        _window_survival,
    )

    t_arr, s_arr = inputs.t_h, inputs.values
    thresholds, windows = inputs.thresholds, inputs.windows
    band = ReliabilityBand(estimator=ESTIMATOR_CRUDE)

    cumulative_time = _cumulative_time_out_of_spec(
        t_arr, s_arr, thresholds, inputs.gap_spans_h
    )
    cumulative_damage = _cumulative_damage(
        t_arr, s_arr, thresholds, inputs.damage_fn, inputs.gap_spans_h
    )
    ttf = _ttf(t_arr, s_arr, thresholds)
    summary = _threshold_summary(t_arr, s_arr, thresholds, inputs.gap_spans_h)
    in_spec_frac = _threshold_in_spec_frac(t_arr, s_arr, thresholds, inputs.gap_spans_h)

    for label, thr_value, _ in thresholds:
        w = windows[windows["threshold_label"] == label]
        uncensored = w[~w["censored"]] if len(w) else w

        band.compliance_state_series[label] = mark_gaps_in_segments(
            _compliance_segments(inputs.reads, label), inputs.gap_spans_h
        )
        band.occupancy[label] = in_spec_frac[label]
        band.survival_curve_min[label] = _window_survival(
            (uncensored["duration_s"].to_numpy(dtype=float) / 60.0).tolist()
            if len(uncensored)
            else []
        )
        band.n_windows[label] = int(len(w))
        band.n_complete[label] = int(len(uncensored))
        band.n_censored[label] = int(w["censored"].sum()) if len(w) else 0
        band.n_endurance_bags[label] = (
            int((w["birth_type"] != BIRTH_UP_CROSSING).sum()) if len(w) else 0
        )
        # The crude estimator drops censored windows; this is how many it discarded.
        band.n_censored_dropped[label] = band.n_censored[label]
        band.n_gap_terminated[label] = (
            int((w["death_type"] == DEATH_GAP_START).sum()) if len(w) else 0
        )

        band.cumulative_time_per_threshold[label] = cumulative_time[label]
        band.cumulative_damage_per_threshold[label] = cumulative_damage[label]
        band.ttf_per_threshold[label] = ttf[label]
        band.threshold_summary[label] = summary[label]
        band.threshold_window_stats[label] = {
            # Nested, not merged: these above/below runs are re-derived from the raw
            # series and are gap-unaware and carve-unaware. Flattening them beside the
            # carve counts invites a reader to take them for the same windows.
            "raw_series_runs": _analyze_threshold_windows(t_arr, s_arr, thr_value),
            "n_windows": band.n_windows[label],
            "n_censored": band.n_censored[label],
            "n_endurance_bags": band.n_endurance_bags[label],
            "n_gaps": band.n_gap_terminated[label],
            "n_dropped_censored": band.n_censored_dropped[label],
        }

    return band
