"""Band 3 of the composed within-calibration panel: reliability.

The 2-state view - the one the carve actually used - and what follows from it: how long
windows last on each side, and how much of the observed record was in spec.

Survival is Kaplan-Meier, read from the job's one Kaplan-Meier node
(`kaplan_meier.KaplanMeierSet`, both sides of every threshold) rather than computed here,
so the panel and any other consumer of that node cannot show two different bands. The set
carries its check outcome; in a job that runs no ledger it reads NOT ASSESSED.

`estimator` is a field, not a comment: the survival title is built from it, so a reader
never has to guess and the label cannot go stale.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from quebra.analyzers.kaplan_meier import KaplanMeierSet
from quebra.analyzers.windows import (
    SIDE_IN_SPEC,
    SIDE_OUT_OF_SPEC,
    DEATH_GAP_START,
    STATE_IN_SPEC,
    STATE_OUT_OF_SPEC,
    mark_gaps_in_segments,
)
from quebra.core._artifact_guard import StaleArtifactGuard

ESTIMATOR_KAPLAN_MEIER = "kaplan_meier"

ESTIMATOR_LABELS = {ESTIMATOR_KAPLAN_MEIER: "Kaplan-Meier"}


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
    # in-spec fraction of observed time; NaN where the record has no observed time
    occupancy: dict[str, float] = field(default_factory=dict)

    estimator: str = ESTIMATOR_KAPLAN_MEIER
    # Both sides of every threshold, curve and log-log band, plus the check outcome.
    # Durations on the curves are MINUTES (`KaplanMeierCurve.time_min`).
    kaplan_meier: KaplanMeierSet | None = None

    # In-spec carve counts, windows whose birth was not observed included.
    n_windows: dict[str, int] = field(default_factory=dict)
    n_censored: dict[str, int] = field(default_factory=dict)
    n_endurance_bags: dict[str, int] = field(default_factory=dict)
    n_gap_terminated: dict[str, int] = field(default_factory=dict)

    # retained from the pre-band panel; all per-threshold, all reliability-side
    cumulative_time_per_threshold: dict[str, np.ndarray] = field(default_factory=dict)
    cumulative_damage_per_threshold: dict[str, np.ndarray] = field(default_factory=dict)
    ttf_per_threshold: dict[str, float | None] = field(default_factory=dict)

    def check_thresholds(self, labels: list[str]) -> None:
        """One entry per threshold in every per-threshold map. See DistinguishBand."""
        maps = {
            "compliance_state_series": self.compliance_state_series,
            "occupancy": self.occupancy,
            "n_windows": self.n_windows,
            "n_censored": self.n_censored,
            "n_endurance_bags": self.n_endurance_bags,
            "n_gap_terminated": self.n_gap_terminated,
            "cumulative_time_per_threshold": self.cumulative_time_per_threshold,
            "cumulative_damage_per_threshold": self.cumulative_damage_per_threshold,
            "ttf_per_threshold": self.ttf_per_threshold,
        }
        for name, mapping in maps.items():
            missing = [label for label in labels if label not in mapping]
            if missing:
                raise ValueError(
                    f"incomplete ReliabilityBand: {name} is missing threshold(s) "
                    f"{missing} - construct via analyzers.reliability_band.run()"
                )
        if not labels:
            return
        if self.kaplan_meier is None:
            raise ValueError("incomplete ReliabilityBand: no Kaplan-Meier set")
        missing_curves = [
            (label, side)
            for label in labels
            for side in (SIDE_IN_SPEC, SIDE_OUT_OF_SPEC)
            if (label, side) not in self.kaplan_meier.curves
        ]
        if missing_curves:
            raise ValueError(
                f"incomplete ReliabilityBand: no Kaplan-Meier curve for {missing_curves}"
            )


@dataclass(slots=True)
class ReliabilityBandInputs:
    t_h: np.ndarray
    values: np.ndarray
    reads: pd.DataFrame
    windows: pd.DataFrame
    thresholds: list[tuple[str, float, bool]]
    gap_spans_h: list[tuple[float, float]]
    kaplan_meier: KaplanMeierSet
    damage_fn: Callable[[np.ndarray], np.ndarray] | None = None


def make_inputs_from_windows(
    *,
    t_h: np.ndarray,
    values: np.ndarray,
    reads: pd.DataFrame,
    windows: pd.DataFrame,
    thresholds: list[tuple[str, float, bool]],
    gap_spans_h: list[tuple[float, float]],
    kaplan_meier: KaplanMeierSet,
    damage_fn: Callable[[np.ndarray], np.ndarray] | None = None,
) -> ReliabilityBandInputs:
    return ReliabilityBandInputs(
        t_h=np.asarray(t_h, dtype=float),
        values=np.asarray(values, dtype=float),
        reads=reads,
        windows=windows,
        thresholds=list(thresholds),
        gap_spans_h=list(gap_spans_h),
        kaplan_meier=kaplan_meier,
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


def _require_same_carve(
    reads: pd.DataFrame,
    windows: pd.DataFrame,
    labels: list[str],
    kaplan_meier: KaplanMeierSet,
) -> None:
    """Raise unless the Kaplan-Meier set was estimated from these tables' carve.

    The band's counts and in-spec fraction come from the tables while the drawn curve
    comes from the set, so a set from another record or another carve would put one
    record's curve under another's numbers. Checked: the tables name one record and every
    curve names it too; every label is on the set's ladder; each in-spec curve counts
    as many carved windows as the table holds under its label.
    """
    off_ladder = [label for label in labels if label not in kaplan_meier.ladder]
    if off_ladder:
        raise ValueError(
            f"the Kaplan-Meier set's ladder {list(kaplan_meier.ladder)} lacks "
            f"threshold(s) {off_ladder}; estimate it from the same carve"
        )
    records = set(reads["dataset_id"].unique()) | set(windows["dataset_id"].unique())
    if len(records) > 1:
        raise ValueError(f"the read and window tables mix records {sorted(records)}")
    if records:
        foreign = sorted(
            {curve.dataset_id for curve in kaplan_meier.curves.values()} - records
        )
        if foreign:
            raise ValueError(
                f"the Kaplan-Meier set was estimated on record(s) {foreign}, the tables "
                f"on {sorted(records)}; estimate it from the same carve"
            )
    for label in labels:
        n_table = int((windows["threshold_label"] == label).sum())
        n_curve = kaplan_meier.curve(label, SIDE_IN_SPEC).n_windows_carved
        if n_table != n_curve:
            raise ValueError(
                f"threshold {label!r}: the in-spec table holds {n_table} window(s) but "
                f"the Kaplan-Meier curve was estimated from {n_curve}; estimate it from "
                f"the same carve"
            )


def run(inputs: ReliabilityBandInputs) -> ReliabilityBand:
    from quebra.analyzers.within_calibration_compute import (
        _cumulative_damage,
        _cumulative_time_out_of_spec,
        _ttf,
        _threshold_in_spec_frac,
    )

    t_arr, s_arr = inputs.t_h, inputs.values
    thresholds, windows = inputs.thresholds, inputs.windows
    # The counts below describe the in-spec side; the out-of-spec table has the same
    # schema and would be read without complaint.
    sides = sorted(set(windows["side"].unique()) - {SIDE_IN_SPEC})
    if sides:
        raise ValueError(
            f"expected the in-spec window table; got windows on side(s) {sides}. "
            f"Pass WindowsResult.windows_in_spec."
        )
    _require_same_carve(
        inputs.reads,
        windows,
        [label for label, _, _ in thresholds],
        inputs.kaplan_meier,
    )
    band = ReliabilityBand(
        estimator=ESTIMATOR_KAPLAN_MEIER, kaplan_meier=inputs.kaplan_meier
    )

    cumulative_time = _cumulative_time_out_of_spec(
        t_arr, s_arr, inputs.reads, thresholds, inputs.gap_spans_h
    )
    cumulative_damage = _cumulative_damage(
        t_arr, s_arr, thresholds, inputs.damage_fn, inputs.gap_spans_h
    )
    ttf = _ttf(t_arr, s_arr, inputs.reads, thresholds)
    in_spec_frac = _threshold_in_spec_frac(
        t_arr, s_arr, inputs.reads, thresholds, inputs.gap_spans_h
    )

    for label, _, _ in thresholds:
        w = windows[windows["threshold_label"] == label]

        band.compliance_state_series[label] = mark_gaps_in_segments(
            _compliance_segments(inputs.reads, label), inputs.gap_spans_h
        )
        band.occupancy[label] = in_spec_frac[label]
        band.n_windows[label] = int(len(w))
        band.n_censored[label] = int(w["censored"].sum()) if len(w) else 0
        band.n_endurance_bags[label] = (
            int((~w["birth_observed"]).sum()) if len(w) else 0
        )
        band.n_gap_terminated[label] = (
            int((w["death_type"] == DEATH_GAP_START).sum()) if len(w) else 0
        )

        band.cumulative_time_per_threshold[label] = cumulative_time[label]
        band.cumulative_damage_per_threshold[label] = cumulative_damage[label]
        band.ttf_per_threshold[label] = ttf[label]

    return band
