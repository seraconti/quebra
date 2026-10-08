"""The typed artifact WithinCalibrationPanel renders: three bands plus what they share.

Composed of one contract per band, each produced by its own step and each complete on
its own:

  signal       analyzers/signal_band.py       what was measured, before any threshold
  distinguish  analyzers/distinguish_band.py  can a reader tell in from out at all
  reliability  analyzers/reliability_band.py  what follows from the 2-state carve

A band is replaceable without touching the others - which is the point. When
Kaplan-Meier lands, `reliability` changes and nothing else does.

The outer class owns only what all three share (the ladder, the axis label, the render
flags) and `meta`. It does NOT own the per-threshold maps any more, so its completeness
check delegates to each band's `check_thresholds`: `StaleArtifactGuard` derives its key
set from `dataclasses.fields(cls)` and would otherwise validate nothing but the four
outer names. Every band inherits the guard for the same reason, and no band may use
`slots=True` - the guard's non-dict branch fires on the slots tuple even for a valid load.

Re-exported from panels.within_calibration, which stays the public import surface.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from quebra.analyzers.distinguish_band import DistinguishBand
from quebra.analyzers.reliability_band import ReliabilityBand
from quebra.analyzers.signal_band import SignalBand
from quebra.core._artifact_guard import StaleArtifactGuard


@dataclass
class WithinCalibrationPanelData(StaleArtifactGuard):
    """Complete typed contract for WithinCalibrationPanel.

    Built by build_within_calibration_panel_data (analyzers/within_calibration_compute.py), the
    sole constructor path. Constructing with thresholds present but a band's
    per-threshold maps unpopulated raises, so an incomplete artifact never materializes.

    thresholds    : (label, value, big_values_good) triples, in the display units of
                    `signal.values`. big_values_good=False -> above threshold is
                    out of spec (e.g. infidelity); True -> below is (e.g. T2*).
    primary_label : y-axis label carrying the unit, e.g. "T2* (µs)"
    traces        : optional extra labeled series overlaid on the signal axis
    meta          : dict with a required non-empty "dataset" entry, which the caption
                    shows (the builder raises without it); other keys are free
    """

    signal: SignalBand
    distinguish: DistinguishBand
    reliability: ReliabilityBand
    meta: dict[str, object]

    thresholds: list[tuple[str, float, bool]] = field(default_factory=list)
    primary_label: str = ""
    traces: list[tuple[str, np.ndarray]] | None = None
    use_log_scale: bool = False
    color: object = None
    include_cumulative_time: bool = True
    include_cumulative_damage: bool = True
    include_ttf: bool = True

    def __post_init__(self) -> None:
        # The ladder lives here; the per-threshold maps live in the bands. Neither side
        # can check completeness alone, so the outer class hands the labels down.
        labels = [label for label, _, _ in self.thresholds]
        self.distinguish.check_thresholds(labels)
        self.reliability.check_thresholds(labels)
