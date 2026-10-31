"""Archetype: across-calibration job - the MTBC of a calibration-event log.

The reliability literature calls this tier repairable-system analysis. This project
names it after the calibration boundary; see panels/across_calibration.py for why, and
for where the literature's terms are deliberately kept.

The only archetype that is NOT a degradation trace. It reads a calibration log through
CalibrationLogSchema (one row per event, no T2*/frequency series), derives inter-event
intervals, and renders AcrossCalibrationPanel. No thresholds, no window carving, no gap policy:
those belong to the within-calibration side, where a metric degrades between
observations.
"""

from __future__ import annotations

import quebra.analyzers.mtbf as mtbf
from quebra.analyzers.mtbf import MtbfResult
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.core.types import Norm
from quebra.panels.across_calibration import (
    AcrossCalibrationPanel,
    AcrossCalibrationPanelData,
    make_mtbf_panel_data,
)
from quebra.schemas.calibration_log import CalibrationLogSchema

# The logical name and category. `include` resolves JOB_ID, so this
# file can move without breaking any composite; recategorising costs one string edit.
JOB_ID = "mtbf_q1"
JOB_FAMILY = "interval"

job = Job("mtbf_q1")

ds = Dataset(
    path="data/real_private/calibration_logs/6D2S_qubit1_freq_log.pickle",
    schema=CalibrationLogSchema,
    qubit=1,
    device="6D2S",
    extra={"run_name": "mtbf_q1"},
)

norm = job.load(ds)


def _run_mtbf(n: Norm) -> MtbfResult:
    return mtbf.run(mtbf.make_inputs_from_norm(n))


def _make_panel_data(
    result: MtbfResult, show_dropped_count: bool
) -> AcrossCalibrationPanelData:
    return make_mtbf_panel_data(result, show_dropped_count=show_dropped_count)


mtbf_result = job.step(_run_mtbf, norm, name="mtbf")
panel_data = job.step(
    _make_panel_data,
    mtbf_result,
    name="mtbf_panel_data",
    # Declared rather than left to the default, so the choice reaches the provenance
    # label: True draws the count of 0-length intervals on the histogram.
    show_dropped_count=False,
)

job.materialize(panel_data, name="mtbf_q1_panel_data")
job.figure(
    AcrossCalibrationPanel, panel_data, targets=["static", "academic"], title="mtbf_q1"
)
