"""Archetype: composite job - overlay two sub-jobs' T2* panel data in one figure.

Composites live in jobs/composite/, not jobs/active/: main.py globs jobs/active for
`run --all`, so a composite would otherwise re-run its sub-jobs on every sweep. Run it
by path instead:

    quebra run jobs/composite/compare_t2star_0704_vs_1004.py
    quebra run jobs/composite/compare_t2star_0704_vs_1004.py --reuse-deps

`include` pulls each sub-job in and runs it through the normal runner, so its datasets
and provenance resolve exactly as a standalone run. Every figure sink's input is always
persisted, so `ref("t2star_panel_data")` reads that artifact back with no edit to the
sub-jobs. The two references become a two-input step; the composite's own provenance
graph click-throughs to each sub-job's graph.

This compares one qubit across two dates rather than two qubits - the mechanism is the
same, and both sub-jobs already exist as archetypes.
"""

from __future__ import annotations

import numpy as np

from quebra.core.job import Job
from quebra.panels.comparison import CompareSeriesData, CompareSeriesPanel
from quebra.panels.within_calibration import WithinCalibrationPanelData

# The logical name and category. `include` resolves JOB_ID, so this
# file can move without breaking any composite; recategorising costs one string edit.
JOB_ID = "compare_t2star_0704_vs_1004"
JOB_FAMILY = "t2star"
# Not swept by a bare `run --all`: it re-runs sub-jobs and/or is long. Selectable
# with `--family t2star` or by path. This is the declaration that replaced the old
# "jobs/composite/ is not swept" directory rule.
JOB_SWEEP = False


def _compare_t2star(*panels: WithinCalibrationPanelData) -> CompareSeriesData:
    """Overlay each sub-job's T2* panel data as a labeled (t_h, values) series.

    Reads the SIGNAL band: after the band split the raw series belongs to band 1, and
    the composite is the adapter that knows how to pull it out.
    """
    series: list[tuple[str, np.ndarray, np.ndarray]] = []
    for pd_ in panels:
        series.append(
            (
                str(pd_.meta["dataset"]),
                np.asarray(pd_.signal.t_h, dtype=float),
                np.asarray(pd_.signal.values, dtype=float),
            )
        )
    return CompareSeriesData(
        series=series,
        x_label="Elapsed time (h)",
        y_label="T2* (µs)",
        title="T2* comparison: qubit 1, 070423 vs 100423",
    )


job = Job("compare_t2star_0704_vs_1004")
d0704 = job.include("t2star_q1_070423", alias="d0704")
d1004 = job.include("t2star_q1_100423", alias="d1004")

_cmp = job.step(
    _compare_t2star,
    d0704.ref("t2star_panel_data"),
    d1004.ref("t2star_panel_data"),
    name="t2star_compare",
)
job.figure(
    CompareSeriesPanel,
    _cmp,
    targets=["static", "academic"],
    title="t2star 0704 vs 1004",
)
