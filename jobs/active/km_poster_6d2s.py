"""Kaplan-Meier in-spec window survival across five 6D2S datasets, at one threshold.

Example job for the poster figure. Five datasets are carved at a SINGLE absolute T2*
threshold and compared; the figure draws the two whose survival curves are furthest
apart, and that pair is chosen by `analyzers.kaplan_meier.compare` in a step, not by eye.

Why one absolute threshold and why 3.0 us. Survival at two different thresholds is
survival of two different events, so a comparison across datasets needs one ladder rung
for all of them. The filtered T2* means of the 6D2S record span 1.62 us (270623 qubit 2)
to 6.19 us (070723 qubit 4), and no single rung covers that whole span: at 3.0 us the
lowest dataset never enters spec at all, and at 2.0 us the highest three sit in spec
essentially forever and produce two or three windows each. 3.0 us is the rung where all
five datasets below carve between 126 and 418 windows, so every curve on the figure is
estimated from real windows rather than from a handful.

The five are spread across the surviving T2* mean range, not picked for their curves:
2.32, 3.08, 4.64, 5.54 and 6.19 us. 270623 qubit 2 is excluded because it has no in-spec
window at this threshold - stated here rather than silently omitted.

Deliberately no interpolate node, for the reason spelled out in
jobs/active/t2star_q1_100423.py: a window must never be carved from a manufactured point.
"""

from __future__ import annotations

import quebra.analyzers.kaplan_meier as kaplan_meier
import quebra.analyzers.t2star as t2star
import quebra.analyzers.windows as windows
from quebra.analyzers.kaplan_meier import KaplanMeierComparison, KaplanMeierCurve
from quebra.analyzers.t2star import T2StarResult
from quebra.analyzers.windows import WindowsResult
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.recipes import RAMSEY_CONFIG, _filter_step, _final_stage
from quebra.plots.km_survival_plot import KMSurvivalPlot
from quebra.schemas.track912 import track912Schema

# The logical name and category. `include` resolves JOB_ID, so this
# file can move without breaking any composite; recategorising costs one string edit.
JOB_ID = "km_poster_6d2s"
JOB_FAMILY = "survival"

PREFIX = "km_poster_6d2s"

# The one rung every curve is carved at. SI seconds, like every threshold in this repo;
# the label is what reaches the figure title and the artifact keys.
THRESHOLD_LABEL = "3.0 µs"
THRESHOLD: list[tuple[str, float, bool]] = [(THRESHOLD_LABEL, 3.0e-6, True)]

# (dataset stem, qubit, display label). The display label is what the legend prints, so
# it is spelled for a reader rather than for the filesystem.
DATASETS: list[tuple[str, int, str]] = [
    ("280623_6D2S_qubit2", 2, "280623 qubit 2"),
    ("040423_6D2S_qubit1", 1, "040423 qubit 1"),
    ("220423_6D2S_qubit1", 1, "220423 qubit 1"),
    ("090623_6D2S_qubit6", 6, "090623 qubit 6"),
    ("070723_6D2S_qubit4", 4, "070723 qubit 4"),
]

job = Job(PREFIX)


def _t2star_run(norm: object) -> T2StarResult:
    return t2star.run(t2star.make_inputs_from_norm(norm))  # type: ignore[arg-type]


def _windows_run(result: T2StarResult, gap_mult: float) -> WindowsResult:
    return windows.run(
        windows.make_inputs_from_frame(
            result.frame,
            time_col="t_rel_s",
            value_col="t2star_s",
            thresholds=THRESHOLD,
            dataset_id=str(result.meta.get("dataset_id", "")),
            gap_mult=gap_mult,
        )
    )


def _km_run(
    window_result: WindowsResult, threshold_label: str, label: str
) -> KaplanMeierCurve:
    return kaplan_meier.run(
        kaplan_meier.make_inputs_from_windows(
            window_result.windows_in_spec,
            threshold_label=threshold_label,
            label=label,
            dataset_id=str(window_result.meta.get("dataset_id", "")),
        )
    )


def _compare(*curves: KaplanMeierCurve, threshold_label: str) -> KaplanMeierComparison:
    """Rank all pairs and select the widest-separated one.

    A step, not a draw-time decision: which two datasets the poster shows is a result,
    and a result belongs in the artifact and on the provenance graph.
    """
    return kaplan_meier.compare(list(curves), threshold_label)


_curves = []
for stem, qubit, display in DATASETS:
    dataset = Dataset(
        path=f"data/real_private/6D2S/{stem}.pickle",
        schema=track912Schema,
        qubit=qubit,
        device="6D2S",
        extra={"run_name": stem},
    )
    _raw = job.load(dataset)
    _filtered = job.step(_filter_step(RAMSEY_CONFIG), _raw, name=f"filter_{stem}")
    _final = job.step(_final_stage, _filtered, name=f"final_filter_stage_{stem}")
    _result = job.step(_t2star_run, _final, name=f"t2star_{stem}")
    _windows = job.step(
        _windows_run,
        _result,
        name=f"windows_{stem}",
        # A step kwarg, not a closure capture, so the gap policy reaches the Mermaid
        # label - the `allan` pattern, not the `filter` pattern.
        gap_mult=10.0,
    )
    _curves.append(
        job.step(
            _km_run,
            _windows,
            name=f"kaplan_meier_{stem}",
            threshold_label=THRESHOLD_LABEL,
            label=display,
        )
    )

_comparison = job.step(
    _compare, *_curves, name="km_compare", threshold_label=THRESHOLD_LABEL
)

# Materialized so the full pairwise ranking - not only the pair that got drawn - is on
# disk beside the figure. The claim "these two differ most" is checkable from it.
job.materialize(_comparison, name=f"{PREFIX}_comparison")
job.figure(
    KMSurvivalPlot,
    _comparison,
    targets=["poster"],
    title=f"{PREFIX}_km_survival",
)
