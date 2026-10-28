"""What the calibration bench established, as five figures.

This job is the one place the bench's numbers reach a figure, and it does so WITHOUT
importing `bench/`. The two result tables are declared as `Dataset`s and loaded with
`job.load_df`, so their sha256 lands in `dataset_hashes` and in the run identity: the
figures cannot silently be built against a different bench, and
`tests/test_bench_isolation.py` still holds because the dependency runs through provenance
rather than through the import graph.

`jobs/bench/results/*.csv` is inside the repository, not under the dataset root, which is why
`core.paths.resolve_dataset_path` grows a repo-root fallback - see its docstring.

The validation figure re-simulates rather than reshaping a table: a P-P plot needs the
statistics themselves and the tables carry only rejection rates. Its `seed` and replicate
count are step kwargs, so both appear on the provenance label and the figure is
reproducible from the record alone.
"""

from __future__ import annotations

from quebra.analyzers.calibration_summary import (
    power_vs_dependence,
    power_vs_trend,
    read_dependence,
    size_vs_n,
    validation_curve,
)
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.plots.calibration_plot import (
    CalibrationPowerPlot,
    CalibrationReadDependencePlot,
    CalibrationSizePlot,
    CalibrationTrendPlot,
    CalibrationValidationPlot,
)

# The logical name and category. `include` resolves JOB_ID, so this
# file can move without breaking any composite; recategorising costs one string edit.
JOB_ID = "check_calibration"
JOB_FAMILY = "validation"

PREFIX = "check_calibration"

# Declared here, not defaulted in the step, so both reach the provenance label.
VALIDATION_SEED = 20260811
VALIDATION_REPLICATES = 4000
VALIDATION_N = (20, 50, 355)

SIZE_TABLE = Dataset(path="jobs/bench/results/size_table.csv", schema=None)
POWER_TABLE = Dataset(path="jobs/bench/results/power_table.csv", schema=None)

job = Job(PREFIX)

_size_table = job.load_df(SIZE_TABLE)
_power_table = job.load_df(POWER_TABLE)

_size = job.step(size_vs_n, _size_table, name="size_vs_n")
_power = job.step(power_vs_dependence, _power_table, name="power_vs_dependence")
_trend = job.step(power_vs_trend, _power_table, _size_table, name="power_vs_trend")
_reads = job.step(read_dependence, _power_table, name="read_dependence")
_validation = job.step(
    validation_curve,
    name="validation_curve",
    n_values=VALIDATION_N,
    n_replicates=VALIDATION_REPLICATES,
    seed=VALIDATION_SEED,
)

job.figure(
    CalibrationSizePlot,
    _size,
    targets=["static", "academic"],
    title=f"{PREFIX} size",
)
job.figure(
    CalibrationPowerPlot,
    _power,
    targets=["static", "academic"],
    title=f"{PREFIX} power",
)
job.figure(
    CalibrationTrendPlot,
    _trend,
    targets=["static", "academic"],
    title=f"{PREFIX} trend",
)
job.figure(
    CalibrationValidationPlot,
    _validation,
    targets=["static", "academic"],
    title=f"{PREFIX} validation",
)
job.figure(
    CalibrationReadDependencePlot,
    _reads,
    targets=["static", "academic"],
    title=f"{PREFIX} read dependence",
)
