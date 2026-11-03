"""Calibration cycles of the Quantum Inspire device Tuna-17, from its status history.

Built by `quebra.cycles_recipe.configure_qi_cycles_job`: the status samples become spells,
the spells cycles and status shares, the cycles get the check battery within each stretch
and, as a diagnostic, over the whole record, and `CalibrationCyclesPanel` draws the
survival, the MTBC, the check outcome and the shares. The rules (a cycle, a gap, a
duplicate timestamp, the observed clock) are in `analyzers/cycles.py`.

No `device` is declared on the Dataset: the schema reads the device from the file's `name`
column and refuses a file that names more than one.
"""

from __future__ import annotations

from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.cycles_recipe import CYCLE_KNOBS, configure_qi_cycles_job
from quebra.schemas.qi_status_history import QiStatusHistorySchema

JOB_ID = "cycles_tuna17"
JOB_FAMILY = "interval"

PREFIX = "cycles_tuna17"
# The poller's spacing between status samples, as the export declares it.
SPACING_S = 30.0
# This job's own permutation seed, so its Monte Carlo error is independent of Tuna-9's.
SEED = 20261101

job = Job(PREFIX)

configure_qi_cycles_job(
    job,
    prefix=PREFIX,
    dataset=Dataset(
        path="data/real_private/qinspire/Tuna-17_status_history1.csv",
        schema=QiStatusHistorySchema,
        extra={"run_name": "Tuna-17"},
    ),
    bench_size_table=Dataset(path="jobs/bench/results/size_table.csv", schema=None),
    spacing_s=SPACING_S,
    seed=SEED,
    **CYCLE_KNOBS,  # type: ignore[arg-type]
)
