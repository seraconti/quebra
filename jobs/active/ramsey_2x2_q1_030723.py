"""Archetype: 2x2-device Ramsey job - Allan and fidelity via configure_ramsey_job.

Shape: Dataset -> configure_ramsey_job, which loads the record and owns the filter, Allan and
fidelity wiring. The flags below are exactly what runs: configure_ramsey_job does not
override them from the device.

Fidelity uses the record's DECLARED Rabi frequency: one row of the declared table, keyed by
the Dataset's device, qubit and `extra["run_date"]`, each of which must match the file
name (checked at graph build). A record with no row of its own raises when the join runs.
"""

from __future__ import annotations

from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.recipes import configure_ramsey_job
from quebra.schemas.track912 import track912Schema

# The logical name and category. `include` resolves JOB_ID, so this
# file can move without breaking any composite; recategorising costs one string edit.
JOB_ID = "ramsey_2x2_q1_030723"
JOB_FAMILY = "ramsey"

PREFIX = "q1_2x2_16h_0307_dataset"

job = Job("ramsey_2x2_q1_030723")
main_ds = Dataset(
    path="data/real_private/2x2/030723_2x2_qubit1.pickle",
    schema=track912Schema,
    qubit=1,
    device="2x2",
    duration_h=16,
    extra={"run_name": PREFIX, "run_date": "2023-07-03"},
)
rabi_table = Dataset(
    path="data/real_private/declared/rabi_frequencies_2x2.csv",
    loader_kwargs={"comment": "#"},
)

configure_ramsey_job(
    job,
    main_ds,
    profile="overnight",
    include_fidelity=True,
    rabi_table=rabi_table,
    figure_prefix=PREFIX,
)
