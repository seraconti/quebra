"""Infidelity of 2x2 qubit 1, record 070723, at its declared Rabi frequency, and its
tellability, with its resolvability tables.

One member of the `fid` family, one job per 2x2 record of the check-ledger jobs, each the
shape of `ramsey_2x2_q1_030723.py` on its own record: the record's row of the declared table
is joined by device, qubit and `run_date`, and a record with no row raises at the join. The
resolvability tables are written as CSVs by the `table` target (`wire_tellability_tables`).
"""

from __future__ import annotations

from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.recipes import configure_ramsey_job
from quebra.record_figure_recipe import wire_tellability_tables
from quebra.schemas.track912 import track912Schema
from quebra.tellability_recipe import wire_infidelity_tellability

JOB_ID = "fid_2x2_q1_070723"
JOB_FAMILY = "fid"
# Not swept by a bare `run --all`; `run --all --family fid` runs every record.
JOB_SWEEP = False

PREFIX = "fid_2x2_q1_070723"

job = Job(JOB_ID)
main_ds = Dataset(
    path="data/real_private/2x2/070723_2x2_qubit1.pickle",
    schema=track912Schema,
    qubit=1,
    device="2x2",
    extra={"run_name": PREFIX, "run_date": "2023-07-07"},
)
rabi_table = Dataset(
    path="data/real_private/declared/rabi_frequencies_2x2.csv",
    loader_kwargs={"comment": "#"},
)

refs = configure_ramsey_job(
    job,
    main_ds,
    profile="overnight",
    include_fidelity=True,
    rabi_table=rabi_table,
    figure_prefix=PREFIX,
)
wire_tellability_tables(job, wire_infidelity_tellability(job, refs, prefix=PREFIX))
