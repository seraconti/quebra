"""Frequency of 2x2 qubit 1, record 200623: Allan, TLF and the fitted frequency's
tellability, with its resolvability tables.

One member of the `freq` family, one job per record of the check-ledger jobs, each the shape
of `ramsey_q1_100423.py` on its own record. The resolvability tables are written as CSVs by
the `table` target (`wire_tellability_tables`); the frequency's tellability has no
thresholds, so the per-threshold table has no row.
"""

from __future__ import annotations

from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.recipes import configure_ramsey_job
from quebra.record_figure_recipe import wire_tellability_tables
from quebra.schemas.track912 import track912Schema
from quebra.tellability_recipe import wire_frequency_tellability

JOB_ID = "freq_2x2_q1_200623"
JOB_FAMILY = "freq"
# Not swept by a bare `run --all`; `run --all --family freq` runs every record.
JOB_SWEEP = False

PREFIX = "freq_2x2_q1_200623"

job = Job(JOB_ID)
main_ds = Dataset(
    path="data/real_private/2x2/200623_2x2_qubit1.pickle",
    schema=track912Schema,
    qubit=1,
    device="2x2",
    extra={"run_name": PREFIX, "run_date": "2023-06-20"},
)

main_node = job.load(main_ds)

refs = configure_ramsey_job(
    job,
    main_node,
    profile="overnight",
    include_fidelity=False,
    include_tlf=True,
    allan_fractional=False,
    figure_prefix=PREFIX,
)
wire_tellability_tables(job, wire_frequency_tellability(job, refs, prefix=PREFIX))
