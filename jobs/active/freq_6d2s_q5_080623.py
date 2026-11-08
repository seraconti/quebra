"""Frequency of 6D2S qubit 5, record 080623: Allan, TLF and the fitted frequency's
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

JOB_ID = "freq_6d2s_q5_080623"
JOB_FAMILY = "freq"
# Not swept by a bare `run --all`; `run --all --family freq` runs every record.
JOB_SWEEP = False

PREFIX = "freq_6d2s_q5_080623"

job = Job(JOB_ID)
main_ds = Dataset(
    path="data/real_private/6D2S/080623_6D2S_qubit5.pickle",
    schema=track912Schema,
    qubit=5,
    device="6D2S",
    extra={"run_name": PREFIX, "run_date": "2023-06-08"},
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
