"""The record figures of every listed record of 6D2S qubit 1, one PDF per block.

Built by `quebra.record_figure_recipe.configure_record_figure_job` on the ledger job
`check_ledger_6d2s_q1`, which it includes and refs once, for its ledger bundle. Each record
the paper's dataset tables list is carved again with that job's own carve kwargs, held to
its ledger, and drawn block by block; a record outside the tables gets no figure, and a
record they shade carries `paper_tables.FLAG_TEXT`.

In `jobs/composite/` and `JOB_SWEEP = False` because a fresh run re-runs the included
ledger job. Run it after the ledger job at a clean commit, with `--reuse-deps`.

Run:  quebra run jobs/composite/record_figures_6d2s_q1.py --reuse-deps
"""

from __future__ import annotations

from quebra.core.job import Job
from quebra.paper_tables import FLAG_TEXT, PAPER_RECORDS, SHADED
from quebra.record_figure_recipe import configure_record_figure_job

JOB_ID = "record_figures_6d2s_q1"
JOB_FAMILY = "record_figure"
JOB_SWEEP = False

job = Job(JOB_ID)

configure_record_figure_job(
    job,
    included=job.include("check_ledger_6d2s_q1"),
    device="6D2S",
    paper_records=PAPER_RECORDS,
    shaded=SHADED,
    flag_text=FLAG_TEXT,
)
