"""The record figures of the 2x2 infidelity record 2x2_q1_030723, one PDF per block.

Built by `quebra.record_figure_recipe.configure_infidelity_record_figure_job` on the Ramsey
job `ramsey_2x2_q1_030723`, which it includes and refs once, for its fidelity result. The
record is carved again with that job's own fidelity carve kwargs (|x| on the nines ladder),
scored by the check battery on those windows against the bench's size table, and drawn
block by block, every block of `record_figure_recipe.BLOCKS`. The paper's tables shade it,
so it carries `paper_tables.FLAG_TEXT`.

In `jobs/composite/` and `JOB_SWEEP = False` because a fresh run re-runs the included
Ramsey job. Run it after that job at a clean commit, with `--reuse-deps`.

Run:  quebra run jobs/composite/record_figures_2x2_q1_030723_infidelity.py --reuse-deps
"""

from __future__ import annotations

from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.paper_tables import FLAG_TEXT, PAPER_RECORDS, SHADED
from quebra.record_figure_recipe import configure_infidelity_record_figure_job

JOB_ID = "record_figures_2x2_q1_030723_infidelity"
JOB_FAMILY = "record_figure"
JOB_SWEEP = False

BENCH_SIZE_TABLE = Dataset(path="jobs/bench/results/size_table.csv", schema=None)

job = Job(JOB_ID)
configure_infidelity_record_figure_job(
    job,
    included=job.include("ramsey_2x2_q1_030723"),
    record="2x2_q1_030723",
    paper_records=PAPER_RECORDS,
    shaded=SHADED,
    flag_text=FLAG_TEXT,
    bench_size_table=BENCH_SIZE_TABLE,
)
