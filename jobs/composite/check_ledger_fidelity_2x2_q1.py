"""Check ledger of every infidelity record of 2x2 qubit 1, and their survey.

One member of the `check_ledger_fidelity` family, built by
`quebra.fidelity_ledger_recipe.configure_fidelity_ledger_job`: per record, the infidelity
carve of its Ramsey job with fidelity (`fid_2x2_q1_<DDMMYY>`), recomputed on the result that
job persists with its own carve kwargs, and the check ledger of those windows with the T2*
ledgers' knobs and bench size table; then one survey over the records. The records are those
of `check_ledger_2x2_q1`, in its order: the 2x2 drift datasets of Capannelli et al.,
Supplementary notes, Supplementary Fig. 11 (p. 11). Each record is named by its Ramsey job's
run name, with its name in the paper's tables beside it.

A ledger that refuses a record's content reads `not computed` with the reason. A failure
inside an included Ramsey job stops this job (the recipe's module docstring). Whether the
bench's size cells transfer to an |x| carve is not measured, and the ties cutoff is not
calibrated on it; every ledger and figure says so.

In `jobs/composite/` and `JOB_SWEEP = False` because a fresh run re-runs every included
Ramsey job. Run it after those jobs at a clean commit, with `--reuse-deps`.

Run:  quebra run jobs/composite/check_ledger_fidelity_2x2_q1.py --reuse-deps
"""

from __future__ import annotations

from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.fidelity_ledger_recipe import configure_fidelity_ledger_job
from quebra.ledger_recipe import LEDGER_KNOBS

JOB_ID = "check_ledger_fidelity_2x2_q1"
JOB_FAMILY = "check_ledger_fidelity"
JOB_SWEEP = False

PREFIX = "fid_2x2_q1"
BENCH_SIZE_TABLE = Dataset(path="jobs/bench/results/size_table.csv", schema=None)

job = Job(JOB_ID)

# One Ramsey job with fidelity per record, in `check_ledger_2x2_q1`'s order. Written out, not
# looped, so `docs/JOBS.md` lists every include.
INCLUDED = [
    job.include("fid_2x2_q1_090623"),
    job.include("fid_2x2_q1_100623"),
    job.include("fid_2x2_q1_200623"),
    job.include("fid_2x2_q1_250623"),
    job.include("fid_2x2_q1_300623"),
    job.include("fid_2x2_q1_030723"),
    job.include("fid_2x2_q1_050723"),
    job.include("fid_2x2_q1_070723"),
    job.include("fid_2x2_q1_100723"),
    job.include("fid_2x2_q1_220823"),
]

configure_fidelity_ledger_job(
    job,
    prefix=PREFIX,
    included=INCLUDED,
    bench_size_table=BENCH_SIZE_TABLE,
    **LEDGER_KNOBS,  # type: ignore[arg-type]
)
