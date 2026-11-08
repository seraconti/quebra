"""Is the independence assumption safe on the 2x2 infidelity windows, or only on one record?

The survey of the infidelity ledgers: one grid per instrument, every 2x2 record down in qubit
order, the infidelity ladder (99%, 99.9%, 99.99%) across, both clocks stacked. One figure per
entry of `independence_survey.SURVEY_KEYS`: the five permutation rows, and one per group of
`check_ledger.GROUP_KEYS`.

It includes the four jobs of the `check_ledger_fidelity` family and refs each job's
`CheckLedgers` bundle, so it carves nothing and scores nothing: every cell is the verdict that
record's own ledger printed, under that job's knobs (`ledger_recipe.LEDGER_KNOBS`). Whether
the bench's size cells transfer to an |x| carve is not measured, and the ties cutoff is not
calibrated on it (`fidelity_ledger_recipe.BENCH_TRANSFER`). The T2* survey of the same
records is `independence_survey`.

In `jobs/composite/` and `JOB_SWEEP = False` because a fresh run re-runs every included
ledger job, once each: one ref per job. Run it with `--reuse-deps`.

Run:  quebra run jobs/composite/independence_survey_fidelity_2x2.py --reuse-deps
"""

from __future__ import annotations

from quebra.core.job import Job
from quebra.ledger_recipe import (
    LEDGER_KNOBS,
    included_ledger_set,
    wire_independence_survey,
)

JOB_ID = "independence_survey_fidelity_2x2"
JOB_FAMILY = "independence"
JOB_SWEEP = False

PREFIX = "independence_survey_fidelity_2x2"

job = Job(JOB_ID)

# Qubit order down the grid, and within a qubit the order of its job's records. Written
# out, not looped, so `docs/JOBS.md` lists them.
_INCLUDED = [
    job.include("check_ledger_fidelity_2x2_q1"),
    job.include("check_ledger_fidelity_2x2_q2"),
    job.include("check_ledger_fidelity_2x2_q3"),
    job.include("check_ledger_fidelity_2x2_q4"),
]

_labels: list[str] = []
_bundles: list[object] = []
for _included in _INCLUDED:
    _job_labels, _bundle = included_ledger_set(_included)
    _labels.extend(_job_labels)
    _bundles.append(_bundle)

wire_independence_survey(
    job,
    _bundles,
    _labels,
    prefix=PREFIX,
    alpha=LEDGER_KNOBS["alpha"],  # type: ignore[arg-type]
)
