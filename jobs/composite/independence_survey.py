"""Is the independence assumption safe on this device, or only on the record we looked at?

The device-wide survey: one grid per instrument, every 6D2S record down, the T2* threshold
ladder across, both clocks stacked. One figure per entry of
`independence_survey.SURVEY_KEYS`, which is `battery.ROW_KEYS` plus C3.

It includes every job of the `check_ledger` family and refs each job's `CheckLedgers`
bundle, so it carves nothing and scores nothing: every cell is the verdict that record's own
ledger printed, under that job's knobs (`ledger_recipe.LEDGER_KNOBS`). A single record
rejecting iid could be that record; many records rejecting it is a property of the
instrument or the device. Each ledger job draws
the same grids for its own qubit; this one stacks the qubits.

**C3 is drawn, and is the one instrument with no bench cell.** A green C3 cell means "did
not reject" and nothing stronger; the figure says what it lacks.

In `jobs/composite/` and `JOB_SWEEP = False` because a fresh run re-runs every included
ledger job, once each: one ref per job.

Run:  quebra run jobs/composite/independence_survey.py
"""

from __future__ import annotations

from quebra.core.job import Job
from quebra.ledger_recipe import (
    LEDGER_KNOBS,
    included_ledger_set,
    wire_independence_survey,
)

# The logical name and category. `include` resolves JOB_ID, so this
# file can move without breaking any composite; recategorising costs one string edit.
JOB_ID = "independence_survey"
JOB_FAMILY = "independence"
JOB_SWEEP = False

PREFIX = "independence_survey"

job = Job(JOB_ID)

# Qubit order down the grid, and within a qubit the calendar order of its job's records.
# Written out, not looped, so `docs/JOBS.md` lists them.
_INCLUDED = [
    job.include("check_ledger_6d2s_q1"),
    job.include("check_ledger_6d2s_q2"),
    job.include("check_ledger_6d2s_q3"),
    job.include("check_ledger_6d2s_q4"),
    job.include("check_ledger_6d2s_q5"),
    job.include("check_ledger_6d2s_q6"),
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
