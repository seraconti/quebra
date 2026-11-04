"""How much of each record's observed time is T2* in spec, at every rung of the ladder?

The compliance survey: for every record of both devices, the time-weighted in-spec fraction
at each rung of the T2* ladder, with its lower and upper bounds
(`analyzers/compliance_survey.py`). It includes every job of the `check_ledger` family and
refs each job's carve bundle (`ledger_recipe.RecordCarves`, the carves the job also
persists one by one as `{record}_windows`), so it carves nothing: every number is read off
the carve that record's ledger scored.

The figure draws the records of the paper's dataset tables, one row per record, 6D2S and
2x2 side by side; a record the paper shades carries `paper_tables.FLAG_TEXT`. A ledger
record outside the tables is computed, kept on the artifact and listed under its panel, not
drawn. A record that failed to load or carve reads `not computed` with its reason.

In `jobs/composite/` and `JOB_SWEEP = False` because a fresh run re-runs every included
ledger job, once each: one ref per job.

Run:  quebra run jobs/composite/compliance_survey.py
"""

from __future__ import annotations

from quebra.analyzers.compliance_survey import build_compliance_survey
from quebra.core.job import Job
from quebra.ledger_recipe import included_carves
from quebra.paper_tables import FLAG_TEXT, PAPER_RECORDS, SHADED, TABLE_SOURCE
from quebra.plots.compliance_survey_plot import ComplianceSurveyPlot

JOB_ID = "compliance_survey"
JOB_FAMILY = "compliance"
JOB_SWEEP = False

PREFIX = "compliance_survey"

job = Job(JOB_ID)

# Device, then qubit order, and within a qubit the calendar order of its job's records.
# Written out, not looped, so `docs/JOBS.md` lists them.
_INCLUDED = [
    job.include("check_ledger_6d2s_q1"),
    job.include("check_ledger_6d2s_q2"),
    job.include("check_ledger_6d2s_q3"),
    job.include("check_ledger_6d2s_q4"),
    job.include("check_ledger_6d2s_q5"),
    job.include("check_ledger_6d2s_q6"),
    job.include("check_ledger_2x2_q1"),
    job.include("check_ledger_2x2_q2"),
    job.include("check_ledger_2x2_q3"),
    job.include("check_ledger_2x2_q4"),
]

_records: list[str] = []
_devices: list[str] = []
_bundles: list[object] = []
for _included in _INCLUDED:
    _job_records, _job_devices, _bundle = included_carves(_included)
    _records.extend(_job_records)
    _devices.extend(_job_devices)
    _bundles.append(_bundle)

_survey = job.step(
    build_compliance_survey,
    *_bundles,  # type: ignore[arg-type]
    name="compliance_survey",
    records=tuple(_records),
    devices=tuple(_devices),
    paper_records=PAPER_RECORDS,
    shaded=SHADED,
    flag_text=FLAG_TEXT,
    table_source=TABLE_SOURCE,
)
job.materialize(_survey, name=f"{PREFIX}_data")
job.figure(
    ComplianceSurveyPlot,
    _survey,
    targets=["static", "academic"],
    title="In-spec fraction on the T2* ladder, every record",
)
