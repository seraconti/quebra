"""F3: the 2x2 records of each condition, their infidelity restricted means pooled per group.

R3's forest on the 2x2 infidelity windows. For each of R3's four 2x2 condition groups, a
DerSimonian-Laird summary of its records' restricted means at every rung of the nines ladder
(99%, 99.9%, 99.99%) and on both sides (`quebra.comparison_recipe.wire_record_forest`). The
restriction time is the comparison tau, `min_largest_age`, taken per group (`per_group`). A
group with fewer than three records with a variance at a threshold and side is listed, not
pooled there: 2x2 200 mK I, a single record, is never pooled, and a larger group can fall
below three at some rung or side.

The groups are R3's: device by temperature period and voltage configuration, the rows of the
paper's 2x2 dataset table (Supp. Fig. 11, 22 records), 11 mK I (4), 200 mK I (1), 200 mK II
(13), 11 mK II (4), with the configurations of Supp. Fig. S13. The 11 mK II group also has the
~2.2 MHz drive, and its evidence line says so. The unit is the record: the between-record
variance mixes qubits and nights.

The windows and the check ledgers are the infidelity ledger jobs' own
(`check_ledger_fidelity_2x2_q1` to `q4`, read through
`quebra.fidelity_comparison_recipe`): this job carves nothing and scores nothing. Each
record's band outcome and its group verdicts at the rung sit beside its numbers in the table.
Whether the bench's size cells transfer to an |x| carve is not measured
(`fidelity_ledger_recipe.BENCH_TRANSFER`).

Records the paper's 2x2 table shades carry `paper_tables.FLAG_TEXT`. Two figures: the forest
at 99% on the in-spec side, and the pooled ladder over every rung and both sides; with more
than one group, the caption says the groups are described, not tested. The ladder is the
ledgers' own, read off the ledger jobs' graphs; this job declares none.

In `jobs/composite/` and `JOB_SWEEP = False`: every record is two refs into a ledger job,
and a composite never reuses its own nested runs, so each ref would re-run its ledger job.
Run the four `check_ledger_fidelity_2x2_q*` jobs first, on their own, at the same clean
commit; then this job with `--reuse-deps`, which reuses those runs.

Run:  quebra run jobs/composite/compare_f3_2x2_conditions.py --reuse-deps
"""

from __future__ import annotations

from quebra.comparison_recipe import wire_record_forest
from quebra.core.job import Job
from quebra.fidelity_comparison_recipe import (
    fidelity_records,
    wire_fidelity_forest_figures,
)
from quebra.paper_tables import FLAG_TEXT, SHADED

JOB_ID = "compare_f3_2x2_conditions"
JOB_FAMILY = "survival"
# Not swept by a bare `run --all`: it re-runs the infidelity ledger jobs it includes.
JOB_SWEEP = False

TAU_RULE = "min_largest_age"
TAU_SCOPE = "per_group"
CONF_LEVEL = 0.95
# The rung and side the forest figure draws; the pooled ladder draws every one.
FOREST_THRESHOLD = "99%"
FOREST_SIDE = "in_spec"
EVIDENCE_2X2 = "temperature and voltage configuration stated (Supp. Fig. S13)"

# R3's 2x2 groups, each record by its run name, in the order of the paper's table.
GROUPS: dict[str, list[str]] = {
    "2x2 11 mK I": [
        "fid_2x2_q1_090623",
        "fid_2x2_q2_090623",
        "fid_2x2_q3_090623",
        "fid_2x2_q1_100623",
    ],
    "2x2 200 mK I": [
        "fid_2x2_q1_200623",
    ],
    "2x2 200 mK II": [
        "fid_2x2_q1_250623",
        "fid_2x2_q3_260623",
        "fid_2x2_q1_300623",
        "fid_2x2_q3_300623",
        "fid_2x2_q1_030723",
        "fid_2x2_q2_030723",
        "fid_2x2_q2_040723",
        "fid_2x2_q3_040723",
        "fid_2x2_q1_050723",
        "fid_2x2_q1_070723",
        "fid_2x2_q1_100723",
        "fid_2x2_q3_110723",
        "fid_2x2_q2_120723",
    ],
    "2x2 11 mK II": [
        "fid_2x2_q1_220823",
        "fid_2x2_q2_220823",
        "fid_2x2_q3_220823",
        "fid_2x2_q4_220823",
    ],
}
EVIDENCE = {
    group: (
        f"{EVIDENCE_2X2}; this group also has the ~2.2 MHz drive"
        if group == "2x2 11 mK II"
        else EVIDENCE_2X2
    )
    for group in GROUPS
}

job = Job(JOB_ID)

# Written out, not looped, so `docs/JOBS.md` lists every include.
RECORDS = fidelity_records(
    [
        job.include("check_ledger_fidelity_2x2_q1"),
        job.include("check_ledger_fidelity_2x2_q2"),
        job.include("check_ledger_fidelity_2x2_q3"),
        job.include("check_ledger_fidelity_2x2_q4"),
    ]
)
FLAGS = {
    record: FLAG_TEXT
    for records in GROUPS.values()
    for record in records
    if RECORDS[record].paper_name in SHADED
}

refs = wire_record_forest(
    job,
    {
        record: RECORDS[record].carve
        for records in GROUPS.values()
        for record in records
    },
    name="f3",
    groups=GROUPS,
    evidence=EVIDENCE,
    tau_scope=TAU_SCOPE,
    tau_rule=TAU_RULE,
    conf_level=CONF_LEVEL,
    flags=FLAGS,
    ledgers={
        record: RECORDS[record].ledger
        for records in GROUPS.values()
        for record in records
    },
)
wire_fidelity_forest_figures(
    job,
    refs["forest"],
    list(RECORDS.values()),
    name="f3",
    evidence=EVIDENCE,
    threshold_label=FOREST_THRESHOLD,
    side=FOREST_SIDE,
    flags=FLAGS,
)
