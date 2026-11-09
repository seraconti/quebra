"""F1: every pair of 2x2 records read in one night, compared on their infidelity windows.

The pairs are a rule, not a choice: every pair of 2x2 records sharing a run date, 12 pairs
over 13 records. Each is compared at every rung of the nines ladder (99%, 99.9%, 99.99%) and
on both sides, as R1 is: the two-sample log-rank test of their windows and the difference of
their restricted means up to the comparison tau (`min_largest_age`, per rung and side over
the two records), by `quebra.comparison_recipe.wire_record_pair`. The pairs were declared
after an exploratory look at these records, and every table and figure says so. Each p is
unadjusted for the other pairs, rungs and sides, and assumes the windows within each record
are independent; any correction is made by hand where the p-values are reported.

The windows and the check ledgers are the infidelity ledger jobs' own
(`check_ledger_fidelity_2x2_q1` to `q4`, read through
`quebra.fidelity_comparison_recipe`): this job carves nothing and scores nothing. Each
record's band outcome (the permutation rows of that rung and bag, both clocks) and its group
verdicts at that rung (2 groups x 2 bags x 2 clocks) sit beside its numbers in the table,
and each panel of the figure prints their tallies. Whether the bench's size cells transfer
to an |x| carve is not measured (`fidelity_ledger_recipe.BENCH_TRANSFER`).

Two records of one night share that night's condition group (temperature and voltage
configuration, Supp. Fig. S13), and the table's evidence line says so. They also share the
night's environment, which the independence the log-rank variance assumes between the two
records does not cover (`analyzers/record_comparison.py`, Validity assumptions).

In `jobs/composite/` and `JOB_SWEEP = False`: every record is two refs into a ledger job,
and a composite never reuses its own nested runs, so each ref would re-run its ledger job.
Run the four `check_ledger_fidelity_2x2_q*` jobs first, on their own, at the same clean
commit; then this job with `--reuse-deps`, which reuses those runs.

Run:  quebra run jobs/composite/compare_f1_2x2_same_night.py --reuse-deps
"""

from __future__ import annotations

from quebra.comparison_recipe import wire_pair_figure, wire_record_pair
from quebra.core.job import Job
from quebra.fidelity_comparison_recipe import fidelity_records
from quebra.paper_tables import FLAG_TEXT, SHADED

JOB_ID = "compare_f1_2x2_same_night"
JOB_FAMILY = "survival"
# Not swept by a bare `run --all`: it re-runs the infidelity ledger jobs it includes.
JOB_SWEEP = False

TAU_RULE = "min_largest_age"
CONF_LEVEL = 0.95
EVIDENCE = (
    "one night of the 2x2: one condition group (temperature and voltage configuration, "
    "Supp. Fig. S13); two qubits"
)
DECLARED = (
    "Pairs declared after an exploratory look at these records: every pair of 2x2 "
    "records sharing a night, 12 pairs at 3 rungs and 2 sides, so up to 72 log-rank p, "
    "each unadjusted for the others"
)

# Every pair of 2x2 records sharing a run date, by run name; record a, then record b, so
# the difference is a - b.
PAIRS = (
    ("fid_2x2_q1_090623", "fid_2x2_q2_090623"),
    ("fid_2x2_q1_090623", "fid_2x2_q3_090623"),
    ("fid_2x2_q2_090623", "fid_2x2_q3_090623"),
    ("fid_2x2_q1_300623", "fid_2x2_q3_300623"),
    ("fid_2x2_q1_030723", "fid_2x2_q2_030723"),
    ("fid_2x2_q2_040723", "fid_2x2_q3_040723"),
    ("fid_2x2_q1_220823", "fid_2x2_q2_220823"),
    ("fid_2x2_q1_220823", "fid_2x2_q3_220823"),
    ("fid_2x2_q1_220823", "fid_2x2_q4_220823"),
    ("fid_2x2_q2_220823", "fid_2x2_q3_220823"),
    ("fid_2x2_q2_220823", "fid_2x2_q4_220823"),
    ("fid_2x2_q3_220823", "fid_2x2_q4_220823"),
)

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

for a, b in PAIRS:
    name = f"f1_{a.removeprefix('fid_2x2_')}_vs_{b.removeprefix('fid_2x2_')}"
    carves = {a: RECORDS[a].carve, b: RECORDS[b].carve}
    flags = {r: FLAG_TEXT for r in (a, b) if RECORDS[r].paper_name in SHADED}
    refs = wire_record_pair(
        job,
        carves,
        name=name,
        evidence=EVIDENCE,
        tau_rule=TAU_RULE,
        conf_level=CONF_LEVEL,
        flags=flags,
        ledgers={a: RECORDS[a].ledger, b: RECORDS[b].ledger},
        declared=DECLARED,
    )
    wire_pair_figure(
        job,
        carves,
        refs["pair"],
        name=name,
        evidence=EVIDENCE,
        flags=flags,
        declared=DECLARED,
    )
