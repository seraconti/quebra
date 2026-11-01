"""R1: two qubits read in one night, compared at every rung and side.

210423_6D2S_qubit1 against 210423_6D2S_qubit2: the two-sample log-rank test of their T2*
windows and the difference of their restricted means up to the comparison tau
(`min_largest_age`, taken per rung and side over the two records), at every rung of the
T2* ladder and on both sides (`quebra.comparison_recipe.wire_record_pair`). Each p is
unadjusted for the other rungs and sides; any correction over them is made by hand where
the p-values are reported, not here. The figure draws both records' window survival at
every rung and side, each panel with its tau, difference and p, and says in its caption
that each p is unadjusted.

Both records are one run of the 6D2S device, in the paper's dataset table (Supp. Fig. 10),
so one voltage setting holds for both, and the table's evidence line says so. The two
qubits also share that night's environment, which the independence the log-rank variance
assumes does not cover (`analyzers/record_comparison.py`, Validity assumptions).

The pair was declared in the approved plan for issue 23 and confirmed at its jobs stop,
before any comparison ran.

The suite builds this graph and never runs it: the records are embargoed.
"""

from __future__ import annotations

from quebra.comparison_recipe import (
    record_dataset,
    wire_pair_figure,
    wire_record_carves,
    wire_record_pair,
)
from quebra.core.job import Job

JOB_ID = "compare_r1_qubits_one_night"
JOB_FAMILY = "survival"
# Not swept by a bare `run --all`: it carves records that other jobs carve too.
JOB_SWEEP = False

TAU_RULE = "min_largest_age"
CONF_LEVEL = 0.95
EVIDENCE = "one run, one voltage setting"

# Record a, then record b: the difference is a - b.
PAIR = ("210423_6D2S_qubit1", "210423_6D2S_qubit2")


job = Job(JOB_ID)
carves = wire_record_carves(job, [record_dataset(stem) for stem in PAIR])
refs = wire_record_pair(
    job,
    carves,
    name="r1",
    evidence=EVIDENCE,
    tau_rule=TAU_RULE,
    conf_level=CONF_LEVEL,
)
wire_pair_figure(job, carves, refs["pair"], name="r1", evidence=EVIDENCE)
