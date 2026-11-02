"""J: three records split at a frequency jump, each one's two halves compared.

For 120423_6D2S_qubit1, 270623_6D2S_qubit2 and 280623_6D2S_qubit2, the `_before` half
against the `_after` half: the two-sample log-rank test of their T2* windows and the
difference of their restricted means up to the comparison tau (`min_largest_age`, taken
per rung and side over the two halves), at every rung of the T2* ladder and on both sides
(`quebra.comparison_recipe.wire_record_pair`), each pair with its figure. Each p is
unadjusted for the other rungs, sides and pairs; any correction over them is made by hand
where the p-values are reported, not here.

The halves are row splits of one night's scan around a frequency jump: contiguous parts of
one run, so they share the qubit, the night and the voltage setting, and the table's
evidence line says so. The independence the log-rank variance assumes is weaker here than
between separate nights (`analyzers/record_comparison.py`, Validity assumptions). A record
and its halves never sit in one comparison (`record_comparison.HALF_SUFFIXES`): the full
records are pooled in `compare_r3_records_one_condition` and `compare_r4_two_devices`, the
halves only here.

280623_6D2S_qubit2 is shaded in the paper's dataset table (Supp. Fig. 10), so both its
halves carry a flag that says they are halves of a shaded record.

The pairs were declared in the approved plan for issue 23 and confirmed at its jobs stop,
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
from quebra.paper_tables import FLAG_TEXT

JOB_ID = "compare_j_across_jump"
JOB_FAMILY = "survival"
# Not swept by a bare `run --all`: it carves records that other jobs carve too.
JOB_SWEEP = False

TAU_RULE = "min_largest_age"
CONF_LEVEL = 0.95
EVIDENCE = "two halves of one run, split at a frequency jump: one voltage setting"
HALF_FLAG = f"half of a record {FLAG_TEXT}"

# Each record split at its jump; its halves are compared before - after.
SPLIT = ("120423_6D2S_qubit1", "270623_6D2S_qubit2", "280623_6D2S_qubit2")
SHADED = ("280623_6D2S_qubit2",)

job = Job(JOB_ID)
for record in SPLIT:
    halves = (f"{record}_before", f"{record}_after")
    carves = wire_record_carves(job, [record_dataset(stem) for stem in halves])
    flags = {stem: HALF_FLAG for stem in halves} if record in SHADED else {}
    name = f"j_{record}"
    refs = wire_record_pair(
        job,
        carves,
        name=name,
        evidence=EVIDENCE,
        tau_rule=TAU_RULE,
        conf_level=CONF_LEVEL,
        flags=flags,
    )
    wire_pair_figure(
        job, carves, refs["pair"], name=name, evidence=EVIDENCE, flags=flags
    )
