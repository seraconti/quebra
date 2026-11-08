"""A2 across the two Quantum Inspire devices: Tuna-17's calibration cycles against Tuna-9's.

It includes `cycles_tuna17` and `cycles_tuna9` and refs each one's cycles and cycles figure
data, so it builds no cycle and scores nothing: each device's cycles and check outcome are
its own job's. `cycle_comparison.compare_cycle_pair` puts the two through
`record_comparison.compare_pair` with no test between the devices: the shared tau, each
device's restricted mean to it with its own interval, and the point difference. The
device that sets the shared tau is restricted where its MTBC is, so its restricted mean is
its MTBC. The other device's is smaller than its MTBC, being restricted short of its own
largest cycle, unless the two largest cycles are equal (`cycle_comparison`). The table is `record_comparison.pair_table`'s, materialized; nothing
is drawn.

In `jobs/composite/` and `JOB_SWEEP = False` because a fresh run re-runs both included jobs.

Run:  quebra run jobs/composite/compare_a2_tuna17_tuna9.py
"""

from __future__ import annotations

from quebra.analyzers.cycle_comparison import compare_cycle_pair
from quebra.analyzers.record_comparison import pair_table
from quebra.core.job import Job

JOB_ID = "compare_a2_tuna17_tuna9"
JOB_FAMILY = "interval"
JOB_SWEEP = False

PREFIX = JOB_ID
TAU_RULE = "min_largest_age"
CONF_LEVEL = 0.95
COMPARISON = "A2 cycles: Tuna-17 vs Tuna-9, restricted means to the shared tau"
EVIDENCE = (
    "two Quantum Inspire devices, each read from its own status export by one cycle rule "
    "(analyzers/cycles.py); nothing else is held fixed between them"
)

job = Job(JOB_ID)

_TUNA17 = job.include("cycles_tuna17")
_TUNA9 = job.include("cycles_tuna9")

pair = job.step(
    compare_cycle_pair,
    _TUNA17.ref("cycles_tuna17_cycles"),
    _TUNA9.ref("cycles_tuna9_cycles"),
    _TUNA17.ref("cycles_tuna17_panel_data"),
    _TUNA9.ref("cycles_tuna9_panel_data"),
    name="a2_pair",
    tau_rule=TAU_RULE,
    conf_level=CONF_LEVEL,
)
table = job.step(
    pair_table, pair, name="a2_pair_table", comparison=COMPARISON, evidence=EVIDENCE
)
job.materialize(table, name=f"{PREFIX}_table")
