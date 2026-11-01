"""R4: the two devices side by side, each a forest of its records at one shared tau.

One DerSimonian-Laird summary per device of its records' T2* restricted means, at every
rung of the T2* ladder and on both sides (`quebra.comparison_recipe.wire_record_forest`).
Tau is shared by both devices (`across_groups`), so their two forests stand on one axis at
one restriction time. No test compares the devices: they are described, not tested.
Under the shared tau, the record with the smallest largest age on either device sets tau
for both, which can leave more records with no death before it, out of the pool
(`analyzers/record_comparison.random_effects_set`, Validity assumptions).

The records are every row of the paper's dataset tables, Supp. Fig. 10 (6D2S, 24) and
Supp. Fig. 11 (2x2, 22): the records `compare_r3_records_one_condition` splits into its
condition groups. Left out: 070723_6D2S_qubit4, 100723_6D2S_qubit4, 070723_6D2S_qubit5
and 070723_6D2S_qubit6, in the ledgers, not in the paper's dataset tables (Supp. Fig. 10).
The `_before` / `_after` files are row splits of a full record and enter no group.

What a device group holds fixed: the device alone. Voltages are constant within a run and
can change from day to day (Capannelli et al., Ramsey section); the 2x2 configurations are
stated in Supp. Fig. S13; the paper draws no conclusion on mixing-chamber temperature for
short drifts and reports nothing on T2* against either. The unit is the record: the
between-record variance mixes qubits, nights and temperature periods, the 2x2 voltage
configurations, and on the 6D2S possibly voltage settings.

Two figures: the forest at one rung and side (`FOREST_THRESHOLD`, `FOREST_SIDE`), and the
pooled ladder over every rung and both sides. Each prints the tau and the record that set it,
the records left out and why, the evidence lines and the check outcome; with more than one
group, the caption says the groups are described, not tested.

The suite builds this graph and never runs it: the records are embargoed.
"""

from __future__ import annotations

from quebra.comparison_recipe import (
    record_dataset,
    wire_forest_figures,
    wire_record_carves,
    wire_record_forest,
)
from quebra.core.job import Job

JOB_ID = "compare_r4_two_devices"
JOB_FAMILY = "survival"
# Not swept by a bare `run --all`: it carves every record of both devices.
JOB_SWEEP = False

TAU_RULE = "min_largest_age"
TAU_SCOPE = "across_groups"
CONF_LEVEL = 0.95
# The rung and side the forest figure draws; the pooled ladder draws every one.
FOREST_THRESHOLD = "3 µs"
FOREST_SIDE = "in_spec"

# Every record by its file stem, in the order of the paper's tables.
GROUPS: dict[str, list[str]] = {
    "6D2S": [
        "040423_6D2S_qubit1",
        "050423_6D2S_qubit1",
        "070423_6D2S_qubit1",
        "100423_6D2S_qubit1",
        "120423_6D2S_qubit1",
        "120423_6D2S_qubit5",
        "120423_6D2S_qubit6",
        "210423_6D2S_qubit1",
        "210423_6D2S_qubit2",
        "210423_6D2S_qubit3",
        "220423_6D2S_qubit1",
        "070623_6D2S_qubit5",
        "070623_6D2S_qubit6",
        "080623_6D2S_qubit5",
        "080623_6D2S_qubit6",
        "090623_6D2S_qubit5",
        "090623_6D2S_qubit6",
        "270623_6D2S_qubit2",
        "280623_6D2S_qubit1",
        "280623_6D2S_qubit2",
        "290623_6D2S_qubit4",
        "290623_6D2S_qubit5",
        "290623_6D2S_qubit6",
        "110723_6D2S_qubit3",
    ],
    "2x2": [
        "090623_2x2_qubit1",
        "090623_2x2_qubit2",
        "090623_2x2_qubit3",
        "100623_2x2_qubit1",
        "200623_2x2_qubit1",
        "250623_2x2_qubit1",
        "260623_2x2_qubit3",
        "300623_2x2_qubit1",
        "300623_2x2_qubit3",
        "030723_2x2_qubit1",
        "030723_2x2_qubit2",
        "040723_2x2_qubit2",
        "040723_2x2_qubit3",
        "050723_2x2_qubit1",
        "070723_2x2_qubit1",
        "100723_2x2_qubit1",
        "110723_2x2_qubit3",
        "120723_2x2_qubit2",
        "220823_2x2_qubit1",
        "220823_2x2_qubit2",
        "220823_2x2_qubit3",
        "220823_2x2_qubit4",
    ],
}
EVIDENCE = {
    "6D2S": (
        "device only: three temperature periods; voltage settings can change between "
        "nights and are not listed per night"
    ),
    "2x2": (
        "device only: four temperature and voltage conditions, each stated "
        "(Supp. Fig. S13)"
    ),
}


job = Job(JOB_ID)
carves = wire_record_carves(
    job, [record_dataset(stem) for stems in GROUPS.values() for stem in stems]
)
refs = wire_record_forest(
    job,
    carves,
    name="r4",
    groups=GROUPS,
    evidence=EVIDENCE,
    tau_scope=TAU_SCOPE,
    tau_rule=TAU_RULE,
    conf_level=CONF_LEVEL,
)
wire_forest_figures(
    job,
    refs["forest"],
    name="r4",
    evidence=EVIDENCE,
    threshold_label=FOREST_THRESHOLD,
    side=FOREST_SIDE,
)
