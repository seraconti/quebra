"""R3: the records of each condition, their restricted means pooled per condition group.

For each condition group, a DerSimonian-Laird summary of its records' T2* restricted means
at every rung of the T2* ladder and on both sides
(`quebra.comparison_recipe.wire_record_forest`). The restriction time is the comparison
tau, `min_largest_age`, taken per group (`per_group`), so a short record lowers tau only in
its own group. A group with fewer than three records with a variance is listed, not
pooled: here that is 2x2 200 mK I, a single record.

The groups are device by temperature period, and on the 2x2 also by voltage configuration.
Their records are the rows of the paper's dataset tables, Supp. Fig. 10 (6D2S, 24) and
Supp. Fig. 11 (2x2, 22):

- 6D2S: 200 mK April (11), 300 mK June (6), 200 mK June-July (7);
- 2x2: 11 mK I (4), 200 mK I (1), 200 mK II (13), 11 mK II (4).

Left out: 070723_6D2S_qubit4, 100723_6D2S_qubit4, 070723_6D2S_qubit5 and
070723_6D2S_qubit6, in the ledgers, not in the paper's dataset tables (Supp. Fig. 10). The
`_before` / `_after` files are row splits of a full record and enter no group.

What a group holds fixed. Voltages are constant within a run and can change from day to day
(Capannelli et al., Ramsey section); the 2x2 configurations are stated in Supp. Fig. S13.
The paper draws no conclusion on mixing-chamber temperature for short drifts and reports
nothing on T2* against either, so a 6D2S group fixes the temperature period only. Each
group's evidence line in the table says which. The unit is the record: the between-record
variance mixes qubits and nights, and on the 6D2S possibly voltage settings.

The suite builds this graph and never runs it: the records are embargoed.
"""

from __future__ import annotations

from quebra.comparison_recipe import (
    record_dataset,
    wire_record_carves,
    wire_record_forest,
)
from quebra.core.job import Job

JOB_ID = "compare_r3_records_one_condition"
JOB_FAMILY = "survival"
# Not swept by a bare `run --all`: it carves every record of both devices.
JOB_SWEEP = False

TAU_RULE = "min_largest_age"
TAU_SCOPE = "per_group"
CONF_LEVEL = 0.95
EVIDENCE_6D2S = (
    "temperature period stated; gate voltages can change between nights and are not "
    "listed per night"
)
EVIDENCE_2X2 = "temperature and voltage configuration stated (Supp. Fig. S13)"

# Every record by its file stem, in the order of the paper's tables.
GROUPS: dict[str, list[str]] = {
    "6D2S 200 mK April": [
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
    ],
    "6D2S 300 mK June": [
        "070623_6D2S_qubit5",
        "070623_6D2S_qubit6",
        "080623_6D2S_qubit5",
        "080623_6D2S_qubit6",
        "090623_6D2S_qubit5",
        "090623_6D2S_qubit6",
    ],
    "6D2S 200 mK June-July": [
        "270623_6D2S_qubit2",
        "280623_6D2S_qubit1",
        "280623_6D2S_qubit2",
        "290623_6D2S_qubit4",
        "290623_6D2S_qubit5",
        "290623_6D2S_qubit6",
        "110723_6D2S_qubit3",
    ],
    "2x2 11 mK I": [
        "090623_2x2_qubit1",
        "090623_2x2_qubit2",
        "090623_2x2_qubit3",
        "100623_2x2_qubit1",
    ],
    "2x2 200 mK I": [
        "200623_2x2_qubit1",
    ],
    "2x2 200 mK II": [
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
    ],
    "2x2 11 mK II": [
        "220823_2x2_qubit1",
        "220823_2x2_qubit2",
        "220823_2x2_qubit3",
        "220823_2x2_qubit4",
    ],
}
EVIDENCE = {
    group: EVIDENCE_6D2S if group.startswith("6D2S ") else EVIDENCE_2X2
    for group in GROUPS
}


job = Job(JOB_ID)
carves = wire_record_carves(
    job, [record_dataset(stem) for stems in GROUPS.values() for stem in stems]
)
wire_record_forest(
    job,
    carves,
    name="r3",
    groups=GROUPS,
    evidence=EVIDENCE,
    tau_scope=TAU_SCOPE,
    tau_rule=TAU_RULE,
    conf_level=CONF_LEVEL,
)
