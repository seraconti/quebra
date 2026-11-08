"""The paper's dataset tables, as the figures and tables that flag their shaded rows quote them.

Capannelli et al., Supplementary notes, Supplementary Fig. 10 (6D2S, p. 10) and
Supplementary Fig. 11 (2x2, p. 11) list the released records and shade some of them. A leaf
that imports nothing, so any layer can import it. Jobs import it and pass the flag text, the
record list and the shaded rows to their steps as kwargs. Records are named as the ledger
jobs name them, each ledger's `dataset_id`. A job that names a record otherwise, as the 2x2
Ramsey jobs do with their run names, maps it to its table name
(`infidelity_recipe.paper_record_name`).
"""

from __future__ import annotations

# The text a record those tables shade carries, the same in every figure and table that
# flags it. Its rule: Supplementary notes, p. 9, the paragraph before Supplementary Fig. 10
# (more than 5% of a dataset's points filtered; unsuitable for the Fourier analysis of
# main-text Fig. 8).
FLAG_TEXT = (
    "shaded in the paper's dataset table (high filtered share; the paper excluded >5% "
    "filtered datasets from its Fourier analysis)"
)

TABLE_SOURCE = (
    "Capannelli et al., Supplementary notes, Supplementary Fig. 10 (6D2S, p. 10) and "
    "Supplementary Fig. 11 (2x2, p. 11)"
)

# The paper's dataset tables as (device, record name), table order. Supplementary Fig. 10 for
# 6D2S, 24 rows; Supplementary Fig. 11 for 2x2, 22 rows. All dates are 2023.
PAPER_RECORDS = (
    ("6D2S", "q1_040423"), ("6D2S", "q1_050423"), ("6D2S", "q1_070423"),
    ("6D2S", "q1_100423"), ("6D2S", "q1_120423"), ("6D2S", "q5_120423"),
    ("6D2S", "q6_120423"), ("6D2S", "q1_210423"), ("6D2S", "q2_210423"),
    ("6D2S", "q3_210423"), ("6D2S", "q1_220423"), ("6D2S", "q5_070623"),
    ("6D2S", "q6_070623"), ("6D2S", "q5_080623"), ("6D2S", "q6_080623"),
    ("6D2S", "q5_090623"), ("6D2S", "q6_090623"), ("6D2S", "q2_270623"),
    ("6D2S", "q1_280623"), ("6D2S", "q2_280623"), ("6D2S", "q4_290623"),
    ("6D2S", "q5_290623"), ("6D2S", "q6_290623"), ("6D2S", "q3_110723"),
    ("2x2", "2x2_q1_090623"), ("2x2", "2x2_q2_090623"), ("2x2", "2x2_q3_090623"),
    ("2x2", "2x2_q1_100623"), ("2x2", "2x2_q1_200623"), ("2x2", "2x2_q1_250623"),
    ("2x2", "2x2_q3_260623"), ("2x2", "2x2_q1_300623"), ("2x2", "2x2_q3_300623"),
    ("2x2", "2x2_q1_030723"), ("2x2", "2x2_q2_030723"), ("2x2", "2x2_q2_040723"),
    ("2x2", "2x2_q3_040723"), ("2x2", "2x2_q1_050723"), ("2x2", "2x2_q1_070723"),
    ("2x2", "2x2_q1_100723"), ("2x2", "2x2_q3_110723"), ("2x2", "2x2_q2_120723"),
    ("2x2", "2x2_q1_220823"), ("2x2", "2x2_q2_220823"), ("2x2", "2x2_q3_220823"),
    ("2x2", "2x2_q4_220823"),
)  # fmt: skip
# The rows those tables shade: 3 of the 24 6D2S rows, 11 of the 22 2x2 rows.
SHADED = (
    "q1_280623", "q2_280623", "q3_110723",
    "2x2_q2_090623", "2x2_q3_090623", "2x2_q1_250623", "2x2_q1_300623",
    "2x2_q1_030723", "2x2_q2_030723", "2x2_q2_040723", "2x2_q1_050723",
    "2x2_q1_070723", "2x2_q1_100723", "2x2_q2_120723",
)  # fmt: skip
