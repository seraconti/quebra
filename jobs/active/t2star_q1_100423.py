"""T2* within-calibration job for 100423_6D2S_qubit1 - carved from observed reads only.

One of the two members of the `t2star` family. The graph is
`quebra.recipes.configure_t2star_job`, plus the T2* record's tellability from
`quebra.tellability_recipe`; what is here is this run's parameter row, which is the part
worth diffing against its sibling. Everything shared lives in the recipe: the two
files are byte-identical once the date and the run duration are normalised.
"""

from __future__ import annotations

from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.recipes import configure_t2star_job
from quebra.schemas.track912 import track912Schema
from quebra.tellability_recipe import wire_t2star_tellability

JOB_ID = "t2star_q1_100423"
JOB_FAMILY = "t2star"

PREFIX = "q1_13h_1004_dataset"

job = Job("t2star_q1_100423")

refs = configure_t2star_job(
    job,
    dataset=Dataset(
        path="data/real_private/6D2S/100423_6D2S_qubit1.pickle",
        schema=track912Schema,
        qubit=1,
        device="6D2S",
        duration_h=13,
        extra={"run_name": PREFIX},
    ),
    prefix=PREFIX,
)
wire_t2star_tellability(job, refs, prefix=PREFIX)
