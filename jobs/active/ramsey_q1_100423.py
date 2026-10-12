"""Archetype: 6D2S Ramsey job - Allan and TLF via configure_ramsey_job.

One of the representative jobs kept in jobs/active. Iterate on this file, then fan it out.

Shape: load -> configure_ramsey_job, which owns the filter, Allan and TLF wiring. Every
analysis reads the filtered reads themselves. The threshold ladder and window carving belong
to the t2star archetype, not here.

No fidelity: an infidelity needs the Rabi frequency of the record, and the 6D2S records carry
no declared one. Fractional Allan is off for the same kind of reason: it needs a reference
frequency this job does not declare.
"""

from __future__ import annotations

from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.recipes import configure_ramsey_job
from quebra.schemas.track912 import track912Schema

# The logical name and category. `include` resolves JOB_ID, so this
# file can move without breaking any composite; recategorising costs one string edit.
JOB_ID = "ramsey_q1_100423"
JOB_FAMILY = "ramsey"

PREFIX = "q1_13h_1004_dataset"

job = Job("ramsey_q1_100423")
main_ds = Dataset(
    path="data/real_private/6D2S/100423_6D2S_qubit1.pickle",
    schema=track912Schema,
    qubit=1,
    device="6D2S",
    duration_h=13,
    extra={"run_name": PREFIX},
)

main_node = job.load(main_ds)

configure_ramsey_job(
    job,
    main_node,
    profile="overnight",
    include_fidelity=False,
    include_tlf=True,
    allan_fractional=False,
    figure_prefix=PREFIX,
)
