"""Check ledger of every record of 2x2 qubit 1, and their survey.

One member of the `check_ledger` family, built by
`quebra.ledger_recipe.configure_check_ledger_job` exactly as the 6D2S members are: per
record, the T2* carve and the check ledger of its windows, then one survey over the records.
What is here is this qubit's record list, the 2x2 drift datasets of Capannelli et al.,
Supplementary notes, Supplementary Fig. 11 (p. 11), in calendar order; each record's
`run_name`, `2x2_q1_<DDMMYY>`, names its nodes, sinks and survey row. A record whose content
fails to load or carve reads `not computed` with the reason.
"""

from __future__ import annotations

from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.ledger_recipe import LEDGER_KNOBS, configure_check_ledger_job
from quebra.schemas.track912 import track912Schema

JOB_ID = "check_ledger_2x2_q1"
JOB_FAMILY = "check_ledger"
# Not swept by a bare `run --all`: it carves and scores every record of the qubit.
JOB_SWEEP = False

BENCH_SIZE_TABLE = Dataset(path="jobs/bench/results/size_table.csv", schema=None)


def _dataset(date: str) -> Dataset:
    return Dataset(
        path=f"data/real_private/2x2/{date}_2x2_qubit1.pickle",
        schema=track912Schema,
        qubit=1,
        device="2x2",
        extra={"run_name": f"2x2_q1_{date}"},
    )


RECORDS = [
    _dataset("090623"),
    _dataset("100623"),
    _dataset("200623"),
    _dataset("250623"),
    _dataset("300623"),
    _dataset("030723"),
    _dataset("050723"),
    _dataset("070723"),
    _dataset("100723"),
    _dataset("220823"),
]

job = Job(JOB_ID)

configure_check_ledger_job(
    job,
    prefix="2x2_q1",
    records=RECORDS,
    bench_size_table=BENCH_SIZE_TABLE,
    **LEDGER_KNOBS,  # type: ignore[arg-type]
)
