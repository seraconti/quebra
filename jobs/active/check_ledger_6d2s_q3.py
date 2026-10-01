"""Check ledger of every record of 6D2S qubit 3, and their survey.

One member of the `check_ledger` family, built by
`quebra.ledger_recipe.configure_check_ledger_job`: per record, the T2* carve and the check
ledger of its windows, then one survey over the records. What is here is this qubit's record
list, in calendar order; each record's `run_name` names its nodes, sinks and survey row.
The `_before` / `_after` files are row splits of a full record, so scoring them beside it
would score the same reads twice; they are left out.
"""

from __future__ import annotations

from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.ledger_recipe import LEDGER_KNOBS, configure_check_ledger_job
from quebra.schemas.track912 import track912Schema

JOB_ID = "check_ledger_6d2s_q3"
JOB_FAMILY = "check_ledger"
# Not swept by a bare `run --all`: it carves and scores every record of the qubit.
JOB_SWEEP = False

BENCH_SIZE_TABLE = Dataset(path="jobs/bench/results/size_table.csv", schema=None)


def _dataset(date: str) -> Dataset:
    return Dataset(
        path=f"data/real_private/6D2S/{date}_6D2S_qubit3.pickle",
        schema=track912Schema,
        qubit=3,
        device="6D2S",
        extra={"run_name": f"q3_{date}"},
    )


RECORDS = [
    _dataset("210423"),
    _dataset("110723"),
]

job = Job(JOB_ID)

configure_check_ledger_job(
    job,
    prefix="q3",
    records=RECORDS,
    bench_size_table=BENCH_SIZE_TABLE,
    **LEDGER_KNOBS,  # type: ignore[arg-type]
)
