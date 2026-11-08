"""The `freq` and `fid` job families: one job per record, on the records the ledgers score.

Oracle: the ten check-ledger jobs' record lists, read off their graphs, and the declared
Rabi table. The table is embargoed, so the fid half resolves it through the data root and
skips when it is not there; the freq half needs no data and always runs.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from quebra.core.discovery import discover
from quebra.core.job import _import_job
from quebra.core.paths import (
    DataRootNotFound,
    DataUnavailable,
    repo_root,
    resolve_data_root,
    resolve_dataset_path,
)

DECLARED = "data/real_private/declared/rabi_frequencies_2x2.csv"


def _records_by_job(family: str) -> dict[str, list[Path]]:
    """Each job of `family`, with the record pickles its graph loads."""
    jobs = discover(repo_root() / "jobs")
    held: dict[str, list[Path]] = {}
    for declared in jobs.values():
        if declared.family != family:
            continue
        dag = _import_job(declared.path).dag
        held[declared.job_id] = [
            Path(node.kwargs["dataset"].path)
            for node in dag.values()
            if "dataset" in node.kwargs
            and Path(node.kwargs["dataset"].path).suffix == ".pickle"
        ]
    return held


def _key(record: Path) -> tuple[str, int, str]:
    """(device, qubit, run_date) of a record named `{DDMMYY}_{device}_qubit{N}.pickle`."""
    date, device, qubit = record.stem.split("_")
    return (
        device,
        int(qubit.removeprefix("qubit")),
        f"20{date[4:]}-{date[2:4]}-{date[:2]}",
    )


@pytest.mark.policy
def test_freq_and_fid_jobs_cover_exactly_the_ledger_records() -> None:
    """Oracle: the ledger jobs' record lists and the declared Rabi table. Each freq and fid
    job loads one record and is named for it; the freq jobs load every ledger record once,
    and the fid jobs every 2x2 ledger record with a declared Rabi row once."""
    ledger = [r for paths in _records_by_job("check_ledger").values() for r in paths]
    assert len(ledger) == len(set(ledger)) and ledger, (
        "premise: one ledger row per record"
    )
    covered: dict[str, list[Path]] = {}
    for family in ("freq", "fid"):
        jobs = _records_by_job(family)
        assert jobs, f"no job declares JOB_FAMILY {family!r}"
        for job_id, records in jobs.items():
            assert len(records) == 1, (job_id, records)
            device, qubit, run_date = _key(records[0])
            date = records[0].stem.split("_")[0]
            assert job_id == f"{family}_{device.lower()}_q{qubit}_{date}", records[0]
        covered[family] = sorted(r for (r,) in jobs.values())
    assert covered["freq"] == sorted(ledger)
    two_by_two = sorted(r for r in ledger if _key(r)[0] == "2x2")
    assert set(covered["fid"]) <= set(two_by_two)
    try:
        table = pd.read_csv(
            resolve_dataset_path(DECLARED, resolve_data_root()), comment="#"
        )
    except (DataRootNotFound, DataUnavailable) as exc:
        pytest.skip(f"freq half checked; the fid half needs the declared table: {exc}")
    declared = {
        (str(d), int(q), str(r))
        for d, q, r in zip(table["device"], table["qubit"], table["run_date"])
    }
    assert covered["fid"] == [r for r in two_by_two if _key(r) in declared]
