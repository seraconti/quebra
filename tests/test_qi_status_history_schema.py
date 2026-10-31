"""The Quantum Inspire status-history schema: declared columns, read as declared.

Every file here is synthetic, written under `tmp_path` in the real format (invented values).
"""

from __future__ import annotations

import os
import time
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import pandas as pd
import pytest

from quebra.core.dataset import Dataset
from quebra.core.job import _load_dataset
from quebra.schemas.qi_status_history import QiStatusHistorySchema

pytestmark = pytest.mark.unit

_ROWS = [
    "2026-01-01 00:00:00+00:00,Tuna-X,0,idle",
    "2026-01-01 00:00:30+00:00,Tuna-X,2,calibrating",
    "2026-01-01 00:01:00+00:00,Tuna-X,2,calibrating",
    "2026-01-01 00:01:00+00:00,Tuna-X,2,calibrating",
    "2026-01-01 00:01:30+00:00,Tuna-X,4,benchmarking",
    "2026-01-01 00:02:00+00:00,Tuna-X,1,executing",
    "2026-01-01 00:08:00+00:00,Tuna-X,3,offline",
]


def _write(
    tmp_path: Path, rows: list[str], header: str = "time,name,value,status"
) -> Path:
    path = tmp_path / "Tuna-X_status_history1.csv"
    path.write_text("\n".join([header, *rows]) + "\n")
    return path


def _load(path: Path, device: str | None = "Tuna-X"):
    return _load_dataset(
        Dataset(
            path=path,
            schema=QiStatusHistorySchema,
            device=device,
            extra={"run_name": "unit"},
        )
    )


@contextmanager
def _host_zone(zone: str) -> Iterator[None]:
    """The process's local time zone set to `zone`, and restored after."""
    before = os.environ.get("TZ")
    os.environ["TZ"] = zone
    time.tzset()
    try:
        yield
    finally:
        if before is None:
            del os.environ["TZ"]
        else:
            os.environ["TZ"] = before
        time.tzset()


def test_qi_status_history_round_trips_a_synthetic_file(tmp_path: Path) -> None:
    """Oracle: the file itself, and one closed-form time. The declared columns rebuilt from
    the Norm, the times written back in their declared format, equal `pd.read_csv` of the
    file, row for row, duplicates kept."""
    path = _write(tmp_path, _ROWS)
    # Read with the local zone away from UTC: a time read as naive local is then wrong on
    # every host, a UTC one included.
    with _host_zone("Europe/Amsterdam"):
        norm = _load(path)
    # 2026-01-01T00:00:00Z: 56 years of 365 days plus 14 leap days, 20454 days of 86400 s.
    assert norm["sample_time_unix_s"][0] == 20454 * 86400.0
    assert norm["sample_value"].dtype.kind == "i"
    rebuilt = pd.DataFrame(
        {
            "time": [
                pd.Timestamp(t, unit="s", tz="UTC").isoformat(sep=" ")
                for t in norm["sample_time_unix_s"]
            ],
            "name": [norm["meta"]["device"]] * norm["meta"]["n_samples"],
            "value": norm["sample_value"],
            "status": norm["sample_status"],
        }
    )
    pd.testing.assert_frame_equal(rebuilt, pd.read_csv(path), check_dtype=False)
    assert norm["meta"]["dataset_id"] == "unit"
    assert norm["meta"]["n_samples"] == len(_ROWS)


def test_qi_status_history_refuses_a_file_that_breaks_its_declaration(
    tmp_path: Path,
) -> None:
    """One case per declared rule: each raises, nothing is coerced, sorted or guessed."""
    good = list(_ROWS)
    cases = [
        (good, "time,name,value", None, KeyError, "lacks columns"),
        (
            [good[0] + ",x"],
            "time,name,value,status,extra",
            None,
            ValueError,
            "undeclared",
        ),
        ([], None, "Tuna-X", ValueError, "no samples"),
        ([good[0].replace("Tuna-X", "")], None, None, ValueError, "blank name"),
        ([good[0].replace("+00:00", "")], None, "Tuna-X", ValueError, "are not"),
        ([good[0].replace("+00:00", "+01:00")], None, "Tuna-X", ValueError, "are not"),
        ([good[0].replace(",0,idle", ",5,idle")], None, "Tuna-X", ValueError, "match"),
        (
            [good[0].replace(",0,idle", ",0,asleep")],
            None,
            "Tuna-X",
            ValueError,
            "unknown",
        ),
        ([good[0].replace(",0,idle", ",1,idle")], None, "Tuna-X", ValueError, "match"),
        (
            [good[0].replace(",0,idle", ",,idle")],
            None,
            "Tuna-X",
            ValueError,
            "integers",
        ),
        ([good[0], good[1].replace("Tuna-X", "Tuna-Y")], None, None, ValueError, "one"),
        ([good[0]], None, "Tuna-Y", ValueError, "declares device"),
        ([good[1], good[0]], None, "Tuna-X", ValueError, "time order"),
    ]
    for rows, header, device, error, match in cases:
        path = _write(tmp_path, rows, header or "time,name,value,status")
        with pytest.raises(error, match=match):
            _load(path, device=device)
