"""Quantum Inspire status history: one status sample per row, one file per device.

The declared columns, exactly these four:

| Column | Type | Meaning |
|---|---|---|
| `time` | str, `YYYY-MM-DD HH:MM:SS+00:00` | the sample time, tz-aware UTC |
| `name` | str | the device; one value in a file |
| `value` | int, 0 to 4 | the status code |
| `status` | str | the status, one to one with `value` |

The codes: idle 0, executing 1, calibrating 2, offline 3, benchmarking 4. The schema
refuses a file that breaks any of this: a missing or an extra column, no rows, a time in
another format or offset, a code outside 0 to 4, a status outside the five, a status that
does not match its code, a blank name, more than one device name (or one other than the
`Dataset`'s `device`), rows out of time order. Rows are kept as the file has them,
duplicate timestamps included; `analyzers/cycles.spells_from_samples` decides what a
duplicate means. The CSV reader types the `value` column, so a cell written `00` or `+0`
reads as 0.

The `Norm` holds `sample_time_unix_s` (float), `sample_status` (str) and `sample_value`
(int), in file order, and `meta` with `dataset_id`, `device` and `n_samples`.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from quebra.core.dataset import Dataset
from quebra.core.types import Norm

COLUMNS = ("time", "name", "value", "status")
STATUS_CODES = {
    "idle": 0,
    "executing": 1,
    "calibrating": 2,
    "offline": 3,
    "benchmarking": 4,
}
# The one time format the files declare: seconds, and the UTC offset written out.
TIME_PATTERN = r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\+00:00"
_EPOCH = pd.Timestamp("1970-01-01", tz="UTC")


class QiStatusHistorySchema:
    """`time, name, value, status` samples -> `Norm`, refusing anything undeclared.

    A `to_norm` schema, not a `pandera.DataFrameModel`, for the reason
    `schemas/ramsey_series.py` gives: it owns the conversion and hands back no frame.
    """

    @classmethod
    def to_norm(cls, frame: pd.DataFrame, dataset: Dataset) -> Norm:
        missing = [c for c in COLUMNS if c not in frame.columns]
        if missing:
            raise KeyError(
                f"{cls.__name__}: {Path(dataset.path).name} lacks columns {missing}"
            )
        # The loader attaches a declared `device` or `qubit` as a column; the file may not.
        attached = {"device": dataset.device, "qubit_id": dataset.qubit}
        extra = [
            c for c in frame.columns if c not in COLUMNS and attached.get(c) is None
        ]
        if extra:
            raise ValueError(
                f"{cls.__name__}: {Path(dataset.path).name} has undeclared columns "
                f"{extra}; declared: {list(COLUMNS)}"
            )
        if len(frame) == 0:
            raise ValueError(f"{Path(dataset.path).name} holds no samples")
        time = frame["time"].astype(str)
        malformed = ~time.str.fullmatch(TIME_PATTERN)
        if malformed.any():
            raise ValueError(
                f"{int(malformed.sum())} time value(s) are not 'YYYY-MM-DD HH:MM:SS+00:00' "
                f"(first: {time[malformed].iloc[0]!r})"
            )
        parsed = pd.to_datetime(time, format="%Y-%m-%d %H:%M:%S%z", utc=True)
        time_unix_s = (parsed - _EPOCH).dt.total_seconds().to_numpy(dtype=float)
        if not pd.api.types.is_integer_dtype(frame["value"]):
            raise ValueError(
                f"the value column must hold integers; got {frame['value'].dtype}"
            )
        value = frame["value"].to_numpy(dtype=int)
        status = frame["status"].astype(str).to_numpy()
        unknown = sorted(set(status) - set(STATUS_CODES))
        if unknown:
            raise ValueError(
                f"unknown status(es) {unknown}; declared: {sorted(STATUS_CODES)}"
            )
        expected = np.array([STATUS_CODES[s] for s in status], dtype=int)
        mismatch = value != expected
        if mismatch.any():
            raise ValueError(
                f"{int(mismatch.sum())} row(s) carry a value that does not match their "
                f"status (codes: {STATUS_CODES})"
            )
        # Checked before `astype(str)`, which would read a blank name as the device "nan".
        if frame["name"].isna().any():
            raise ValueError(
                f"{int(frame['name'].isna().sum())} row(s) have a blank name"
            )
        names = sorted(set(frame["name"].astype(str)))
        if len(names) != 1:
            raise ValueError(f"a file holds one device; the name column has {names}")
        device = names[0]
        if dataset.device is not None and dataset.device != device:
            raise ValueError(
                f"the Dataset declares device {dataset.device!r}; the file is {device!r}"
            )
        if np.any(np.diff(time_unix_s) < 0.0):
            raise ValueError(
                f"{int(np.sum(np.diff(time_unix_s) < 0.0))} sample(s) are earlier than "
                "the one before; the file must be in time order"
            )
        return Norm(
            {
                "sample_time_unix_s": time_unix_s,
                "sample_status": status,
                "sample_value": value,
                "meta": {
                    "dataset_id": str(
                        dataset.extra.get("run_name", Path(dataset.path).stem)
                    ),
                    "device": device,
                    "n_samples": int(len(time_unix_s)),
                },
            }
        )
