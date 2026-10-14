"""The default normaliser: a timestamp/frequency series into a `Norm`.

Keeping this out of `quebra.core.job._load_dataset` is what keeps the generic loading path
generic. Inside it, `job.load` would require `timestamp` and `frequency` columns and a
resolvable run-start whether or not a schema was involved, so a job with nothing to do with
frequency tracking could not use `job.load` at all.

It is a schema, which is what it always was. `job.load` dispatches to
`schema.to_norm(frame, dataset)` and this class is only the DEFAULT - pass your own schema
and none of the requirements below apply to you. See `docs/WRITING_A_SCHEMA.md`.

Nothing about the computation changed in the move. The column names, the sort, the
`delta_hz` convention (relative to the mean, no Rabi drive synthesised) and the three-level
run-start resolution are as they were. The optional passthrough columns are listed in
`docs/WRITING_A_SCHEMA.md`.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from quebra.core.dataset import Dataset
from quebra.core.types import Norm
from quebra.transforms.lookup_prior import check_unix_s


class RamseySeriesSchema:
    """`timestamp` + `frequency` -> `Norm`, with run-start resolution.

    Deliberately NOT a `pandera.DataFrameModel`. The two schema shapes are different
    things: a `validate` schema checks and cleans a frame and lets the default normaliser
    run afterwards; a `to_norm` schema takes the frame and owns the conversion outright.
    This is the second kind, and there is no frame to hand back at the end of it.
    """

    @classmethod
    def to_norm(cls, frame: pd.DataFrame, dataset: Dataset) -> Norm:
        if "timestamp" not in frame.columns or "frequency" not in frame.columns:
            raise KeyError(
                f"{cls.__name__} requires 'timestamp' and 'frequency' columns; "
                f"{Path(dataset.path).name} has: {list(frame.columns)[:12]}. "
                f"Pass a schema with a to_norm() method to load data of another shape, "
                f"or use job.load_df() for a raw DataFrame."
            )

        timestamp = pd.to_numeric(frame["timestamp"], errors="coerce")
        frequency = pd.to_numeric(frame["frequency"], errors="coerce")
        valid = np.isfinite(timestamp.to_numpy(dtype=float)) & np.isfinite(
            frequency.to_numpy(dtype=float)
        )
        if not np.any(valid):
            raise ValueError(
                f"Dataset {dataset.path} has no valid numeric timestamp/frequency rows"
            )

        t_raw = timestamp.to_numpy(dtype=float)[valid]
        f_hz = frequency.to_numpy(dtype=float)[valid]
        order = np.argsort(t_raw)
        t_raw = t_raw[order]
        f_hz = f_hz[order]
        t_rel_s = t_raw - float(t_raw[0])

        meta = dict(dataset.extra)
        meta.update(
            {
                "dataset_id": str(meta.get("run_name", Path(dataset.path).stem)),
                "run_name": str(meta.get("run_name", Path(dataset.path).stem)),
                "qubit": dataset.qubit,
                "device": dataset.device,
                "duration_h": dataset.duration_h,
                "n_points": int(len(t_rel_s)),
            }
        )
        # Determine run_start_unix_s via three resolution levels (see docs/TIME_SEMANTICS.md):
        #   1. Explicit: Dataset.extra['run_start_unix_s'] already in meta - validate and use.
        #   2. date_only_midnight: DDMMYY_ filename prefix - midnight of that date (local naive).
        #   3. No valid source → raise; do NOT fall back to t_raw[0] (yields ~1970 epoch).
        if meta.get("run_start_unix_s") is not None:
            meta["run_start_unix_s"] = check_unix_s(
                meta["run_start_unix_s"], label="Dataset.extra['run_start_unix_s']"
            )
            meta["run_start_resolution"] = "explicit"
        else:
            date_match = re.match(r"^(\d{2})(\d{2})(\d{2})_", Path(dataset.path).stem)
            if date_match is not None:
                try:
                    day = int(date_match.group(1))
                    month = int(date_match.group(2))
                    year = 2000 + int(date_match.group(3))
                    run_day = pd.Timestamp(year=year, month=month, day=day)
                    meta["run_start_unix_s"] = float(run_day.value / 1e9)
                    meta["run_start_resolution"] = "date_only_midnight"
                    warnings.warn(
                        f"[{meta['dataset_id']}] run_start_unix_s derived from DDMMYY filename prefix "
                        f"as midnight local time (resolution: date_only_midnight). Intra-day "
                        f"precision is lost, so nothing that needs the time of day may rely "
                        f"on it. "
                        f"Provide Dataset.extra['run_start_unix_s'] for precision.",
                        UserWarning,
                        # 3, not 2: after the move out of `_load_dataset` there is one more
                        # frame between here and the caller, and `stacklevel=2` attributed
                        # every one of these warnings to `core/job.py`, a library line.
                        stacklevel=3,
                    )
                except Exception as exc:
                    raise ValueError(
                        f"Cannot determine run start time for '{meta['dataset_id']}': DDMMYY prefix "
                        f"found but parsing failed ({exc}). "
                        f"Provide Dataset.extra['run_start_unix_s'] explicitly."
                    ) from exc
            else:
                raise ValueError(
                    f"Cannot determine run start time for '{meta['dataset_id']}'. "
                    f"Provide Dataset.extra['run_start_unix_s'] explicitly, "
                    f"or use a filename with a DDMMYY_ prefix encoding the run date."
                )

        # default delta (relative to mean); rabi drive is not synthesized here
        delta_hz = f_hz - float(np.mean(f_hz))
        meta.setdefault("rabi_source", "missing")

        norm = Norm(
            {
                "t_rel_s": t_rel_s,
                "delta_hz": delta_hz,
                "raw_frequency_hz": f_hz,
                "meta": meta,
            }
        )
        if "normalised chi-square" in frame.columns:
            chi = pd.to_numeric(
                frame["normalised chi-square"], errors="coerce"
            ).to_numpy(dtype=float)[valid][order]
            norm["chi_squared"] = chi
        if "T2star" in frame.columns:
            norm["T2star_s"] = pd.to_numeric(frame["T2star"], errors="coerce").to_numpy(
                dtype=float
            )[valid][order]
        if "T2star error" in frame.columns:
            norm["T2star_error_s"] = pd.to_numeric(
                frame["T2star error"], errors="coerce"
            ).to_numpy(dtype=float)[valid][order]
        if "frequency error" in frame.columns:
            norm["frequency_error_hz"] = pd.to_numeric(
                frame["frequency error"], errors="coerce"
            ).to_numpy(dtype=float)[valid][order]
        return norm
