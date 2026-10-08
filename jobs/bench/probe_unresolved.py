"""How affordable is a resolved/unresolved split on the T2* ladder?

Standalone. Wired into nothing, imported by nothing. It answers one question with
numbers so the decision about a 4-state timeline is made on evidence rather than taste.

A read is `unresolved` at threshold `a` when abs(value - a) < k * sigma_v, with sigma_v
the read's `t2star_error_s`. A read whose sigma is not finite is `unknown`: counted
separately, never silently treated as resolved.

THREE THINGS TO KNOW BEFORE READING THE OUTPUT
1. `frac_reads_unknown_error` is a FLOOR, not a measurement. schemas/track912.py's
   `validate` drops every row with a null "T2star error" before the Norm is built, so
   this count sees only surviving inf/NaN. The true pre-validation unknown fraction is
   not observable from a single load, and this script loads once by design.
2. No chi-squared or sigma filter is applied here, so window counts will not match what
   a job produces. This measures the ladder against the raw record.
3. Every job file declares a ladder byte-identical to T2STAR_DEFAULT_LADDER, so that
   ladder is used for all datasets.

Usage:  python jobs/bench/probe_unresolved.py
Writes: bench/results/probe_unresolved_out.csv
"""

from __future__ import annotations

import dataclasses
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from quebra.analyzers import windows
from quebra.analyzers.t2star import T2STAR_DEFAULT_LADDER
from quebra.core.dataset import Dataset
from quebra.core.job import _load_dataset
from quebra.core.paths import default_dataset_root, repo_root, resolve_dataset_path
from quebra.schemas.track912 import track912Schema

K_VALUES = (0.5, 1.0, 2.0)
LADDER_SPACING_S = 1e-6  # the ladder steps in 1 µs
DATASET_SUBDIR = "data/real_private/6D2S"
# Beside this file, like every other bench output. Writing through repo_root()
# must point at `jobs/bench/results/`.
OUT_CSV = Path(__file__).resolve().parent / "results" / "probe_unresolved_out.csv"


def _dataset_files() -> list[Path]:
    root = default_dataset_root() / DATASET_SUBDIR
    if not root.is_dir():
        raise FileNotFoundError(
            f"dataset directory not found: {root}. Datasets live outside the repo, "
            f"one level up (see core/paths.py)."
        )
    return sorted(root.glob("*.pickle"))


def _qubit_from_name(stem: str) -> int:
    # Some stems carry a suffix after the id, e.g. 120423_6D2S_qubit1_before.
    match = re.search(r"qubit(\d+)", stem)
    if match is None:
        raise ValueError(f"cannot read a qubit id out of {stem!r}")
    return int(match.group(1))


def _load(path: Path) -> dict[str, object]:
    ds = Dataset(
        path=str(path),
        schema=track912Schema,
        qubit=_qubit_from_name(path.stem),
        device="6D2S",
        extra={"run_name": path.stem},
    )
    ds = dataclasses.replace(
        ds, path=resolve_dataset_path(ds.path, default_dataset_root())
    )
    return _load_dataset(ds)


def _probe_one(name: str, norm: dict[str, object]) -> tuple[list[dict], dict]:
    t_s = np.asarray(norm["t_rel_s"], dtype=float)
    v_s = np.asarray(norm["T2star_s"], dtype=float)
    sigma = (
        np.asarray(norm["T2star_error_s"], dtype=float)
        if "T2star_error_s" in norm
        else np.full(len(v_s), np.nan)
    )

    finite = np.isfinite(t_s) & np.isfinite(v_s)
    t_f, v_f, sig_f = t_s[finite], v_s[finite], sigma[finite]
    sigma_known = np.isfinite(sig_f)
    n_reads = int(len(v_f))

    carved = windows.run(
        windows.WindowsInputs(
            t_rel_s=t_s,
            values=v_s,
            thresholds=T2STAR_DEFAULT_LADDER,
            dataset_id=name,
        )
    )
    reads = carved.reads
    windows_df = carved.windows_in_spec

    rows: list[dict] = []
    for k in K_VALUES:
        for label, thr_s, bvg in T2STAR_DEFAULT_LADDER:
            unresolved = np.zeros(n_reads, dtype=bool)
            unresolved[sigma_known] = (
                np.abs(windows.margin(v_f[sigma_known], thr_s, bvg))
                < k * sig_f[sigma_known]
            )

            sub = reads[(reads["threshold_label"] == label) & reads["in_spec"]]
            n_windows = int((windows_df["threshold_label"] == label).sum())
            purities: list[float] = []
            if len(sub) and n_windows:
                # The read table's t_read_s comes from the same finite array in the
                # same order, so a search gives each in-spec read its index back.
                order = np.searchsorted(t_f, sub["t_read_s"].to_numpy(dtype=float))
                resolved_flag = ~unresolved[order]
                grouped = pd.DataFrame(
                    {
                        "window_index": sub["window_index"].to_numpy(),
                        "resolved": resolved_flag,
                    }
                ).groupby("window_index")["resolved"]
                purities = grouped.mean().tolist()

            rows.append(
                {
                    "dataset": name,
                    "threshold": label,
                    "k": k,
                    "frac_reads_unresolved": float(unresolved.mean())
                    if n_reads
                    else np.nan,
                    "frac_reads_unknown_error": float((~sigma_known).mean())
                    if n_reads
                    else np.nan,
                    "n_windows": n_windows,
                    "frac_windows_fully_resolved": float(
                        np.mean([p == 1.0 for p in purities])
                    )
                    if purities
                    else np.nan,
                    "median_window_purity": float(np.median(purities))
                    if purities
                    else np.nan,
                }
            )

    median_sigma_us = (
        float(np.median(sig_f[sigma_known])) * 1e6 if sigma_known.any() else np.nan
    )
    per_dataset = {
        "dataset": name,
        "n_reads": n_reads,
        "median_sigma_us": median_sigma_us,
        "sigma_over_ladder_spacing": median_sigma_us / (LADDER_SPACING_S * 1e6),
        "frac_reads_unknown_error": float((~sigma_known).mean()) if n_reads else np.nan,
        "n_gaps": carved.diagnostics["n_gaps"],
    }
    return rows, per_dataset


def main() -> None:
    warnings.filterwarnings("ignore")
    files = _dataset_files()
    print(
        f"probing {len(files)} datasets from {default_dataset_root() / DATASET_SUBDIR}\n"
    )

    all_rows: list[dict] = []
    per_dataset: list[dict] = []
    for path in files:
        name = path.stem
        try:
            norm = _load(path)
        except (
            Exception
        ) as exc:  # a dataset that will not load is a finding, not a stop
            print(f"  SKIP {name}: {type(exc).__name__}: {exc}")
            continue
        rows, summary = _probe_one(name, norm)
        all_rows.extend(rows)
        per_dataset.append(summary)

    df = pd.DataFrame(all_rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)

    pd.set_option("display.width", 200)
    fmt = lambda x: f"{x:.4g}"  # noqa: E731

    print("\n=== PER DATASET (once) ===")
    print(pd.DataFrame(per_dataset).to_string(index=False, float_format=fmt))

    print("\n=== PER DATASET x k  (mean over the ladder) ===")
    by_ds_k = (
        df.groupby(["dataset", "k"])[
            [
                "frac_reads_unresolved",
                "frac_windows_fully_resolved",
                "median_window_purity",
            ]
        ]
        .mean()
        .reset_index()
    )
    print(by_ds_k.to_string(index=False, float_format=fmt))

    print("\n=== PER THRESHOLD x k  (across all datasets) ===")
    by_thr_k = (
        df.groupby(["threshold", "k"])[
            ["frac_reads_unresolved", "n_windows", "frac_windows_fully_resolved"]
        ]
        .mean()
        .reset_index()
    )
    print(by_thr_k.to_string(index=False, float_format=fmt))

    print(f"\nwrote {OUT_CSV.relative_to(repo_root())} ({len(df)} rows)")


if __name__ == "__main__":
    main()
