"""Execute the grid and write the size and power tables.

One cell -> many replicates -> one row per (check, calibration, variant). Every row carries
what is needed to read it honestly:

- `rejection_rate` with `mc_se` beside it. At 2000 replicates the Monte Carlo standard
  error at a true 0.05 is `sqrt(.05*.95/2000) = 0.0049`, so a +/-0.01 band around nominal
  is two SE wide and a correct test fails it by chance about 5% of the time. A size table
  without its own error bar invites exactly the over-reading this bench exists to prevent.
- `n_failed` and `fail_reason`. A replicate whose statistic is undefined is not a
  non-rejection; counting it as one would drag every degenerate cell toward "conservative".
- `mean_n_distinct_durations`. A cell where the quantised generator produced 7 distinct
  values out of 20 is a cell where the permutation null is nearly atomic, and its size
  should be read as such rather than compared with a continuous cell's.
- `mean_induced_lag1`. The swept parameter is a latent correlation; this is the
  duration-level dependence it actually produced, which is the quantity a check sees.

Parallelism is over CELLS, via joblib. A cell's cost grows with its n, so the work is
handed out cell by cell and the scheduler balances it; each worker
seeds its own generator from `(MASTER_SEED, cell index, replicate)` so the run reproduces
exactly regardless of how the cells were distributed.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from jobs.bench.arms import (
    ARM_A,
    ARM_B,
    ARM_C,
    ARM_D,
    ARM_E,
    arm_a,
    arm_b,
    arm_c,
    arm_d,
    arm_e,
)
from jobs.bench.grid import (
    ALPHA,
    KIND_SIZE,
    MASTER_SEED,
    N_PERM,
    Cell,
    all_cells,
)
from quebra.analyzers.checks._permutation import block_permutations
from quebra.analyzers.checks.battery import row_key, run_battery
from quebra.analyzers.checks.result import Segment

RESULTS_DIR = Path(__file__).resolve().parent / "results"

TABLE_COLUMNS = [
    "kind",
    "arm",
    "clock",
    "n_target",
    "shape",
    "quantised",
    "censoring_target",
    "rho",
    "b",
    "check",
    "calibration",
    "variant",
    "alpha",
    "n_replicates",
    "n_used",
    "n_failed",
    "n_failed_tau_only",
    "rejection_rate",
    "mc_se",
    "mean_n_events",
    "mean_n_segments",
    "censoring_realised",
    "mean_n_distinct_durations",
    "mean_induced_lag1",
    "fail_reason",
]


@dataclass(frozen=True)
class Replicate:
    """One generated record, before any check touches it."""

    segments: list[Segment]
    n_events: int
    n_segments: int
    n_distinct: int


def generate(cell: Cell, rng: np.random.Generator) -> Replicate:
    """Build one record for `cell`. Raises if the arm cannot serve this cell."""
    if cell.arm == ARM_A:
        segments = arm_a(
            cell.n,
            rng,
            shape=cell.shape,
            quantised=cell.quantised,
            censoring=cell.censoring,
        )
    elif cell.arm == ARM_D:
        segments = arm_d(
            cell.n,
            rng,
            shape=cell.shape,
            b=cell.b,
            quantised=cell.quantised,
            censoring=cell.censoring,
        )
    elif cell.arm == ARM_E:
        segments = arm_e(
            cell.n,
            rng,
            shape=cell.shape,
            rho=cell.rho,
            quantised=cell.quantised,
            censoring=cell.censoring,
        )
    elif cell.arm == ARM_B:
        segments, _dropped = arm_b(cell.n, rng, clock=cell.clock)
    elif cell.arm == ARM_C:
        segments, _dropped = arm_c(cell.n, rng, rho=cell.rho, clock=cell.clock)
    else:
        raise ValueError(f"unknown arm {cell.arm!r}")
    if not segments:
        raise ValueError("the arm produced no usable segment")
    durations = np.concatenate([s.x for s in segments])
    return Replicate(
        segments=segments,
        n_events=int(len(durations)),
        n_segments=len(segments),
        n_distinct=int(len(np.unique(durations))),
    )


def run_cell(cell: Cell, cell_index: int) -> list[dict[str, object]]:
    """Every row for one cell. Never raises: a failed replicate is counted, not fatal."""
    hits: dict[tuple[str, str, str], int] = {}
    seen: dict[tuple[str, str, str], int] = {}
    n_events: list[int] = []
    n_segments: list[int] = []
    n_distinct: list[int] = []
    induced_lag1: list[float] = []
    n_failed = 0
    n_failed_tau_only = 0
    first_failure = ""

    def note_failure(exc: Exception) -> None:
        nonlocal first_failure
        if not first_failure:
            first_failure = f"{type(exc).__name__}: {exc}"[:160]

    for replicate_index in range(cell.n_replicates):
        rng = np.random.default_rng([MASTER_SEED, cell_index, replicate_index])
        try:
            replicate = generate(cell, rng)
            perm = block_permutations(
                [s.n_events for s in replicate.segments], N_PERM, rng
            )
        except (ValueError, KeyError) as exc:
            n_failed += 1
            note_failure(exc)
            continue

        try:
            results = run_battery(
                replicate.segments,
                clock=cell.clock,
                perm=perm,
                include_tau_checks=cell.include_tau_checks,
            )
        except (ValueError, KeyError) as exc:
            note_failure(exc)
            # C1 and C2 are undefined on this record - typically a segment whose
            # durations are all one quantum, so gamma_hat is zero. C5 and C6 never touch
            # gamma_hat or tau and are perfectly well defined, so they are retried rather
            # than discarded: dropping them too would shrink their sample on exactly the
            # degenerate cells, and a rank check's size there is worth knowing. Each row's
            # own `n_used` is the denominator, so the surviving rows stay unbiased.
            if not cell.include_tau_checks:
                n_failed += 1
                continue
            try:
                results = run_battery(
                    replicate.segments,
                    clock=cell.clock,
                    perm=perm,
                    include_tau_checks=False,
                )
                n_failed_tau_only += 1
            except (ValueError, KeyError):
                n_failed += 1
                continue

        n_events.append(replicate.n_events)
        n_segments.append(replicate.n_segments)
        n_distinct.append(replicate.n_distinct)
        for result in results:
            key = row_key(result)
            seen[key] = seen.get(key, 0) + 1
            if result.p_value is not None and result.p_value < ALPHA:
                hits[key] = hits.get(key, 0) + 1
            if "r_lag1" in result.extra and key[2] == "studentized":
                induced_lag1.append(float(result.extra["r_lag1"]))

    rows: list[dict[str, object]] = []
    mean_events = float(np.mean(n_events)) if n_events else float("nan")
    mean_segments = float(np.mean(n_segments)) if n_segments else float("nan")
    realised_censoring = (
        mean_segments / (mean_events + mean_segments) if n_events else float("nan")
    )
    for key, used in sorted(seen.items()):
        rate = hits.get(key, 0) / used if used else float("nan")
        rows.append(
            {
                "kind": cell.kind,
                "arm": cell.arm,
                "clock": cell.clock,
                "n_target": cell.n,
                "shape": cell.shape,
                "quantised": cell.quantised,
                "censoring_target": cell.censoring,
                "rho": cell.rho,
                "b": cell.b,
                "check": key[0],
                "calibration": key[1],
                "variant": key[2],
                "alpha": ALPHA,
                "n_replicates": cell.n_replicates,
                "n_used": used,
                "n_failed": n_failed,
                "n_failed_tau_only": n_failed_tau_only,
                "rejection_rate": rate,
                "mc_se": float(np.sqrt(max(rate * (1.0 - rate), 0.0) / used))
                if used
                else float("nan"),
                "mean_n_events": mean_events,
                "mean_n_segments": mean_segments,
                "censoring_realised": realised_censoring,
                "mean_n_distinct_durations": float(np.mean(n_distinct))
                if n_distinct
                else float("nan"),
                "mean_induced_lag1": float(np.mean(induced_lag1))
                if induced_lag1
                else float("nan"),
                "fail_reason": first_failure,
            }
        )
    if not rows:
        # Every replicate failed. Emit one row saying so rather than silently omitting
        # the cell, which would leave a hole the report could read as "not run".
        rows.append(
            {
                **{column: None for column in TABLE_COLUMNS},
                "kind": cell.kind,
                "arm": cell.arm,
                "clock": cell.clock,
                "n_target": cell.n,
                "shape": cell.shape,
                "quantised": cell.quantised,
                "censoring_target": cell.censoring,
                "rho": cell.rho,
                "b": cell.b,
                "alpha": ALPHA,
                "n_replicates": cell.n_replicates,
                "n_used": 0,
                "n_failed": n_failed,
                "n_failed_tau_only": n_failed_tau_only,
                "fail_reason": first_failure or "no rows produced",
            }
        )
    return rows


def main(n_jobs: int = -1, limit: int | None = None) -> pd.DataFrame:
    cells = all_cells()
    if limit is not None:
        cells = cells[:limit]
    total_replicates = sum(c.n_replicates for c in cells)
    print(
        f"[bench] {len(cells)} cells, {total_replicates} replicates, "
        f"B={N_PERM}, alpha={ALPHA}, n_jobs={n_jobs}",
        flush=True,
    )
    started = time.perf_counter()
    batches = Parallel(n_jobs=n_jobs, verbose=5)(
        delayed(run_cell)(cell, index) for index, cell in enumerate(cells)
    )
    elapsed_s = time.perf_counter() - started

    frame = pd.DataFrame([row for batch in batches for row in batch])
    frame = frame.reindex(columns=TABLE_COLUMNS)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    size = frame[frame["kind"] == KIND_SIZE]
    power = frame[frame["kind"] != KIND_SIZE]
    size.to_csv(RESULTS_DIR / "size_table.csv", index=False)
    power.to_csv(RESULTS_DIR / "power_table.csv", index=False)
    (RESULTS_DIR / "runtime.txt").write_text(
        f"cells={len(cells)}\n"
        f"replicates={total_replicates}\n"
        f"n_perm={N_PERM}\n"
        f"n_jobs={n_jobs}\n"
        f"wall_clock_s={elapsed_s:.1f}\n"
        f"wall_clock_min={elapsed_s / 60.0:.1f}\n"
    )
    print(
        f"[bench] done in {elapsed_s / 60.0:.1f} min; "
        f"{len(size)} size rows, {len(power)} power rows -> {RESULTS_DIR}",
        flush=True,
    )
    return frame


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run the check calibration bench.")
    parser.add_argument("--jobs", type=int, default=-1)
    parser.add_argument("--limit", type=int, default=None, help="first N cells only")
    parsed = parser.parse_args()
    main(n_jobs=parsed.jobs, limit=parsed.limit)
