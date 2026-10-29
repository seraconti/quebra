"""Execute the grid and write the size, power and serial tables.

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
  duration-level dependence it actually produced, which is the quantity a check sees. It
  averages C5's uncentred `r_1` (`extra["r_lag1"]`), so on a single-segment cell it
  carries that statistic's null offset: a cell with no dependence reads about -1/N, not 0,
  N being its event count.
  On a censored cell, which is several segments, the offset also moves with the layout, as
  in the serial cells below; there it is not the dependence the cell carries either.

The serial cells get a table of their own, `serial_table.csv`, with `TABLE_COLUMNS` plus
`delta`, the realised segment counts (min and max beside the mean), the quantiles of the
segment sizes over every segment of the computed replicates, and `n_not_computed`, the
replicates a row has no p-value for. Rejection rates are over computed replicates only. It
leaves out `mean_induced_lag1`: across several segments the uncentred `r_1` also moves with
the layout and with any difference in level between segments, so there it does not measure
the dependence a cell carries.

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
    ARM_E2,
    ARM_EL,
    arm_a,
    arm_b,
    arm_c,
    arm_d,
    arm_e,
    arm_e2,
    arm_el,
)
from jobs.bench.grid import (
    ALPHA,
    KIND_POWER,
    KIND_SERIAL,
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

SEGMENT_SIZE_QUANTILES = {
    "segment_size_min": 0.0,
    "segment_size_q25": 0.25,
    "segment_size_median": 0.5,
    "segment_size_q75": 0.75,
    "segment_size_max": 1.0,
}


def _serial_columns() -> list[str]:
    """`TABLE_COLUMNS` with the serial table's additions, each beside its relative, and
    without `mean_induced_lag1`."""
    columns: list[str] = []
    for column in TABLE_COLUMNS:
        if column == "mean_induced_lag1":
            continue
        if column == "mean_n_segments":
            columns += [
                "n_segments_min",
                column,
                "n_segments_max",
                *SEGMENT_SIZE_QUANTILES,
            ]
            continue
        columns.append(column)
        if column == "b":
            columns.append("delta")
        elif column == "n_failed_tau_only":
            columns.append("n_not_computed")
    return columns


SERIAL_COLUMNS = _serial_columns()


@dataclass(frozen=True)
class Replicate:
    """One generated record, before any check touches it."""

    segments: list[Segment]
    n_events: int
    n_segments: int
    n_distinct: int


def _refuse_unservable(cell: Cell) -> None:
    """Raise for a cell its arm would serve without what the cell asks for."""
    if cell.arm == ARM_EL and (cell.quantised or cell.censoring):
        raise ValueError(
            "Arm EL draws two continuous segments of its own; a cell asking for "
            f"quantised={cell.quantised} or censoring={cell.censoring} would not get it"
        )


def generate(cell: Cell, rng: np.random.Generator) -> Replicate:
    """Build one record for `cell`. Raises if the arm cannot serve this cell."""
    _refuse_unservable(cell)
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
    elif cell.arm == ARM_E2:
        segments = arm_e2(
            cell.n,
            rng,
            shape=cell.shape,
            rho=cell.rho,
            quantised=cell.quantised,
            censoring=cell.censoring,
        )
    elif cell.arm == ARM_EL:
        segments = arm_el(
            rng,
            events_per_segment=cell.n // 2,
            shape=cell.shape,
            phi=cell.rho,
            delta=cell.delta,
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
    """Every row for one cell. A failed replicate is counted, not fatal; a cell its arm
    cannot serve raises before the first replicate, so a misconfigured grid stops the run."""
    _refuse_unservable(cell)
    hits: dict[tuple[str, str, str], int] = {}
    seen: dict[tuple[str, str, str], int] = {}
    n_events: list[int] = []
    n_segments: list[int] = []
    n_distinct: list[int] = []
    segment_events: list[int] = []
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

        # One set of arguments for both calls below, so the retry cannot drift from the
        # first call: C5-raw is a bench row and both must ask for it.
        battery_kwargs = {"clock": cell.clock, "perm": perm, "include_c5_raw": True}
        try:
            results = run_battery(
                replicate.segments,
                include_tau_checks=cell.include_tau_checks,
                **battery_kwargs,
            )
        except (ValueError, KeyError) as exc:
            note_failure(exc)
            # Something in the battery is undefined on this record. The retry drops C1
            # and C2 and keeps C5, C6 and CvM by permutation, which stay defined when tau
            # is singular (`tau == T_N`). It does not rescue a record whose `gamma_hat` is
            # zero, a segment whose durations are all one quantum: CvM's statistic divides
            # by `gamma_hat` too, so the retry raises and the replicate counts as failed
            # for every row. Each row's own `n_used` is the denominator, so the surviving
            # rows stay unbiased.
            if not cell.include_tau_checks:
                n_failed += 1
                continue
            try:
                results = run_battery(
                    replicate.segments,
                    include_tau_checks=False,
                    **battery_kwargs,
                )
                n_failed_tau_only += 1
            except (ValueError, KeyError):
                n_failed += 1
                continue

        n_events.append(replicate.n_events)
        n_segments.append(replicate.n_segments)
        n_distinct.append(replicate.n_distinct)
        segment_events.extend(s.n_events for s in replicate.segments)
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
    layout = {
        "delta": cell.delta,
        "n_segments_min": int(min(n_segments)) if n_segments else None,
        "n_segments_max": int(max(n_segments)) if n_segments else None,
        **{
            column: float(np.quantile(segment_events, q)) if segment_events else None
            for column, q in SEGMENT_SIZE_QUANTILES.items()
        },
    }
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
                "n_not_computed": cell.n_replicates - used,
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
                **layout,
            }
        )
    if not rows:
        # Every replicate failed. Emit one row saying so rather than silently omitting
        # the cell, which would leave a hole the report could read as "not run".
        rows.append(
            {
                **{column: None for column in SERIAL_COLUMNS},
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
                "n_not_computed": cell.n_replicates,
                "fail_reason": first_failure or "no rows produced",
                **layout,
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
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    size = frame[frame["kind"] == KIND_SIZE].reindex(columns=TABLE_COLUMNS)
    power = frame[frame["kind"] == KIND_POWER].reindex(columns=TABLE_COLUMNS)
    serial = frame[frame["kind"] == KIND_SERIAL].reindex(columns=SERIAL_COLUMNS)
    if len(size) + len(power) + len(serial) != len(frame):
        unknown = sorted(set(frame["kind"]) - {KIND_SIZE, KIND_POWER, KIND_SERIAL})
        raise ValueError(f"rows of unknown kind {unknown} would be written to no table")
    size.to_csv(RESULTS_DIR / "size_table.csv", index=False)
    power.to_csv(RESULTS_DIR / "power_table.csv", index=False)
    serial.to_csv(RESULTS_DIR / "serial_table.csv", index=False)
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
        f"{len(size)} size rows, {len(power)} power rows, {len(serial)} serial rows "
        f"-> {RESULTS_DIR}",
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
