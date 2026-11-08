"""The check-ledger jobs of the 2x2 infidelity windows: every record of one qubit, and a survey.

`configure_fidelity_ledger_job` is the infidelity twin of
`ledger_recipe.configure_check_ledger_job`. Per record it takes the record's included Ramsey
job with fidelity (a `fid_2x2_*` job) and wires its carve through
`infidelity_recipe.wire_infidelity_carve`: one ref, to the fidelity result that job persists,
and the |x| windows recomputed with that job's own carve kwargs, on the nines ladder. It runs
the check ledger on both bags of those windows with the T2* ledgers' knobs
(`ledger_recipe.LEDGER_KNOBS`, passed by the job) and bench size table, draws one figure per
bag, and materializes the carve as `{record}_windows` and the ledger as
`{record}_check_ledger_data`. Over the job's records it bundles the ledgers as
`{prefix}_check_ledgers` and draws one independence survey; `ledger_recipe.included_ledger_set`
reads that bundle as it reads a T2* ledger job's.

A record is named by its carve's `dataset_id`, the run name of its Ramsey job, and its ledger
carries the record's name in the paper's dataset tables beside it
(`InfidelityCheckLedger.paper_name`).

Failures. A ledger that raises ValueError or KeyError, which is how the steps refuse a
record's content (`recipes._carve_stage`), is recorded: that record's ledger reads `not
computed` with the reason at every
threshold and clock, and the job's other records and its survey are still scored, as the T2*
ledgers record a failed carve. The record wraps the whole ledger run, so a fault common to
every record, such as a bench table the ledger cannot read, reads as every record `not
computed` with the same reason, and the reason names it. The carve needs no such path. It is
the included job's own carve node, the same step on the same persisted input with the same
kwargs, which that job's run computed, fresh or reused, before the ref resolves; at the same
commit and input it cannot fail here unless the included job failed first. A reused run made
under other library versions is the exception, and it raises here, loudly. A failure inside
an included job stops this job, because the runner resolves the ref by running that job,
outside any step of this one.

The bench. Whether the bench's size cells, measured on synthetic records and keyed by event
count and clock, transfer to an |x| carve is not measured, and the ties cutoff is not
calibrated on it; every ledger and every figure says so (`BENCH_TRANSFER`).

A module of its own, apart from `ledger_recipe`, so the T2* ledger jobs' identity is
unchanged: whatever a module imports joins the identity closure of every step function it
defines.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, fields
from typing import Any

import pandas as pd

from quebra.analyzers.check_ledger import CheckLedger
from quebra.analyzers.event_table import SIDES
from quebra.analyzers.windows import WindowsResult
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.core.types import RecordFailure
from quebra.infidelity_recipe import wire_infidelity_carve
from quebra.ledger_recipe import (
    CARVE_SUFFIX,
    LEDGER_DATA_SUFFIX,
    TARGETS,
    collect_ledgers,
    ledger_bag,
    ledger_step,
    record_label,
    wire_independence_survey,
)
from quebra.panels.infidelity_check_ledger import InfidelityCheckLedgerPanel

# What every infidelity ledger and its figure say about the bench they read.
BENCH_TRANSFER = (
    "Whether the bench's size cells, measured on synthetic records (some generated as "
    "durations, some as reads put through the windows carve) and keyed by event count and "
    "clock, transfer to this |x| carve is not measured, and the ties cutoff "
    "(tie_cutoff_distinct) is not calibrated on it (docs/iid_checks/LIMITATIONS.md, "
    "sections 6 and 8)."
)
# The stage a recorded ledger failure names.
LEDGER_STAGE = "check_ledger"


@dataclass
class InfidelityCheckLedger(CheckLedger):
    """The check ledger of one record's infidelity windows.

    Every field of `CheckLedger`, so every reader of a ledger reads this one unchanged, plus
    the record's name in the paper's dataset tables and the bench statement its figure
    prints (`BENCH_TRANSFER`).
    """

    paper_name: str = ""
    bench_transfer: str = ""


def fidelity_ledger_step(
    window_result: WindowsResult,
    bench_size_table: pd.DataFrame,
    *,
    paper_name: str,
    thresholds: list[tuple[str, float, bool]],
    alpha: float,
    min_events_pass: int,
    tie_cutoff_distinct: int,
    lag_max: int,
    n_permutations: int,
    seed: int,
    include_c3: bool,
    c3_n_null_sim: int,
) -> InfidelityCheckLedger:
    """The check ledger of one record's infidelity windows, by `ledger_recipe.ledger_step`.

    A ValueError or KeyError from it is recorded as a `RecordFailure` of stage
    `LEDGER_STAGE` for the carve's record, whose ledger reads `not computed` with the reason
    at every threshold and clock (module docstring). Every knob is a required kwarg, so each
    reaches the provenance label.
    """
    record = str(window_result.meta["dataset_id"])
    knobs: dict[str, Any] = dict(
        thresholds=thresholds,
        alpha=alpha,
        min_events_pass=min_events_pass,
        tie_cutoff_distinct=tie_cutoff_distinct,
        lag_max=lag_max,
        n_permutations=n_permutations,
        seed=seed,
        include_c3=include_c3,
        c3_n_null_sim=c3_n_null_sim,
    )
    try:
        ledger = ledger_step(window_result, bench_size_table, **knobs)
    except (ValueError, KeyError) as exc:
        failure = RecordFailure(
            dataset_id=record,
            stage=LEDGER_STAGE,
            reason=f"{type(exc).__name__}: {exc}"[:300],
        )
        ledger = ledger_step(failure, bench_size_table, **knobs)
    return InfidelityCheckLedger(
        **{f.name: getattr(ledger, f.name) for f in fields(CheckLedger)},
        paper_name=paper_name,
        bench_transfer=BENCH_TRANSFER,
    )


def configure_fidelity_ledger_job(
    job: Job,
    *,
    prefix: str,
    included: Sequence[Any],
    bench_size_table: Dataset,
    alpha: float,
    min_events_pass: int,
    tie_cutoff_distinct: int,
    lag_max: int,
    n_permutations: int,
    seed: int,
    include_c3: bool,
    c3_n_null_sim: int,
) -> None:
    """Per record: the carve, its check ledger and one figure per bag. Then the bundle and
    one survey (module docstring).

    `included` holds the returns of `job.include`, one Ramsey job with fidelity per record,
    in the job's order; the job file calls `job.include` itself, so `docs/JOBS.md` lists
    every include. `prefix` names the job-level sinks: `{prefix}_check_ledgers` and the
    `{prefix}_independence_survey_*` grids and figures. Raises at graph build when no record
    is given, or when two records were carved on different ladders.
    """
    if not included:
        raise ValueError(f"{prefix}: a ledger job needs at least one record")
    bench = job.load_df(bench_size_table)
    ledgers: list[object] = []
    labels: list[str] = []
    ladders: dict[str, tuple[tuple[str, float, bool], ...]] = {}
    for sub in included:
        carve = wire_infidelity_carve(job, sub)
        record = carve.run_name
        thresholds = carve.carve_kwargs["thresholds"]
        ladders[record] = tuple(tuple(rung) for rung in thresholds)
        job.materialize(carve.windows, name=f"{record}{CARVE_SUFFIX}")
        ledger = job.step(
            fidelity_ledger_step,
            carve.windows,
            bench,
            name=f"check_ledger_{record}",
            paper_name=carve.paper_name,
            thresholds=thresholds,
            alpha=alpha,
            min_events_pass=min_events_pass,
            tie_cutoff_distinct=tie_cutoff_distinct,
            lag_max=lag_max,
            n_permutations=n_permutations,
            seed=seed,
            include_c3=include_c3,
            c3_n_null_sim=c3_n_null_sim,
        )
        for side in SIDES:
            bag = job.step(
                ledger_bag, ledger, name=f"check_ledger_{record}_{side}", side=side
            )
            job.figure(
                InfidelityCheckLedgerPanel,
                bag,
                targets=TARGETS,
                title=f"{record} check ledger {side} bag",
            )
        job.materialize(ledger, name=f"{record}{LEDGER_DATA_SUFFIX}")
        ledgers.append(ledger)
        labels.append(record_label(record))
    if len(set(ladders.values())) != 1:
        raise ValueError(f"{prefix}: records carved on different ladders: {ladders}")
    bundle = job.step(
        collect_ledgers,
        *ledgers,  # type: ignore[arg-type]
        name="check_ledgers",
        labels=tuple(labels),
    )
    job.materialize(bundle, name=f"{prefix}_check_ledgers")
    wire_independence_survey(
        job, [bundle], labels, prefix=f"{prefix}_independence_survey", alpha=alpha
    )
