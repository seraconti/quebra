"""The check-ledger jobs: the checks of every record of one qubit, and their survey.

A ledger job runs the checks and nothing else. Per record it wires the T2* carve
(`recipes.wire_t2star_carve`, the definition the T2* panel is drawn from), the check ledger
of those windows, its `CheckLedgerPanel` and a materialized `{prefix}_check_ledger_data`.
Over all of the job's records it then builds one independence survey: one grid per
instrument, records down, thresholds across.

The job also bundles its ledgers into one `CheckLedgers` artifact, `{prefix}_check_ledgers`.
The device-wide survey (`jobs/composite/independence_survey.py`) includes the ledger jobs
and refs that bundle, once per job, through `included_ledger_set`. It neither carves nor
scores, so every verdict it draws is the one the per-record ledger printed. One ref per job
matters: a composite run without reuse runs an included job once for every distinct ref.

A module of its own, not part of `quebra.recipes`, because of the identity closure
(`quebra.core.closure`): whatever a module imports joins the closure of every step function
it defines. Kept here, the checks reach only the ledger jobs; inside `recipes` they would
re-key every T2* job whenever a check changed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd

from quebra.analyzers.check_ledger import CheckLedger, make_inputs_from_windows
from quebra.analyzers.check_ledger import run as run_ledger
from quebra.analyzers.independence_survey import (
    IndependenceSurveyData,
    build_independence_survey,
)
from quebra.analyzers.windows import WindowsResult
from quebra.core._artifact_guard import StaleArtifactGuard
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.panels.check_ledger import CheckLedgerPanel
from quebra.plots.independence_survey_plot import (
    SURVEY_PLOTS,
    IndependenceSurveyOverviewPlot,
)
from quebra.recipes import T2STAR_THRESHOLDS, wire_t2star_carve

# The knobs that turn a p-value into a verdict, defined once so every qubit's ledger is
# scored on the same rule. Each still reaches the provenance label as a step kwarg.
LEDGER_KNOBS: dict[str, object] = {
    # The level the bench measured size at (`jobs/bench/grid.py:ALPHA`), which
    # `bench_acceptance_at_n` assumes.
    "alpha": 0.05,
    # The smallest event count at which the promotion report scores a check
    # (`jobs/bench/report.py:scored_cells`, `min_n=35`). A policy value: it is not
    # derived from power.
    "min_events_pass": 35,
    # The bench sweeps quantisation only as a boolean, so it does not calibrate this value
    # (`docs/iid_checks/LIMITATIONS.md`, section 6).
    "tie_cutoff_distinct": 5,
    "lag_max": 5,
    # The bench's `N_PERM`, so a ledger p-value has the resolution the bench measured.
    "n_permutations": 999,
    "seed": 20260812,
    # C3 has no bench cell, so its cells can never read `pass`; it is run for its p-value.
    # Its null simulation grows steeply in cost with n. At 200 the p-value floor is
    # 1/201 = 0.005, an order below alpha, and the count is recorded on the artifact as
    # `CheckLedger.c3_n_null_sim`.
    "include_c3": True,
    "c3_n_null_sim": 200,
}

# The materialize name of a record's ledger is its prefix plus this.
LEDGER_DATA_SUFFIX = "_check_ledger_data"

TARGETS = ["static", "academic"]


def record_label(prefix: str) -> str:
    """`q1_070423` -> `q1 070423`: the row label of a record in every survey grid."""
    return prefix.replace("_", " ")


@dataclass
class CheckLedgers(StaleArtifactGuard):
    """Every record ledger of one ledger job, in the job's order, each with its row label."""

    labels: tuple[str, ...]
    ledgers: tuple[CheckLedger, ...]


def ledger_step(
    window_result: WindowsResult,
    bench_size_table: pd.DataFrame,
    *,
    thresholds: list[tuple[str, float, bool]],
    alpha: float,
    min_events_pass: int,
    tie_cutoff_distinct: int,
    lag_max: int,
    n_permutations: int,
    seed: int,
    include_c3: bool,
    c3_n_null_sim: int,
) -> CheckLedger:
    """The check ledger of one record's windows.

    Every knob is a required kwarg, so each one reaches the provenance label.
    """
    return run_ledger(
        make_inputs_from_windows(
            window_result,
            bench_size_table,
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
    )


def collect_ledgers(*ledgers: CheckLedger, labels: tuple[str, ...]) -> CheckLedgers:
    """Bundle a job's ledgers, so one ref carries all of them.

    Each label must be its own ledger's record label, `record_label(ledger.dataset_id)`, so
    a label list out of step with the ledgers raises here instead of mislabelling a row.
    """
    if len(labels) != len(ledgers):
        raise ValueError(f"{len(labels)} labels for {len(ledgers)} ledgers")
    for label, ledger in zip(labels, ledgers):
        if record_label(ledger.dataset_id) != label:
            raise ValueError(
                f"label {label!r} is attached to the ledger of {ledger.dataset_id!r}"
            )
    return CheckLedgers(labels=tuple(labels), ledgers=tuple(ledgers))


def survey_ledgers(
    *sets: CheckLedgers, dataset_labels: tuple[str, ...], alpha: float
) -> IndependenceSurveyData:
    """The independence survey of every ledger in `sets`, rows in `dataset_labels` order.

    `dataset_labels` is declared when the graph is built, so the row order reaches the
    provenance label. It must equal the sets' own labels in order, or this raises.
    """
    labels = tuple(label for s in sets for label in s.labels)
    if labels != tuple(dataset_labels):
        raise ValueError(
            f"declared rows {dataset_labels} differ from the ledgers' own {labels}"
        )
    return build_independence_survey(
        *(ledger for s in sets for ledger in s.ledgers),
        dataset_labels=tuple(dataset_labels),
        alpha=alpha,
    )


def wire_independence_survey(
    job: Job,
    ledger_sets: list[object],
    labels: list[str],
    *,
    prefix: str,
    alpha: float,
) -> None:
    """One survey over `ledger_sets`: the grids, the overview and one figure per check.

    `ledger_sets` are `CheckLedgers` nodes of this job or refs into included ledger jobs;
    the survey reshapes their ledgers and scores nothing.
    """
    survey = job.step(
        survey_ledgers,
        *ledger_sets,  # type: ignore[arg-type]
        name="independence_survey",
        dataset_labels=tuple(labels),
        alpha=alpha,
    )
    job.materialize(survey, name=f"{prefix}_grids")
    # The overview first: it is the one a reader opens to orient, and the per-instrument
    # figures are where a p-value is read.
    job.figure(
        IndependenceSurveyOverviewPlot,
        survey,
        targets=TARGETS,
        title=f"{prefix} overview all checks",
    )
    for plot_class in SURVEY_PLOTS:
        job.figure(
            plot_class, survey, targets=TARGETS, title=f"{prefix} {plot_class.__name__}"
        )


def configure_check_ledger_job(
    job: Job,
    *,
    prefix: str,
    records: list[Dataset],
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
    """Per record: the carve, then its check ledger. Then the bundle, and one survey.

    A record's prefix is its dataset's `extra["run_name"]`, which also becomes the ledger's
    `dataset_id`, so the label a row carries is derived from the record and not retyped.
    `prefix` names the job-level sinks: `{prefix}_check_ledgers` and the
    `{prefix}_independence_survey_*` grids and figures. No T2* panel and no materialized
    windows: those belong to the T2* jobs.
    """
    bench = job.load_df(bench_size_table)
    ledgers: list[object] = []
    labels: list[str] = []
    for dataset in records:
        record = str(dataset.extra["run_name"])
        _, window_node = wire_t2star_carve(
            job, dataset=dataset, thresholds=T2STAR_THRESHOLDS, node_suffix=f"_{record}"
        )
        ledger = job.step(
            ledger_step,
            window_node,
            bench,
            name=f"check_ledger_{record}",
            thresholds=T2STAR_THRESHOLDS,
            alpha=alpha,
            min_events_pass=min_events_pass,
            tie_cutoff_distinct=tie_cutoff_distinct,
            lag_max=lag_max,
            n_permutations=n_permutations,
            seed=seed,
            include_c3=include_c3,
            c3_n_null_sim=c3_n_null_sim,
        )
        job.figure(
            CheckLedgerPanel, ledger, targets=TARGETS, title=f"{record} check ledger"
        )
        # Not `{record}_check_ledger`: that is the figure's safe-named title, and two sinks
        # may not share a provenance name.
        job.materialize(ledger, name=f"{record}{LEDGER_DATA_SUFFIX}")
        ledgers.append(ledger)
        labels.append(record_label(record))
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


def included_ledger_set(included: Any) -> tuple[list[str], object]:
    """The row labels and a ref to the `CheckLedgers` bundle of an included ledger job.

    Read off the included job's graph, so the device-wide survey retypes no record list.
    Raises unless the job has exactly one bundle node and exactly one sink persisting it: a
    survey silently missing a job's rows is the failure this guards.
    """
    nodes = [n for n in included.job.dag.values() if n.fn is collect_ledgers]
    if len(nodes) != 1:
        raise ValueError(
            f"included job {included.alias!r} has {len(nodes)} `collect_ledgers` nodes, "
            "not one; it is not a ledger job"
        )
    (node,) = nodes
    # A materialize sink carries `node`; a figure sink carries `input` instead.
    names = [
        sink.name
        for sink in included.job.sinks
        if getattr(sink, "node", None) is not None and sink.node.node_id == node.node_id
    ]
    if len(names) != 1:
        raise ValueError(
            f"included job {included.alias!r} persists its ledger bundle under "
            f"{len(names)} sinks, not one"
        )
    return list(node.kwargs["labels"]), included.ref(names[0])
