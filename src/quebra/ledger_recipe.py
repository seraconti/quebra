"""The check-ledger jobs: the checks of every record of one qubit, and their survey.

A ledger job runs the checks and nothing else. Per record it wires the T2* carve
(`recipes.wire_t2star_carve`, the definition the T2* panel is drawn from), the check ledger
of both bags of those windows, one `CheckLedgerPanel` per bag and a materialized
`{record}_check_ledger_data`.
Over all of the job's records it then builds one independence survey: one grid per
instrument, records down, thresholds across.

The job also bundles its ledgers into one `CheckLedgers` artifact, `{prefix}_check_ledgers`.
The device-wide survey (`jobs/composite/independence_survey.py`) includes the ledger jobs
and refs that bundle, once per job, through `included_ledger_set`. It neither carves nor
scores, so every verdict it draws is the one the per-record ledger printed. One ref per job
matters: a composite run without reuse runs an included job once for every distinct ref.
The compliance survey (`jobs/composite/compliance_survey.py`) refs each job's carve bundle
instead (`RecordCarves`, `{prefix}_carves`), once per job, through `included_carves`.

A module of its own, not part of `quebra.recipes`, because of the identity closure
(`quebra.core.closure`): whatever a module imports joins the closure of every step function
it defines. Kept here, the checks reach only the ledger jobs; inside `recipes` they would
re-key every T2* job whenever a check changed.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import pandas as pd

from quebra.analyzers.check_ledger import (
    CheckLedger,
    CheckLedgerInputs,
    make_inputs_from_windows,
)
from quebra.analyzers.check_ledger import run as run_ledger
from quebra.analyzers.event_table import SIDES
from quebra.analyzers.independence_survey import (
    IndependenceSurveyData,
    build_independence_survey,
)
from quebra.analyzers.windows import WindowsResult
from quebra.core._artifact_guard import StaleArtifactGuard
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.core.types import RecordFailure
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
    # C3 is not run: the ledgers carry the permutation rows and the group rows built from
    # them only, and C3 has no bench cell,
    # so its cells could never read `pass`. `c3_n_null_sim` stays a declared knob, recorded
    # on the artifact as `CheckLedger.c3_n_null_sim`, for a caller that asks for C3.
    "include_c3": False,
    "c3_n_null_sim": 200,
}

# The materialize name of a record's ledger is its prefix plus this.
LEDGER_DATA_SUFFIX = "_check_ledger_data"
# The materialize name of a record's carve, a `WindowsResult` or a `RecordFailure`. Other
# jobs read the carve by this name; keep it stable.
CARVE_SUFFIX = "_windows"
# The materialize name of a ledger job's carve bundle is its prefix plus this.
CARVES_DATA_SUFFIX = "_carves"

TARGETS = ["static", "academic"]


def record_label(prefix: str) -> str:
    """`q1_070423` -> `q1 070423`: the row label of a record in every survey grid."""
    return prefix.replace("_", " ")


@dataclass
class CheckLedgers(StaleArtifactGuard):
    """Every record ledger of one ledger job, in the job's order, each with its row label."""

    labels: tuple[str, ...]
    ledgers: tuple[CheckLedger, ...]


@dataclass
class RecordCarves(StaleArtifactGuard):
    """Every record carve of one ledger job, in the job's order, each under its run name.

    The same carves the job persists one by one as `{record}_windows`, bundled so a
    composite reading every record takes one ref per job.
    """

    records: tuple[str, ...]
    carves: tuple[WindowsResult | RecordFailure, ...]
    # The ladder the carves were cut on, so a reader of the bundle declares none of its own.
    thresholds: tuple[tuple[str, float, bool], ...]


def ledger_step(
    window_result: WindowsResult | RecordFailure,
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

    Every knob is a required kwarg, so each one reaches the provenance label. A record that
    failed to load or carve gets a ledger whose every threshold and clock reads `not
    computed` with the failure's note.
    """
    if isinstance(window_result, RecordFailure):
        return run_ledger(
            CheckLedgerInputs(
                windows=pd.DataFrame({"threshold_label": []}),
                bench_size_table=bench_size_table,
                thresholds=list(thresholds),
                dataset_id=window_result.dataset_id,
                alpha=alpha,
                min_events_pass=min_events_pass,
                tie_cutoff_distinct=tie_cutoff_distinct,
                lag_max=lag_max,
                n_permutations=n_permutations,
                seed=seed,
                include_c3=include_c3,
                c3_n_null_sim=c3_n_null_sim,
                record_failure=window_result.note(),
            )
        )
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


def ledger_bag(ledger: CheckLedger, *, side: str) -> CheckLedger:
    """One bag's rows of a record's ledger, for the figure that draws that bag."""
    if side not in SIDES:
        raise ValueError(f"side must be one of {SIDES}; got {side!r}")
    rows = ledger.rows[ledger.rows["side"] == side].reset_index(drop=True)
    if not len(rows):
        raise ValueError(f"the ledger of {ledger.dataset_id!r} has no {side} row")
    return replace(ledger, rows=rows)


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


def collect_carves(
    *carves: WindowsResult | RecordFailure,
    records: tuple[str, ...],
    thresholds: list[tuple[str, float, bool]],
) -> RecordCarves:
    """Bundle a job's carves under their run names, checked against each carve's own id,
    with the ladder they were cut on."""
    if len(records) != len(carves):
        raise ValueError(f"{len(records)} run names for {len(carves)} carves")
    for record, carve in zip(records, carves):
        carved = (
            carve.dataset_id
            if isinstance(carve, RecordFailure)
            else str(carve.meta["dataset_id"])
        )
        if carved != record:
            raise ValueError(
                f"run name {record!r} is attached to the carve of {carved!r}"
            )
    return RecordCarves(
        records=tuple(records),
        carves=tuple(carves),
        thresholds=tuple((str(lb), float(v), bool(b)) for lb, v, b in thresholds),
    )


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
    """Per record: the carve, then its check ledger and one figure per bag. Then the
    bundle, and one survey.

    A record's prefix is its dataset's `extra["run_name"]`, which also becomes the ledger's
    `dataset_id`, so the label a row carries is derived from the record and not retyped.
    `prefix` names the job-level sinks: `{prefix}_check_ledgers`, the carve bundle
    `{prefix}_carves` the compliance survey reads, and the `{prefix}_independence_survey_*`
    grids and figures. No T2* panel: that belongs to the T2* jobs. Each record's carve is
    also persisted on its own, as `{record}_windows` (`CARVE_SUFFIX`).

    The carve records a failure (`wire_t2star_carve(on_failure="record")`): a record whose
    content fails to load or carve is persisted as a `RecordFailure` under the same name,
    and its ledger reads `not computed` with the reason. It is never dropped.
    """
    bench = job.load_df(bench_size_table)
    ledgers: list[object] = []
    labels: list[str] = []
    carves: list[object] = []
    records_seen: list[str] = []
    for dataset in records:
        record = str(dataset.extra["run_name"])
        _, window_node = wire_t2star_carve(
            job,
            dataset=dataset,
            thresholds=T2STAR_THRESHOLDS,
            node_suffix=f"_{record}",
            on_failure="record",
        )
        job.materialize(window_node, name=f"{record}{CARVE_SUFFIX}")
        carves.append(window_node)
        records_seen.append(record)
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
        # One figure per bag. The titles name the bag; `ledger_bag` hands each its rows.
        for side in SIDES:
            bag = job.step(
                ledger_bag, ledger, name=f"check_ledger_{record}_{side}", side=side
            )
            job.figure(
                CheckLedgerPanel,
                bag,
                targets=TARGETS,
                title=f"{record} check ledger {side} bag",
            )
        # Not `{record}_check_ledger_{side}_bag`: that is a figure's safe-named title, and
        # two sinks may not share a provenance name.
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
    carve_bundle = job.step(
        collect_carves,
        *carves,  # type: ignore[arg-type]
        name="record_carves",
        records=tuple(records_seen),
        thresholds=T2STAR_THRESHOLDS,
    )
    job.materialize(carve_bundle, name=f"{prefix}{CARVES_DATA_SUFFIX}")
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


def included_carves(included: Any) -> tuple[list[str], list[str], object]:
    """The run names and devices of an included ledger job's records, and a ref to its
    `RecordCarves` bundle.

    Read off the included job's graph, so the compliance survey retypes no record list: the
    run names are the bundle's own `records`, each device its record's dataset's. Raises
    unless the job has exactly one bundle node and exactly one sink persisting it, as
    `included_ledger_set` does.
    """
    job = included.job
    nodes = [n for n in job.dag.values() if n.fn is collect_carves]
    if len(nodes) != 1:
        raise ValueError(
            f"included job {included.alias!r} has {len(nodes)} `collect_carves` nodes, "
            "not one; it is not a ledger job"
        )
    (node,) = nodes
    names = [
        sink.name
        for sink in job.sinks
        if getattr(sink, "node", None) is not None and sink.node.node_id == node.node_id
    ]
    if len(names) != 1:
        raise ValueError(
            f"included job {included.alias!r} persists its carve bundle under "
            f"{len(names)} sinks, not one"
        )
    devices = []
    for parent in node.inputs:
        load = job.dag[parent.node_id]
        while "dataset" not in load.kwargs:
            (up,) = load.inputs
            load = job.dag[up.node_id]
        dataset = load.kwargs["dataset"]
        if dataset.device is None:
            raise ValueError(
                f"included job {included.alias!r}: {dataset.path!r} declares no device, "
                "so the survey cannot place it in a panel"
            )
        devices.append(str(dataset.device))
    return list(node.kwargs["records"]), devices, included.ref(names[0])
