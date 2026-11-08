"""Record-figure jobs: the figures of every record of one ledger job, one PDF per block,
and the figures of the infidelity record of one Ramsey job.

A record-figure job includes exactly one ledger job and takes one ref into it, the ledger
bundle (`ledger_recipe.included_ledger_set`): a composite run without reuse runs an
included job once for every distinct ref. Every other input is read off the included job's
graph when the job is built: each record's `Dataset`, the kwargs of its carve, and the
run-set its ledger ran. The carve is recomputed here with those kwargs (`wire_t2star_carve`)
and held to the ledger by `analyzers/record_figure.cross_check_carve`, so every figure is
drawn on the windows its ledger scored.

The records drawn are the ledger job's records that the paper's dataset tables list
(`paper_tables.PAPER_RECORDS`). A ledger record outside them gets no figure. A record the
tables shade (`paper_tables.SHADED`) is drawn and carries the flag text in its caption. A
listed record of this ledger job's qubits that the job does not hold raises.

**Blocks.** `BLOCKS` is the one registration point: each entry wires one block of one
record, its builder step from `analyzers/record_figure.py` and its figure from
`plots/record_figure_plot.py`, off the record's `RecordNodes`, and returns its block's node.
A new block is one function appended there, and given its row in `PAGE_ORDER`: after the
blocks, one page per record stacks them in that reading order (`RecordFigurePlot`). The
job's tables are wired after the records, each written by the `table` target: the curve
summaries, one row per record, rung and bag, with the records outside the tables as
`excluded=`, and the two resolvability tables, every drawn record stacked
(`wire_resolvability_tables`). A record names its metric (`RecordNodes.metric`, one of
`METRICS`), and every figure prints `RecordNodes.display`, its name in the paper's tables
and its metric. A single-record job writes the same two tables for its own tellability with
`wire_tellability_tables`, and no figure; it names its record and metric the same way, read
off its own graph.

**Infidelity.** `configure_infidelity_record_figure_job` draws every block of the one
infidelity record of a Ramsey job with fidelity. It includes that job and wires its carve
through `infidelity_recipe.wire_infidelity_carve`: one ref, to the fidelity result, and the
|x| windows recomputed with the job's own carve kwargs. The check rows are the battery run
here on those windows (`ledger_step` with `LEDGER_KNOBS`), materialized: the Ramsey job runs
no ledger. The record's name in the paper's tables, read off the job's declared Rabi row,
must equal the name the figure job declares. Whether the bench's size cells, measured on
synthetic duration arms keyed by event count and clock, transfer to an |x| carve is not
measured; the rows say which cell they read, as for T2*.

**Tellability.** Each record's tellability is wired once, on the carve's own fit and
windows, with the carve's kwargs and the standalone jobs' permutation count and seed
(`tellability_recipe`), and every block that reads it takes that one node. A record whose
tellability refuses its content (a ValueError) reads `not computed` with the reason, as a
carve stage's refusal does, so the job's other records and blocks are still drawn; the
kwargs are checked when the job is built, so a wiring error still raises.

A module of its own, not part of `quebra.recipes`, for the reason `ledger_recipe.py` gives:
what a module imports joins the identity closure of every step function it defines.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from datetime import date
from typing import Any

from quebra.analyzers.check_ledger import CheckLedger
from quebra.analyzers.compliance_survey import record_compliance
from quebra.analyzers.check_outcome import RowKey
from quebra.analyzers.fidelity import FidelityResult
from quebra.analyzers.kaplan_meier import CONF_LEVEL, KaplanMeierSet
from quebra.analyzers.record_figure import (
    RESOLVABILITY_TABLES,
    assemble_record_figure,
    build_fraction_block,
    build_hazard_block,
    build_record_summaries,
    build_resolvability_block,
    build_survival_block,
    build_timeline_block,
    cross_check_carve,
    curve_summaries_of,
    kaplan_meier_of,
    ledger_run_set,
    nelson_aalen_of,
    panel_data_of,
    resolvability_table,
    summaries_table,
)
from quebra.analyzers.t2star import T2StarResult
from quebra.analyzers.tellability import TellabilityResult
from quebra.analyzers.within_calibration_data import WithinCalibrationPanelData
from quebra.analyzers.windows import WindowsResult
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.core.types import RecordFailure
from quebra.ledger_recipe import (
    LEDGER_DATA_SUFFIX,
    LEDGER_KNOBS,
    CheckLedgers,
    collect_carves,
    included_ledger_set,
    ledger_step,
)
from quebra.plots.record_figure_plot import (
    RecordFigurePlot,
    FractionBlockPlot,
    HazardBlockPlot,
    ResolvabilityBlockPlot,
    ResolvabilityTablePlot,
    SummariesTablePlot,
    SurvivalBlockPlot,
    TimelineBlockPlot,
)
from quebra.infidelity_recipe import (
    CARVE_KWARGS,
    paper_record_name,
    wire_infidelity_carve,
)
from quebra.recipes import (
    XI_SEED,
    _fidelity_panel_data,
    wire_t2star_carve,
)
from quebra.tellability_recipe import (
    TELLABILITY_PERMUTATIONS,
    TELLABILITY_SEED,
    _frequency_tellability,
    _infidelity_tellability,
    _t2star_tellability,
)

TARGETS = ["static", "academic"]
TABLE_TARGETS = ["table"]
# The T2* panel data's shape statistic needs this many reads per window, the value
# `recipes.configure_t2star_job` uses, so the timelines read the T2* panels' own states.
SHAPE_MIN_READS = 5
# The metrics a record figure draws, as `tellability` names them.
METRICS = ("t2star", "infidelity")
# How a figure names each metric in its titles and captions (`RecordNodes.display`).
METRIC_LABELS = {"t2star": "T2*", "frequency": "frequency", "infidelity": "infidelity"}


@dataclass(frozen=True)
class RecordNodes:
    """What every block of one record is wired from: the record's name in the paper's
    tables, which names its nodes and figures; its run name, the `dataset_id` its artifacts
    carry (the same string for a 6D2S record); its carve held to its ledger, its ledger, the
    run-set the ledger ran, and its flag text (empty unless the paper's tables shade it).

    `fit` and `windows` are the carve's own fit and windows nodes, and `carve_kwargs` the
    kwargs it was carved with. A step that must trace to observed reads alone, as
    `tests/test_windows_not_interpolated.py` holds a tellability step to, is wired on these:
    `carve` descends from the ledger bundle too. The cross-check still stands behind it,
    since the runner computes every node before any sink. `metric` names what the record
    measures (`METRICS`) and `tellability` is its one tellability node. `display` is what
    every figure prints: the record's name and its metric's label (`METRIC_LABELS`)."""

    record: str
    run_name: str
    carve: object
    ledger: object
    run_set: tuple[RowKey, ...]
    flag: str
    fit: object
    windows: object
    carve_kwargs: Mapping[str, Any]
    metric: str
    tellability: object

    def __post_init__(self) -> None:
        if self.metric not in METRICS:
            raise ValueError(
                f"{self.record}: metric {self.metric!r} is not one of {METRICS}"
            )

    @property
    def display(self) -> str:
        return f"{self.record}, {METRIC_LABELS[self.metric]}"


def ledger_of(ledgers: CheckLedgers, *, record: str) -> CheckLedger:
    """The ledger of `record` in a ledger job's bundle; raises unless exactly one."""
    found = [ledger for ledger in ledgers.ledgers if ledger.dataset_id == record]
    if len(found) != 1:
        raise ValueError(f"{len(found)} ledgers of {record!r} in the bundle, not one")
    return found[0]


def checked_carve(
    carve: WindowsResult | RecordFailure, ledger: CheckLedger
) -> WindowsResult | RecordFailure:
    """The carve, once `cross_check_carve` holds it to the ledger; every block reads the
    carve through this node."""
    cross_check_carve(carve, ledger)
    return carve


def wire_survival_block(job: Job, nodes: RecordNodes) -> object:
    """The survival block: Kaplan-Meier per rung and bag, each band with its check
    outcome."""
    km = job.step(
        kaplan_meier_of,
        nodes.carve,
        name=f"kaplan_meier_{nodes.record}",
        conf_level=CONF_LEVEL,
    )
    block = job.step(
        build_survival_block,
        km,
        nodes.ledger,
        name=f"survival_block_{nodes.record}",
        run_set=nodes.run_set,
        conf_level=CONF_LEVEL,
        flag=nodes.flag,
        display=nodes.display,
    )
    job.figure(
        SurvivalBlockPlot, block, targets=TARGETS, title=f"{nodes.record} survival"
    )
    return block


def wire_hazard_block(job: Job, nodes: RecordNodes) -> object:
    """The hazard block: Nelson-Aalen per rung and bag, each item with its band's check
    outcome in its `checks` field."""
    hazards = job.step(
        nelson_aalen_of,
        nodes.carve,
        name=f"nelson_aalen_{nodes.record}",
        conf_level=CONF_LEVEL,
    )
    block = job.step(
        build_hazard_block,
        hazards,
        nodes.ledger,
        name=f"hazard_block_{nodes.record}",
        run_set=nodes.run_set,
        conf_level=CONF_LEVEL,
        flag=nodes.flag,
        display=nodes.display,
    )
    job.figure(
        HazardBlockPlot,
        block,
        targets=TARGETS,
        title=f"{nodes.record} cumulative hazard",
    )
    return block


def wire_fraction_block(job: Job, nodes: RecordNodes) -> object:
    """The fraction block: `record_compliance` on the checked carve, the compliance
    survey's numbers for this record, drawn rung by rung."""
    compliance = job.step(
        record_compliance,
        nodes.carve,
        name=f"compliance_{nodes.record}",
        thresholds=tuple(
            (str(lb), float(v), bool(b))
            for lb, v, b in nodes.carve_kwargs["thresholds"]
        ),
    )
    block = job.step(
        build_fraction_block,
        compliance,
        nodes.carve,
        name=f"fraction_block_{nodes.record}",
        flag=nodes.flag,
        display=nodes.display,
    )
    job.figure(
        FractionBlockPlot,
        block,
        targets=TARGETS,
        title=f"{nodes.record} in-spec fraction",
    )
    return block


def t2star_tellability_of(
    fit: T2StarResult | RecordFailure,
    windows: WindowsResult | RecordFailure,
    *,
    thresholds: Sequence[Sequence[object]],
    gap_mult: float,
    unsure_k: float,
    n_permutations: int,
    seed: int,
) -> TellabilityResult | RecordFailure:
    """The tellability of a T2* record (`tellability_recipe._t2star_tellability`). A failed
    carve passes through, the windows' failure first: when only the windows stage failed,
    the fit is a result and the windows carry the failure. A ValueError, how the step
    refuses a record's content, becomes a `RecordFailure` of the stage `tellability`."""
    for upstream in (windows, fit):
        if isinstance(upstream, RecordFailure):
            return upstream
    try:
        return _t2star_tellability(
            fit,
            windows,
            thresholds=thresholds,
            gap_mult=gap_mult,
            unsure_k=unsure_k,
            n_permutations=n_permutations,
            seed=seed,
        )
    except ValueError as exc:
        return RecordFailure(
            dataset_id=str(windows.meta["dataset_id"]),
            stage="tellability",
            reason=f"{type(exc).__name__}: {exc}"[:300],
        )


def infidelity_tellability_of(
    result: FidelityResult | RecordFailure,
    windows: WindowsResult | RecordFailure,
    *,
    thresholds: Sequence[Sequence[object]],
    gap_mult: float,
    n_permutations: int,
    seed: int,
) -> TellabilityResult | RecordFailure:
    """The tellability of an infidelity record (`tellability_recipe._infidelity_tellability`,
    which reads k off the windows, so it takes no `unsure_k`). Failures pass through and a
    ValueError becomes a `RecordFailure` of the stage `tellability`, as for T2*."""
    for upstream in (windows, result):
        if isinstance(upstream, RecordFailure):
            return upstream
    try:
        return _infidelity_tellability(
            result,
            windows,
            thresholds=thresholds,
            gap_mult=gap_mult,
            n_permutations=n_permutations,
            seed=seed,
        )
    except ValueError as exc:
        return RecordFailure(
            dataset_id=str(windows.meta["dataset_id"]),
            stage="tellability",
            reason=f"{type(exc).__name__}: {exc}"[:300],
        )


def check_tellability_kwargs(
    record: str, carve_kwargs: Mapping[str, Any], *, n_permutations: int, seed: int
) -> None:
    """Raise on the kwargs a tellability step would refuse whatever the record held, so a
    wiring error raises when the job is built instead of reading `not computed`."""
    k, gap_mult = float(carve_kwargs["k"]), float(carve_kwargs["gap_mult"])
    if not (math.isfinite(k) and k > 0.0):
        raise ValueError(
            f"{record}: the carve's k must be finite and positive, got {k}"
        )
    if not (math.isfinite(gap_mult) and gap_mult > 0.0):
        raise ValueError(
            f"{record}: the carve's gap_mult must be finite and positive, got {gap_mult}"
        )
    labels = [str(row[0]) for row in carve_kwargs["thresholds"]]
    if not labels or len(set(labels)) != len(labels):
        raise ValueError(f"{record}: the carve's rungs {labels} are empty or repeated")
    if int(n_permutations) < 1 or int(seed) < 0:
        raise ValueError(
            f"{record}: tellability needs n_permutations >= 1 and seed >= 0, got "
            f"{n_permutations} and {seed}"
        )


def wire_resolvability_block(job: Job, nodes: RecordNodes) -> object:
    """The resolvability block: the record's tellability read for whether its reads resolve
    the events the windows count. Draws either metric; returns the block node, which the
    job's resolvability tables read."""
    block = job.step(
        build_resolvability_block,
        nodes.tellability,
        name=f"resolvability_block_{nodes.record}",
        flag=nodes.flag,
        display=nodes.display,
        paper_name=nodes.record,
    )
    job.figure(
        ResolvabilityBlockPlot,
        block,
        targets=TARGETS,
        title=f"{nodes.record} resolvability",
    )
    return block


def fidelity_panel_data_of(
    result: FidelityResult | RecordFailure,
    carve: WindowsResult | RecordFailure,
    km: KaplanMeierSet | RecordFailure,
    *,
    xi_seed: int,
    thresholds: tuple[tuple[str, float, bool], ...],
) -> WithinCalibrationPanelData | RecordFailure:
    """The record's infidelity panel data, built as the Ramsey jobs build it
    (`recipes._fidelity_panel_data`, which reads k and use_uncertainty off the carve); the
    first failure of its inputs passes through."""
    for upstream in (result, carve, km):
        if isinstance(upstream, RecordFailure):
            return upstream
    if not (
        isinstance(result, FidelityResult)
        and isinstance(carve, WindowsResult)
        and isinstance(km, KaplanMeierSet)
    ):
        raise TypeError(
            "infidelity panel data needs a fidelity result, its carve and its "
            "Kaplan-Meier set"
        )
    return _fidelity_panel_data(
        result, carve, km, xi_seed=xi_seed, thresholds=list(thresholds)
    )


def wire_timeline_block(job: Job, nodes: RecordNodes) -> object:
    """The timeline block: read state by rung from the record's panel data, built as the
    record's own metric job builds it (`panel_data_of` for T2*, `fidelity_panel_data_of`
    for the infidelity) on the carve held to the ledger, the culled rungs named."""
    km = job.step(
        kaplan_meier_of,
        nodes.carve,
        name=f"kaplan_meier_panel_{nodes.record}",
        conf_level=CONF_LEVEL,
    )
    thresholds = tuple(
        (str(lb), float(v), bool(b)) for lb, v, b in nodes.carve_kwargs["thresholds"]
    )
    if nodes.metric == "t2star":
        panel = job.step(
            panel_data_of,
            nodes.fit,
            nodes.carve,
            km,
            name=f"panel_data_{nodes.record}",
            shape_min_reads=SHAPE_MIN_READS,
            xi_seed=XI_SEED,
            k=nodes.carve_kwargs["k"],
            use_uncertainty=nodes.carve_kwargs["use_uncertainty"],
            thresholds=thresholds,
        )
    elif nodes.metric == "infidelity":
        panel = job.step(
            fidelity_panel_data_of,
            nodes.fit,
            nodes.carve,
            km,
            name=f"panel_data_{nodes.record}",
            xi_seed=XI_SEED,
            thresholds=thresholds,
        )
    else:
        raise ValueError(f"{nodes.record}: no panel data step for {nodes.metric!r}")
    block = job.step(
        build_timeline_block,
        panel,
        name=f"timeline_block_{nodes.record}",
        record=nodes.run_name,
        flag=nodes.flag,
        display=nodes.display,
    )
    job.figure(
        TimelineBlockPlot, block, targets=TARGETS, title=f"{nodes.record} timelines"
    )
    return block


def wire_resolvability_tables(job: Job, blocks: Sequence[object]) -> None:
    """The job's two resolvability tables (`RESOLVABILITY_TABLES`), every drawn record's
    block stacked, written by the `table` target beside the curve summaries. Nothing is
    wired when no block was."""
    if not blocks:
        return
    for part in RESOLVABILITY_TABLES:
        table = job.step(
            resolvability_table,
            *blocks,  # type: ignore[arg-type]
            name=f"resolvability_{part}_table",
            part=part,
        )
        job.figure(
            ResolvabilityTablePlot,
            table,
            targets=TABLE_TARGETS,
            title=f"{job.name} resolvability {part.replace('_', ' ')}",
        )


# The metric of each tellability step a single-record job writes tables for.
_TELLABILITY_METRIC = {
    _t2star_tellability: "t2star",
    _frequency_tellability: "frequency",
    _infidelity_tellability: "infidelity",
}


def _single_record_names(job: Job, tellability: object) -> tuple[str, str]:
    """The `display` and `paper_name` of a single-record job's tables, read off its graph:
    the metric of its tellability step, and the paper's name for the record its Dataset
    declares (device, qubit and `run_date`, by `infidelity_recipe.paper_record_name`)."""
    node = job.dag[tellability.node_id]  # type: ignore[attr-defined]
    if node.fn not in _TELLABILITY_METRIC:
        raise ValueError(f"{node.node_id!r} is not a tellability step")
    load = node
    while "dataset" not in load.kwargs:
        load = job.dag[load.inputs[0].node_id]
    dataset = load.kwargs["dataset"]
    if "run_date" not in dataset.extra:
        raise ValueError(f"{job.name}: its record's Dataset declares no run_date")
    name = paper_record_name(
        str(dataset.device),
        int(dataset.qubit),
        date.fromisoformat(str(dataset.extra["run_date"])),
    )
    return f"{name}, {METRIC_LABELS[_TELLABILITY_METRIC[node.fn]]}", name


def wire_tellability_tables(job: Job, tellability: object) -> object:
    """The resolvability tables of a single-record job's tellability, as the figure jobs
    write them: the block built on `tellability` (a `wire_*_tellability` node), then its two
    tables by the `table` target. No figure is drawn. The block's `display` and
    `paper_name` are read off the job's graph (`_single_record_names`). Returns the block
    node."""
    display, paper_name = _single_record_names(job, tellability)
    block = job.step(
        build_resolvability_block,
        tellability,
        name="resolvability_block",
        flag="",
        display=display,
        paper_name=paper_name,
    )
    wire_resolvability_tables(job, [block])
    return block


BLOCKS: tuple[Callable[[Job, RecordNodes], object], ...] = (
    wire_survival_block,
    wire_hazard_block,
    wire_fraction_block,
    wire_resolvability_block,
    wire_timeline_block,
)
# The page's rows, top to bottom, in the method's reading order: can the reads resolve the
# crossings, what the reads did, how much time was in spec, how long windows last, how the
# hazard grows. Every block of `BLOCKS` has its place here.
PAGE_ORDER: tuple[Callable[[Job, RecordNodes], object], ...] = (
    wire_resolvability_block,
    wire_timeline_block,
    wire_fraction_block,
    wire_survival_block,
    wire_hazard_block,
)


def wire_record_summaries(job: Job, nodes: RecordNodes) -> object:
    """One record's curve summaries with their check outcomes and intervals; returns the
    node the job's summaries table reads."""
    summaries = job.step(
        curve_summaries_of,
        nodes.carve,
        name=f"curve_summaries_{nodes.record}",
        conf_level=CONF_LEVEL,
    )
    return job.step(
        build_record_summaries,
        summaries,
        nodes.ledger,
        name=f"record_summaries_{nodes.record}",
        run_set=nodes.run_set,
        conf_level=CONF_LEVEL,
        flag=nodes.flag,
        paper_name=nodes.record,
    )


def _wire_record(job: Job, nodes: RecordNodes) -> tuple[object, object | None]:
    """Every block of one record, its page, then its curve summaries; returns the summaries
    node and the resolvability block node, which the job's tables read."""
    wired = {block: block(job, nodes) for block in BLOCKS}
    # `BLOCKS` and `PAGE_ORDER` hold the same blocks (tests/test_record_figure.py); a
    # registry swapped for a test draws no page.
    if all(block in wired for block in PAGE_ORDER):
        page = job.step(
            assemble_record_figure,
            *(wired[block] for block in PAGE_ORDER),  # type: ignore[arg-type]
            name=f"record_figure_{nodes.record}",
            record=nodes.run_name,
        )
        job.figure(
            RecordFigurePlot,
            page,
            targets=TARGETS,
            title=f"{nodes.record} record figure",
        )
    return wire_record_summaries(job, nodes), wired.get(wire_resolvability_block)


def _wire_job_tables(
    job: Job,
    summaries: Sequence[object],
    resolvability: Sequence[object | None],
    *,
    excluded: tuple[str, ...],
) -> None:
    """The job's tables: the curve summaries of its records, then the records outside the
    paper's tables as `excluded`, and the two resolvability tables."""
    table = job.step(
        summaries_table,
        *summaries,  # type: ignore[arg-type]
        name="curve_summaries_table",
        excluded=excluded,
    )
    job.figure(
        SummariesTablePlot,
        table,
        targets=TABLE_TARGETS,
        title=f"{job.name} curve summaries",
    )
    wire_resolvability_tables(job, [node for node in resolvability if node is not None])


@dataclass(frozen=True)
class _LedgerRecord:
    dataset: Dataset
    carve_kwargs: dict[str, Any]


def _ledger_records(
    included: Any,
) -> tuple[dict[str, _LedgerRecord], tuple[RowKey, ...]]:
    """Each record of an included ledger job, in the job's order, with its `Dataset` and
    carve kwargs, and the run-set its ledgers ran, all read off the job's graph."""
    job = included.job
    bundles = [n for n in job.dag.values() if n.fn is collect_carves]
    if len(bundles) != 1:
        raise ValueError(
            f"included job {included.alias!r} has {len(bundles)} `collect_carves` nodes, "
            "not one; it is not a ledger job"
        )
    (bundle,) = bundles
    records: dict[str, _LedgerRecord] = {}
    for parent, record in zip(bundle.inputs, bundle.kwargs["records"], strict=True):
        windows = job.dag[parent.node_id]
        if (
            windows.kwargs.get("stage") != "windows"
            or windows.kwargs.get("record") != record
        ):
            raise ValueError(
                f"included job {included.alias!r}: {record!r} is not carved by the T2* "
                "carve's windows stage"
            )
        load = windows
        while "dataset" not in load.kwargs:
            (up,) = load.inputs
            load = job.dag[up.node_id]
        records[str(record)] = _LedgerRecord(
            dataset=load.kwargs["dataset"],
            carve_kwargs={key: windows.kwargs[key] for key in CARVE_KWARGS},
        )
    switches = {n.kwargs["include_c3"] for n in job.dag.values() if n.fn is ledger_step}
    if len(switches) != 1:
        raise ValueError(
            f"included job {included.alias!r}: its ledgers ran with include_c3 "
            f"{sorted(switches)}, not one value"
        )
    return records, ledger_run_set(include_c3=switches.pop())


def configure_record_figure_job(
    job: Job,
    *,
    included: Any,
    device: str,
    paper_records: Sequence[tuple[str, str]],
    shaded: Sequence[str],
    flag_text: str,
    tellability_permutations: int = TELLABILITY_PERMUTATIONS,
) -> None:
    """Every block of every listed record of the ledger job `included` (the return of
    `job.include`, written in the job file so the job manifest lists it), one figure per
    block per record (module docstring). `device` picks the rows of `paper_records` that are
    this job's. `tellability_permutations` is the standalone jobs' count unless a test
    lowers it."""
    ledger_job = included.alias
    records, run_set = _ledger_records(included)
    labels, bundle = included_ledger_set(included)
    if labels != [record.replace("_", " ") for record in records]:
        raise ValueError(
            f"{ledger_job}: its ledger bundle lists {labels}, its carves {list(records)}"
        )
    prefixes = {record.rsplit("_", 1)[0] for record in records}
    listed = [
        record
        for table_device, record in paper_records
        if table_device == device and record.rsplit("_", 1)[0] in prefixes
    ]
    missing = [record for record in listed if record not in records]
    if missing:
        raise ValueError(
            f"{ledger_job}: the paper's tables list {missing}, which it does not hold"
        )
    summaries: list[object] = []
    resolvability: list[object] = []
    for record, held in records.items():
        if record not in listed:
            continue
        fit, windows = wire_t2star_carve(
            job,
            dataset=held.dataset,
            node_suffix=f"_{record}",
            on_failure="record",
            **held.carve_kwargs,
        )
        ledger = job.step(ledger_of, bundle, name=f"ledger_{record}", record=record)
        carve = job.step(checked_carve, windows, ledger, name=f"checked_carve_{record}")
        check_tellability_kwargs(
            record,
            held.carve_kwargs,
            n_permutations=tellability_permutations,
            seed=TELLABILITY_SEED,
        )
        tellability = job.step(
            t2star_tellability_of,
            fit,
            windows,
            name=f"tellability_t2star_{record}",
            thresholds=held.carve_kwargs["thresholds"],
            gap_mult=held.carve_kwargs["gap_mult"],
            unsure_k=held.carve_kwargs["k"],
            n_permutations=tellability_permutations,
            seed=TELLABILITY_SEED,
        )
        nodes = RecordNodes(
            record=record,
            run_name=record,
            carve=carve,
            ledger=ledger,
            run_set=run_set,
            flag=flag_text if record in shaded else "",
            fit=fit,
            windows=windows,
            carve_kwargs=held.carve_kwargs,
            metric="t2star",
            tellability=tellability,
        )
        summary, block = _wire_record(job, nodes)
        summaries.append(summary)
        resolvability.append(block)
    _wire_job_tables(
        job,
        summaries,
        resolvability,
        excluded=tuple(record for record in records if record not in listed),
    )


def configure_infidelity_record_figure_job(
    job: Job,
    *,
    included: Any,
    record: str,
    paper_records: Sequence[tuple[str, str]],
    shaded: Sequence[str],
    flag_text: str,
    bench_size_table: Dataset,
    ledger_knobs: Mapping[str, object] = LEDGER_KNOBS,
    tellability_permutations: int = TELLABILITY_PERMUTATIONS,
) -> None:
    """Every block of the infidelity record of the Ramsey job `included` (the return of
    `job.include`), one figure per block, and the job's tables (module docstring). `record`
    is the record's name in the paper's tables; it must be the name the job's declared Rabi
    row gives and a row of `paper_records`. `ledger_knobs` and `tellability_permutations`
    are the ledger jobs' and the standalone jobs' unless a test lowers them."""
    if set(ledger_knobs) != set(LEDGER_KNOBS):
        raise ValueError(
            f"ledger_knobs must name {sorted(LEDGER_KNOBS)}, got {sorted(ledger_knobs)}"
        )
    infidelity = wire_infidelity_carve(job, included)
    if infidelity.paper_name != record:
        raise ValueError(
            f"{included.alias}: its declared Rabi row names the record "
            f"{infidelity.paper_name!r}, not {record!r}"
        )
    if (infidelity.device, record) not in set(paper_records):
        raise ValueError(
            f"{record!r} is not a {infidelity.device} record of the paper's tables"
        )
    result, windows = infidelity.result, infidelity.windows
    carve_kwargs = infidelity.carve_kwargs
    run_name = infidelity.run_name
    ledger = job.step(
        ledger_step,
        windows,
        job.load_df(bench_size_table),
        name=f"check_ledger_{record}",
        thresholds=carve_kwargs["thresholds"],
        **ledger_knobs,
    )
    job.materialize(ledger, name=f"{record}{LEDGER_DATA_SUFFIX}")
    carve = job.step(checked_carve, windows, ledger, name=f"checked_carve_{record}")
    check_tellability_kwargs(
        record,
        carve_kwargs,
        n_permutations=tellability_permutations,
        seed=TELLABILITY_SEED,
    )
    tellability = job.step(
        infidelity_tellability_of,
        result,
        windows,
        name=f"tellability_infidelity_{record}",
        thresholds=carve_kwargs["thresholds"],
        gap_mult=carve_kwargs["gap_mult"],
        n_permutations=tellability_permutations,
        seed=TELLABILITY_SEED,
    )
    nodes = RecordNodes(
        record=record,
        run_name=run_name,
        carve=carve,
        ledger=ledger,
        run_set=ledger_run_set(include_c3=bool(ledger_knobs["include_c3"])),
        flag=flag_text if record in shaded else "",
        fit=result,
        windows=windows,
        carve_kwargs=carve_kwargs,
        metric="infidelity",
        tellability=tellability,
    )
    summary, block = _wire_record(job, nodes)
    _wire_job_tables(job, [summary], [block], excluded=())
