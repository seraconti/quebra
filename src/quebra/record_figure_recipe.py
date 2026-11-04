"""Record-figure jobs: the figures of every record of one ledger job, one PDF per block.

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
`plots/record_figure_plot.py`, off the record's `RecordNodes`. A new block is one function
appended there. The job's tables are wired after the records: the curve summaries, one row
per record, rung and bag, with the records outside the tables as `excluded=`, written by
the `table` target. A record names its metric (`RecordNodes.metric`, one of
`METRICS`).

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
from typing import Any

from quebra.analyzers.check_ledger import CheckLedger
from quebra.analyzers.compliance_survey import record_compliance
from quebra.analyzers.check_outcome import RowKey
from quebra.analyzers.kaplan_meier import CONF_LEVEL
from quebra.analyzers.record_figure import (
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
    summaries_table,
)
from quebra.analyzers.t2star import T2StarResult
from quebra.analyzers.tellability import TellabilityResult
from quebra.analyzers.windows import WindowsResult
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.core.types import RecordFailure
from quebra.ledger_recipe import (
    CheckLedgers,
    collect_carves,
    included_ledger_set,
    ledger_step,
)
from quebra.plots.record_figure_plot import (
    FractionBlockPlot,
    HazardBlockPlot,
    ResolvabilityBlockPlot,
    SummariesTablePlot,
    SurvivalBlockPlot,
    TimelineBlockPlot,
)
from quebra.recipes import XI_SEED, wire_t2star_carve
from quebra.tellability_recipe import (
    TELLABILITY_PERMUTATIONS,
    TELLABILITY_SEED,
    _t2star_tellability,
)

TARGETS = ["static", "academic"]
TABLE_TARGETS = ["table"]
# The T2* panel data's shape statistic needs this many reads per window, the value
# `recipes.configure_t2star_job` uses, so the timelines read the T2* panels' own states.
SHAPE_MIN_READS = 5
# The carve kwargs a ledger job's windows node carries, which the figure job carves with.
CARVE_KWARGS = ("gap_mult", "k", "use_uncertainty", "thresholds")
# The metrics a record figure draws, as `tellability` names them.
METRICS = ("t2star", "infidelity")


@dataclass(frozen=True)
class RecordNodes:
    """What every block of one record is wired from: the record's run name, its carve held
    to its ledger, its ledger, the run-set the ledger ran, and its flag text (empty unless
    the paper's tables shade it).

    `fit` and `windows` are the carve's own fit and windows nodes, and `carve_kwargs` the
    kwargs it was carved with. A step that must trace to observed reads alone, as
    `tests/test_windows_not_interpolated.py` holds a tellability step to, is wired on these:
    `carve` descends from the ledger bundle too. The cross-check still stands behind it,
    since the runner computes every node before any sink. `metric` names what the record
    measures (`METRICS`) and `tellability` is its one tellability node."""

    record: str
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


def wire_survival_block(job: Job, nodes: RecordNodes) -> None:
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
    )
    job.figure(
        SurvivalBlockPlot, block, targets=TARGETS, title=f"{nodes.record} survival"
    )


def wire_hazard_block(job: Job, nodes: RecordNodes) -> None:
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
    )
    job.figure(
        HazardBlockPlot,
        block,
        targets=TARGETS,
        title=f"{nodes.record} cumulative hazard",
    )


def wire_fraction_block(job: Job, nodes: RecordNodes) -> None:
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
    )
    job.figure(
        FractionBlockPlot,
        block,
        targets=TARGETS,
        title=f"{nodes.record} in-spec fraction",
    )


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


def wire_resolvability_block(job: Job, nodes: RecordNodes) -> None:
    """The resolvability block: the record's tellability read for whether its reads resolve
    the events the windows count. Draws either metric."""
    block = job.step(
        build_resolvability_block,
        nodes.tellability,
        name=f"resolvability_block_{nodes.record}",
        flag=nodes.flag,
    )
    job.figure(
        ResolvabilityBlockPlot,
        block,
        targets=TARGETS,
        title=f"{nodes.record} resolvability",
    )


def wire_timeline_block(job: Job, nodes: RecordNodes) -> None:
    """The timeline block: read state by rung from the record's T2* panel data, built as
    the T2* jobs build it on the carve held to the ledger, the culled rungs named."""
    km = job.step(
        kaplan_meier_of,
        nodes.carve,
        name=f"kaplan_meier_panel_{nodes.record}",
        conf_level=CONF_LEVEL,
    )
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
        thresholds=tuple(
            (str(lb), float(v), bool(b))
            for lb, v, b in nodes.carve_kwargs["thresholds"]
        ),
    )
    block = job.step(
        build_timeline_block,
        panel,
        name=f"timeline_block_{nodes.record}",
        record=nodes.record,
        flag=nodes.flag,
    )
    job.figure(
        TimelineBlockPlot, block, targets=TARGETS, title=f"{nodes.record} timelines"
    )


BLOCKS: tuple[Callable[[Job, RecordNodes], None], ...] = (
    wire_survival_block,
    wire_hazard_block,
    wire_fraction_block,
    wire_resolvability_block,
    wire_timeline_block,
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
    )


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
        for block in BLOCKS:
            block(job, nodes)
        summaries.append(wire_record_summaries(job, nodes))
    table = job.step(
        summaries_table,
        *summaries,  # type: ignore[arg-type]
        name="curve_summaries_table",
        excluded=tuple(record for record in records if record not in listed),
    )
    job.figure(
        SummariesTablePlot,
        table,
        targets=TABLE_TARGETS,
        title=f"{job.name} curve summaries",
    )
