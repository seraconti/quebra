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
appended there.

A module of its own, not part of `quebra.recipes`, for the reason `ledger_recipe.py` gives:
what a module imports joins the identity closure of every step function it defines.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from quebra.analyzers.check_ledger import CheckLedger
from quebra.analyzers.check_outcome import RowKey
from quebra.analyzers.kaplan_meier import CONF_LEVEL
from quebra.analyzers.record_figure import (
    build_hazard_block,
    build_survival_block,
    cross_check_carve,
    kaplan_meier_of,
    ledger_run_set,
    nelson_aalen_of,
)
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
from quebra.plots.record_figure_plot import HazardBlockPlot, SurvivalBlockPlot
from quebra.recipes import wire_t2star_carve

TARGETS = ["static", "academic"]
# The carve kwargs a ledger job's windows node carries, which the figure job carves with.
CARVE_KWARGS = ("gap_mult", "k", "use_uncertainty", "thresholds")


@dataclass(frozen=True)
class RecordNodes:
    """What every block of one record is wired from: the record's run name, its carve held
    to its ledger, its ledger, the run-set the ledger ran, and its flag text (empty unless
    the paper's tables shade it).

    `fit` and `windows` are the carve's own fit and windows nodes, and `carve_kwargs` the
    kwargs it was carved with. A step that must trace to observed reads alone, as
    `tests/test_windows_not_interpolated.py` holds a tellability step to, is wired on these:
    `carve` descends from the ledger bundle too. The cross-check still stands behind it,
    since the runner computes every node before any sink."""

    record: str
    carve: object
    ledger: object
    run_set: tuple[RowKey, ...]
    flag: str
    fit: object
    windows: object
    carve_kwargs: Mapping[str, Any]


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


BLOCKS: tuple[Callable[[Job, RecordNodes], None], ...] = (
    wire_survival_block,
    wire_hazard_block,
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
) -> None:
    """Every block of every listed record of the ledger job `included` (the return of
    `job.include`, written in the job file so the job manifest lists it), one figure per
    block per record (module docstring). `device` picks the rows of `paper_records` that are
    this job's."""
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
        nodes = RecordNodes(
            record=record,
            carve=carve,
            ledger=ledger,
            run_set=run_set,
            flag=flag_text if record in shaded else "",
            fit=fit,
            windows=windows,
            carve_kwargs=held.carve_kwargs,
        )
        for block in BLOCKS:
            block(job, nodes)
