"""The infidelity carve of a Ramsey record, wired in a job that includes the Ramsey job.

`wire_infidelity_carve` is the one wiring of it. It reads the included job's graph when the
job is built: the job's one fidelity carve node (`recipes._fidelity_windows`), the fidelity
result that node reads, persisted under the job's one sink on it, and the job's one declared
Rabi row. It takes one ref, to that result, and recomputes the carve on it with the node's
own kwargs (|x| on the nines ladder, the carve's k and sigma switch, the gap rule), so the
windows are the Ramsey job's own. The record is named twice: by its run name, the
`dataset_id` every artifact of it carries, and by its name in the paper's dataset tables
(`paper_record_name`), built from the declared Rabi row's device, qubit and run date. The
name follows the tables' rule; whether the tables list the record is the caller's check.
The Rabi row is keyed on device, qubit and date only, so a record's condition suffix, if its
file name carries one, is not part of the name.

A job that is not a Ramsey job with fidelity (no fidelity carve, or more than one), or that
persists the result under no sink or several, raises when the job is built.

A module of its own so every job that carves an infidelity record wires it one way, without
importing the record figures' recipe. Its one step is `recipes._fidelity_windows`, so
`quebra.recipes` joins the identity closure of every job that wires it.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date
from typing import Any

from quebra.core.job import Job, record_id
from quebra.recipes import _fidelity_windows
from quebra.transforms.declared_rabi import join_declared_rabi

# The kwargs a carve node carries, the fidelity carve's and the T2* windows stage's alike;
# the recomputed carve takes them unchanged.
CARVE_KWARGS = ("gap_mult", "k", "use_uncertainty", "thresholds")
# The devices of the paper's dataset tables, and the prefix each puts before a record's
# qubit and date (`paper_tables.PAPER_RECORDS`: `q1_070423`, `2x2_q1_030723`).
TABLE_PREFIX = {"6D2S": "", "2x2": "2x2_"}


def paper_record_name(device: str, qubit: int, run_date: date) -> str:
    """A record's name in the paper's dataset tables: `q{qubit}_{DDMMYY}`, prefixed `2x2_`
    for the 2x2 device. Raises for a device the tables do not list."""
    if device not in TABLE_PREFIX:
        raise ValueError(f"the paper's tables list no device {device!r}")
    return f"{TABLE_PREFIX[device]}q{int(qubit)}_{run_date:%d%m%y}"


@dataclass(frozen=True)
class InfidelityCarve:
    """The nodes and names `wire_infidelity_carve` returns.

    `result` is the ref to the included job's fidelity result and `windows` the carve
    recomputed on it with `carve_kwargs`, the included carve node's own. `run_name` is the
    record's `dataset_id`; `device` and `paper_name` come from the declared Rabi row.
    """

    result: object
    windows: object
    carve_kwargs: dict[str, Any]
    run_name: str
    device: str
    paper_name: str


def _one(nodes: list[Any], what: str, alias: str) -> Any:
    if len(nodes) != 1:
        raise ValueError(
            f"included job {alias!r} has {len(nodes)} {what}, not one; it is not a Ramsey "
            "job with fidelity"
        )
    return nodes[0]


def wire_infidelity_carve(job: Job, included: Any) -> InfidelityCarve:
    """Wire the infidelity carve of the one record of the Ramsey job `included` (the return
    of `job.include`) onto `job` (module docstring). The carve node is named
    `fidelity_windows_{paper name}`."""
    graph = included.job
    carve = _one(
        [n for n in graph.dag.values() if n.fn is _fidelity_windows],
        "fidelity carves",
        included.alias,
    )
    (source,) = carve.inputs
    names = [
        sink.name
        for sink in graph.sinks
        if getattr(sink, "node", None) is not None
        and sink.node.node_id == source.node_id
    ]
    if len(names) != 1:
        raise ValueError(
            f"included job {included.alias!r} persists its fidelity result under "
            f"{len(names)} sinks, not one"
        )
    load = graph.dag[source.node_id]
    while "dataset" not in load.kwargs:
        load = graph.dag[load.inputs[0].node_id]
    rabi = _one(
        [n for n in graph.dag.values() if n.fn is join_declared_rabi],
        "declared Rabi rows",
        included.alias,
    ).kwargs
    device = str(rabi["device"])
    paper_name = paper_record_name(
        device, int(rabi["qubit"]), date.fromisoformat(str(rabi["run_date"]))
    )
    carve_kwargs = {key: carve.kwargs[key] for key in CARVE_KWARGS}
    result = included.ref(names[0])
    windows = job.step(
        _fidelity_windows, result, name=f"fidelity_windows_{paper_name}", **carve_kwargs
    )
    return InfidelityCarve(
        result=result,
        windows=windows,
        carve_kwargs=carve_kwargs,
        run_name=record_id(load.kwargs["dataset"]),
        device=device,
        paper_name=paper_name,
    )
