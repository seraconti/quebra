"""The 2x2 infidelity records as the comparisons read them, off the infidelity ledger jobs.

`included_fidelity_records(included)` reads one included infidelity ledger job
(`fidelity_ledger_recipe.configure_fidelity_ledger_job`) off its graph. Per record, in the
job's order, it gives the record's run name, its name in the paper's dataset tables, and refs
to the carve and the ledger the job persists under `{record}_windows` and
`{record}_check_ledger_data`. A comparison job includes the ledger jobs and hands these
refs to `comparison_recipe.wire_record_pair` or `wire_record_forest`, with `ledgers=`, so
it carves and scores nothing: every window and every check cell is the ledger job's own.
Records are keyed by run name, the carve's `dataset_id`; the paper name travels beside it.

Each ref resolves by running its ledger job, unless a run of that job at the same identity
and commit on a clean tree can be reused (`--reuse-deps`). A comparison's own nested runs are
never reused within that run (`runner._locate_artifact`), and a comparison job holds two refs
per record, so each ref would re-run its ledger job. Run the ledger jobs first, on their own,
at the same clean commit; then the comparison job with `--reuse-deps`.

The records' ladder is the ledgers' own, read off the ledger jobs' graphs, so a comparison
job declares none: `wire_fidelity_forest_figures` names the forest's rung against it.

A module of its own so `ledger_recipe` and `comparison_recipe` keep their imports.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from quebra.comparison_recipe import wire_forest_figures
from quebra.core.job import Job
from quebra.fidelity_ledger_recipe import fidelity_ledger_step
from quebra.ledger_recipe import CARVE_SUFFIX, LEDGER_DATA_SUFFIX


@dataclass(frozen=True)
class FidelityRecord:
    """One record of an infidelity ledger job: its run name, its name in the paper's
    tables, refs to its persisted carve and ledger, and the ladder its ledger scored."""

    record: str
    paper_name: str
    carve: Any
    ledger: Any
    thresholds: tuple[tuple[str, float, bool], ...]


def included_fidelity_records(included: Any) -> list[FidelityRecord]:
    """Every record of the included infidelity ledger job `included` (the return of
    `job.include`), in the job's order (module docstring).

    Raises when the job scores no infidelity ledger, or when a record's carve or ledger is
    persisted under no sink or several, or under a name its partner does not share.
    """
    graph = included.job
    ledgers = [n for n in graph.dag.values() if n.fn is fidelity_ledger_step]
    if not ledgers:
        raise ValueError(
            f"included job {included.alias!r} scores no infidelity ledger; it is not an "
            "infidelity ledger job"
        )
    sinks: dict[str, list[str]] = {}
    for sink in graph.sinks:
        node = getattr(sink, "node", None)
        if node is not None:
            sinks.setdefault(node.node_id, []).append(sink.name)
    records = []
    for ledger in ledgers:
        carve_id = ledger.inputs[0].node_id
        ledger_names = [
            n for n in sinks.get(ledger.node_id, []) if n.endswith(LEDGER_DATA_SUFFIX)
        ]
        carve_names = [n for n in sinks.get(carve_id, []) if n.endswith(CARVE_SUFFIX)]
        if len(ledger_names) != 1 or len(carve_names) != 1:
            raise ValueError(
                f"included job {included.alias!r} persists the ledger {ledger.node_id!r} "
                f"under {ledger_names} and its carve under {carve_names}, not one each"
            )
        record = ledger_names[0].removesuffix(LEDGER_DATA_SUFFIX)
        if carve_names[0] != f"{record}{CARVE_SUFFIX}":
            raise ValueError(
                f"included job {included.alias!r} persists the ledger of {record!r} beside "
                f"the carve {carve_names[0]!r}"
            )
        records.append(
            FidelityRecord(
                record=record,
                paper_name=str(ledger.kwargs["paper_name"]),
                carve=included.ref(carve_names[0]),
                ledger=included.ref(ledger_names[0]),
                thresholds=tuple(tuple(rung) for rung in ledger.kwargs["thresholds"]),
            )
        )
    return records


def fidelity_records(included: Sequence[Any]) -> dict[str, FidelityRecord]:
    """Every record of every included infidelity ledger job, by run name, in job order.
    Raises when two jobs hold one record."""
    out: dict[str, FidelityRecord] = {}
    for sub in included:
        for record in included_fidelity_records(sub):
            if record.record in out:
                raise ValueError(f"{record.record} is held by two included jobs")
            out[record.record] = record
    return out


def wire_fidelity_forest_figures(
    job: Job,
    forest: object,
    records: Sequence[FidelityRecord],
    *,
    name: str,
    evidence: dict[str, str],
    threshold_label: str,
    side: str,
    flags: dict[str, str] | None = None,
) -> dict[str, object]:
    """`comparison_recipe.wire_forest_figures` on the ladder `records` were scored on, read
    off the ledger jobs' graphs. Raises when the records carry two ladders."""
    ladders = {record.thresholds for record in records}
    if len(ladders) != 1:
        raise ValueError(f"records scored on {len(ladders)} ladders, not one")
    (ladder,) = ladders
    return wire_forest_figures(
        job,
        forest,
        name=name,
        evidence=evidence,
        threshold_label=threshold_label,
        side=side,
        thresholds=[tuple(rung) for rung in ladder],  # type: ignore[misc]
        flags=flags,
    )
