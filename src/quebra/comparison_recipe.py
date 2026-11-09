"""Comparison wiring: records carved on one ladder, then a pair or a forest of them.

Apart from `recipes.py`, which never imports it, so `analyzers/record_comparison.py` stays
out of the identity closure of every job that compares no records, the reason
`tellability_recipe.py` gives. A comparison job carves its records with
`wire_record_carves`, then calls `wire_record_pair` or `wire_record_forest` on the carves.

A record is named by its file stem, which is also its carve's `dataset_id` and the name the
halves rule reads (`record_comparison.HALF_SUFFIXES`). Every setting that changes a number
is a step kwarg, so it reaches the identity and the provenance label: the ladder, the tau
rule, the tau scope, the confidence level and a forest's groups. The evidence line of a
pair or of each group, the flag a record may carry and how a pair was declared change no
number, so they are kwargs of the table and figure-data steps alone.

Given `ledgers=`, a ledger per compared record, `wire_record_pair` and `wire_record_forest`
put one more step between the comparison and its table and figures:
`comparison_checks.pairs_with_checks` or `forest_with_checks`, which fills each record's
check outcome and group-row verdicts from its own ledger and holds its carve to it. Without
it the graph is the one it was.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

from quebra.analyzers import comparison_checks, record_comparison
from quebra.analyzers.check_outcome import PERMUTATION_KEYS
from quebra.analyzers.event_table import SIDES
from quebra.analyzers.kaplan_meier import CONF_LEVEL
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.plots.record_comparison_plot import (
    ForestPlot,
    PairLadderPlot,
    PooledLadderPlot,
)
from quebra.recipes import T2STAR_THRESHOLDS, wire_t2star_carve
from quebra.schemas.track912 import track912Schema

FIGURE_TARGETS = ("static", "academic")


def record_dataset(stem: str) -> Dataset:
    """One embargoed record as a Dataset, named by its file stem `DDMMYY_<device>_qubit<N>`.

    The path is `data/real_private/<device>/<stem>.pickle` (AGENTS.md, Data) and the schema
    `track912Schema`. No `run_name` is set, so the stem names the record, as
    `wire_record_carves` requires. A stem of another shape raises.
    """
    _, device, qubit, *_ = stem.split("_")
    if not qubit.startswith("qubit"):
        raise ValueError(f"{stem!r} is not a DDMMYY_<device>_qubit<N> file stem")
    return Dataset(
        path=f"data/real_private/{device}/{stem}.pickle",
        schema=track912Schema,
        qubit=int(qubit.removeprefix("qubit")),
        device=device,
    )


def wire_record_carves(
    job: Job,
    datasets: Sequence[Dataset],
    *,
    thresholds: list[tuple[str, float, bool]] | None = None,
) -> dict[str, object]:
    """The T2* carve of every record on one ladder; return the carve nodes by record name.

    `thresholds` is the one ladder, `recipes.T2STAR_THRESHOLDS` when None. A `run_name`
    other than the file stem would rename the record and hide a half from the halves rule,
    so it raises at graph build, as a record given twice does.
    """
    carves: dict[str, object] = {}
    for dataset in datasets:
        name = Path(dataset.path).stem
        run_name = dataset.extra.get("run_name", name)
        if run_name != name:
            raise ValueError(
                f"{name}: a compared record is named by its file stem; got run_name "
                f"{run_name!r}"
            )
        if name in carves:
            raise ValueError(f"{name} is given twice")
        _, carves[name] = wire_t2star_carve(
            job, dataset=dataset, thresholds=thresholds, node_suffix=f"_{name}"
        )
    return carves


def wire_record_pair(
    job: Job,
    carves: Mapping[str, object],
    *,
    name: str,
    evidence: str,
    tau_rule: str,
    conf_level: float = CONF_LEVEL,
    flags: Mapping[str, str] | None = None,
    ledgers: Mapping[str, object] | None = None,
    declared: str = "",
) -> dict[str, object]:
    """Two records compared at every threshold and side, and their table; both materialised.

    `carves` holds the two records, a first. `evidence` is the table's evidence line and
    `flags` its per-record flags (`record_comparison.pair_table`); a blank line, or a flag
    on a record not in the pair, raises at graph build. `ledgers`, when given, holds each
    record's check ledger under the same names as `carves`; the pair is then filled from
    them (module docstring), and the filled pair is what the table reads and what is
    materialised. `declared`, when given, says how the pair was chosen. Returns the `pair`
    and `table` nodes.
    """
    if len(carves) != 2:
        raise ValueError(f"a pair is two records; got {list(carves)}")
    record_comparison.require_no_record_beside_its_halves(list(carves))
    record_comparison.require_evidence({name: evidence}, [name])
    record_comparison.require_flags(dict(flags or {}), list(carves))
    carve_a, carve_b = carves.values()
    pair = job.step(
        record_comparison.compare_pair_set,
        carve_a,  # type: ignore[arg-type]
        carve_b,  # type: ignore[arg-type]
        name=f"{name}_pair",
        tau_rule=tau_rule,
        conf_level=conf_level,
    )
    if ledgers is not None:
        if set(ledgers) != set(carves):
            raise ValueError(
                f"ledgers of {sorted(ledgers)} given for the pair {list(carves)}"
            )
        pair = job.step(
            comparison_checks.pairs_with_checks,
            pair,
            carve_a,
            carve_b,
            *(ledgers[record] for record in carves),
            name=f"{name}_pair_checks",
            run_set=PERMUTATION_KEYS,
        )
    table = job.step(
        record_comparison.pair_table,
        pair,
        name=f"{name}_table",
        comparison=name,
        evidence=evidence,
        flags=dict(flags or {}),
        **({"declared": declared} if declared else {}),  # type: ignore[arg-type]
    )
    job.materialize(pair, name=f"{name}_pair")
    job.materialize(table, name=f"{name}_table")
    return {"pair": pair, "table": table}


def wire_record_forest(
    job: Job,
    carves: Mapping[str, object],
    *,
    name: str,
    groups: Mapping[str, Sequence[str]],
    evidence: Mapping[str, str],
    tau_scope: str,
    tau_rule: str,
    conf_level: float = CONF_LEVEL,
    flags: Mapping[str, str] | None = None,
    ledgers: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """A random-effects summary per group at every threshold and side, and their table.

    Both are materialised. Only the carves `groups` names are wired. Groups that
    `record_comparison.require_disjoint_groups` refuses, or that name a record with no
    carve, raise at graph build, as do `evidence` without one line per group and a flag
    on a record no group names. `ledgers`, when given, holds the check ledger of every
    record the groups name; the forest is then filled from them (module docstring), and
    the filled forest is what the table reads and what is materialised. Returns the
    `forest` and `table` nodes.
    """
    everyone = record_comparison.require_disjoint_groups(groups)
    record_comparison.require_evidence(evidence, list(groups))
    record_comparison.require_flags(dict(flags or {}), everyone)
    unknown = sorted(set(everyone) - set(carves))
    if unknown:
        raise KeyError(f"groups name records with no carve: {unknown}")
    used = [carves[record] for record in everyone]
    forest = job.step(
        record_comparison.random_effects_set,
        *used,  # type: ignore[arg-type]
        name=f"{name}_forest",
        groups={group: list(names) for group, names in groups.items()},
        tau_scope=tau_scope,
        tau_rule=tau_rule,
        conf_level=conf_level,
    )
    if ledgers is not None:
        missing = sorted(set(everyone) - set(ledgers))
        extra = sorted(set(ledgers) - set(everyone))
        if missing or extra:
            raise KeyError(
                f"ledgers must be the groups' records: missing {missing}, extra {extra}"
            )
        forest = job.step(
            comparison_checks.forest_with_checks,
            forest,
            *used,
            *(ledgers[record] for record in everyone),
            name=f"{name}_forest_checks",
            records=tuple(everyone),
            run_set=PERMUTATION_KEYS,
        )
    table = job.step(
        record_comparison.random_effects_table,
        forest,
        name=f"{name}_table",
        comparison=name,
        evidence=dict(evidence),
        flags=dict(flags or {}),
    )
    job.materialize(forest, name=f"{name}_forest")
    job.materialize(table, name=f"{name}_table")
    return {"forest": forest, "table": table}


def wire_pair_figure(
    job: Job,
    carves: Mapping[str, object],
    pair: object,
    *,
    name: str,
    evidence: str,
    flags: Mapping[str, str] | None = None,
    rungs: Sequence[str] | None = None,
    thresholds: list[tuple[str, float, bool]] | None = None,
    targets: Sequence[str] = FIGURE_TARGETS,
    declared: str = "",
) -> object:
    """The pair figure: both sides of each drawn threshold of `pair`, beside the curves.

    `carves` and `pair` are what `wire_record_carves` and `wire_record_pair` gave for this
    pair, a first, and `evidence` and `flags` what the pair's table was given. `rungs`
    names the thresholds the figure draws, every rung of `thresholds` (the carves' ladder,
    `recipes.T2STAR_THRESHOLDS` when None) by default; one the ladder lacks raises at graph
    build. It is a kwarg of the figure-data step alone, passed only when given, so the
    table keeps every rung. `declared` is what the pair's table was given, printed in the
    caption. Returns the figure-data node.
    """
    if len(carves) != 2:
        raise ValueError(f"a pair is two records; got {list(carves)}")
    record_comparison.require_evidence({name: evidence}, [name])
    record_comparison.require_flags(dict(flags or {}), list(carves))
    drawn: dict[str, object] = {}
    if rungs is not None:
        ladder = [label for label, _, _ in (thresholds or T2STAR_THRESHOLDS)]
        record_comparison.select_rungs(ladder, rungs)
        drawn["rungs"] = list(rungs)
    if declared:
        drawn["declared"] = declared
    carve_a, carve_b = carves.values()
    data = job.step(
        record_comparison.build_pair_ladder_data,
        pair,  # type: ignore[arg-type]
        carve_a,  # type: ignore[arg-type]
        carve_b,  # type: ignore[arg-type]
        name=f"{name}_pair_data",
        comparison=name,
        evidence=evidence,
        flags=dict(flags or {}),
        **drawn,  # type: ignore[arg-type]
    )
    job.figure(PairLadderPlot, data, targets=list(targets), title=f"{name}_pair_figure")
    return data


def wire_forest_figures(
    job: Job,
    forest: object,
    *,
    name: str,
    evidence: Mapping[str, str],
    threshold_label: str,
    side: str,
    thresholds: list[tuple[str, float, bool]] | None = None,
    flags: Mapping[str, str] | None = None,
    targets: Sequence[str] = FIGURE_TARGETS,
) -> dict[str, object]:
    """The forest at one threshold and side, and the pooled ladder, of a `forest` node.

    `threshold_label` must name a rung of `thresholds`, the carves' ladder
    (`recipes.T2STAR_THRESHOLDS` when None), and `side` one of `event_table.SIDES`;
    either raises at graph build otherwise. Both are kwargs of the forest's figure-data
    step alone, so drawing another rung leaves the comparison step and its artifact unchanged. `evidence` and
    `flags` are what the forest's table was given; the figure-data step checks them
    against the groups. Returns the `forest_data` and `ladder_data` nodes.
    """
    ladder = [label for label, _, _ in (thresholds or T2STAR_THRESHOLDS)]
    if threshold_label not in ladder:
        raise ValueError(f"{threshold_label!r} is not a rung of the ladder {ladder}")
    if side not in SIDES:
        raise ValueError(f"side must be one of {SIDES}; got {side!r}")
    forest_data = job.step(
        record_comparison.build_forest_data,
        forest,  # type: ignore[arg-type]
        name=f"{name}_forest_data",
        comparison=name,
        evidence=dict(evidence),
        threshold_label=threshold_label,
        side=side,
        flags=dict(flags or {}),
    )
    ladder_data = job.step(
        record_comparison.build_pooled_ladder_data,
        forest,  # type: ignore[arg-type]
        name=f"{name}_ladder_data",
        comparison=name,
        evidence=dict(evidence),
    )
    job.figure(
        ForestPlot, forest_data, targets=list(targets), title=f"{name}_forest_figure"
    )
    job.figure(
        PooledLadderPlot,
        ladder_data,
        targets=list(targets),
        title=f"{name}_ladder_figure",
    )
    return {"forest_data": forest_data, "ladder_data": ladder_data}
