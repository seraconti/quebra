"""Comparison wiring: records carved on one ladder, then a pair or a forest of them.

Apart from `recipes.py`, which never imports it, so `analyzers/record_comparison.py` stays
out of the identity closure of every job that compares no records, the reason
`tellability_recipe.py` gives. A comparison job carves its records with
`wire_record_carves`, then calls `wire_record_pair` or `wire_record_forest` on the carves.

A record is named by its file stem, which is also its carve's `dataset_id` and the name the
halves rule reads (`record_comparison.HALF_SUFFIXES`). Every setting that changes a number
is a step kwarg, so it reaches the identity and the provenance label: the ladder, the tau
rule, the tau scope, the confidence level and a forest's groups. The evidence line of a
pair or of each group changes no number, so it is a kwarg of the table step alone.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

from quebra.analyzers import record_comparison
from quebra.analyzers.kaplan_meier import CONF_LEVEL
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.recipes import wire_t2star_carve
from quebra.schemas.track912 import track912Schema


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
) -> dict[str, object]:
    """Two records compared at every threshold and side, and their table; both materialised.

    `carves` holds the two records, a first. `evidence` is the table's evidence line
    (`record_comparison.pair_table`); a blank one raises at graph build. Returns the `pair`
    and `table` nodes.
    """
    if len(carves) != 2:
        raise ValueError(f"a pair is two records; got {list(carves)}")
    record_comparison.require_no_record_beside_its_halves(list(carves))
    record_comparison.require_evidence({name: evidence}, [name])
    carve_a, carve_b = carves.values()
    pair = job.step(
        record_comparison.compare_pair_set,
        carve_a,  # type: ignore[arg-type]
        carve_b,  # type: ignore[arg-type]
        name=f"{name}_pair",
        tau_rule=tau_rule,
        conf_level=conf_level,
    )
    table = job.step(
        record_comparison.pair_table,
        pair,
        name=f"{name}_table",
        comparison=name,
        evidence=evidence,
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
) -> dict[str, object]:
    """A random-effects summary per group at every threshold and side, and their table.

    Both are materialised. Only the carves `groups` names are wired. Groups that
    `record_comparison.require_disjoint_groups` refuses, or that name a record with no
    carve, raise at graph build, as does `evidence` without one line per group. Returns
    the `forest` and `table` nodes.
    """
    everyone = record_comparison.require_disjoint_groups(groups)
    record_comparison.require_evidence(evidence, list(groups))
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
    table = job.step(
        record_comparison.random_effects_table,
        forest,
        name=f"{name}_table",
        comparison=name,
        evidence=dict(evidence),
    )
    job.materialize(forest, name=f"{name}_forest")
    job.materialize(table, name=f"{name}_table")
    return {"forest": forest, "table": table}
