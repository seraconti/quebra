"""Tellability wiring: where it reaches, and where it must not.

Oracles:
- the identity closure, computed by the same `code_closure` that keys every run: a check
  ledger job asks for no tellability, so none of its modules may be in that closure, or
  every ledger would re-run whenever tellability changes;
- the job graphs: every tellability node reads the windows its own job carved, at the same
  thresholds and under the same gap rule, is told the floor its metric was clipped at, and
  is materialised under its job's prefix;
- the step functions: the record's id and the load's counts reach the artifact, and a
  missing id is refused.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from quebra import tellability_recipe
from quebra.analyzers import fidelity
from quebra.cli import _module_from_path
from quebra.core.closure import code_closure
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.core.paths import repo_root
from quebra.recipes import configure_ramsey_job

pytestmark = pytest.mark.policy

TELLABILITY_MODULES = {
    "quebra.tellability_recipe",
    "quebra.analyzers.tellability",
    "quebra.analyzers.spectrum",
}


def _module(stem: str):
    return _module_from_path(Path(repo_root() / "jobs" / "active" / f"{stem}.py"))


def _job(stem: str):
    return _module(stem).job


@pytest.mark.parametrize("stem", ["check_ledger_6d2s_q1", "km_with_checks_6d2s"])
def test_jobs_without_tellability_do_not_carry_it_in_their_identity(stem) -> None:
    """Oracle: the identity closure itself. `recipes.py` must never import the tellability
    wiring, or these jobs' identities would move with it."""
    job = _job(stem)
    closure = set(code_closure([node.fn for node in job.dag.values()]))
    assert "quebra.recipes" in closure, "control: the closure reaches the recipes"
    assert not (closure & TELLABILITY_MODULES), sorted(closure & TELLABILITY_MODULES)


@pytest.mark.parametrize(
    "stem, node, windows_node",
    [
        ("t2star_q1_070423", "tellability_t2star", "windows"),
        ("ramsey_2x2_q1_030723", "tellability_infidelity", "fidelity_windows"),
    ],
)
def test_tellability_reads_its_own_jobs_windows(stem, node, windows_node) -> None:
    """Oracle: the job graph. The tellability step consumes the carve node of its own job,
    and passes the carve's thresholds and gap rule, so the run's own check that the two
    agree is exercised on every run."""
    job = _job(stem)
    told = job.dag[node]
    carve = job.dag[windows_node]
    assert any(getattr(ref, "node_id", None) == windows_node for ref in told.inputs)
    assert told.kwargs["thresholds"] == carve.kwargs["thresholds"]
    assert told.kwargs["gap_mult"] == carve.kwargs["gap_mult"]


@pytest.mark.parametrize(
    "stem, node",
    [
        ("t2star_q1_070423", "tellability_t2star"),
        ("t2star_q1_100423", "tellability_t2star"),
        ("ramsey_q1_100423", "tellability_frequency"),
        ("ramsey_2x2_q1_030723", "tellability_infidelity"),
    ],
)
def test_each_tellability_node_is_materialised_under_its_jobs_prefix(
    stem, node
) -> None:
    """Oracle: the job file's own PREFIX. The artifact is found by that name."""
    module = _module(stem)
    names = {
        sink.name
        for sink in module.job.sinks
        if getattr(getattr(sink, "node", None), "node_id", None) == node
    }
    assert names == {f"{module.PREFIX}_tellability"}


def test_the_infidelity_step_is_told_the_floor_fidelity_clips_at() -> None:
    """Oracle: `fidelity.INFIDELITY_FLOOR`, the one constant both clip sites use."""
    told = _job("ramsey_2x2_q1_030723").dag["tellability_infidelity"]
    assert told.kwargs["clip_floor"] == fidelity.INFIDELITY_FLOOR


def test_a_value_key_other_than_the_default_reaches_the_step() -> None:
    """Oracle: the caller's argument. A frequency tellability asked to read `delta_hz` is
    labelled with it, under the prefix it was given."""
    job = Job("synthetic_tellability_wiring")
    refs = configure_ramsey_job(
        job,
        Dataset(path="data/none.pickle", qubit=1, device="6D2S"),
        profile="overnight",
        include_fidelity=False,
    )
    node = tellability_recipe.wire_frequency_tellability(
        job, refs, prefix="unit", value_key="delta_hz"
    )
    assert job.dag[node.node_id].kwargs["value_key"] == "delta_hz"
    materialised = [
        s.name
        for s in job.sinks
        if getattr(getattr(s, "node", None), "node_id", None) == node.node_id
    ]
    assert materialised == ["unit_tellability"]


def _frequency_norm(meta):
    rng = np.random.default_rng(3)
    t = np.r_[0.0, np.cumsum(rng.uniform(10.0, 40.0, size=299))]
    return {
        "t_rel_s": t,
        "raw_frequency_hz": 2.1e6 + rng.normal(scale=1e3, size=t.size),
        "delta_hz": rng.normal(scale=1e3, size=t.size),
        "frequency_error_hz": np.full(t.size, 500.0),
        "meta": meta,
    }


def test_the_frequency_step_names_its_record_and_repeats_the_load_counts() -> None:
    """Oracle: the norm's own meta. The record id and the load's counts reach the
    artifact; the value read is the one asked for."""
    meta = {
        "dataset_id": "unit",
        "n_rows_dropped_by_schema": 2,
        "n_rows_dropped_invalid": 0,
    }
    result = tellability_recipe._frequency_tellability(
        _frequency_norm(meta),
        value_key="delta_hz",
        sigma_key="frequency_error_hz",
        gap_mult=10.0,
        n_permutations=19,
        seed=1,
    )
    assert result.dataset_id == "unit"
    assert result.metric.value_key == "delta_hz"
    metric = result.metric
    assert (metric.n_rows_dropped_by_schema, metric.n_rows_dropped_invalid) == (2, 0)


@pytest.mark.parametrize(
    "meta", [{}, {"dataset_id": ""}, None], ids=["absent", "empty", "no-meta"]
)
def test_a_frequency_record_with_no_id_is_refused(meta) -> None:
    """Oracle: specification, as fidelity refuses it: an artifact must name its record."""
    with pytest.raises(KeyError, match="dataset_id"):
        tellability_recipe._frequency_tellability(
            _frequency_norm(meta),
            value_key="raw_frequency_hz",
            sigma_key="frequency_error_hz",
            gap_mult=10.0,
            n_permutations=19,
            seed=1,
        )


@pytest.mark.parametrize("metric", ["t2star", "infidelity"])
def test_the_t2star_and_infidelity_steps_repeat_the_load_counts(metric) -> None:
    """Oracle: the result meta the load counts came through. Both steps hand them to the
    artifact, so it says how many of the file's rows never reached it."""
    from quebra.analyzers import t2star, windows

    rng = np.random.default_rng(5)
    n = 300
    t = np.r_[0.0, np.cumsum(rng.uniform(10.0, 40.0, size=n - 1))]
    counts = {"n_rows_dropped_by_schema": 4, "n_rows_dropped_invalid": 1}
    norm = {
        "t_rel_s": t,
        "T2star_s": 3e-6 + 5e-7 * rng.normal(size=n),
        "T2star_error_s": np.full(n, 2e-7),
        "delta_hz": rng.normal(scale=2e4, size=n),
        "rabi_hz": np.full(n, 704e3),
        "meta": {"dataset_id": "unit", "rabi_source": "declared", **counts},
    }
    if metric == "t2star":
        result = t2star.run(t2star.make_inputs_from_norm(norm))
        thresholds = [("u", 3e-6, True)]
        carved = windows.run(
            windows.WindowsInputs(
                t_rel_s=t,
                values=result.frame["t2star_s"].to_numpy(),
                thresholds=thresholds,
                gap_mult=10.0,
                dataset_id="unit",
            )
        )
        told = tellability_recipe._t2star_tellability(
            result,
            carved,
            thresholds=thresholds,
            gap_mult=10.0,
            unsure_k=1.0,
            n_permutations=19,
            seed=1,
        )
    else:
        result = fidelity.run(fidelity.make_inputs_from_norm(norm, {}))
        thresholds = [("u", 1e-3, False)]
        carved = windows.run(
            windows.WindowsInputs(
                t_rel_s=result.frame["t_rel_s"].to_numpy(),
                values=fidelity.panel_series(result),
                thresholds=thresholds,
                gap_mult=10.0,
                dataset_id="unit",
            )
        )
        told = tellability_recipe._infidelity_tellability(
            result,
            carved,
            thresholds=thresholds,
            gap_mult=10.0,
            clip_floor=fidelity.INFIDELITY_FLOOR,
            n_permutations=19,
            seed=1,
        )
    metric_row = told.metric
    assert (metric_row.n_rows_dropped_by_schema, metric_row.n_rows_dropped_invalid) == (
        4,
        1,
    )
