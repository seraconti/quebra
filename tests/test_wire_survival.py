"""The survival estimators wired as nodes, resolved on a synthetic carve."""

from __future__ import annotations

import re

import numpy as np
import pytest

from quebra import recipes
from quebra.analyzers import event_table, windows
from quebra.core.closure import parameter_row
from quebra.core.job import Job
from quebra.core.reference import LocalRef
from quebra.core.runner import _toposort

pytestmark = pytest.mark.integration

LADDER = [("4", 4.0, True), ("6", 6.0, True)]


def _synthetic_carve() -> windows.WindowsResult:
    t = np.arange(60.0)
    v = 5.0 + 3.0 * np.sin(t / 3.0)
    return windows.run(windows.WindowsInputs(t_rel_s=t, values=v, thresholds=LADDER))


def _wired():
    job = Job("survival_wiring")
    carve = job.step(_synthetic_carve, name="carve")
    tables, km = recipes.wire_kaplan_meier(job, carve, conf_level=0.9)
    refs = recipes.wire_survival(
        job, tables, rmst_tau_s=7.0, conf_level=0.9, turnbull_tol=1e-9
    )
    return job, {"kaplan_meier": km, **refs}


def _resolve(job: Job, refs: dict[str, LocalRef]) -> dict[str, object]:
    results: dict[str, object] = {}
    for node_id in _toposort(job, [ref.node_id for ref in refs.values()]):
        node = job.dag[node_id]
        inputs = [results[ref.node_id] for ref in node.inputs]
        results[node_id] = node.fn(*inputs, **node.kwargs)
    return {name: results[ref.node_id] for name, ref in refs.items()}


def test_every_estimator_node_resolves_for_every_threshold_and_side() -> None:
    """Oracle: specification. One result per (label, side), carrying the wired parameters."""
    job, refs = _wired()
    out = _resolve(job, refs)
    keys = {(label, side) for label, _, _ in LADDER for side in event_table.SIDES}

    assert set(out["kaplan_meier"].curves) == keys
    assert out["kaplan_meier"].conf_level == 0.9
    for name in ("nelson_aalen", "curve_summaries", "turnbull", "placement_bracket"):
        assert set(out[name].items) == keys, name
    for key in keys:
        summary = out["curve_summaries"].get(*key)
        bracket = out["placement_bracket"].get(*key)
        # The sine crosses both levels several times each way, so no cell is empty and
        # the parameter asserts below never read a typed-empty result.
        assert summary.n_windows > 0, key
        assert summary.conf_level == 0.9, key
        assert summary.rmst_tau_s == 7.0 and bracket.tau_s == 7.0, key
        assert out["nelson_aalen"].get(*key).conf_level == 0.9
        turnbull = out["turnbull"].get(*key)
        assert turnbull.converged and turnbull.tol == 1e-9, key


def test_every_parameter_that_changes_a_number_reaches_the_identity() -> None:
    """Oracle: `core.closure.parameter_row`, the identity's record of each step's kwargs.

    A parameter captured in a closure would change the numbers without changing the run's
    identity or its provenance label.
    """
    job, _ = _wired()
    row = "\n".join(parameter_row(list(job.dag.values())))
    for node, kwarg in (
        ("kaplan_meier", "conf_level=0.9"),
        ("nelson_aalen", "conf_level=0.9"),
        ("curve_summaries", "conf_level=0.9"),
        ("curve_summaries", "rmst_tau_s=7.0"),
        ("curve_summaries", "quantiles=(0.25, 0.5, 0.75)"),
        ("placement_bracket", "rmst_tau_s=7.0"),
        ("turnbull", "tol=1e-09"),
        ("turnbull", "max_iter=100000"),
    ):
        line = next(r for r in row.splitlines() if r.startswith(f"{node}("))
        # Bounded on both sides, so "conf_level=0.9" does not match "conf_level=0.95".
        assert re.search(rf"[(,]{re.escape(kwarg)}[,)]", line), (node, kwarg, line)
