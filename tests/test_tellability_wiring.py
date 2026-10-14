"""Tellability wiring: where it reaches, and where it must not.

Oracles:
- the identity closure, computed by the same `code_closure` that keys every run: a check
  ledger job asks for no tellability, so none of its modules may be in that closure, or
  every ledger would re-run whenever tellability changes;
- the job graphs: every tellability node reads the windows its own job carved, at the same
  thresholds and under the same gap rule.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from quebra.cli import _module_from_path
from quebra.core.closure import code_closure
from quebra.core.paths import repo_root

pytestmark = pytest.mark.policy

TELLABILITY_MODULES = {
    "quebra.tellability_recipe",
    "quebra.analyzers.tellability",
    "quebra.analyzers.spectrum",
}


def _job(stem: str):
    path = repo_root() / "jobs" / "active" / f"{stem}.py"
    return _module_from_path(Path(path)).job


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
