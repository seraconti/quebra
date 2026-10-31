"""The T2* ladder: declared once, the same for every record.

Oracles:
- issue #22: "T2* declared once, 2 to 5 µs in steps of 0.5 µs", built here from its start,
  stop and step rather than retyped, and used as the reference by every test below;
- the job graphs: every node of every declared job whose `thresholds` kwarg holds a
  big-good rung, the T2* side, carries exactly that ladder;
- the job files' source: a job-local carve hides its rung from the graph scan, so the jobs
  that could carve one are named, and they are exactly the two one-rung Kaplan-Meier
  figures, each with one rung.
"""

from __future__ import annotations

import ast
from pathlib import Path

import numpy as np
import pytest

from quebra.analyzers.t2star import T2STAR_DEFAULT_LADDER
from quebra.cli import _module_from_path
from quebra.core.discovery import discover
from quebra.core.paths import repo_root
from quebra.recipes import T2STAR_THRESHOLDS

# One-rung Kaplan-Meier figures at their own 3 µs rung: not ladder carves.
OWN_RUNG_JOBS = {"km_poster_6d2s", "km_with_checks_6d2s"}
_WINDOWS_MODULE = "quebra.analyzers.windows"


def _job_files() -> list[Path]:
    """Every declared job under jobs/, found as `quebra run --all` finds it."""
    jobs = discover(Path(repo_root()) / "jobs")
    return sorted(declared.path for declared in jobs.values())


def _issue_ladder() -> list[tuple[str, float, bool]]:
    """The ladder issue #22 declares, from its start, stop and step."""
    return [
        (f"{rung_us:g} µs", rung_us * 1e-6, True)
        for rung_us in np.arange(2.0, 5.0 + 0.25, 0.5)
    ]


def _is_issue_ladder(ladder) -> bool:
    expected = _issue_ladder()
    return len(ladder) == len(expected) and all(
        label == want_label
        and value == pytest.approx(want_value, rel=1e-15)
        and big_good is want_big_good
        for (label, value, big_good), (want_label, want_value, want_big_good) in zip(
            ladder, expected
        )
    )


@pytest.mark.unit
def test_the_t2star_ladder_is_2_to_5_us_in_steps_of_half_a_us() -> None:
    """Oracle: issue #22. Seven rungs from 2 to 5 µs, 0.5 µs apart, big values good,
    labelled in µs without a trailing zero; the recipes' ladder is the same list."""
    assert len(_issue_ladder()) == 7
    assert _is_issue_ladder(T2STAR_DEFAULT_LADDER), T2STAR_DEFAULT_LADDER
    assert T2STAR_THRESHOLDS == list(T2STAR_DEFAULT_LADDER)


@pytest.mark.policy
def test_every_t2star_node_of_every_job_carves_on_the_declared_ladder() -> None:
    """Oracle: issue #22, read off every job graph. A node whose `thresholds` holds a
    big-good rung is on the T2* side (fidelity's ladder is small-good) and must carry
    the issue's ladder. The scan must see the T2* jobs' carve, panel and tellability
    nodes and the ledgers' carves, so an empty scan cannot pass."""
    seen: dict[tuple[str, str], list] = {}
    for path in _job_files():
        job = _module_from_path(path).job
        for node_id, node in job.dag.items():
            ladder = node.kwargs.get("thresholds")
            if ladder is not None and any(bool(rung[2]) for rung in ladder):
                seen[(path.stem, node_id)] = list(ladder)
    expected_nodes = {
        (stem, node)
        for stem in ("t2star_q1_070423", "t2star_q1_100423")
        for node in ("windows", "t2star_panel_data", "tellability_t2star")
    }
    assert expected_nodes <= set(seen), sorted(expected_nodes - set(seen))
    assert any(
        stem.startswith("check_ledger_6d2s") and node.startswith("check_ledger_")
        for stem, node in seen
    ), "control: the scan reaches the ledger jobs"
    wrong = sorted(key for key, ladder in seen.items() if not _is_issue_ladder(ladder))
    assert not wrong, wrong


def _imports_the_carve(tree: ast.AST) -> bool:
    """Whether a module imports `analyzers.windows`, which a job-local carve needs."""
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            if any(alias.name == _WINDOWS_MODULE for alias in node.names):
                return True
        elif isinstance(node, ast.ImportFrom):
            if node.module == _WINDOWS_MODULE:
                return True
            if node.module == "quebra.analyzers" and any(
                alias.name == "windows" for alias in node.names
            ):
                return True
    return False


@pytest.mark.policy
def test_only_the_two_km_figures_carve_a_rung_of_their_own() -> None:
    """Oracle: the job files' source. A job can carve on a rung of its own through a
    `thresholds=` keyword, a literal list of (label, value, big_values_good) rungs, or a
    carve of its own, which needs `analyzers.windows`; such a carve hides its rung from the
    graph scan above. This scan names every such job: exactly the two Kaplan-Meier
    figures, and each declares one big-good rung."""

    def is_rung(element: ast.expr) -> bool:
        return (
            isinstance(element, ast.Tuple)
            and len(element.elts) == 3
            and isinstance(element.elts[2], ast.Constant)
            and isinstance(element.elts[2].value, bool)
        )

    declaring: dict[str, Path] = {}
    for path in _job_files():
        tree = ast.parse(path.read_text(encoding="utf-8"))
        if _imports_the_carve(tree):
            declaring[path.stem] = path
        for node in ast.walk(tree):
            if isinstance(node, ast.keyword) and node.arg == "thresholds":
                declaring[path.stem] = path
            if isinstance(node, (ast.List, ast.Tuple)) and any(map(is_rung, node.elts)):
                declaring[path.stem] = path
    assert set(declaring) == OWN_RUNG_JOBS, sorted(set(declaring) ^ OWN_RUNG_JOBS)
    for stem, path in declaring.items():
        rungs = _module_from_path(path).THRESHOLD
        assert len(rungs) == 1 and rungs[0][2] is True, (stem, rungs)
