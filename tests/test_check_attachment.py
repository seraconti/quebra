"""The leaf module that names what an estimate rests on."""

from __future__ import annotations

import ast
import subprocess
import sys
from pathlib import Path

import pytest

from quebra.analyzers import check_attachment
from quebra.analyzers.assumptions import A1_RENEWAL_DURATIONS
from quebra.analyzers.kaplan_meier import KaplanMeierComparison

pytestmark = pytest.mark.unit


def test_the_assumption_id_literal_names_the_real_record() -> None:
    """Oracle: `analyzers/assumptions.py`. The literal keeps quebra out of the leaf."""
    assert check_attachment.A1_RENEWAL_DURATIONS_ID == A1_RENEWAL_DURATIONS.id


def test_an_empty_attachment_reads_not_assessed_and_the_comparison_says_the_same() -> (
    None
):
    """Oracle: AGENTS.md section 5, empty means NOT ASSESSED, said one way everywhere.

    Three checks asked and two verdicts back, so the count asked cannot be read off the
    verdicts by mistake.
    """
    empty = check_attachment.CheckAttachment()
    assert "NOT ASSESSED" in empty.summary()
    bare = KaplanMeierComparison(curves=[], threshold_label="t")
    assert bare.check_summary() == empty.summary()
    asked = check_attachment.CheckAttachment(
        checks_asked=("c1", "c2", "c3"),
        checks_unanswered=("c3",),
        check_verdicts=(("c1", "q1", "pass"), ("c2", "q1", "fail")),
    )
    assert "3 checks asked" in asked.summary()
    comparison = KaplanMeierComparison(
        curves=[],
        threshold_label="t",
        checks_asked=asked.checks_asked,
        checks_unanswered=asked.checks_unanswered,
        check_verdicts=asked.check_verdicts,
    )
    assert comparison.check_summary() == asked.summary()


# The leaf imports nothing from quebra. Allowed: the two standard-library modules it uses.
_ALLOWED = {"__future__", "dataclasses"}


def _imports(source: str) -> set[str]:
    """Every module the source imports: plain, relative (leading dots), or dynamic."""
    found: set[str] = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            found.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            found.add("." * node.level + (node.module or ""))
        elif isinstance(node, ast.Call):
            func = node.func
            name = (
                func.attr
                if isinstance(func, ast.Attribute)
                else getattr(func, "id", "")
            )
            if name in {"import_module", "__import__"}:
                found.add(f"<{name}>")
    return found


@pytest.mark.parametrize(
    "source",
    [
        "from quebra.analyzers.assumptions import A1_RENEWAL_DURATIONS",
        "from . import assumptions",
        "import importlib\nimportlib.import_module('quebra.analyzers.assumptions')",
        "__import__('quebra.analyzers.assumptions')",
        "def f():\n    import quebra.analyzers.assumptions",
    ],
    ids=["absolute", "relative", "importlib", "dunder-import", "inside-a-function"],
)
def test_the_import_scan_sees_every_form_of_import(source: str) -> None:
    """Oracle: each snippet imports a quebra module, so the scan must report one."""
    assert not _imports(source) <= _ALLOWED


def test_the_leaf_imports_nothing_from_quebra() -> None:
    """Oracle: the module docstring's contract, read off the source and off a fresh process.

    The scan covers imports inside functions; the process covers whatever runs at import.
    """
    assert _imports(Path(check_attachment.__file__).read_text()) <= _ALLOWED
    probe = (
        "import sys, quebra.analyzers.check_attachment; "
        "print(' '.join(sorted(m for m in sys.modules if m.split('.')[0] == 'quebra')))"
    )
    loaded = subprocess.run(
        [sys.executable, "-c", probe], capture_output=True, text=True, check=True
    ).stdout.split()
    assert loaded == ["quebra", "quebra.analyzers", "quebra.analyzers.check_attachment"]
