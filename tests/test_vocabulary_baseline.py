"""Ratchet on the word "gate" in src/ and tests/.

The project avoids "gate" because it reads a check as something that stops or licenses a
computation, and no check here does (AGENTS.md section 10). The physics term is exempt:
gate fidelity, gate infidelity, gate error, a quantum gate, a single- or two-qubit gate.

Scope: every `.py` and `.md` file under `src/` and `tests/`, this file excepted. `spec/` is
out of scope; it keeps its older wording as history.

A line is first split at case changes, so `ReuseGate` reads as `Reuse_Gate`. A match is
then "gate", "gates", "gated", "gating", "ungated" or a "gatekeep" word, in any case and
bounded by non-letters: `reuse_gate`, `ReuseGate`, `isGated` and `REUSE_GATE` count, while
"aggregate", "navigate" and "gateway" do not. Only "gate" and "gates" can be exempt: when
followed by fidelity, infidelity or error, singular or plural, or preceded by quantum,
single-qubit or two-qubit, joined by spaces, hyphens or underscores. That is the list in
section 10 and no wider, so "a per-qubit gate", "gating error" and "gated fidelity" count.

The count is zero and must stay zero. A ratchet that scans nothing also counts zero, so
the first test checks two things before it reads the count: the walk reached a known file
under each scope directory, and the scan finds at least one hit per `COUNTED` phrase in
this file, which the walk skips because it holds them on purpose. The second test pins the
matcher on phrases whose answer is known, so a pattern that stopped matching would fail
there rather than pass here.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.policy

REPO_ROOT = Path(__file__).resolve().parents[1]
SCOPE = ("src", "tests")
SUFFIXES = (".py", ".md")
THIS_FILE = Path(__file__).resolve()
# One file the walk must reach under each scope directory.
KNOWN = ("src/quebra/core/runner.py", "tests/test_reuse_rule.py")

# A case change inside an identifier is a word boundary: ReuseGate reads as Reuse_Gate.
CASE_CHANGE = re.compile(r"(?<=[a-z])(?=[A-Z])|(?<=[A-Z])(?=[A-Z][a-z])")
WORD = re.compile(
    r"(?<![A-Za-z])(?:un)?gat(?:e|es|ed|ing|ekeep[a-z]*)(?![A-Za-z])", re.IGNORECASE
)
PHYSICS_AFTER = re.compile(
    r"^[\s_-]+(?:fidelit(?:y|ies)|infidelit(?:y|ies)|errors?)", re.IGNORECASE
)
PHYSICS_BEFORE = re.compile(
    r"(?:quantum|(?:single|two)[\s_-]+qubit)[\s_-]+$", re.IGNORECASE
)


def _forbidden_in(text: str) -> int:
    text = CASE_CHANGE.sub("_", text)
    n = 0
    for match in WORD.finditer(text):
        after = text[match.end() : match.end() + 20]
        before = text[max(0, match.start() - 20) : match.start()]
        physics = PHYSICS_AFTER.search(after) or PHYSICS_BEFORE.search(before)
        if physics and match.group().lower() in ("gate", "gates"):
            continue
        n += 1
    return n


def _scanned_files() -> list[Path]:
    files: list[Path] = []
    for directory in SCOPE:
        for path in sorted((REPO_ROOT / directory).rglob("*")):
            if path.suffix not in SUFFIXES or "__pycache__" in path.parts:
                continue
            if path.resolve() == THIS_FILE:
                continue
            files.append(path)
    return files


def _occurrences(paths: list[Path]) -> list[str]:
    found: list[str] = []
    for path in paths:
        for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if _forbidden_in(line):
                found.append(
                    f"  {path.relative_to(REPO_ROOT)}:{lineno}: {line.strip()}"
                )
    return found


def test_no_gate_outside_the_physics_term() -> None:
    """Oracle: AGENTS.md section 10 sets the count at zero; the controls are the file tree.

    The walk must reach `KNOWN`, and the scan must find a hit for every `COUNTED` phrase
    in this file, lines past the first included, before a zero count means anything.
    """
    files = _scanned_files()
    seen = {path.relative_to(REPO_ROOT).as_posix() for path in files}
    missed = [known for known in KNOWN if known not in seen]
    assert not missed, f"the walk did not reach {missed}"
    assert len(_occurrences([THIS_FILE])) >= len(COUNTED), (
        "the scan missed hits in a file that holds one per COUNTED phrase"
    )
    found = _occurrences(files)
    assert not found, (
        'The word "gate" is back in src/ or tests/. Say "reuse rule" for the run-reuse '
        'condition, "test", "step" or "enforced" for the suite and CI, "condition" for a '
        'check row, and "decides", "stops" or "changes" for a statistic (AGENTS.md section '
        "10).\n" + "\n".join(found)
    )


# Hand-classified phrases with one forbidden match each.
COUNTED = [
    "this check gates the band",
    "the reuse gate admits it",
    "a gated row",
    "gating on the verdict",
    "def test_the_gate_fires():",
    "reuse_gate = True",
    "class ReuseGate:",
    "GateResult = 1",
    "isGated",
    "REUSE_GATE = 1",
    "ungated rows",
    "the gatekeeper check",
    "a per-qubit gate on the ledger",
    "gating error",
    "a gated fidelity",
]


def test_matcher_counts_the_forbidden_sense_and_spares_physics() -> None:
    """Oracle: hand-classified phrases, each with a known count."""
    exempt = [
        "Compute gate fidelity from typed FidelityInputs.",
        "gate_fidelity = 1 - infidelity",
        "a gate-error metric",
        "the two-qubit gate on the 2x2 device",
        "a single-qubit gate",
        "TwoQubitGate",
        "the gate infidelity",
        "a quantum gate",
        "GateFidelity",
        "aggregate",
        "AggregateResult",
        "navigate the gateway",
        "gate fidelities",
        "gate infidelities",
        "gate errors",
        "the two-qubit gates",
    ]
    assert [_forbidden_in(p) for p in COUNTED] == [1] * len(COUNTED)
    assert [_forbidden_in(p) for p in exempt] == [0] * len(exempt)
