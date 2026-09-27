"""The tier axis is total: every collected test answers exactly one kind of question.

Oracle: the marker list in `pyproject.toml`. This is a policy test, not a statistical one.

`spec/quebraplan.md` 5.1 proposed six tier DIRECTORIES. `spec/spectests06.md` declined them
and put the tier on a marker instead, because directories cross-cut the oracle-and-subject
index in `AGENTS.md` section 7 and several files here legitimately hold more than one tier.
A marker scheme decays silently unless something makes it total, which is what this file is.

The check is decidable, which is why it exists at all: "does this item carry one of these
seven names" is a property of the item. The neighbouring oracle rule ("does this docstring
name a source of truth") is a judgement, and `spec/spectests06.md` R6.1 declines to guard it
rather than enforce a regex that would check something else.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]

# The pre-deselection probe runs a CHILD pytest over `tests`, which re-collects this file.
# Without a sentinel the child re-enters the probe and spawns a grandchild, unbounded.
_CHILD_PROBE_ENV = "QUEBRA_PRE_DESELECTION_CHILD"

# The tier axis. Exactly one per test.
TIER_MARKERS = frozenset(
    {
        "unit",
        "properties",
        "statistical",
        "integration",
        "validation",
        "regression",
        "policy",
    }
)

# The cost axis. Orthogonal: zero or more per test, and never a substitute for a tier.
COST_MARKERS = frozenset({"slow", "heavy", "real", "r"})


def tier_markers_of(item) -> set[str]:
    """The tier markers on one collected item. The detector both tests below share."""
    return {mark.name for mark in item.iter_markers()} & TIER_MARKERS


def _violations(items) -> list[tuple[str, list[str]]]:
    """(nodeid, tiers) for every item not carrying exactly one tier marker."""
    out = []
    for item in items:
        tiers = sorted(tier_markers_of(item))
        if len(tiers) != 1:
            out.append((item.nodeid, tiers))
    return out


class _FakeMark:
    def __init__(self, name: str) -> None:
        self.name = name


class _FakeItem:
    """A planted item, for the positive control. Mirrors the `iter_markers` surface."""

    def __init__(self, nodeid: str, names: tuple[str, ...]) -> None:
        self.nodeid = nodeid
        self._names = names

    def iter_markers(self):
        return [_FakeMark(name) for name in self._names]


@pytest.mark.policy
def test_every_collected_test_carries_exactly_one_tier_marker(
    collected_pre_deselection,
) -> None:
    """Reads the real session, so it cannot pass by checking a copy of itself.

    `conftest.collected_pre_deselection` and not `request.session.items`, because the
    latter is what survived the selector. `scripts/fast-selector.txt` excludes `r`, and it
    is the only selector CI runs, so a guard reading the post-deselection list would go
    green on an `r`-marked test carrying no tier and go red only on a full local run. The
    collection record removes that dependence: whenever this test is itself selected, it
    sees every item the run collected, not only the ones its own `-m` expression kept.

    A path restriction (`pytest tests/test_foo.py`) still narrows what is collected at all,
    which no hook can widen. The full run is what makes this total.
    """
    assert collected_pre_deselection, (
        "the pre-deselection record is empty, so this guard would pass vacuously; "
        "conftest.pytest_itemcollected did not run"
    )
    bad = _violations(collected_pre_deselection)
    assert not bad, (
        f"{len(bad)} collected test(s) do not carry exactly one tier marker "
        f"({'/'.join(sorted(TIER_MARKERS))}):\n"
        + "\n".join(f"  {nodeid}  ->  {tiers or 'none'}" for nodeid, tiers in bad[:20])
    )


@pytest.mark.policy
def test_the_record_is_pre_deselection_not_post(tmp_path) -> None:
    """Oracle: pytest's own deselection, observed in a child run under a narrowing selector.

    The guard above is only stronger than the old `request.session.items` if the record it
    reads is built BEFORE deselection. Nothing else pins that: a replacement returning the
    post-deselection list would still be non-empty, would still make the guard pass, and only
    a hand-run mutation would find it.

    The child exercises the FIXTURE, not the stash behind it. Probing the stash directly
    leaves the swap this test exists to catch undetected, because the guard consumes
    `collected_pre_deselection` and a replacement fixture body is exactly what would be
    edited. The child writes its two numbers to a file rather than to stdout, so a parse
    failure cannot be mistaken for a pass.
    """
    if os.environ.get(_CHILD_PROBE_ENV):
        pytest.skip(
            "child of the pre-deselection probe; spawning another would recurse"
        )

    # The probe lives under `tests/`, not in tmp_path, for the same reason the race test
    # below plants its victim there: a hook and a fixture declared in the root conftest are
    # dispatched through `node.ihook`, which skips conftests that do not apply to the item's
    # path. A probe outside the rootdir cannot see `collected_pre_deselection` at all. That
    # resolves differently under an editable install than under the clean non-editable
    # resolve `make check-ci` builds, so a tmp_path probe passes here and fails there.
    probe = REPO_ROOT / "tests" / "test_zz_record_probe_tmp.py"
    result = tmp_path / "counts.txt"
    try:
        probe.write_text(
            "import pytest\n"
            "@pytest.mark.policy\n"
            "def test_probe(collected_pre_deselection, request):\n"
            f"    open({str(result)!r}, 'w').write(\n"
            "        f'{len(collected_pre_deselection)} {len(request.session.items)}'\n"
            "    )\n",
            encoding="utf-8",
        )
        completed = subprocess.run(
            [sys.executable, "-m", "pytest", "-q", "-m", "policy", "tests"],
            cwd=REPO_ROOT,
            env={**os.environ, _CHILD_PROBE_ENV: "1"},
            capture_output=True,
            text=True,
        )
    finally:
        probe.unlink(missing_ok=True)

    # The child's exit code is deliberately NOT asserted: it couples this test to every other
    # policy test in the child run, and a failure elsewhere would surface here as someone
    # else's traceback. The result file is the real condition.
    assert result.exists(), (
        f"the probe test did not run.\nstdout:\n{completed.stdout}\n"
        f"stderr:\n{completed.stderr}"
    )
    record, selected = (int(n) for n in result.read_text().split())
    assert record > selected, (
        f"the fixture handed the guard {record} items while the `-m policy` selector kept "
        f"{selected}, so it is not a pre-deselection record"
    )


@pytest.mark.policy
def test_the_record_survives_a_competing_tryfirst_deselector(tmp_path) -> None:
    """Oracle: pytest's hook ordering, observed against a deliberately hostile plugin.

    This pins the property that CHOSE `pytest_itemcollected` over
    `pytest_collection_modifyitems(tryfirst=True)`. The modifyitems form loses a race: a
    competing `tryfirst` implementation can remove an item before the recorder sees it.
    `pytest_itemcollected` fires per item inside collection, before any modifyitems runs at
    all, so it cannot lose.

    Without this test the rejected form can be restored with the suite green, because under
    an ordinary selector the two behave identically. The competitor is what separates them.
    """
    if os.environ.get(_CHILD_PROBE_ENV):
        pytest.skip(
            "child of the pre-deselection probe; spawning another would recurse"
        )

    # The competitor is a conftest in a SUBDIRECTORY of tests/, not a `-p` plugin, and that
    # is load-bearing twice over. A `-p` plugin registers BEFORE the root conftest, so under
    # pytest's LIFO ordering among `tryfirst` implementations the root conftest would win the
    # race anyway and the test would pass whichever hook form ships. A subdirectory conftest
    # registers LATER, so its `tryfirst` runs FIRST and genuinely beats a modifyitems
    # recorder. The victim must also live under tests/: a hook declared in the root conftest
    # is dispatched through `node.ihook`, which skips conftests that do not apply to the
    # item's path, so an item collected from outside the rootdir never enters the record.
    race = REPO_ROOT / "tests" / "_race_tmp"
    probe = REPO_ROOT / "tests" / "test_zz_probe_tmp.py"
    result = tmp_path / "seen.txt"
    try:
        race.mkdir(exist_ok=True)
        (race / "conftest.py").write_text(
            "import pytest\n"
            "@pytest.hookimpl(tryfirst=True)\n"
            "def pytest_collection_modifyitems(session, config, items):\n"
            "    items[:] = [i for i in items if 'test_zz_victim_tmp' not in i.nodeid]\n",
            encoding="utf-8",
        )
        (race / "test_zz_victim_tmp.py").write_text(
            "import pytest\n@pytest.mark.policy\ndef test_victim():\n    pass\n",
            encoding="utf-8",
        )
        probe.write_text(
            "import pytest\n"
            "@pytest.mark.policy\n"
            "def test_probe(collected_pre_deselection, request):\n"
            "    recorded = any(\n"
            "        'test_zz_victim_tmp' in i.nodeid for i in collected_pre_deselection\n"
            "    )\n"
            "    selected = any(\n"
            "        'test_zz_victim_tmp' in i.nodeid for i in request.session.items\n"
            "    )\n"
            f"    open({str(result)!r}, 'w').write(f'{{recorded}} {{selected}}')\n",
            encoding="utf-8",
        )
        completed = subprocess.run(
            [sys.executable, "-m", "pytest", "-q", "-m", "policy", "tests"],
            cwd=REPO_ROOT,
            env={**os.environ, _CHILD_PROBE_ENV: "1"},
            capture_output=True,
            text=True,
        )
    finally:
        probe.unlink(missing_ok=True)
        # rmtree, not unlink: the child leaves a __pycache__ directory behind and a
        # cleanup that raises inside `finally` would mask the real assertion.
        shutil.rmtree(race, ignore_errors=True)

    assert result.exists(), (
        f"the probe test did not run.\nstdout:\n{completed.stdout}\n"
        f"stderr:\n{completed.stderr}"
    )
    recorded, selected = result.read_text().split()
    assert selected == "False", (
        "the competing deselector did not actually remove the victim, so this test would "
        "pass without exercising the race it exists to guard"
    )
    assert recorded == "True", (
        "a competing tryfirst deselector removed the victim before the record saw it: the "
        "collection record is not immune to hook ordering"
    )


@pytest.mark.policy
def test_the_detector_catches_a_missing_and_a_doubled_tier() -> None:
    """Positive control. Without it the guard above passes whenever collection is empty."""
    planted = [
        _FakeItem("planted::no_tier", ()),
        _FakeItem("planted::cost_only", ("slow", "real")),
        _FakeItem("planted::two_tiers", ("unit", "statistical")),
        _FakeItem("planted::ok", ("unit",)),
        _FakeItem("planted::ok_with_cost", ("statistical", "slow")),
    ]
    bad = dict(_violations(planted))
    assert set(bad) == {
        "planted::no_tier",
        "planted::cost_only",
        "planted::two_tiers",
    }
    assert bad["planted::two_tiers"] == ["statistical", "unit"]


@pytest.mark.policy
def test_the_two_axes_are_exactly_what_pyproject_declares(request) -> None:
    """Oracle: `[tool.pytest.ini_options] markers`, read rather than re-typed.

    The sets above are literals, so without this they are a second copy of the
    declaration and could drift from it silently - deleting a marker from `pyproject.toml`
    would leave this file green while the declared and asserted axes disagreed. That is the
    defect class this whole checkpoint exists to remove, so the guard must not carry it.
    """
    import tomllib

    # `config.getini("markers")` also returns builtin and plugin markers (parametrize,
    # skipif, hypothesis, ...), so it cannot say what THIS project declares. Read the file.
    # `rootpath` is derived from the ini file's own location, not the cwd, so this keeps
    # working from an unrelated working directory.
    pyproject = request.config.rootpath / "pyproject.toml"
    entries = tomllib.loads(pyproject.read_text())["tool"]["pytest"]["ini_options"][
        "markers"
    ]
    declared = {entry.split(":", 1)[0].strip() for entry in entries}
    assert TIER_MARKERS | COST_MARKERS == declared, (
        "the two axes here and the markers declared in pyproject.toml have drifted:\n"
        f"  declared not classified: {sorted(declared - (TIER_MARKERS | COST_MARKERS))}\n"
        f"  classified not declared: {sorted((TIER_MARKERS | COST_MARKERS) - declared)}"
    )
    assert not (TIER_MARKERS & COST_MARKERS), (
        "a name on both axes would make 'exactly one tier' ambiguous"
    )
