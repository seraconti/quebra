"""Unbounded composite nesting.

A composite may include a composite to any depth. Cycles fail loudly at import
time; artifact reuse is depth-complete (the locator searches one shared pool
recursively) under the unchanged reuse rule (identity + commit + clean tree). These
build nested composites as tmp job files and run them through run_job.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from quebra.core import runner
from quebra.core.job import Job
from quebra.core.runner import run_job
from quebra.cli import _module_from_path

pytestmark = pytest.mark.integration


def _load(job_py: Path) -> Job:
    return _module_from_path(job_py).job


def _leaf(d: Path, name: str, csv: str = "a,b\n1,2\n") -> Path:
    (d / f"{name}.csv").write_text(csv)
    p = d / f"{name}.py"
    p.write_text(
        "from quebra.core.job import Job\n"
        "from quebra.core.dataset import Dataset\n"
        f'job = Job(name="{name}")\n'
        f'node = job.load_df(Dataset(path="{name}.csv", schema=None))\n'
        'job.materialize(node, name="out")\n'
    )
    return p


def _composite(d: Path, name: str, child: Path, alias: str = "c") -> Path:
    p = d / f"{name}.py"
    p.write_text(
        "from quebra.core.job import Job\n"
        f'job = Job(name="{name}")\n'
        f'inc = job.include(r"{child}", alias="{alias}")\n'
        f'm = job.step(lambda x: x, inc.ref("out"), name="{name}_step")\n'
        'job.materialize(m, name="out")\n'
    )
    return p


def _n(out: Path, pattern: str) -> int:
    return len(list(out.glob(pattern)))


# --- cycle detection (import/definition time) -------------------------------


def test_include_cycle_raises_at_import(tmp_path: Path) -> None:
    a = tmp_path / "a.py"
    b = tmp_path / "b.py"
    a.write_text(
        f'from quebra.core.job import Job\njob = Job("a")\njob.include(r"{b}")\n'
    )
    b.write_text(
        f'from quebra.core.job import Job\njob = Job("b")\njob.include(r"{a}")\n'
    )
    with pytest.raises(ValueError, match="include cycle:"):
        Job("top").include(a)


def test_self_include_raises_at_import(tmp_path: Path) -> None:
    s = tmp_path / "s.py"
    s.write_text(
        f'from quebra.core.job import Job\njob = Job("s")\njob.include(r"{s}")\n'
    )
    with pytest.raises(ValueError, match="include cycle:"):
        Job("top").include(s)


def test_import_stack_unwinds_on_mid_import_error(tmp_path: Path) -> None:
    # a sub-job that raises during import must not leak a stale stack entry -
    # a subsequent unrelated include of the same file must still work
    from quebra.core.job import _IMPORT_STACK

    boom = tmp_path / "boom.py"
    boom.write_text('raise RuntimeError("boom during import")\n')
    with pytest.raises(RuntimeError, match="boom"):
        Job("top").include(boom)
    assert _IMPORT_STACK == []  # popped despite the exception

    boom.write_text('from quebra.core.job import Job\njob = Job("fixed")\n')
    assert Job("top2").include(boom).job.name == "fixed"


# --- depth-2 / depth-3 end to end -------------------------------------------


def test_depth2_fresh_end_to_end(tmp_path: Path) -> None:
    leaf = _leaf(tmp_path, "leaf")
    c1 = _composite(tmp_path, "c1", leaf, alias="l")
    c2 = _composite(tmp_path, "c2", c1, alias="c1")
    out = tmp_path / "output"

    run_job(_load(c2), out, force=True, data_root=tmp_path)

    # leaf artifact materialized two levels down
    assert _n(out, "c2_*/subjobs_output/c1_*/subjobs_output/leaf_*/out.pkl") == 1
    assert _n(out, "c2_*/out.pkl") == 1  # top composite produced its own


def test_depth3_fresh_then_reuse(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(runner, "is_tree_clean", lambda: True)
    leaf = _leaf(tmp_path, "leaf")
    c1 = _composite(tmp_path, "c1", leaf, alias="l")
    c2 = _composite(tmp_path, "c2", c1, alias="c1")
    c3 = _composite(tmp_path, "c3", c2, alias="c2")
    out = tmp_path / "output"

    run_job(_load(c3), out, force=True, data_root=tmp_path)  # fresh, depth 3
    assert _n(out, "**/leaf_*/") == 1

    # second run reuses the whole c2 subtree (matching identity + commit, clean tree)
    run_job(_load(c3), out, force=True, data_root=tmp_path, reuse_deps=True)
    assert _n(out, "**/leaf_*/") == 1  # leaf NOT re-run
    assert _n(out, "**/c2_*/") == 1  # c2 subtree reused, not rebuilt


# --- diamond dedup within one invocation (the accepted 3b semantics) --------


def test_same_invocation_diamond_dedup(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # g includes c1 and c2; both include the SAME leaf. Under --reuse-deps on a
    # clean tree, c2's locator reuses the leaf c1's branch just produced (the
    # sibling artifact is under the grandparent g, not under c2's own dir).
    monkeypatch.setattr(runner, "is_tree_clean", lambda: True)
    leaf = _leaf(tmp_path, "leaf")
    c1 = _composite(tmp_path, "c1", leaf, alias="l")
    c2 = _composite(tmp_path, "c2", leaf, alias="l")
    g = tmp_path / "g.py"
    g.write_text(
        "from quebra.core.job import Job\n"
        'job = Job(name="g")\n'
        f'a = job.include(r"{c1}", alias="c1")\n'
        f'b = job.include(r"{c2}", alias="c2")\n'
        'm = job.step(lambda x, y: x, a.ref("out"), b.ref("out"), name="merge")\n'
        'job.materialize(m, name="out")\n'
    )
    out = tmp_path / "output"

    run_job(_load(g), out, force=True, data_root=tmp_path, reuse_deps=True)
    assert _n(out, "**/leaf_*/") == 1  # located once across the two branches


# --- identity iff, one level deeper -----------------------------------------


def test_depth2_identity_folds_leaf_through_both_levels(tmp_path: Path) -> None:
    # The iff property one level deeper. Rather than mutate a dataset (content_hash
    # is process-cached by path) or compare across roots (the composite source
    # embeds the child path, so identical trees in different dirs legitimately
    # differ), verify the recursive fold STRUCTURALLY: the top's children hold c1's
    # digest, c1's children hold the leaf's digest. Combined with fold injectivity
    # (test_identity), that IS "top identity changes iff the leaf's does".
    r = tmp_path / "r"
    r.mkdir()
    leaf = _leaf(r, "leaf")
    c1 = _composite(r, "c1", leaf, alias="l")
    c2 = _composite(r, "c2", c1, alias="c1")
    top = _load(c2)

    # deterministic: the same tree re-imported yields the same top identity
    assert top.build_identity(r).digest == _load(c2).build_identity(r).digest

    c1_job = top.includes[0].job
    leaf_job = c1_job.includes[0].job
    top_id = top.build_identity(r)
    c1_id = c1_job.build_identity(r)
    leaf_id = leaf_job.build_identity(r)
    assert top_id.children == (c1_id.digest,)  # top folds c1's identity
    assert c1_id.children == (leaf_id.digest,)  # c1 folds the leaf's identity


# --- reuse across depths, both directions -----------------------------------


def test_top_standalone_leaf_reused_by_nested_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(runner, "is_tree_clean", lambda: True)
    leaf = _leaf(tmp_path, "leaf")
    c1 = _composite(tmp_path, "c1", leaf, alias="l")
    c2 = _composite(tmp_path, "c2", c1, alias="c1")
    out = tmp_path / "output"

    run_job(_load(leaf), out, force=True, data_root=tmp_path)  # top-level standalone
    assert _n(out, "leaf_*/") == 1

    run_job(_load(c2), out, force=True, data_root=tmp_path, reuse_deps=True)
    # the depth-2 leaf locator reused the top-level standalone artifact: no nested leaf
    assert _n(out, "**/leaf_*/") == 1
    assert _n(out, "*/subjobs_output/**/leaf_*/") == 0


def test_nested_leaf_reused_by_separate_top_composite(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(runner, "is_tree_clean", lambda: True)
    leaf = _leaf(tmp_path, "leaf")
    c1 = _composite(tmp_path, "c1", leaf, alias="l")
    g1 = _composite(tmp_path, "g1", c1, alias="c1")
    g2 = _composite(tmp_path, "g2", c1, alias="c1")  # separate top, same child c1
    out = tmp_path / "output"

    run_job(_load(g1), out, force=True, data_root=tmp_path)  # produces c1 deep in g1
    assert _n(out, "**/c1_*/") == 1

    run_job(_load(g2), out, force=True, data_root=tmp_path, reuse_deps=True)
    # g2 reuses c1's artifact buried in g1's tree (pool-wide recursive search)
    assert _n(out, "**/c1_*/") == 1
    assert _n(out, "**/leaf_*/") == 1


# --- rider: step() ownership ------------------------------------------------


def test_step_rejects_foreign_artifactref(tmp_path: Path) -> None:
    leaf = _leaf(tmp_path, "leaf")
    x = Job("x")
    x.job_file = (tmp_path / "x.py").resolve()
    ref = x.include(leaf, alias="l").ref("out")  # ref owned by composite X
    y = Job("y")
    with pytest.raises(ValueError, match="different composite"):
        y.step(lambda a: a, ref, name="s")  # used in composite Y


# --- provenance click-through paths at depth 2 ------------------------------


def test_prov_subjob_dirs_pool_relative_at_depth2(tmp_path: Path) -> None:
    import json

    leaf = _leaf(tmp_path, "leaf")
    c1 = _composite(tmp_path, "c1", leaf, alias="l")
    c2 = _composite(tmp_path, "c2", c1, alias="c1")
    out = tmp_path / "output"

    run_job(_load(c2), out, force=True, data_root=tmp_path)

    # the nested c1 run's includes entry points at the leaf's provenance dir,
    # recorded relative to the top pool root and actually present there
    c1_prov = next(out.glob("c2_*/subjobs_output/c1_*/provenance/out.prov.json"))
    rec = json.loads(c1_prov.read_text())
    sub = rec["includes"][0]["subjob_prov_dir"]
    assert not Path(sub).is_absolute()
    assert sub.startswith("c2_")  # pool-root-relative, spanning the full depth
    assert (out / sub).is_dir()
