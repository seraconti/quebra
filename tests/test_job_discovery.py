"""Jobs declare their own identity and category, and discovery never imports them.

Two things live in the file rather than the filesystem: `JOB_ID`
(so `include` stops naming a path) and `JOB_FAMILY` (so recategorising stops meaning moving a
file).

The no-import rule is the load-bearing one. Importing a job file BUILDS its graph, so a
discovery pass that imported would execute every job in the tree to read two strings, and one
broken job would take the whole listing down.
"""

from __future__ import annotations

import pytest

from quebra.core.discovery import (
    DuplicateJobId,
    by_family,
    discover,
    read_declaration,
    resolve,
    swept,
)

pytestmark = pytest.mark.unit

REPO_JOBS = (
    pytest.importorskip("pathlib").Path(__file__).resolve().parent.parent / "jobs"
)


def _write(tmp_path, name, body):
    path = tmp_path / name
    path.write_text(body)
    return path


def test_discovery_does_not_import_the_job(tmp_path):
    """The rule, tested against a file whose import would blow up.

    If discovery imported, this would raise instead of returning a declaration.
    """
    _write(
        tmp_path,
        "explosive.py",
        'JOB_ID = "explosive"\nJOB_FAMILY = "demo"\n'
        'raise RuntimeError("importing this job must never happen during discovery")\n',
    )
    found = discover(tmp_path)
    assert found["explosive"].family == "demo"


def test_a_file_with_no_job_id_is_skipped(tmp_path):
    """`jobs/bench/` holds study modules that are not jobs; they must not appear."""
    _write(tmp_path, "helper.py", "X = 1\n")
    assert discover(tmp_path) == {}


def test_a_duplicate_job_id_raises_and_names_both_files(tmp_path):
    """Not last-one-wins: the ID is what `include` resolves, so a silent winner would make a
    composite depend on whichever file the walk reached first."""
    _write(tmp_path, "a.py", 'JOB_ID = "dup"\n')
    _write(tmp_path, "b.py", 'JOB_ID = "dup"\n')
    with pytest.raises(DuplicateJobId) as excinfo:
        discover(tmp_path)
    assert "a.py" in str(excinfo.value) and "b.py" in str(excinfo.value)


def test_a_computed_job_id_is_not_read(tmp_path):
    """Only a plain literal counts. A computed value would have to be executed to be read."""
    _write(tmp_path, "computed.py", 'JOB_ID = "pre" + "fix"\n')
    assert discover(tmp_path) == {}


def test_an_unparseable_file_is_not_silently_skipped(tmp_path):
    """INVERTED from its first form, which asserted the walk continued past a broken file.

    That was wrong: the old directory glob handed every file to the importer, so a syntax
    error surfaced as a failed job in the batch. Skipping it here would make a broken job
    VANISH from `--all` at exit 0 — quieter than the thing it replaced.
    """
    from quebra.core.discovery import UnreadableJobFile

    _write(tmp_path, "broken.py", "def (\n")
    _write(tmp_path, "good.py", 'JOB_ID = "good"\n')
    with pytest.raises(UnreadableJobFile):
        discover(tmp_path)


def test_resolve_names_what_was_available_when_it_misses(tmp_path):
    _write(tmp_path, "a.py", 'JOB_ID = "alpha"\n')
    with pytest.raises(ValueError, match="alpha"):
        resolve("beta", tmp_path)


def test_by_family_selects_and_sorts(tmp_path):
    _write(tmp_path, "b.py", 'JOB_ID = "b"\nJOB_FAMILY = "fam"\n')
    _write(tmp_path, "a.py", 'JOB_ID = "a"\nJOB_FAMILY = "fam"\n')
    _write(tmp_path, "c.py", 'JOB_ID = "c"\nJOB_FAMILY = "other"\n')
    assert [j.job_id for j in by_family(tmp_path, "fam")] == ["a", "b"]


def test_an_annotated_assignment_is_read(tmp_path):
    path = _write(tmp_path, "ann.py", 'JOB_ID: str = "annotated"\n')
    declared = read_declaration(path)
    assert declared is not None and declared.job_id == "annotated"


# ------------------------------------------------------------------ the real repository


def test_every_in_repo_job_declares_an_id_and_a_family():
    """A job with no `JOB_FAMILY` is reported, not silently skipped."""
    found = discover(REPO_JOBS)
    assert found, "no jobs discovered in this repository"
    missing = sorted(k for k, v in found.items() if not v.family)
    assert not missing, f"jobs declaring no JOB_FAMILY: {missing}"


def test_the_declared_id_matches_the_filename_today():
    """Not a rule the tool enforces — the whole point is that it need not — but true now, and
    a divergence should be deliberate rather than a typo."""
    for job_id, declared in discover(REPO_JOBS).items():
        assert job_id == declared.path.stem, (
            f"{declared.path} declares JOB_ID {job_id!r}; if that divergence is intended, "
            f"delete this assertion and say why"
        )


# ------------------------------------------------------------------ the sweep declaration


def test_sweep_defaults_to_true_and_can_be_declined(tmp_path):
    _write(tmp_path, "in.py", 'JOB_ID = "in"\n')
    _write(tmp_path, "out.py", 'JOB_ID = "out"\nJOB_SWEEP = False\n')
    found = discover(tmp_path)
    assert found["in"].sweep is True
    assert found["out"].sweep is False
    assert [j.job_id for j in swept(tmp_path)] == ["in"]


def test_a_bare_sweep_excludes_the_composites_as_the_docs_promise():
    """The regression this guards is measured, not hypothetical.

    Replacing the declaration with "every discovered job" silently widens `run --all`
    from the 9 swept jobs to all 18, pulling in the nine that decline the sweep: the
    composites, `jobs/active/km_with_checks_6d2s.py` and the six `check_ledger` jobs.
    `AGENTS.md` and `docs/WRITING_A_JOB.md` promise the composites are not swept. One is the
    independence survey, which someone typing `run --all` would start without asking for
    it, and each `check_ledger` job runs the R bridge over every record of a qubit.
    """
    all_jobs = discover(REPO_JOBS)
    swept_ids = {j.job_id for j in swept(REPO_JOBS)}
    for composite in (
        "compare_t2star_0704_vs_1004",
        "independence_survey",
        *(f"check_ledger_6d2s_q{q}" for q in range(1, 7)),
    ):
        assert composite in all_jobs, f"{composite} should still be discoverable"
        assert composite not in swept_ids, (
            f"{composite} is swept by a bare --all; AGENTS.md and docs/WRITING_A_JOB.md "
            f"both promise composites are not"
        )
    assert len(swept_ids) == 9, f"expected 9 swept jobs, got {sorted(swept_ids)}"


def test_an_unparseable_job_file_raises_rather_than_vanishing(tmp_path):
    """A broken job must fail the batch loudly rather than becoming invisible."""
    from quebra.core.discovery import UnreadableJobFile

    _write(tmp_path, "broken.py", "def (\n")
    with pytest.raises(UnreadableJobFile, match="could not be parsed"):
        discover(tmp_path)


def test_a_job_file_that_forgets_JOB_ID_is_detectable():
    """Discovery skips it — it is indistinguishable from a helper module — so the gap has to
    be closed by comparing against the FILES on disk, which is what this does. Without it a
    job that forgets the declaration is excluded from `--all` forever and silently."""

    declared = {d.path for d in discover(REPO_JOBS).values()}
    undeclared = []
    for path in sorted(REPO_JOBS.rglob("*.py")):
        if path.name == "__init__.py" or "bench" in path.parts:
            continue
        if path in declared:
            continue
        # A file that builds a job graph but declares no JOB_ID is the failure case.
        source = path.read_text(encoding="utf-8")
        if "Job(" in source and "JOB_ID" not in source:
            undeclared.append(str(path.relative_to(REPO_JOBS)))
    assert not undeclared, f"job files with no JOB_ID, invisible to --all: {undeclared}"
