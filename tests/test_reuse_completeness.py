"""The reuse rule, exercised against real git and real run directories.

Reuse requires ALL of: a matching content identity, a matching commit, a clean tree at both
ends, and a run that actually finished. Every other reuse test in the suite monkeypatches
`is_tree_clean` to `lambda: True`, so nothing here may do that when the point is the function
itself.

Two properties carry most of the weight:

- **Both git helpers must describe the CWD's repository**, not the one this module is
  installed into, because they answer the commit and clean-tree conditions of one rule, and
  two conditions read off different trees assert nothing. `"nogit"` therefore never
  satisfies the commit condition, not even against itself.

- **A candidate run must be complete.** The sink loop writes artifacts one at a time, so a run
  that raises partway through still leaves a readable provenance record whose identity, commit
  and tree-clean all match. Matching on directory name alone would read that as reusable and
  skip the job forever, leaving the unwritten sinks unproduced.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from quebra.core import runner
from quebra.core.paths import DataRootNotFound
from quebra.core.runner import _expected_sink_pkls, run_job
from quebra.cli import _module_from_path
from quebra.provenance import get_git_commit, is_tree_clean

pytestmark = pytest.mark.integration


def _git(repo: Path, *args: str) -> None:
    subprocess.run(
        ["git", *args],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    )


@pytest.fixture
def tiny_repo(tmp_path: Path) -> Path:
    """A real git repository with one committed file, so a clean tree is reachable."""
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q")
    _git(repo, "config", "user.email", "t@example.invalid")
    _git(repo, "config", "user.name", "t")
    (repo / "tracked.txt").write_text("one\n")
    _git(repo, "add", "tracked.txt")
    _git(repo, "commit", "-qm", "initial")
    return repo


# --- is_tree_clean / get_git_commit describe the CWD's repository ------------


def test_both_git_helpers_describe_the_cwd_repository(tiny_repo, monkeypatch):
    """Not the repository this module happens to be installed in.

    This is the whole point of the cwd anchor: under a wheel install `__file__` is
    site-packages, and under an editable install inside another checkout it is that other
    checkout - either way a `__file__`-anchored helper answers about a repository the run has
    nothing to do with.
    """
    monkeypatch.chdir(tiny_repo)
    commit = get_git_commit()
    assert commit != "nogit"
    assert is_tree_clean() is True

    head = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"],
        cwd=tiny_repo,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    assert commit == head, "the commit reported is not this repository's HEAD"


def test_a_dirty_tracked_file_makes_the_tree_dirty(tiny_repo, monkeypatch):
    monkeypatch.chdir(tiny_repo)
    assert is_tree_clean() is True
    (tiny_repo / "tracked.txt").write_text("two\n")
    assert is_tree_clean() is False


def test_an_untracked_file_does_not_make_the_tree_dirty(tiny_repo, monkeypatch):
    """`--untracked-files=no`: a scratch file must not disable artifact reuse.

    The claim the reuse rule makes is only "no uncommitted TRACKED changes", because an untracked
    file cannot have contributed to a result that was produced from committed code.
    """
    monkeypatch.chdir(tiny_repo)
    (tiny_repo / "scratch.log").write_text("noise\n")
    assert is_tree_clean() is True


def test_outside_a_repository_reports_nogit_and_dirty(tmp_path, monkeypatch):
    """ "Do not know" must read as dirty, so reuse never fires on an unvouched tree."""
    plain = tmp_path / "not_a_repo"
    plain.mkdir()
    monkeypatch.chdir(plain)
    assert get_git_commit() == "nogit"
    assert is_tree_clean() is False


# --- a partial run is not reusable ------------------------------------------


def _two_sink_job(d: Path, name: str) -> Path:
    (d / f"{name}.csv").write_text("a,b\n1,2\n")
    p = d / f"{name}.py"
    p.write_text(
        "from quebra.core.job import Job\n"
        "from quebra.core.dataset import Dataset\n"
        f'job = Job(name="{name}")\n'
        f'node = job.load_df(Dataset(path="{name}.csv", schema=None))\n'
        'first = job.step(lambda x: x, node, name="first")\n'
        'second = job.step(lambda x: x, node, name="second")\n'
        'job.materialize(first, name="one")\n'
        'job.materialize(second, name="two")\n'
    )
    return p


def test_expected_sink_pkls_names_every_sink(tmp_path):
    job = _module_from_path(_two_sink_job(tmp_path, "twosink")).job
    assert _expected_sink_pkls(job) == {"one.pkl", "two.pkl"}


def test_expected_sink_pkls_uses_the_source_node_for_a_figure(tmp_path):
    """The non-obvious branch: a figure's artifact is named for the node it reads, not for
    the figure's own title, and that holds whether or not a PDF is rendered."""
    (tmp_path / "f.csv").write_text("a,b\n1,2\n")
    p = tmp_path / "figjob.py"
    p.write_text(
        "from quebra.core.job import Job\n"
        "from quebra.core.dataset import Dataset\n"
        "from quebra.plots.km_survival_plot import KMSurvivalPlot\n"
        'job = Job(name="figjob")\n'
        'node = job.load_df(Dataset(path="f.csv", schema=None))\n'
        'shaped = job.step(lambda x: x, node, name="shaped")\n'
        'job.figure(KMSurvivalPlot, shaped, targets=["static"], title="a title")\n'
    )
    job = _module_from_path(p).job
    assert _expected_sink_pkls(job) == {"shaped.pkl"}


# --- "nogit" never satisfies the commit condition of the reuse rule ---------


def test_nogit_is_not_reuse_eligible_even_against_itself(tmp_path):
    """A repository with no commits answers `status` cleanly while `rev-parse HEAD` fails, so
    the pair ("nogit", clean) is reachable. Admitting it would reduce the reuse rule to an identity
    match, and identity does not cover a plot class - those reach a run through `job.figure`
    rather than through a step function, so they are absent from the import closure.
    """
    run_dir = tmp_path / "run"
    (run_dir / "provenance").mkdir(parents=True)
    (run_dir / "provenance" / "out.prov.json").write_text(
        '{"identity": "abc123", "git_commit": "nogit", "tree_clean": true}'
    )
    assert (
        runner._reuse_eligible_dir([run_dir], "abc123", "nogit", tree_clean=True)
        is None
    )


def test_a_real_commit_is_reuse_eligible(tmp_path):
    """Negative control for the line above."""
    run_dir = tmp_path / "run"
    (run_dir / "provenance").mkdir(parents=True)
    (run_dir / "provenance" / "out.prov.json").write_text(
        '{"identity": "abc123", "git_commit": "deadbee", "tree_clean": true}'
    )
    assert (
        runner._reuse_eligible_dir([run_dir], "abc123", "deadbee", tree_clean=True)
        == run_dir
    )


# --- a figure and a materialize may not share a name ------------------------


def test_a_figure_and_a_materialize_sharing_a_name_is_rejected(tmp_path):
    """Even on the SAME node. Both write `provenance/{name}.prov.json`, and the two records
    differ - only the figure's carries `targets_rendered` and a figure label - so the second
    declared silently overwrites the first, and a run that rendered PDFs reports none.
    """
    (tmp_path / "c.csv").write_text("a,b\n1,2\n")
    p = tmp_path / "clash.py"
    p.write_text(
        "from quebra.core.job import Job\n"
        "from quebra.core.dataset import Dataset\n"
        "from quebra.plots.km_survival_plot import KMSurvivalPlot\n"
        'job = Job(name="clash")\n'
        'node = job.load_df(Dataset(path="c.csv", schema=None))\n'
        'shaped = job.step(lambda x: x, node, name="shaped")\n'
        'job.figure(KMSurvivalPlot, shaped, targets=["static"], title="dup name")\n'
        'job.materialize(shaped, name="dup_name")\n'
    )
    job = _module_from_path(p).job
    with pytest.raises(ValueError, match="sink artifact name collision"):
        run_job(job, tmp_path / "output", force=True, data_root=tmp_path)


def test_two_materializes_of_one_node_under_one_name_are_allowed(tmp_path):
    """Negative control: same kind, same node, same name is one artifact asked for twice."""
    (tmp_path / "d.csv").write_text("a,b\n1,2\n")
    p = tmp_path / "dup.py"
    p.write_text(
        "from quebra.core.job import Job\n"
        "from quebra.core.dataset import Dataset\n"
        'job = Job(name="dup")\n'
        'node = job.load_df(Dataset(path="d.csv", schema=None))\n'
        'job.materialize(node, name="same")\n'
        'job.materialize(node, name="same")\n'
    )
    job = _module_from_path(p).job
    runner._check_sink_artifact_names(job)  # must not raise


# --- a named-but-missing data root fails the RUN, not just the resolver ------


def test_run_job_refuses_a_data_root_that_does_not_exist(tmp_path):
    """The demand check has to bite where the tool actually passes a root.

    `tests/test_data_root_resolution.py` pins `resolve_data_root` directly, which is not
    enough: the CLI resolves `--data-root` itself and hands the result to `run_job`, so a
    hardened resolver nothing calls with an explicit argument changes nothing about the
    shipped behaviour.

    And the failure mode is not "no data root" - it is a SILENT SUBSTITUTION. Every job in
    `jobs/` declares a repo-relative `Dataset.path`, and `resolve_dataset_path` tries the repo
    root as its second candidate, so with a bogus root the run finds those files anyway,
    completes, and records the repository tree's dataset hashes as if they were the requested
    tree's.
    """
    job_py = _two_sink_job(tmp_path, "rooted")
    with pytest.raises(DataRootNotFound, match="demand, not a candidate"):
        run_job(
            _module_from_path(job_py).job,
            tmp_path / "output",
            force=True,
            data_root=tmp_path / "not_a_real_root",
        )


def test_run_job_still_accepts_a_real_data_root(tmp_path, monkeypatch):
    """Negative control: the demand check must not reject a root that exists."""
    monkeypatch.setattr(runner, "is_tree_clean", lambda: True)
    job_py = _two_sink_job(tmp_path, "rooted")
    out = tmp_path / "output"
    run_job(_module_from_path(job_py).job, out, force=True, data_root=tmp_path)
    assert sorted(p.name for p in out.glob("rooted_*/*.pkl")) == ["one.pkl", "two.pkl"]


def test_a_run_missing_one_sink_is_not_reused(tmp_path, monkeypatch):
    """The positive control for the completeness filter.

    Deleting one sink's pkl reproduces exactly what a mid-loop failure leaves behind: the
    directory, the other sink's artifact, and a provenance record whose identity, commit and
    tree-clean all match. Before the filter this was indistinguishable from a finished run.
    """
    monkeypatch.setattr(runner, "is_tree_clean", lambda: True)
    job_py = _two_sink_job(tmp_path, "twosink")
    out = tmp_path / "output"

    run_job(_module_from_path(job_py).job, out, force=True, data_root=tmp_path)
    run_dirs = sorted(out.glob("twosink_*"))
    assert len(run_dirs) == 1
    assert (run_dirs[0] / "one.pkl").is_file()
    assert (run_dirs[0] / "two.pkl").is_file()

    # Amputate the second sink, leaving the first sink's record readable and matching.
    (run_dirs[0] / "two.pkl").unlink()

    run_job(_module_from_path(job_py).job, out, data_root=tmp_path)

    # Globbed across run dirs rather than checked inside `run_dirs[0]`. The run-dir name
    # carries a one-second timestamp, so whether the re-run lands back in the first directory
    # or creates a second one depends on which side of a second boundary the two calls fall.
    # Either outcome satisfies the claim being made, which is only that the sink got produced.
    assert list(out.glob("twosink_*/two.pkl")), (
        "the missing sink was not produced anywhere: the partial run was treated as "
        "reusable and the job was skipped"
    )


def test_a_complete_run_is_still_reused(tmp_path, monkeypatch, capsys):
    """The negative control. A filter that rejected every candidate would also pass the test
    above, so the reuse path has to be shown still working.

    Asserted on the runner's own "Skipping" line, for the same reason the test above avoids a
    directory count: within one second the two runs are indistinguishable on disk.
    """
    monkeypatch.setattr(runner, "is_tree_clean", lambda: True)
    job_py = _two_sink_job(tmp_path, "twosink")
    out = tmp_path / "output"

    run_job(_module_from_path(job_py).job, out, force=True, data_root=tmp_path)
    capsys.readouterr()
    run_job(_module_from_path(job_py).job, out, data_root=tmp_path)
    assert "Skipping twosink" in capsys.readouterr().out, (
        "a complete run was not reused; the completeness filter is too strict"
    )
