"""Reuse gate + artifact completeness.

The gate: an artifact is reusable only when its content identity AND git commit
match the current run AND the working tree is clean; any mismatch re-runs fresh.
These unit-test the pure gate helper (so the clean-tree/commit logic is pinned
without needing a real clean repo), plus the panel-data completeness validators.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.core.runner import _reuse_eligible_dir, run_job

pytestmark = pytest.mark.unit


def _passthrough(x: object) -> object:
    return x


def test_sink_artifact_name_collision_raises(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # two sinks writing the same {basename}.pkl from DIFFERENT nodes would clobber
    # each other and let a composite reuse a wrong-but-plausible artifact - reject
    # it at run start (before any output dir exists)
    (tmp_path / "data.csv").write_text("a,b\n1,2\n")
    job_py = tmp_path / "collide.py"
    job_py.write_text("# job\n")
    job = Job(name="collide")
    a = job.load_df(Dataset(path="data.csv", schema=None))
    b = job.step(_passthrough, a, name="b")
    job.materialize(a, "dup")
    job.materialize(b, "dup")  # same basename, different source node
    job.job_file = job_py.resolve()

    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match="collision"):
        run_job(job, tmp_path / "out", force=True, data_root=tmp_path)
    assert not list((tmp_path / "out").glob("collide_*"))  # failed before mkdir


def _seed_run(
    root: Path, name: str, identity: str, commit: str, tree_clean: bool = True
) -> Path:
    run_dir = root / name
    (run_dir / "provenance").mkdir(parents=True)
    (run_dir / "provenance" / "out.prov.json").write_text(
        json.dumps(
            {"identity": identity, "git_commit": commit, "tree_clean": tree_clean}
        )
    )
    return run_dir


def test_reuse_requires_identity_commit_and_clean_tree(tmp_path: Path) -> None:
    d = _seed_run(tmp_path, "job_aaaaaa_20240101_000000", "ID", "COMMIT")
    cands = [d]

    # all conditions met => reusable
    assert _reuse_eligible_dir(cands, "ID", "COMMIT", tree_clean=True) == d
    # consumer's dirty tree => never reusable, regardless of match
    assert _reuse_eligible_dir(cands, "ID", "COMMIT", tree_clean=False) is None
    # identity mismatch => not reusable
    assert _reuse_eligible_dir(cands, "OTHER", "COMMIT", tree_clean=True) is None
    # commit mismatch => not reusable (commit-safe, not content-safe)
    assert _reuse_eligible_dir(cands, "ID", "OTHER", tree_clean=True) is None


def test_reuse_rejects_artifact_produced_on_dirty_tree(tmp_path: Path) -> None:
    # an artifact built on a dirty tree records tree_clean=False and is never
    # reusable even at a matching commit on a now-clean tree (it does not
    # correspond to any commit's code)
    d = _seed_run(
        tmp_path, "job_aaaaaa_20240101_000000", "ID", "COMMIT", tree_clean=False
    )
    assert _reuse_eligible_dir([d], "ID", "COMMIT", tree_clean=True) is None


def test_reuse_picks_newest_matching(tmp_path: Path) -> None:
    old = _seed_run(tmp_path, "job_aaaaaa_20240101_000000", "ID", "C")
    new = _seed_run(tmp_path, "job_aaaaaa_20240102_000000", "ID", "C")
    # candidates sorted oldest->newest; the newest eligible one wins
    assert _reuse_eligible_dir([old, new], "ID", "C", tree_clean=True) == new


def test_reuse_skips_provless_dir(tmp_path: Path) -> None:
    bare = tmp_path / "job_aaaaaa_20240101_000000"
    bare.mkdir()
    assert _reuse_eligible_dir([bare], "ID", "C", tree_clean=True) is None


# --- completeness validators ------------------------------------------------


def test_within_calibration_incomplete_construction_raises() -> None:
    from quebra.panels.within_calibration import WithinCalibrationPanelData

    from quebra.analyzers.distinguish_band import DistinguishBand
    from quebra.analyzers.reliability_band import ReliabilityBand
    from quebra.analyzers.signal_band import SignalBand

    # A threshold with no per-threshold entries in any band is incomplete. The outer
    # class owns the ladder and delegates the sweep, so the error comes from the band.
    with pytest.raises(ValueError, match="incomplete DistinguishBand"):
        WithinCalibrationPanelData(
            signal=SignalBand(
                t_h=np.array([0.0, 1.0]), values=np.array([1.0, 2.0]), gap_spans_h=[]
            ),
            distinguish=DistinguishBand(),
            reliability=ReliabilityBand(),
            meta={},
            thresholds=[("thr", 1.5, False)],
            primary_label="x",
        )


def test_across_calibration_incomplete_construction_raises() -> None:
    from quebra.panels.across_calibration import AcrossCalibrationPanelData

    with pytest.raises(ValueError, match="incomplete AcrossCalibrationPanelData"):
        # non-empty raw inputs but default (empty) derived arrays => length mismatch
        AcrossCalibrationPanelData(
            intervals_s=np.array([1.0, 2.0]),
            event_times_unix_s=np.array([10.0, 20.0]),
            stats={},
            meta={},
        )
