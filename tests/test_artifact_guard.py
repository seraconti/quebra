"""Staleness guard for materialized panel-data artifacts.

A pre-split pickle restores __dict__ without the derived fields; the guard must
fail loudly at the pickle boundary instead of crashing mid-render or silently
drawing class defaults. Synthetic tests drive __setstate__ exactly as pickle
protocol 2 does (``cls.__new__(cls)`` then ``__setstate__(state)``); the
real-artifact test loads the pickles on disk through `load_artifact`, end to end.
"""

from __future__ import annotations

import dataclasses
import json
import pickle
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from quebra.analyzers import windows
from quebra.analyzers.within_calibration_compute import (
    build_within_calibration_panel_data,
)
from quebra.panels._across_calibration_compute import (
    build_across_calibration_panel_data,
)
from quebra.core._artifact_guard import StaleArtifactGuard
from quebra.panels.within_calibration import WithinCalibrationPanelData
from quebra.panels.across_calibration import AcrossCalibrationPanelData

pytestmark = pytest.mark.integration


REPO_ROOT = Path(__file__).resolve().parents[1]


def _carved(t_h, series, thresholds):
    """Carve through the real analyzer, as a job does.

    The builder does not carve: it consumes the window and read tables, so a test
    that constructs panel data has to produce them the same way production does.
    """
    result = windows.run(
        windows.WindowsInputs(
            t_rel_s=np.asarray(t_h, dtype=float) * 3600.0,
            values=np.asarray(series, dtype=float),
            thresholds=thresholds,
            dataset_id="unit",
        )
    )
    return result.windows, result.reads


_GUARD_THRESHOLDS = [("1e-4", 1e-4, False)]


def _valid_within_calibration() -> WithinCalibrationPanelData:
    rng = np.random.default_rng(7)
    t_h = np.linspace(0.0, 12.0, 200)
    series = 1e-5 + 1e-5 * rng.random(200)
    return build_within_calibration_panel_data(
        t_h=t_h,
        primary_series=series,
        primary_label="Infidelity",
        thresholds=_GUARD_THRESHOLDS,
        meta={"qubit": "1"},
        windows=_carved(t_h, series, _GUARD_THRESHOLDS)[0],
        reads=_carved(t_h, series, _GUARD_THRESHOLDS)[1],
        gap_spans_s=[],  # uniform spacing: the carve records no gap
    )


def _valid_across_calibration() -> AcrossCalibrationPanelData:
    rng = np.random.default_rng(8)
    ev = 1.6e9 + np.cumsum(rng.gamma(2.0, 3600.0, 60))
    iv = np.diff(ev)
    ev = ev[1:]
    return build_across_calibration_panel_data(
        intervals_s=iv,
        event_times_unix_s=ev,
        stats={
            "count": int(len(iv)),
            "mean_s": float(np.mean(iv)),
            "std_s": float(np.std(iv)),
            "min_s": float(np.min(iv)),
            "max_s": float(np.max(iv)),
        },
        meta={"qubit": "2", "device": "6D2S", "dataset_id": "unit"},
    )


@pytest.mark.parametrize(
    ("cls", "make", "derived_field"),
    [
        (WithinCalibrationPanelData, _valid_within_calibration, "distinguish"),
        (AcrossCalibrationPanelData, _valid_across_calibration, "elapsed_days"),
    ],
)
def test_stale_state_raises(cls, make, derived_field) -> None:
    state = dict(make().__dict__)
    del state[derived_field]
    obj = cls.__new__(cls)
    with pytest.raises(ValueError, match="stale .* artifact"):
        obj.__setstate__(state)


def test_missing_cv_alone_raises() -> None:
    # As a plain class default, cv lets a stale instance silently draw CV=nan.
    # With default_factory + the guard, its absence must raise like any field. It now
    # lives on the signal band, so this exercises the NESTED guard.
    band = _valid_within_calibration().signal
    state = dict(band.__dict__)
    del state["cv"]
    obj = type(band).__new__(type(band))
    with pytest.raises(ValueError, match="cv"):
        obj.__setstate__(state)


def test_nested_band_guard_fires_on_a_stale_band() -> None:
    """The whole point of making every band inherit the guard.

    StaleArtifactGuard derives its key set from dataclasses.fields(cls), so the OUTER
    class only ever validates {signal, distinguish, reliability, meta, ...}. A band
    missing half its fields would load clean unless the band guards itself.
    """
    obj = _valid_within_calibration()
    blob = pickle.dumps(obj)
    # rename `occupancy` inside the reliability band's pickled state
    tampered = blob.replace(b"\x8c\toccupancy", b"\x8c\toccupanZy")
    assert tampered != blob, "fixture did not contain the expected pickled key"
    with pytest.raises(ValueError, match="stale ReliabilityBand artifact"):
        pickle.loads(tampered)


def test_outer_guard_alone_would_not_catch_it() -> None:
    """Documents WHY the per-band guards exist: the outer key set is only four names."""
    import dataclasses

    outer_fields = {f.name for f in dataclasses.fields(WithinCalibrationPanelData)}
    assert "occupancy" not in outer_fields
    assert {"signal", "distinguish", "reliability"} <= outer_fields


@pytest.mark.parametrize("make", [_valid_within_calibration, _valid_across_calibration])
def test_builder_pickle_round_trip(make) -> None:
    obj = make()
    loaded = pickle.loads(pickle.dumps(obj))
    assert type(loaded) is type(obj)
    assert set(loaded.__dict__) == set(obj.__dict__)


def test_non_dict_state_raises() -> None:
    # a slots=True dataclass would pickle state as (dict, slots_dict); the guard
    # must answer with its own error, not AttributeError on .keys()
    obj = WithinCalibrationPanelData.__new__(WithinCalibrationPanelData)
    with pytest.raises(ValueError, match="unexpected pickle state"):
        obj.__setstate__(({}, {"x": 1}))


class _StalePickle:
    """Pickles into <cls> + a partial state dict - byte-wise what a pre-split
    artifact looks like on load (object.__new__ then __setstate__)."""

    def __init__(self, cls: type, state: dict) -> None:
        self.cls, self.state = cls, state

    def __reduce__(self):
        return (object.__new__, (self.cls,), self.state)


def _identity(x: object) -> object:
    return x


def _import_job_module(job_py: Path):
    # delegate to the real CLI loader so this test tracks its behavior
    from quebra.cli import _module_from_path

    return _module_from_path(job_py).job


def test_composite_reuse_of_reuse_eligible_stale_artifact_aborts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Defense-in-depth: even a cached sub-job artifact that PASSES the reuse gate
    (matching identity + commit, clean tree) must abort at load if its pickle is
    schema-stale - the guard, not a mid-render crash. Requires seeding the cache
    dir by the sub-job's IDENTITY and forcing a clean tree so the gate admits it."""
    from quebra.core.dataset import Dataset  # noqa: F401  (sub-job module imports it)
    from quebra.core.job import Job
    from quebra.core import runner
    from quebra.core.runner import run_job

    sub_py = tmp_path / "tiny_sub.py"
    sub_py.write_text(
        "from quebra.core.job import Job\n"
        "from quebra.core.dataset import Dataset\n"
        'job = Job(name="tiny_sub")\n'
        'node = job.load_df(Dataset(path="data.csv", schema=None))\n'
        'job.materialize(node, name="panel_data")\n'
    )
    (tmp_path / "data.csv").write_text("a,b\n1,2\n")
    out = tmp_path / "out"

    # chdir BEFORE seeding: both git helpers are anchored on the working directory, so
    # seeding from the repository and running from tmp_path would record two different
    # commits, the gate would reject the artifact, and the stale pickle would never be
    # loaded - leaving the guard this test exists for unexercised.
    monkeypatch.chdir(tmp_path)

    # Both halves of the gate are forced, because tmp_path is not a repository: `git status`
    # there fails (→ dirty) and `rev-parse HEAD` fails (→ "nogit"), and "nogit" is not a
    # commit match even against itself. A fixed stand-in commit is what lets the gate ADMIT
    # the artifact, which is the precondition for testing what happens at load.
    fake_commit = "cafe123"
    monkeypatch.setattr(runner, "get_git_commit", lambda: fake_commit)
    monkeypatch.setattr(runner, "is_tree_clean", lambda: True)

    # the sub-job's real identity names its cache dir; seed a prov record whose
    # identity + commit MATCH this run so the gate admits it, plus a stale pickle
    sub_identity = _import_job_module(sub_py).build_identity(tmp_path).digest
    cached = out / f"tiny_sub_{sub_identity[:6]}_20200101_000000"
    (cached / "provenance").mkdir(parents=True)
    with (cached / "panel_data.pkl").open("wb") as fh:
        pickle.dump(_StalePickle(WithinCalibrationPanelData, {"t_h": None}), fh)
    (cached / "provenance" / "panel_data.prov.json").write_text(
        json.dumps(
            {
                "identity": sub_identity,
                "git_commit": fake_commit,
                "tree_clean": True,
            }
        )
    )

    comp_py = tmp_path / "tiny_comp.py"
    comp_py.write_text("# synthetic composite job file\n")
    comp = Job(name="tiny_comp")
    inc = comp.include(sub_py, alias="s1")
    node = comp.step(_identity, inc.ref("panel_data"), name="reused")
    comp.materialize(node, name="reused_out")
    comp.job_file = comp_py.resolve()

    with pytest.raises(ValueError, match="stale WithinCalibrationPanelData"):
        run_job(comp, out, force=True, data_root=tmp_path, reuse_deps=True)


def test_real_pre_split_artifact_raises() -> None:
    """The actual pre-split artifacts on disk must fail to load through `load_artifact`.

    Two failure modes, and the second one is a deliberate, accepted loss.

    `ValueError` is a stale artifact: it loaded, and a dataclass in it lacks a field its class
    declares now, which `StaleArtifactGuard` or the loader's completeness check reports.

    `ModuleNotFoundError` or `AttributeError` is a class that was renamed or restructured, not
    moved: the vocabulary rename (`panels.non_repairable` -> `panels.within_calibration`) and
    the move of the within-calibration compute out of the render package. The src-layout move
    itself is resolved by the alias, so what remains here is what an alias must refuse.
    Recovering these means re-running the jobs that wrote them.

    This is recorded rather than hidden because `output/` is append-only and those
    artifacts are evidence. They are not recoverable by editing code; they are recovered by
    re-running the jobs that produced them, which is the accepted plan.

    Post-rename artifacts in the same locations must still load cleanly. Skips if no
    pre-split artifact exists (a fresh checkout has none).
    """
    import io

    from quebra.core._artifact_guard import load_artifact

    candidates = [
        p
        for root in ("output", "output_backup2")
        for p in (REPO_ROOT / root).glob("*/subjobs_output/*/t2star_panel_data.pkl")
    ]
    stale_errors: list[str] = []
    for pkl in candidates:
        try:
            obj = load_artifact(io.BytesIO(pkl.read_bytes()))
        except ValueError as exc:
            assert "stale" in str(exc) and "--reuse-deps" in str(exc)
            stale_errors.append(str(exc))
        except (ModuleNotFoundError, AttributeError) as exc:
            # The loader says "renamed or restructured, not moved" only after the package's
            # new home imported, so a wheel that failed to ship `quebra/panels/` raises its own
            # error and fails here instead of passing as an accepted loss.
            missing = str(exc)
            assert "renamed or restructured, not moved" in missing, exc
            stale_errors.append(missing)
        else:
            assert isinstance(obj, WithinCalibrationPanelData)
    if not stale_errors:
        pytest.skip("no pre-split artifact present on this machine")


# ------------------------------------------------------------------ moved-module aliasing


@dataclasses.dataclass(frozen=True)
class _HashableRecord:
    """A hashable record with a defaulted field, so a stale one can be a key or set member."""

    name: str
    note: str = ""


@dataclasses.dataclass
class _WithDerivedField:
    """A current class whose `init=False` fields are never in a fresh instance's state."""

    name: str
    derived: int = dataclasses.field(default=0, init=False)
    never_set: int = dataclasses.field(init=False)


@dataclasses.dataclass
class _GuardedWithDerivedField(StaleArtifactGuard):
    """The guarded twin: `StaleArtifactGuard.__setstate__` must exempt `init=False` too."""

    name: str
    derived: int = dataclasses.field(default=0, init=False)


def _as_if_written_before_the_move(obj: object, package: str = "analyzers") -> bytes:
    """Bytes of `obj` as a pickle written when its package sat at the repository root.

    Protocol 2 records each class as a text GLOBAL, `c<module>\\n<name>\\n`, so the module path
    can be rewritten without disturbing any length prefix. That is a faithful stand-in for the
    artifacts under `output/` that predate `src/quebra/`: same class, old path.
    """
    raw = pickle.dumps(obj, protocol=2)
    new_path = f"cquebra.{package}.".encode()
    assert new_path in raw, "fixture no longer exercises a quebra class"
    return raw.replace(new_path, f"c{package}.".encode())


def _check_result():
    from quebra.analyzers.checks.result import CheckResult

    return CheckResult(
        check="c1_lewis_robinson",
        statistic=1.25,
        p_value=0.5,
        calibration="permutation",
        clock="in_spec",
        n_events=40,
        n_segments=1,
        n_censored_dropped=0,
        notes="",
    )


def test_an_artifact_written_before_its_package_moved_loads_through_the_alias():
    """Oracle: the object written, compared field by field with the object read back.

    Checks the negative first. If a plain unpickle could read these bytes, the alias would not
    be what is being tested and the positive half would pass for the wrong reason.
    """
    import io

    from quebra.core._artifact_guard import load_artifact

    written = _check_result()
    old_bytes = _as_if_written_before_the_move(written)

    with pytest.raises(ModuleNotFoundError):
        pickle.loads(old_bytes)

    loaded = load_artifact(io.BytesIO(old_bytes))
    assert type(loaded) is type(written), "the alias resolved to a different class"
    assert loaded == written, (
        "the alias resolved the class but the fields came back wrong"
    )


def test_a_stale_artifact_is_refused_whatever_its_class() -> None:
    """Oracle: hand-built pickles, each with one field its class declares removed.

    The alias makes pickles written before the move loadable, and some of those classes have
    gained fields since. Only `StaleArtifactGuard` subclasses check themselves, so
    `load_artifact` checks every dataclass it can reach, and a class default never counts as
    the missing field: that default is what a stale pickle would silently read. Three shapes,
    each the shape of a real cached artifact: a plain dataclass with a default (the
    `KaplanMeierComparison` case), a slots dataclass inside a container (the TLF case), and a
    guarded class that moved and went stale.
    """
    import dataclasses
    import io

    from quebra.analyzers.tlf import TLFResult
    from quebra.core._artifact_guard import load_artifact

    def complete_tlf():
        fields = [
            f.name for f in dataclasses.fields(TLFResult) if f.name != "diagnostics"
        ]
        return TLFResult(**dict.fromkeys(fields))

    plain = _check_result()
    del plain.__dict__["notes"]  # frozen blocks setattr, not the state dict
    holder = complete_tlf()
    holder.gmm1 = plain  # a complete dataclass whose FIELD holds the stale one
    for where, root in (
        ("a dict value", {"result": plain}),
        ("a list", [plain]),
        ("a tuple", (plain,)),
        ("a dataclass field", holder),
    ):
        with pytest.raises(
            ValueError, match=r"stale CheckResult artifact.*\['notes'\]"
        ):
            load_artifact(io.BytesIO(_as_if_written_before_the_move(root)))
            pytest.fail(f"a stale dataclass in {where} loaded")

    # A key or a set member must hash while it unpickles, which a stale frozen record does
    # only when its missing field has a class default: the silent case. `CheckResult` holds
    # a dict and cannot hash, so a hashable record stands in.
    record = _HashableRecord(name="grid")
    del record.__dict__["note"]
    for where, root in (("a dict key", {record: "grid"}), ("a set", {record})):
        with pytest.raises(
            ValueError, match=r"stale _HashableRecord artifact.*\['note'\]"
        ):
            load_artifact(io.BytesIO(pickle.dumps(root, protocol=2)))
            pytest.fail(f"a stale dataclass in {where} loaded")

    tlf = TLFResult(
        **{
            f.name: None
            for f in dataclasses.fields(TLFResult)
            if f.name != "diagnostics"
        }
    )
    raw_complete = _as_if_written_before_the_move({"tlf": [tlf]})
    assert load_artifact(io.BytesIO(raw_complete))["tlf"][0] == tlf, (
        "control: complete loads"
    )
    del tlf.fit_failed
    with pytest.raises(ValueError, match=r"stale TLFResult artifact.*\['fit_failed'\]"):
        load_artifact(io.BytesIO(_as_if_written_before_the_move({"tlf": [tlf]})))

    # Control: `__init__` never writes an `init=False` field, so a FRESH instance lacks it
    # too. Refusing it would refuse a current artifact.
    fresh = _WithDerivedField(name="x")
    assert "derived" not in fresh.__dict__, "control no longer exercises init=False"
    assert "never_set" not in fresh.__dict__, (
        "control no longer exercises an unset field"
    )
    loaded = load_artifact(io.BytesIO(pickle.dumps(fresh, protocol=2)))
    assert (type(loaded), loaded.name) == (_WithDerivedField, "x")
    guarded_fresh = _GuardedWithDerivedField(name="y")
    assert "derived" not in guarded_fresh.__dict__, (
        "control no longer exercises init=False"
    )
    assert (
        load_artifact(io.BytesIO(pickle.dumps(guarded_fresh, protocol=2))).name == "y"
    )

    guarded = _valid_across_calibration()
    del guarded.__dict__["histogram_edges"]
    with pytest.raises(ValueError, match="stale AcrossCalibrationPanelData artifact"):
        load_artifact(io.BytesIO(_as_if_written_before_the_move(guarded, "panels")))


def test_an_alias_failure_names_its_real_cause(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Oracle: nine hand-built failure causes, each checked for the error it must raise.

    The alias never renames a class. A moved module whose class is not there under the same
    name, and a module that did not arrive at all, are the case of the three retired panel
    classes (`NonRepairablePanelData` and its siblings): renamed and restructured, not moved.
    Mapping them onto their successors would load old field values into a new shape, the
    wrong-but-plausible result AGENTS.md section 3 refuses, so both stay unreadable and say why.

    Everything else is a real error and must reach the reader as itself: a broken import
    inside the moved module, a missing dependency, and a missing new home (a broken install).
    None of them may be reported as a rename.
    """
    import io

    from quebra.core import _artifact_guard
    from quebra.core._artifact_guard import load_artifact

    renamed = _as_if_written_before_the_move(_check_result()).replace(
        b"\nCheckResult\n", b"\nNoSuchRenamedClass\n"
    )
    with pytest.raises(AttributeError, match="renamed or restructured, not moved"):
        load_artifact(io.BytesIO(renamed))

    gone = _as_if_written_before_the_move(_check_result()).replace(
        b"canalyzers.checks.result\n", b"canalyzers.checks.no_such_module\n"
    )
    with pytest.raises(ModuleNotFoundError, match="renamed or restructured, not moved"):
        load_artifact(io.BytesIO(gone))
    # A whole moved subpackage gone: the missing name is a dotted prefix of the module.
    gone_pkg = b"canalyzers.no_such_pkg.mod\nX\n."
    with pytest.raises(ModuleNotFoundError, match="renamed or restructured, not moved"):
        load_artifact(io.BytesIO(gone_pkg))

    # Protocol 4 names a nested class by its dotted qualname; a missing inner name is a
    # rename too, not a bare lookup failure.
    module, qualname = b"analyzers.checks.result", b"CheckResult.Missing"
    nested = (
        b"\x80\x04\x8c"
        + bytes([len(module)])
        + module
        + b"\x8c"
        + bytes([len(qualname)])
        + qualname
        + b"\x93."
    )
    with pytest.raises(AttributeError, match="renamed or restructured, not moved"):
        load_artifact(io.BytesIO(nested))

    pkg = tmp_path / "qre_moved_probe"
    pkg.mkdir()
    (pkg / "__init__.py").write_text("")
    (pkg / "broken.py").write_text(
        'raise ImportError("boom: a bug inside the moved module")\n'
    )
    (pkg / "needs_dep.py").write_text("import qre_no_such_dependency_probe\n")
    # A missing module whose name is a string prefix, not a dotted prefix, of the moved one.
    (pkg / "needs_sibling.py").write_text("import qre_moved_probe.needs_sib\n")
    (pkg / "raises_bare.py").write_text(
        'raise ModuleNotFoundError("bare, with no name")\n'
    )
    monkeypatch.syspath_prepend(str(tmp_path))
    monkeypatch.setitem(_artifact_guard.MODULE_ALIASES, "oldprobe", "qre_moved_probe")
    monkeypatch.setitem(
        _artifact_guard.MODULE_ALIASES, "lostprobe", "qre_no_such_home_probe"
    )
    try:
        with pytest.raises(ImportError, match="boom") as bug:
            load_artifact(io.BytesIO(b"coldprobe.broken\nX\n."))
        assert bug.type is ImportError, (
            f"a broken import surfaced as {bug.type.__name__}"
        )

        with pytest.raises(ModuleNotFoundError) as dep:
            load_artifact(io.BytesIO(b"coldprobe.needs_dep\nX\n."))
        assert dep.value.name == "qre_no_such_dependency_probe", dep.value
        assert "renamed" not in str(dep.value), "a missing dependency read as a rename"

        with pytest.raises(ModuleNotFoundError) as home:
            load_artifact(io.BytesIO(b"clostprobe.mod\nX\n."))
        assert home.value.name == "qre_no_such_home_probe", home.value
        assert "renamed" not in str(home.value), "a broken install read as a rename"

        with pytest.raises(ModuleNotFoundError) as sibling:
            load_artifact(io.BytesIO(b"coldprobe.needs_sibling\nX\n."))
        assert sibling.value.name == "qre_moved_probe.needs_sib", sibling.value
        assert "renamed" not in str(sibling.value), (
            "a sibling-prefix name read as a rename"
        )

        with pytest.raises(ModuleNotFoundError, match="bare, with no name") as bare:
            load_artifact(io.BytesIO(b"coldprobe.raises_bare\nX\n."))
        assert "renamed" not in str(bare.value), (
            "a nameless import error read as a rename"
        )
    finally:
        sys.modules.pop("qre_moved_probe", None)


def test_a_composite_reuses_a_sub_job_artifact_written_before_its_package_moved(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Oracle: the object seeded into the cache, compared with what the composite received.

    Proves the WIRING only: the runner's composite transport, the pipeline's one artifact read,
    goes through `load_artifact`, both halves of it. The alias resolves an old path, and the
    completeness check refuses a stale artifact at that read. The state built here cannot arise
    in production. The reuse gate admits only a run at the current commit on a clean tree, and
    such a run is written under current paths, so this test fakes a provenance record claiming
    the current commit over an old-path pickle, set up as
    `test_composite_reuse_of_reuse_eligible_stale_artifact_aborts` sets up a stale one.
    """
    from quebra.core.dataset import Dataset  # noqa: F401  (sub-job module imports it)
    from quebra.core.job import Job
    from quebra.core import runner
    from quebra.core.runner import run_job

    (tmp_path / "data.csv").write_text("a,b\n1,2\n")
    out = tmp_path / "out"
    monkeypatch.chdir(tmp_path)
    fake_commit = "cafe123"
    monkeypatch.setattr(runner, "get_git_commit", lambda: fake_commit)
    monkeypatch.setattr(runner, "is_tree_clean", lambda: True)

    def compose(tag: str, cached_bytes: bytes) -> Path:
        """Seed a reuse-eligible sub-job run holding `cached_bytes`, then run a composite."""
        sub_py = tmp_path / f"{tag}_sub.py"
        sub_py.write_text(
            "from quebra.core.job import Job\n"
            "from quebra.core.dataset import Dataset\n"
            f'job = Job(name="{tag}_sub")\n'
            'node = job.load_df(Dataset(path="data.csv", schema=None))\n'
            'job.materialize(node, name="panel_data")\n'
        )
        sub_identity = _import_job_module(sub_py).build_identity(tmp_path).digest
        cached = out / f"{tag}_sub_{sub_identity[:6]}_20200101_000000"
        (cached / "provenance").mkdir(parents=True)
        (cached / "panel_data.pkl").write_bytes(cached_bytes)
        (cached / "provenance" / "panel_data.prov.json").write_text(
            json.dumps(
                {
                    "identity": sub_identity,
                    "git_commit": fake_commit,
                    "tree_clean": True,
                }
            )
        )
        comp_py = tmp_path / f"{tag}_comp.py"
        comp_py.write_text("# synthetic composite job file\n")
        comp = Job(name=f"{tag}_comp")
        inc = comp.include(sub_py, alias="s1")
        node = comp.step(_identity, inc.ref("panel_data"), name="reused")
        comp.materialize(node, name="reused_out")
        comp.job_file = comp_py.resolve()
        run_job(comp, out, force=True, data_root=tmp_path, reuse_deps=True)
        produced = sorted(out.glob(f"{tag}_comp_*/reused_out.pkl"))
        assert produced, (
            "the composite wrote nothing, so it never reached the transport"
        )
        return produced[-1]

    written = _check_result()
    with compose("moved", _as_if_written_before_the_move(written)).open("rb") as fh:
        received = pickle.load(fh)
    assert received == written, (
        "the composite received a different object than was seeded"
    )

    stale = _check_result()
    del stale.__dict__["notes"]
    with pytest.raises(ValueError, match=r"stale CheckResult artifact.*\['notes'\]"):
        compose("stale", _as_if_written_before_the_move(stale))
