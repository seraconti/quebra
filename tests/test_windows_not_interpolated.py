"""The carve must never consume a resampled series.

A gap policy is meaningless on a uniform grid, and an interpolated point is not an
observation - a window built from one reports lifetime the instrument never measured.

This is a WHITELIST, not a blacklist. Asserting "no node called interpolate" would pass
for any future decimation or resampling step under a different name, and would be
vacuous on jobs that have no interpolate node at all. Asserting the windows node's
ancestors are a subset of a known-observed set fails on anything new until someone
deliberately adds it here.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from quebra.core.dataset import Dataset
from quebra.core.discovery import discover
from quebra.core.job import Job
from quebra.core.reference import LocalRef
from quebra.core.runner import _ancestors

pytestmark = pytest.mark.policy


REPO_ROOT = Path(__file__).resolve().parents[1]
JOBS_DIR = REPO_ROOT / "jobs"

# Nodes whose output is observed reads, filtered but never resampled.
OBSERVED_ANCESTORS = {
    "load",
    "load_df",
    "lookup_prior",
    "t2star_filter",
    "t2star_final_filter_stage",
    "t2star",
    "filter",
    "final_filter_stage",
    "fidelity_raw",
    "windows",
    "fidelity_windows",
}


def _load_job(path: Path) -> Job:
    spec = importlib.util.spec_from_file_location(f"_jobmod_{path.stem}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    job = getattr(module, "job", None)
    if not isinstance(job, Job):
        raise AssertionError(f"{path} does not define a module-level `job`")
    return job


def _job_files() -> list[Path]:
    """Every DECLARED job, composites included.

    Composites included because a composite may carve in-graph, and a guard scoped to
    `jobs/active` would never examine it. `discover` is what makes widening safe: it reads
    `JOB_ID` statically without importing, so the eight `jobs/bench/*.py` study modules and the three `__init__.py` files,
    none of which define a module-level `job`, are skipped rather than crashing the loader.
    """
    return sorted(declared.path for declared in discover(JOBS_DIR).values())


# Names of operations that put a value on the series the instrument did not measure. Checked
# before the whitelist and independently of the suffix rule, so neither can cover for the
# other: a name carrying one of these is rejected whatever its prefix and whatever label
# vouches for it.
RESAMPLING_TOKENS = (
    "interp",
    "sample",  # resample, upsample, downsample, subsample
    "grid",  # regrid, gridded, uniform_grid
    "bin",  # rebin, binned
    "uniform",
    "decimate",
    "smooth",
    "spline",
    "pchip",
    "lerp",
    "fill",  # backfill, ffill
    "impute",
)


def _carve_label(job: Job, node_id: str) -> str | None:
    """The per-dataset label carried by ONE carve node, or None if it carries none.

    A job carving several datasets needs one node per dataset, and the discriminator it picks
    is its own business: `km_poster_6d2s` uses the file stem, the `check_ledger` jobs use a
    readable label (`q1_040423`) that appears nowhere in the Dataset. Requiring the suffix to
    BE a stem only fits the first convention, and this test has no business dictating node
    names that end up in a provenance graph. A per-dataset label appears on every step of that
    dataset's chain, so `windows_<label>` vouches for `filter_<label>`.

    Taken from the node UNDER TEST, never from the union over the job's carve nodes. Under the
    union a resampler named `windows_pchip`, sitting between `filter` and `windows`, matches
    `_windows_nodes`, contributes `pchip` to the label set, and is then vouched for by the
    label it injected. Scoped to one node it can only vouch for its own chain, and a carve
    node that IS the resampler is caught by the token list instead.

    Longest prefix wins, so `final_filter_stage_x` yields `x` rather than the `filter_stage_x`
    that the shorter `final` would give.
    """
    name = job.dag[node_id].fn_name
    for known in sorted(OBSERVED_ANCESTORS, key=len, reverse=True):
        if name.startswith(f"{known}_"):
            return name[len(known) + 1 :]
    return None


def _unexpected(job: Job, node_id: str) -> set[str]:
    """Ancestor step names that are not known-observed. ONE definition.

    The positive control below asserts against this function rather than re-deriving the
    set expression, so a control cannot keep passing while the real check drifts.
    """
    label = _carve_label(job, node_id)
    labels = {label} if label else set()
    unexpected = set()
    for ancestor in _ancestors(job, node_id):
        name = job.dag[ancestor].fn_name
        if any(token in name.lower() for token in RESAMPLING_TOKENS):
            unexpected.add(name)
            continue
        if name in OBSERVED_ANCESTORS:
            continue
        # `{whitelisted}_{label}`, where the label is one a carve node in this job carries.
        # Both halves must hold: the prefix names a step kind known to leave reads alone, and
        # the suffix names a dataset this job actually carves.
        base, _, suffix = name.rpartition("_")
        while suffix and base:
            if base in OBSERVED_ANCESTORS and suffix in labels:
                break
            base, _, tail = base.rpartition("_")
            suffix = f"{tail}_{suffix}"
        else:
            unexpected.add(name)
    return unexpected


def _windows_nodes(job: Job) -> list[str]:
    # Key on node_id / fn_name, never fn.__name__: every closure-factory step in
    # quebra/recipes.py names its inner function `step`, so filter, interpolate, allan
    # and fidelity all report the same fn.__name__.
    # Substring, not equality: the fidelity carve is named `fidelity_windows`, and an
    # exact match silently skipped it - the test passed while that node was in fact
    # carving interpolated data. A matcher that can miss a node makes this whole file
    # vacuous, so it must over-match rather than under-match.
    return [
        node_id
        for node_id, node in job.dag.items()
        if "windows" in node.fn_name or "windows" in node_id
    ]


def test_there_are_jobs_to_check() -> None:
    # Guards against the whole suite passing vacuously if jobs/active empties out.
    assert _job_files(), f"no job files found in {JOBS_DIR}"


@pytest.mark.parametrize("job_file", _job_files(), ids=lambda p: p.stem)
def test_windows_ancestors_are_observed_reads(job_file: Path) -> None:
    job = _load_job(job_file)
    for node_id in _windows_nodes(job):
        unexpected = _unexpected(job, node_id)
        assert not unexpected, (
            f"{job_file.name}: the windows node {node_id!r} consumes {sorted(unexpected)}, "
            f"which is not a known-observed step. If it does not resample the series, add "
            f"its step kind to OBSERVED_ANCESTORS deliberately - unless the name carries a "
            f"RESAMPLING_TOKENS substring, which no whitelist entry overrides."
        )


def test_the_check_catches_a_resampled_input() -> None:
    """Positive control: without it, the whitelist could be vacuously true."""

    def _passthrough(*xs: object) -> object:
        return xs[0]

    job = Job("synthetic_bad_job")
    job.job_file = REPO_ROOT / "jobs" / "active" / "__init__.py"
    raw = job.step(_passthrough, name="filter")
    resampled = job.step(_passthrough, raw, name="interpolate")
    carved = job.step(_passthrough, resampled, name="windows")

    assert isinstance(carved, LocalRef)
    assert _unexpected(job, "windows") == {"interpolate"}


def test_the_guard_actually_matches_something() -> None:
    """The guard can go silently dead, and 16 green tests would not say so.

    `_windows_nodes` matches on the substring `windows`. Half the discovered jobs carve
    nothing, so if the carve step were ever renamed - to `segments`, say - every parametrised
    case would iterate an empty list and pass, the three synthetic controls would keep passing
    on their own hardcoded names, and the suite would report full health over zero coverage.
    That is the failure the module docstring records `fidelity_windows` causing once already.

    So: assert the matcher still finds carves, in more than one job, including the
    `check_ledger` family, which carves the most: one carve per record, 28 records over six
    jobs.
    """
    per_job = {path.stem: len(_windows_nodes(_load_job(path))) for path in _job_files()}
    carving = {name: n for name, n in per_job.items() if n}
    assert len(carving) >= 2, f"the carve matcher found nothing to check: {per_job}"
    ledger = {
        name: n for name, n in per_job.items() if name.startswith("check_ledger_")
    }
    assert len(ledger) == 6 and sum(ledger.values()) >= 28, (
        f"the check_ledger family carves 28 records; the matcher sees {ledger}"
    )


def test_the_check_catches_a_resampled_input_behind_a_dataset_suffix() -> None:
    """Positive control for the suffix rule: the widening must not open a hole.

    A resampling step named after a dataset the job really loads is the exact case the
    suffix rule could have let through, so it is asserted rather than assumed.
    """

    def _passthrough(*xs: object) -> object:
        return xs[0]

    stem = "100423_6D2S_qubit1"
    job = Job("synthetic_suffixed_bad_job")
    job.job_file = REPO_ROOT / "jobs" / "active" / "__init__.py"
    # Loaded, so `stem` really is an allowed suffix for this job - without this the
    # control would pass for the trivial reason that no suffix is allowed at all.
    loaded = job.load(Dataset(path=f"data/real_private/6D2S/{stem}.pickle", qubit=1))
    filtered = job.step(_passthrough, loaded, name=f"filter_{stem}")
    resampled = job.step(_passthrough, filtered, name=f"interpolate_{stem}")
    job.step(_passthrough, resampled, name=f"windows_{stem}")

    # `filter_<stem>` is accepted; `interpolate_<stem>` is not. Naming a resampling step
    # after a dataset the job really loads is the case the suffix rule could have let
    # through, so it is asserted rather than assumed.
    assert _unexpected(job, f"windows_{stem}") == {f"interpolate_{stem}"}


def test_the_check_catches_a_resample_hidden_behind_a_whitelisted_prefix() -> None:
    """The control for the suffix half of the rule.

    `filter` is whitelisted, so a prefix-only rule would accept `filter_then_interpolate` -
    and any of `filter_smoothed`, `filter_pchip`, `filter_upsampled`, which no token list
    would reliably enumerate. Requiring the suffix to be a label a carve node in this job
    actually carries is what closes that, so it is asserted on a name whose prefix is
    whitelisted and whose suffix is not a carve label.
    """

    def _passthrough(*xs: object) -> object:
        return xs[0]

    job = Job("synthetic_prefixed_bad_job")
    job.job_file = REPO_ROOT / "jobs" / "active" / "__init__.py"
    loaded = job.load(
        Dataset(path="data/real_private/6D2S/100423_6D2S_qubit1.pickle", qubit=1)
    )
    sneaky = job.step(_passthrough, loaded, name="filter_then_interpolate")
    job.step(_passthrough, sneaky, name="windows_x")

    assert _unexpected(job, "windows_x") == {"filter_then_interpolate"}

    # Rejected on the suffix rule (the label `x` does not match), and independently on the
    # token list. Both layers are asserted, because either alone would leave a hole: the
    # suffix rule is defeated by naming the carve node to match, and the token list can never
    # enumerate every way to invent a value.
    for name in ("filter_smoothed", "filter_pchip", "filter_imputed"):
        by_suffix = Job("synthetic_mismatched_label")
        by_suffix.job_file = REPO_ROOT / "jobs" / "active" / "__init__.py"
        src = by_suffix.load(
            Dataset(path="data/real_private/6D2S/100423_6D2S_qubit1.pickle", qubit=1)
        )
        by_suffix.step(
            _passthrough, by_suffix.step(_passthrough, src, name=name), name="windows_x"
        )
        assert _unexpected(by_suffix, "windows_x") == {name}, name

        # And with the carve node renamed to match, which defeats the suffix rule outright.
        matched = Job("synthetic_matched_label")
        matched.job_file = REPO_ROOT / "jobs" / "active" / "__init__.py"
        src = matched.load(
            Dataset(path="data/real_private/6D2S/100423_6D2S_qubit1.pickle", qubit=1)
        )
        carve = f"windows_{name.split('_', 1)[1]}"
        matched.step(
            _passthrough, matched.step(_passthrough, src, name=name), name=carve
        )
        assert _unexpected(matched, carve) == {name, carve}, name
