"""A record that fails to load or carve, in a job that scores many records.

Oracles:
- the loader's and the carve stages' own errors: a recorded failure names the error the plain
  path raises, at the stage that raised it;
- the failure's contract: it passes through every later stage unchanged, and its ledger reads
  `not computed` with the failure's note at every threshold and both clocks, so the record is
  never dropped;
- a ledger job on two synthetic records, one sound and one the schema refuses, run by the
  runner into a temporary output folder: the sound record is scored, the other is persisted
  as a `RecordFailure` with a ledger that says so.
"""

from __future__ import annotations

import pickle
from itertools import product
from pathlib import Path

import numpy as np
import pytest

from quebra.analyzers.check_ledger import VERDICT_NOT_COMPUTED, CheckLedger
from quebra.analyzers.checks.result import CLOCK_CALENDAR, CLOCK_IN_SPEC
from quebra.analyzers.windows import WindowsResult
from quebra.cli import _module_from_path
from quebra.core.dataset import Dataset
from quebra.core.job import _load_dataset, _load_dataset_or_failure
from quebra.core.runner import run_job
from quebra.core.types import Norm, RecordFailure
from quebra.ledger_recipe import CARVE_SUFFIX, LEDGER_DATA_SUFFIX, ledger_step
from quebra.recipes import T2STAR_THRESHOLDS, _carve_stage
from quebra.schemas.track912 import track912Schema

REPO = Path(__file__).resolve().parents[1]
SIZE_TABLE = REPO / "jobs" / "bench" / "results" / "size_table.csv"
_COLUMNS = [
    "timestamp",
    "amp",
    "amp error",
    "offset",
    "offset error",
    "T2star",
    "T2star error",
    "frequency",
    "frequency error",
    "phase",
    "phase_error",
    "normalised chi-square",
]
_KNOBS = {
    "alpha": 0.05,
    "min_events_pass": 35,
    "tie_cutoff_distinct": 5,
    "lag_max": 5,
    "n_permutations": 19,
    "seed": 1,
    "include_c3": False,
    "c3_n_null_sim": 200,
}


def _track912_csv(path: Path, *, drop: str | None = None) -> Path:
    """A synthetic record whose T2* wanders across the 2 to 5 us ladder."""
    rows = []
    for i in range(240):
        t2star = 3.5e-6 + 1.2e-6 * np.sin(i / 9.0)
        rows.append(
            [i * 30.0, 1.0, 0.02, 0.5, 0.005, t2star, 2e-7, 2.1e6, 6e3, 0.1, 0.03, 0.1]
        )
    columns = [c for c in _COLUMNS if c != drop]
    keep = [i for i, c in enumerate(_COLUMNS) if c != drop]
    lines = [",".join(columns)] + [",".join(str(r[i]) for i in keep) for r in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def _dataset(path: Path, name: str) -> Dataset:
    return Dataset(
        path=str(path),
        schema=track912Schema,
        qubit=1,
        device="probe",
        extra={"run_name": name},
    )


@pytest.mark.unit
def test_a_record_the_schema_refuses_is_recorded_with_the_error_it_raised(tmp_path):
    """Oracle: the plain loader's own error. The recording loader returns a failure at the
    load stage naming that error, and loads a sound record exactly as the plain one does."""
    bad = _dataset(_track912_csv(tmp_path / "070423_bad.csv", drop="T2star"), "bad")
    with pytest.raises(Exception) as raised:
        _load_dataset(bad)
    failure = _load_dataset_or_failure(bad)
    assert isinstance(failure, RecordFailure)
    assert (failure.dataset_id, failure.stage) == ("bad", "load")
    assert failure.reason.startswith(type(raised.value).__name__ + ":")

    good = _dataset(_track912_csv(tmp_path / "070423_ok.csv"), "ok")
    loaded = _load_dataset_or_failure(good)
    assert isinstance(loaded, Norm)
    assert np.array_equal(loaded["T2star_s"], _load_dataset(good)["T2star_s"])


@pytest.mark.unit
def test_a_missing_file_still_raises(tmp_path):
    """Oracle: the declaration. A file that is not there is the job's error, not the
    record's content, and the recording loader raises it."""
    with pytest.raises(FileNotFoundError):
        _load_dataset_or_failure(_dataset(tmp_path / "070423_none.csv", "none"))


@pytest.mark.unit
def test_a_carve_stage_records_its_refusal_and_passes_a_failure_through():
    """Oracle: the stage's own error. A failure from upstream comes out of every stage
    unchanged; a stage that refuses its input with a KeyError records it under its own
    name; an error the stages do not use to refuse a record still raises."""
    upstream = RecordFailure(dataset_id="r", stage="load", reason="X: y")
    windows_kwargs = {
        "gap_mult": 10.0,
        "k": 1.0,
        "use_uncertainty": True,
        "thresholds": list(T2STAR_THRESHOLDS),
    }
    for stage in ("t2star_filter", "t2star_final_filter_stage", "t2star"):
        assert _carve_stage(upstream, stage=stage, record="r") is upstream
    assert _carve_stage(upstream, stage="windows", record="r", **windows_kwargs) is (
        upstream
    )
    refused = _carve_stage({"meta": {"dataset_id": "r"}}, stage="t2star", record="r")
    assert isinstance(refused, RecordFailure)
    assert (refused.stage, refused.reason.split(":")[0]) == ("t2star", "KeyError")
    with pytest.raises(AttributeError):
        _carve_stage(None, stage="windows", record="r", **windows_kwargs)


@pytest.mark.unit
def test_a_failed_record_reads_not_computed_at_every_threshold_and_clock():
    """Oracle: the failure's contract. Its ledger has one `(all)` row per threshold and
    clock, each `not computed` with the failure's note, and the ledger panel draws it."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    from quebra.panels.check_ledger import CheckLedgerPanel

    failure = RecordFailure(dataset_id="2x2_q1_090623", stage="t2star", reason="E: z")
    ledger = ledger_step(
        failure,
        pd.read_csv(SIZE_TABLE),
        thresholds=list(T2STAR_THRESHOLDS),
        **_KNOBS,
    )
    rows = ledger.rows
    labels = [label for label, _v, _b in T2STAR_THRESHOLDS]
    assert ledger.dataset_id == "2x2_q1_090623"
    assert sorted(zip(rows["threshold_label"], rows["clock"])) == sorted(
        product(labels, (CLOCK_IN_SPEC, CLOCK_CALENDAR))
    )
    assert set(rows["check_id"]) == {"(all)"}
    assert set(rows["verdict"]) == {VERDICT_NOT_COMPUTED}
    assert set(rows["notes"]) == {"t2star failed: E: z"}
    fig = CheckLedgerPanel(name="failed_probe").build_matplotlib(ledger)
    plt.close(fig)


@pytest.mark.integration
def test_a_ledger_job_scores_the_sound_record_and_records_the_broken_one(tmp_path):
    """Oracle: the runner's own artifacts. Of two records, the one the schema refuses is
    persisted as a `RecordFailure` with a ledger reading `load failed` everywhere; the
    sound one is carved and scored; the job's bundle names both."""
    _track912_csv(tmp_path / "070423_ok.csv")
    _track912_csv(tmp_path / "070423_bad.csv", drop="T2star")
    job_py = tmp_path / "ledger_probe.py"
    job_py.write_text(
        "from quebra.core.dataset import Dataset\n"
        "from quebra.core.job import Job\n"
        "from quebra.ledger_recipe import configure_check_ledger_job\n"
        "from quebra.schemas.track912 import track912Schema\n"
        'job = Job(name="ledger_probe")\n'
        "def _ds(name):\n"
        "    return Dataset(path=f'070423_{name}.csv', schema=track912Schema, qubit=1,\n"
        "                   device='probe', extra={'run_name': f'probe_{name}'})\n"
        "configure_check_ledger_job(\n"
        "    job, prefix='probe', records=[_ds('ok'), _ds('bad')],\n"
        f"    bench_size_table=Dataset(path=r'{SIZE_TABLE}', schema=None),\n"
        f"    **{_KNOBS!r},\n"
        ")\n"
    )
    out = tmp_path / "output"
    run_job(
        _module_from_path(job_py).job,
        out,
        force=True,
        data_root=tmp_path,
        render_figures=False,
    )

    def artifact(name: str) -> object:
        (path,) = list(out.glob(f"ledger_probe_*/{name}.pkl"))
        return pickle.loads(path.read_bytes())

    broken = artifact(f"probe_bad{CARVE_SUFFIX}")
    assert isinstance(broken, RecordFailure)
    assert (broken.dataset_id, broken.stage) == ("probe_bad", "load")
    broken_ledger = artifact(f"probe_bad{LEDGER_DATA_SUFFIX}")
    assert isinstance(broken_ledger, CheckLedger)
    assert set(broken_ledger.rows["verdict"]) == {VERDICT_NOT_COMPUTED}
    assert all(n.startswith("load failed: ") for n in broken_ledger.rows["notes"])

    assert isinstance(artifact(f"probe_ok{CARVE_SUFFIX}"), WindowsResult)
    sound_ledger = artifact(f"probe_ok{LEDGER_DATA_SUFFIX}")
    assert (sound_ledger.rows["verdict"] != VERDICT_NOT_COMPUTED).any()
    assert artifact("probe_check_ledgers").labels == ("probe ok", "probe bad")
