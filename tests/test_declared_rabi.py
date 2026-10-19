"""The declared Rabi frequency: one row per record, joined exactly, refused otherwise.

Fidelity needs the record's Rabi frequency. It is a DECLARED value keyed by (device,
run_date, qubit), taken from the record's own `Dataset`, never searched by time. These tests
pin the join (exact match, kHz to Hz, its raises), the provenance it hands the fidelity
artifact, the checks `configure_ramsey_job` makes at graph build, the 2x2 archetype's wiring
and its declared ladder, and, where the private tree is present, that the real table names
every 2x2 record exactly once.
"""

from __future__ import annotations

import importlib.util
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from quebra import recipes
from quebra.analyzers import fidelity
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.core.runner import _ancestors
from quebra.recipes import configure_ramsey_job
from quebra.transforms.declared_rabi import join_declared_rabi

REPO_ROOT = Path(__file__).resolve().parents[1]
TABLE_PATH = (
    REPO_ROOT / "data" / "real_private" / "declared" / "rabi_frequencies_2x2.csv"
)
MANIFEST_PATH = REPO_ROOT / "data" / "real_private" / "MANIFEST.toml"
RECORD_2X2 = "data/real_private/2x2/030723_2x2_qubit1.pickle"
ARCHETYPE_2X2 = REPO_ROOT / "jobs" / "active" / "ramsey_2x2_q1_030723.py"


def _table(rows: list[tuple[str, str, int, float]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["device", "run_date", "qubit", "rabi_khz"])


# Qubit 1 on two dates with different values: a join that ignores the date sees two rows.
_TABLE = _table(
    [
        ("2x2", "2023-07-03", 1, 704.0),
        ("2x2", "2023-07-05", 1, 674.0),
        ("2x2", "2023-07-03", 2, 566.0),
    ]
)


def _norm(n: int = 5) -> dict[str, object]:
    return {
        "t_rel_s": np.arange(n, dtype=float) * 60.0,
        "delta_hz": np.zeros(n),
        "meta": {"dataset_id": "unit"},
    }


def _join(table: pd.DataFrame = _TABLE, run_date: str = "2023-07-03", qubit: int = 1):
    return join_declared_rabi(
        _norm(), table, device="2x2", run_date=run_date, qubit=qubit
    )


# --- the join ------------------------------------------------------------------------


@pytest.mark.unit
def test_the_join_writes_the_declared_row_in_hz_on_every_read() -> None:
    """Oracle: the table row. 704 kHz is 704000 Hz on every read, with its source named."""
    norm = _norm()
    out = join_declared_rabi(norm, _TABLE, device="2x2", run_date="2023-07-03", qubit=1)
    assert np.array_equal(out["rabi_hz"], np.full(5, 704_000.0))
    assert out["meta"]["rabi_source"] == "declared"
    assert out["meta"]["rabi_declared_key"] == "2x2/2023-07-03/qubit1"
    assert "rabi_source" not in norm["meta"], "the input norm must not be mutated"


@pytest.mark.unit
def test_the_join_matches_the_date_not_only_the_qubit() -> None:
    """Oracle: the table row. Qubit 1 has two dates; the second date's value is joined."""
    assert np.all(_join(run_date="2023-07-05")["rabi_hz"] == 674_000.0)


@pytest.mark.unit
def test_a_record_without_its_own_row_raises() -> None:
    """Oracle: specification. No row is borrowed from another date of the same qubit."""
    with pytest.raises(KeyError, match="no declared Rabi frequency for 2x2/2023-07-04"):
        _join(run_date="2023-07-04")


@pytest.mark.unit
def test_a_key_naming_two_rows_raises() -> None:
    """Oracle: specification. A duplicated key is refused, never resolved by picking one."""
    doubled = pd.concat([_TABLE, _TABLE.iloc[[0]]], ignore_index=True)
    with pytest.raises(
        ValueError, match="2 declared Rabi rows for 2x2/2023-07-03/qubit1"
    ):
        _join(table=doubled)


@pytest.mark.unit
@pytest.mark.parametrize("value", [0.0, -5.0, float("nan")])
def test_a_non_positive_or_non_finite_value_raises(value: float) -> None:
    """Oracle: specification. A Rabi frequency is finite and positive."""
    with pytest.raises(ValueError, match="must be finite and positive"):
        _join(table=_table([("2x2", "2023-07-03", 1, value)]))


@pytest.mark.unit
def test_a_table_missing_a_key_column_raises() -> None:
    """Oracle: specification. The join refuses a table it cannot key."""
    with pytest.raises(KeyError, match="no column"):
        _join(table=_TABLE.drop(columns=["run_date"]))


@pytest.mark.unit
def test_the_join_will_not_overwrite_an_existing_rabi_frequency() -> None:
    """Oracle: specification. A record's Rabi frequency has one source."""
    norm = {**_norm(), "rabi_hz": np.full(5, 1e6)}
    with pytest.raises(KeyError, match="already carries 'rabi_hz'"):
        join_declared_rabi(norm, _TABLE, device="2x2", run_date="2023-07-03", qubit=1)


@pytest.mark.unit
@pytest.mark.parametrize("source", [None, "missing"], ids=["absent", "schema_default"])
def test_fidelity_refuses_a_rabi_frequency_with_no_named_source(source) -> None:
    """Oracle: specification. The artifact must say where its Rabi frequency came from."""
    norm = {**_norm(), "rabi_hz": np.full(5, 704e3)}
    if source is not None:
        norm["meta"] = {**norm["meta"], "rabi_source": source}
    with pytest.raises(KeyError, match="requires meta\\['rabi_source'\\]"):
        fidelity.make_inputs_from_norm(norm, {})


@pytest.mark.unit
def test_the_fidelity_artifact_names_where_its_rabi_frequency_came_from() -> None:
    """Oracle: the join's metadata. The fidelity result carries the source and the key."""
    norm = _join()
    result = fidelity.run(fidelity.make_inputs_from_norm(norm, {}))
    assert result.meta["rabi_source"] == "declared"
    assert result.meta["rabi_declared_key"] == "2x2/2023-07-03/qubit1"
    assert np.all(result.frame["rabi_frequency_hz"] == 704_000.0)


# --- graph build -----------------------------------------------------------------------


def _record(qubit: object = 1, device: object = "2x2", **extra: object) -> Dataset:
    return Dataset(path=RECORD_2X2, qubit=qubit, device=device, extra=dict(extra))


_TABLE_DS = Dataset(path=str(TABLE_PATH), loader_kwargs={"comment": "#"})


def _configure(
    dataset: object, rabi_table: object, *, include_fidelity: bool = True
) -> Job:
    job = Job("synthetic_rabi_wiring")
    job.job_file = REPO_ROOT / "jobs" / "active" / "__init__.py"
    configure_ramsey_job(
        job,
        dataset,
        profile="overnight",
        include_fidelity=include_fidelity,
        include_tlf=True,
        rabi_table=rabi_table,  # type: ignore[arg-type]
    )
    return job


def _loaded_node() -> object:
    job = Job("synthetic_loaded")
    job.job_file = REPO_ROOT / "jobs" / "active" / "__init__.py"
    return job.load(_record(run_date="2023-07-03"))


@pytest.mark.integration
@pytest.mark.parametrize(
    ("dataset", "table", "error", "match"),
    [
        (_loaded_node, _TABLE_DS, TypeError, "needs the record's Dataset"),
        (lambda: _record(run_date="2023-07-03"), None, ValueError, "needs rabi_table"),
        (_record, _TABLE_DS, ValueError, r"join needs Dataset extra\['run_date'\]"),
        (
            lambda: _record(device=None, run_date="2023-07-03"),
            _TABLE_DS,
            ValueError,
            "join needs Dataset device",
        ),
        (
            lambda: _record(qubit=None, run_date="2023-07-03"),
            _TABLE_DS,
            ValueError,
            "join needs Dataset qubit",
        ),
        (
            lambda: _record(run_date="2023-07-04"),
            _TABLE_DS,
            ValueError,
            "declares run_date='2023-07-04' but the file",
        ),
        (
            lambda: _record(qubit=2, run_date="2023-07-03"),
            _TABLE_DS,
            ValueError,
            "declares qubit=2 but the file",
        ),
        (
            lambda: _record(device="6D2S", run_date="2023-07-03"),
            _TABLE_DS,
            ValueError,
            "declares device='6D2S' but the file",
        ),
        (
            lambda: _record(run_date="2023-07-03"),
            "rabi_frequencies_2x2.csv",
            TypeError,
            "rabi_table must be a Dataset",
        ),
        (
            lambda: Dataset(
                path="data/real_private/2x2/no_date.pickle",
                qubit=1,
                device="2x2",
                extra={"run_date": "2023-07-03"},
            ),
            _TABLE_DS,
            ValueError,
            "is not DDMMYY_<device>_qubit<N>",
        ),
    ],
    ids=[
        "loaded_node",
        "no_table",
        "no_run_date",
        "no_device",
        "no_qubit",
        "wrong_date",
        "wrong_qubit",
        "wrong_device",
        "table_not_a_dataset",
        "undated_name",
    ],
)
def test_a_mis_declared_fidelity_job_raises_at_graph_build(
    dataset, table, error, match
) -> None:
    """Oracle: specification. Each way a job can mis-declare the Rabi key fails at build."""
    with pytest.raises(error, match=match):
        _configure(dataset(), table)


@pytest.mark.integration
def test_a_rabi_table_without_fidelity_raises() -> None:
    """Oracle: specification. A table nothing would read is refused, not ignored."""
    with pytest.raises(
        ValueError, match="include_fidelity=False; nothing would read it"
    ):
        _configure(_record(run_date="2023-07-03"), _TABLE_DS, include_fidelity=False)


@pytest.mark.integration
def test_the_join_feeds_fidelity_and_nothing_else() -> None:
    """Oracle: specification. The join sits on the filtered reads, upstream of fidelity
    only; neither Allan nor TLF depends on the Rabi table."""
    job = _configure(_record(run_date="2023-07-03"), _TABLE_DS)
    assert "declared_rabi" in _ancestors(job, "fidelity_raw")
    assert "final_filter_stage" in _ancestors(job, "declared_rabi")
    assert "declared_rabi" not in _ancestors(job, "allan")
    assert "declared_rabi" not in _ancestors(job, "tlf")
    assert job.dag["declared_rabi"].kwargs == {
        "device": "2x2",
        "run_date": "2023-07-03",
        "qubit": 1,
    }


def _load_job(path: Path) -> Job:
    spec = importlib.util.spec_from_file_location(f"_jobmod_{path.stem}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.job


@pytest.mark.integration
def test_the_2x2_archetype_computes_fidelity() -> None:
    """Oracle: the job file. The flags it declares are what runs; nothing turns fidelity off
    behind it."""
    job = _load_job(ARCHETYPE_2X2)
    assert {"declared_rabi", "fidelity_raw", "fidelity_windows"} <= set(job.dag)


# --- the declared ladder ---------------------------------------------------------------


@pytest.mark.unit
def test_the_fidelity_ladder_is_the_eight_nines_rungs() -> None:
    """Oracle: specification. 99% to 99.9999999%, i.e. infidelity 1e-2 to 1e-9, small good."""
    assert [value for _, value, _ in fidelity.FIDELITY_LADDER] == [
        10.0**-k for k in range(2, 10)
    ]
    assert [label for label, _, _ in fidelity.FIDELITY_LADDER][:2] == ["99%", "99.9%"]
    assert not any(big for _, _, big in fidelity.FIDELITY_LADDER)


@pytest.mark.integration
def test_both_fidelity_steps_receive_the_whole_declared_ladder() -> None:
    """Oracle: the declared ladder. The carve and the panel get the same eight rungs."""
    job = _load_job(ARCHETYPE_2X2)
    ladder = list(fidelity.FIDELITY_LADDER)
    assert job.dag["fidelity_windows"].kwargs["thresholds"] == ladder
    assert job.dag["fidelity_panel_data"].kwargs["thresholds"] == ladder


@pytest.mark.unit
def test_the_carve_keeps_rungs_the_record_never_crosses() -> None:
    """Oracle: the declared ladder. A record spanning 1e-3 to 9e-3 is still carved on all
    eight rungs; the ones it never reaches give one-sided windows instead of vanishing."""
    n = 40
    t_rel_s = np.arange(n, dtype=float) * 60.0
    frame = pd.DataFrame(
        {"t_rel_s": t_rel_s, "infidelity": np.geomspace(1e-3, 9e-3, n)}
    )
    result = fidelity.FidelityResult(frame=frame, meta={"dataset_id": "unit"})
    carved = recipes._fidelity_windows(result, 10.0, list(fidelity.FIDELITY_LADDER))
    labels = [label for label, _, _ in fidelity.FIDELITY_LADDER]
    assert sorted(set(carved.reads["threshold_label"])) == sorted(labels)
    windowed = set(carved.windows_in_spec["threshold_label"]) | set(
        carved.windows_out_of_spec["threshold_label"]
    )
    assert windowed == set(labels), "every rung has windows, on one side or both"


# --- the real table against the manifest -------------------------------------------------


def _unmatched(
    table: pd.DataFrame, records: list[tuple[str, int]]
) -> tuple[list[tuple[str, int]], list[tuple[str, int]]]:
    """(records without exactly one row, rows naming no record), keyed by (run_date, qubit)."""
    rows = [(str(d), int(q)) for d, q in zip(table["run_date"], table["qubit"])]
    without = [r for r in records if rows.count(r) != 1]
    orphans = [r for r in rows if r not in records]
    return without, orphans


@pytest.mark.unit
def test_the_table_check_sees_a_missing_row() -> None:
    """Oracle: construction. A table one row short is reported, so the real check can fail."""
    records = [("2023-07-03", 1), ("2023-07-05", 1)]
    short = _table([("2x2", "2023-07-03", 1, 704.0)])
    assert _unmatched(short, records) == ([("2023-07-05", 1)], [])


@pytest.mark.policy
@pytest.mark.real
def test_every_2x2_record_has_exactly_one_declared_row() -> None:
    """Oracle: data/real_private/MANIFEST.toml. One row per 2x2 record, and no orphan row."""
    if not TABLE_PATH.exists():
        pytest.skip(f"no declared Rabi table at {TABLE_PATH}")
    table = pd.read_csv(TABLE_PATH, comment="#")
    with MANIFEST_PATH.open("rb") as handle:
        manifest = tomllib.load(handle)
    records = [
        (r["run_date"], int(r["qubit"]))
        for r in manifest["record"]
        if r.get("device") == "2x2"
    ]
    assert records, "the manifest lists no 2x2 record"
    assert set(table["device"]) == {"2x2"}
    assert _unmatched(table, records) == ([], [])
