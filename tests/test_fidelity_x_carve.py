"""The infidelity record carved on x = |Delta| / Omega: the reference, the units, the carve
and the panel's relative-error flag.

Subject: `fidelity.run` (the first kept read as reference, x and sigma_x), the unit refusals
at load, `recipes._fidelity_windows` and its wiring, and the panel builder's
`include_relative_error`. Each test names its oracle: the exact interval on the
infidelity, the definition of x and its sigma, the declaration, or the flag's contract.
"""

from __future__ import annotations

import importlib.util
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from quebra import recipes, tellability_recipe
from quebra.analyzers import fidelity, windows
from quebra.analyzers.fidelity import DETUNING_LADDER, FIDELITY_LADDER, infidelity_at
from quebra.analyzers.within_calibration_compute import (
    build_within_calibration_panel_data,
)
from quebra.core.dataset import Dataset
from quebra.core.job import _load_dataset
from quebra.schemas.track912 import track912Schema
from quebra.transforms.declared_rabi import join_declared_rabi

REPO_ROOT = Path(__file__).resolve().parents[1]
ARCHETYPE_2X2 = REPO_ROOT / "jobs" / "active" / "ramsey_2x2_q1_030723.py"
RABI_HZ = 1.0e6
_RUNG_OF = {label: value for label, value, _ in FIDELITY_LADDER}


def _result(
    x_signed: np.ndarray,
    s: np.ndarray,
    *,
    t_step_s: float = 60.0,
    rabi_hz: float = RABI_HZ,
):
    """A fidelity result whose reads have x = x_signed and sigma_x = s, read 0 first."""
    n = len(x_signed)
    norm = {
        "t_rel_s": np.arange(n, dtype=float) * t_step_s,
        "raw_frequency_hz": np.asarray(x_signed, dtype=float) * rabi_hz,
        "frequency_error_hz": np.asarray(s, dtype=float) * rabi_hz,
        "rabi_hz": np.full(n, rabi_hz),
        "meta": {"dataset_id": "unit", "rabi_source": "declared"},
    }
    return fidelity.run(fidelity.make_inputs_from_norm(norm, {}))


def _exact_interval_states(value, sigma, k, rung):
    """The state each read gets from the exact infidelity interval
    [I(max(x - k s, 0)), I(x + k s)] against the rung T: in spec when I(x) <= T, unsure
    unless the whole interval lies on one side."""
    lower = infidelity_at(np.maximum(value - k * sigma, 0.0))
    upper = infidelity_at(value + k * sigma)
    inside = infidelity_at(value) <= rung
    sure = (upper <= rung) | (lower >= rung)
    return np.where(
        inside,
        np.where(sure, windows.STATE_IN_SPEC, windows.STATE_IN_SPEC_UNCERTAIN),
        np.where(sure, windows.STATE_OUT_OF_SPEC, windows.STATE_OUT_OF_SPEC_UNCERTAIN),
    )


# --- the carve on |x| ----------------------------------------------------------------


@pytest.mark.statistical
def test_x_carve_labels_equal_the_exact_infidelity_interval() -> None:
    """Oracle: the exact interval on the infidelity, [I(max(x - k s, 0)), I(x + k s)],
    which is exact because I rises strictly with |x| below PEAK_X. On 10 000 random
    (x, s, k), 100 carves of 100 reads each with its own k, every read's state at every
    rung equals the interval's, built from the generated x and s, never from the carve's
    own read table. Omega is a power of two, so x * Omega / Omega is x exactly and the
    carve must read the generated values bit for bit."""
    rng = np.random.default_rng(20261031)
    x_t = np.array([value for _, value, _ in DETUNING_LADDER])
    rabi_hz = 2.0**20
    t_step_s = 60.0
    n_compared = 0
    for _ in range(100):
        k = float(rng.uniform(0.5, 3.0))
        near = x_t[rng.integers(0, len(x_t), size=100)]
        x_hat = near * rng.uniform(0.0, 2.5, size=100)
        s = near * 10.0 ** rng.uniform(-2.0, 0.5, size=100)
        x_hat[0] = 0.0  # read 0 is the reference
        assert np.all(x_hat + k * s < fidelity.PEAK_X)
        sign = np.where(rng.random(100) < 0.5, -1.0, 1.0)
        carved = recipes._fidelity_windows(
            _result(sign * x_hat, s, t_step_s=t_step_s, rabi_hz=rabi_hz),
            10.0,
            k,
            True,
            list(DETUNING_LADDER),
        )
        reads = carved.reads
        value = np.maximum(x_hat, fidelity.ABS_X_FLOOR)
        for label, _, _ in DETUNING_LADDER:
            rows = reads[reads["threshold_label"] == label]
            assert np.array_equal(
                rows["t_read_s"].to_numpy(), np.arange(100) * t_step_s
            )
            assert np.array_equal(rows["value"].to_numpy(dtype=float), value), label
            assert np.array_equal(rows["sigma_v"].to_numpy(dtype=float), s), label
            expected = _exact_interval_states(value, s, k, _RUNG_OF[label])
            assert np.array_equal(rows["state"].to_numpy(dtype=object), expected), label
            n_compared += len(rows)
    assert n_compared == 100 * 100 * len(DETUNING_LADDER)


@pytest.mark.unit
def test_abs_x_and_sigma_x_do_not_depend_on_k() -> None:
    """Oracle: the definition, x = |Delta| / Omega and s = sigma_f / Omega carry no k.
    Carved at k = 1 and k = 2, the read table's value and sigma are identical and the
    in/out side too; only the unsure states move."""
    rng = np.random.default_rng(7)
    x_signed = rng.normal(scale=0.08, size=200)
    x_signed[0] = 0.0
    result = _result(x_signed, np.full(200, 0.01))
    one, two = (
        recipes._fidelity_windows(result, 10.0, k, True, list(DETUNING_LADDER)).reads
        for k in (1.0, 2.0)
    )
    for column in ("value", "sigma_v", "in_spec", "margin"):
        assert np.array_equal(one[column].to_numpy(), two[column].to_numpy()), column
    assert not np.array_equal(one["state"].to_numpy(), two["state"].to_numpy())


@pytest.mark.unit
def test_the_reference_is_the_first_read() -> None:
    """Oracle: the definition, Delta_i = f_i - f_0. Read 0 has Delta = 0, |x| = 0 and
    infidelity 0 exactly; every other read's Delta is its frequency minus read 0's."""
    f = np.array([2.10e6, 2.13e6, 2.07e6, 2.16e6, 2.10e6])
    norm = {
        "t_rel_s": np.arange(5, dtype=float) * 60.0,
        "raw_frequency_hz": f,
        "frequency_error_hz": np.full(5, 3e3),
        "rabi_hz": np.full(5, 7.0e5),
        "meta": {"dataset_id": "unit", "rabi_source": "declared"},
    }
    frame = fidelity.run(fidelity.make_inputs_from_norm(norm, {})).frame
    assert frame["delta_from_first_hz"].iloc[0] == 0.0
    assert frame["abs_x"].iloc[0] == 0.0
    assert frame["infidelity"].iloc[0] == 0.0
    assert np.array_equal(frame["delta_from_first_hz"].to_numpy(), f - f[0])
    assert np.array_equal(frame["x"].to_numpy(), (f - f[0]) / 7.0e5)
    assert np.array_equal(frame["sigma_x"].to_numpy(), np.full(5, 3e3 / 7.0e5))


_TRACK912_HEADER = (
    "timestamp,amp,amp error,offset,offset error,T2star,T2star error,frequency,"
    "frequency error,phase,phase_error,normalised chi-square\n"
)


def _row(t: float, frequency: str, t2star: str = "3e-6", chi: str = "0.1") -> str:
    return f"{t:.1f},1.0,0.02,0.5,0.005,{t2star},4e-7,{frequency},6e3,0.1,0.03,{chi}\n"


def _through_the_load(tmp_path, rows: list[str]) -> fidelity.FidelityResult:
    """Schema, normaliser, filter, declared join and `fidelity.run`, as the recipe wires
    them."""
    path = tmp_path / "record.csv"
    path.write_text(_TRACK912_HEADER + "".join(rows))
    norm = _load_dataset(
        Dataset(
            path=path,
            schema=track912Schema,
            qubit=1,
            device="2x2",
            extra={"run_start_unix_s": 1.7e9},
        )
    )
    filtered = recipes._final_stage(recipes._filter_step(recipes.RAMSEY_CONFIG)(norm))
    table = pd.DataFrame(
        [("2x2", "2023-07-03", 1, 700.0)],
        columns=["device", "run_date", "qubit", "rabi_khz"],
    )
    joined = join_declared_rabi(
        filtered, table, device="2x2", run_date="2023-07-03", qubit=1
    )
    return fidelity.run(fidelity.make_inputs_from_norm(joined, recipes.RAMSEY_CONFIG))


@pytest.mark.integration
@pytest.mark.parametrize(
    "first_row",
    [_row(0, "inf"), _row(0, "2.0e6", t2star="")],
    ids=["non_finite_frequency", "missing_t2star"],
)
def test_the_reference_is_the_first_read_the_load_keeps(tmp_path, first_row) -> None:
    """Oracle: the definition, the reference is the first read that survives the load.
    A first raw row the load drops (a non-finite frequency, by the normaliser; a missing
    T2*, by the schema) moves the reference to the next kept read."""
    kept = [_row(20, "2.11e6"), _row(40, "2.14e6"), _row(60, "2.09e6")]
    frame = _through_the_load(tmp_path, [first_row, *kept]).frame
    assert len(frame) == 3
    assert np.array_equal(
        frame["delta_from_first_hz"].to_numpy(), np.array([0.0, 0.03e6, -0.02e6])
    )
    assert frame["infidelity"].iloc[0] == 0.0


@pytest.mark.integration
def test_a_first_read_over_the_chi_square_line_refuses_the_record(tmp_path) -> None:
    """Oracle: the filter's contract, a guard and not a cut. A read at or above the
    chi-square threshold refuses the record, so the reference never moves on that
    account."""
    rows = [_row(0, "2.0e6", chi="0.9"), _row(20, "2.11e6"), _row(40, "2.14e6")]
    with pytest.raises(ValueError, match="normalised chi-square"):
        _through_the_load(tmp_path, rows)


# --- units ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("key", ["frequency_error_rad_s", "rabi_rad_s"])
def test_fidelity_refuses_a_frequency_declared_in_radians(key) -> None:
    """Oracle: the cyclic-hertz declaration. x divides one frequency by another, so a
    frequency declared in rad/s beside the hertz ones is refused, naming its key."""
    norm = {
        "t_rel_s": np.arange(3, dtype=float),
        "raw_frequency_hz": np.full(3, 2e6),
        "frequency_error_hz": np.full(3, 1e3),
        "rabi_hz": np.full(3, 7e5),
        key: np.full(3, 1.0),
        "meta": {"dataset_id": "unit", "rabi_source": "declared"},
    }
    with pytest.raises(ValueError, match=key):
        fidelity.make_inputs_from_norm(norm, {})


@pytest.mark.unit
def test_the_declared_join_refuses_an_angular_rabi_column() -> None:
    """Oracle: the cyclic-hertz declaration. A Rabi table that gives an angular column is
    refused, naming the column, whether or not it also carries `rabi_khz`."""
    norm = {"t_rel_s": np.arange(3, dtype=float), "meta": {"dataset_id": "unit"}}
    for columns, row in (
        (
            ["device", "run_date", "qubit", "rabi_rad_s"],
            ("2x2", "2023-07-03", 1, 4.4e6),
        ),
        (
            ["device", "run_date", "qubit", "rabi_khz", "rabi_rad_s"],
            ("2x2", "2023-07-03", 1, 700.0, 4.4e6),
        ),
    ):
        table = pd.DataFrame([row], columns=columns)
        with pytest.raises(ValueError, match="rabi_rad_s"):
            join_declared_rabi(
                norm, table, device="2x2", run_date="2023-07-03", qubit=1
            )
    for name in ("Omega (rad/s)", "rabi_angular_hz", "omega_rad_per_s"):
        table = pd.DataFrame(
            [("2x2", "2023-07-03", 1, 700.0, 4.4e6)],
            columns=["device", "run_date", "qubit", "rabi_khz", name],
        )
        with pytest.raises(ValueError, match="angular unit"):
            join_declared_rabi(
                norm, table, device="2x2", run_date="2023-07-03", qubit=1
            )


@pytest.mark.unit
def test_the_declared_join_reads_columns_that_only_contain_the_letters_rad() -> None:
    """Oracle: the token rule. `rad` inside a longer word, or an `omega` column with no
    unit token, declares no angular unit, so the join reads `rabi_khz` as usual."""
    norm = {"t_rel_s": np.arange(3, dtype=float), "meta": {"dataset_id": "unit"}}
    table = pd.DataFrame(
        [("2x2", "2023-07-03", 1, 700.0, 712.0, "set by hand")],
        columns=[
            "device",
            "run_date",
            "qubit",
            "rabi_khz",
            "rabi_khz_upgraded",
            "omega_notes",
        ],
    )
    joined = join_declared_rabi(
        norm, table, device="2x2", run_date="2023-07-03", qubit=1
    )
    assert np.array_equal(joined["rabi_hz"], np.full(3, 700.0e3))


# --- wiring --------------------------------------------------------------------------


def _load_job(path: Path):
    spec = importlib.util.spec_from_file_location(f"_jobmod_{path.stem}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.job


@pytest.mark.integration
def test_the_fidelity_carve_declares_k_and_uncertainty_and_tellability_reads_them() -> (
    None
):
    """Oracle: the declaration. The 2x2 job's carve gets FIDELITY_K and
    FIDELITY_USE_UNCERTAINTY as step kwargs, so both reach provenance, and the infidelity
    tellability reads its unsure level from the carve, not from a default."""
    job = _load_job(ARCHETYPE_2X2)
    kwargs = job.dag["fidelity_windows"].kwargs
    assert kwargs["k"] == recipes.FIDELITY_K
    assert kwargs["use_uncertainty"] is recipes.FIDELITY_USE_UNCERTAINTY is True

    rng = np.random.default_rng(3)
    x_signed = rng.normal(scale=0.05, size=60)
    x_signed[0] = 0.0
    result = _result(x_signed, np.full(60, 0.01))
    carved = recipes._fidelity_windows(result, 10.0, 2.0, True, list(DETUNING_LADDER))
    told = tellability_recipe._infidelity_tellability(
        result,
        carved,
        thresholds=list(DETUNING_LADDER),
        gap_mult=10.0,
        n_permutations=19,
        seed=1,
    )
    assert told.settings["unsure_k"] == 2.0


# --- the panel's relative-error flag -------------------------------------------------


@pytest.mark.unit
def test_the_relative_error_flag_changes_what_is_drawn_and_no_data() -> None:
    """Oracle: the flag's contract, a render flag. Built with `include_relative_error` on
    and off from one carve, the three bands are byte-identical when pickled, the flag is
    the only field that differs, and the fidelity adapter sets it off."""
    rng = np.random.default_rng(11)
    x_signed = rng.normal(scale=0.08, size=120)
    x_signed[0] = 0.0
    result = _result(x_signed, np.full(120, 0.01))
    carved = recipes._fidelity_windows(result, 10.0, 1.0, True, list(DETUNING_LADDER))
    km = recipes._kaplan_meier_set(recipes._event_tables(carved), 0.95)
    common = dict(
        t_h=result.frame["t_rel_s"].to_numpy(dtype=float) / 3600.0,
        primary_series=fidelity.panel_series(result),
        primary_label="|x|",
        thresholds=list(DETUNING_LADDER),
        meta={"dataset": "unit"},
        windows=carved.windows_in_spec,
        reads=carved.reads,
        gap_spans_s=carved.diagnostics["gap_spans_s"],
        kaplan_meier=km,
        primary_sigma=result.frame["sigma_x"].to_numpy(dtype=float),
        use_log_scale=True,
        use_uncertainty=True,
    )
    on, off = (
        build_within_calibration_panel_data(**common, include_relative_error=flag)
        for flag in (True, False)
    )
    for band in ("signal", "distinguish", "reliability"):
        assert pickle.dumps(getattr(on, band)) == pickle.dumps(getattr(off, band)), band
    assert (on.include_relative_error, off.include_relative_error) == (True, False)
    panel = recipes._fidelity_panel_data(result, carved, km, 0, list(DETUNING_LADDER))
    assert panel.include_relative_error is False
    assert panel.meta["caption_note"] == fidelity.DAMAGE_NOTE


# --- tellability on x ----------------------------------------------------------------


@pytest.mark.unit
def test_the_infidelity_threshold_rows_use_sigma_x_on_a_hand_record() -> None:
    """Oracle: a hand-computed table. Six reads 60 s apart at the 99% rung (x_T = 0.1737),
    k = 1: |x| = 0, 0.16, 0.10, 0.30, 0.25, 0.05 with sigma_x 0.01 except read 1 (0.02,
    so |x_T - x| = 0.0137 < k s: unsure). Intervals, each held by its left read: 0-1
    settled in; 1-2 unsure; 2-3 and 4-5 cross the rung; 3-4 settled out. The time-weighted
    in-spec fraction lies in [60/300, 240/300] = [0.2, 0.8], against a read fraction of
    4/6; sigma is declared, and |x| has no floor, so the clipped row is not computed."""
    x_signed = np.array([0.0, 0.16, -0.10, 0.30, -0.25, 0.05])
    s = np.array([0.01, 0.02, 0.01, 0.01, 0.01, 0.01])
    result = _result(x_signed, s)
    rung = [DETUNING_LADDER[0]]
    carved = recipes._fidelity_windows(result, 10.0, 1.0, True, rung)
    told = tellability_recipe._infidelity_tellability(
        result, carved, thresholds=rung, gap_mult=10.0, n_permutations=19, seed=1
    )
    (row,) = told.thresholds
    assert told.metric.sigma_status == "declared"
    assert row.status["bounds"] == row.status["crossing_counts"] == "computed"
    assert row.in_spec_fraction_lower == pytest.approx(0.2, abs=1e-12)
    assert row.in_spec_fraction_upper == pytest.approx(0.8, abs=1e-12)
    assert told.metric.status["clipped"] == "not computed: no floor for this metric"
    assert told.settings["threshold_value_key"] == "abs_x"
    assert told.metric.value_key == "x"


def _two_by_two_shaped_norm(seed: int = 9, n: int = 400) -> dict[str, object]:
    """Synthetic, shaped like a 2x2 record: irregular reads with one gap, a drifting
    fitted frequency with per-read errors, and one declared Rabi frequency."""
    rng = np.random.default_rng(seed)
    steps = rng.uniform(10.0, 40.0, size=n - 1)
    steps[n // 2] += 5000.0
    t = np.r_[0.0, np.cumsum(steps)]
    sigma_f = rng.uniform(2e3, 6e3, size=n)
    walk = np.cumsum(rng.normal(scale=1.5e3, size=n))
    return {
        "t_rel_s": t,
        "raw_frequency_hz": 2.1e6 + walk + rng.normal(scale=sigma_f),
        "frequency_error_hz": sigma_f,
        "rabi_hz": np.full(n, 6.5e5),
        "meta": {"dataset_id": "unit", "rabi_source": "declared"},
    }


@pytest.mark.unit
def test_the_infidelity_metric_rows_are_the_frequency_rows_rescaled() -> None:
    """Oracle: the identity x = (f - f_0) / Omega with sigma_x = sigma_f / Omega, and a
    floating-mean spectrum. On the same reads, every scale-free metric row of the infidelity
    record (G10 and its halves, the second-difference ratio, the noise-floor ratio, the
    self-consistency and sigma-bound rows) equals the frequency record's, and the spectrum,
    edge and noise-floor levels are the frequency's over Omega^2."""
    norm = _two_by_two_shaped_norm()
    omega = 6.5e5
    frequency = tellability_recipe._frequency_tellability(
        norm,
        value_key="raw_frequency_hz",
        sigma_key="frequency_error_hz",
        gap_mult=10.0,
        n_permutations=19,
        seed=1,
    ).metric
    result = fidelity.run(fidelity.make_inputs_from_norm(norm, {}))
    carved = recipes._fidelity_windows(result, 10.0, 1.0, True, list(DETUNING_LADDER))
    infidelity = tellability_recipe._infidelity_tellability(
        result,
        carved,
        thresholds=list(DETUNING_LADDER),
        gap_mult=10.0,
        n_permutations=19,
        seed=1,
    ).metric
    for group in (
        "spectrum",
        "noise_floor",
        "second_difference",
        "g10",
        "g10_halves",
        "self_consistency",
        "sigma_bound",
    ):
        assert frequency.status[group] == infidelity.status[group] == "computed", group
    for name in (
        "g10",
        "g10_first_half",
        "g10_second_half",
        "second_difference_ratio",
        "edge_noise_ratio",
        "self_consistency",
        "sigma_bound_ratio",
    ):
        assert getattr(infidelity, name) == pytest.approx(
            getattr(frequency, name), rel=1e-9
        ), name
    np.testing.assert_allclose(
        infidelity.spectrum_per_hz, frequency.spectrum_per_hz / omega**2, rtol=1e-9
    )
    for name in ("edge_level_per_hz", "noise_floor_per_hz"):
        assert getattr(infidelity, name) == pytest.approx(
            getattr(frequency, name) / omega**2, rel=1e-9
        ), name
