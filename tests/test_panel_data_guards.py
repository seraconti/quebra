"""Panel data refuses inputs from another carve and states NaN for what it cannot measure.

Every case carves one record through `windows.run` and estimates the Kaplan-Meier set
from that carve, as the recipe does, then hands the builder or one band a single input
that does not belong to it. Oracle: construction. The test makes the mismatch, so a raise
naming it is the expected outcome and a silent result is the defect.

The record is T2*-like: carved in seconds, drawn in µs, so every accepted case runs at the
display scale the T2* adapter uses.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from quebra import recipes
from quebra.analyzers import (
    distinguish_band,
    event_table,
    fidelity,
    kaplan_meier,
    reliability_band,
    t2star,
    windows,
)
from quebra.analyzers.within_calibration_compute import (
    _threshold_in_spec_frac,
    build_within_calibration_panel_data,
)

pytestmark = pytest.mark.unit

_N = 400
_T_H = np.arange(_N) * 0.25 + 20000.0
_SCALE = 1e6
# Rung "z" sits above every read, so it has no in-spec window at all.
_LADDER_S = [("x", 25e-6, True), ("y", 30e-6, True), ("z", 100e-6, True)]


def _series_s(phase: float = 0.0) -> np.ndarray:
    rng = np.random.default_rng(1)
    s = 25e-6 + 8e-6 * np.sin(np.arange(_N) / 9.0 + phase) + rng.normal(0, 1e-6, _N)
    s[[5, 77]] = np.nan  # failed fits: the carve and the panel drop the same reads
    return s


def _carve(series_s, ladder=_LADDER_S, dataset_id="A"):
    return windows.run(
        windows.WindowsInputs(
            t_rel_s=_T_H * 3600.0,
            values=series_s,
            thresholds=ladder,
            dataset_id=dataset_id,
        )
    )


def _km(carved):
    return kaplan_meier.kaplan_meier_set(event_table.event_tables_from_carve(carved))


def _display(ladder):
    return [(label, value * _SCALE, good) for label, value, good in ladder]


def _build(carved, *, thresholds, series, km=None):
    return build_within_calibration_panel_data(
        t_h=_T_H,
        primary_series=series,
        primary_label="T2* (µs)",
        thresholds=thresholds,
        meta={"dataset": "A"},
        windows=carved.windows_in_spec,
        reads=carved.reads,
        gap_spans_s=carved.diagnostics["gap_spans_s"],
        kaplan_meier=_km(carved) if km is None else km,
    )


def test_a_carve_in_seconds_drawn_in_microseconds_is_accepted() -> None:
    """Oracle: the T2* adapter's scale. 1e6 on the series and the ladder is one record.

    Positive control for every refusal below: the same carve builds drawn in seconds and
    drawn in µs, and the in-spec fractions agree, since classification is the carve's.
    """
    s = _series_s()
    carved = _carve(s)
    assert not (carved.windows_in_spec["threshold_label"] == "z").any()
    in_s = _build(carved, thresholds=_LADDER_S, series=s)
    in_us = _build(carved, thresholds=_display(_LADDER_S), series=s * _SCALE)
    assert in_us.reliability.occupancy == in_s.reliability.occupancy


@pytest.mark.parametrize(
    ("carve_ladder", "panel_ladder", "series_scale", "named"),
    [
        # one rung: only the series can pin the scale, a ladder of one cannot
        ([("x", 25e-6, True)], [("x", 30.0, True)], _SCALE, ["threshold 'x'"]),
        (
            _LADDER_S,
            [("x", 25.0, True), ("y", 30.0, False), ("z", 100.0, True)],
            _SCALE,
            ["threshold 'y'"],
        ),
        # a rung with no in-spec window is checked too
        (
            _LADDER_S,
            [("x", 25.0, True), ("y", 30.0, True), ("z", 90.0, True)],
            _SCALE,
            ["threshold 'z'"],
        ),
        # one scale on the ladder, another on the series
        (
            _LADDER_S,
            _display(_LADDER_S),
            1e3,
            ["threshold 'x'", "threshold 'y'", "threshold 'z'"],
        ),
    ],
    ids=["one_rung_value", "direction", "rung_without_windows", "two_scales"],
)
def test_a_rung_drawn_at_another_value_or_direction_is_refused(
    carve_ladder, panel_ladder, series_scale, named
) -> None:
    """Oracle: construction. The panel draws a rung the carve did not classify against."""
    s = _series_s()
    carved = _carve(s, ladder=carve_ladder)
    with pytest.raises(ValueError, match="the carve and the panel disagree") as err:
        _build(carved, thresholds=panel_ladder, series=s * series_scale)
    assert f"disagree on {named!r}" in str(err.value)


def test_a_negative_display_scale_is_refused() -> None:
    """Oracle: construction. Negated series, flipped ladder: every margin but its sign."""
    s = _series_s()
    carved = _carve(s)
    flipped = [(label, -value * _SCALE, not good) for label, value, good in _LADDER_S]
    with pytest.raises(ValueError, match="no positive multiple"):
        _build(carved, thresholds=flipped, series=-s * _SCALE)


@pytest.mark.parametrize(
    ("other", "message"),
    [
        (lambda: _carve(_series_s(phase=2.0), dataset_id="B"), r"record\(s\) \['B'\]"),
        (lambda: _carve(_series_s(), ladder=_LADDER_S[:2]), r"lacks threshold\(s\)"),
        (lambda: _carve(_series_s(phase=2.0)), "the in-spec table holds"),
    ],
    ids=["another_record", "another_ladder", "same_id_other_windows"],
)
def test_a_kaplan_meier_set_from_another_carve_is_refused(other, message) -> None:
    """Oracle: construction. The drawn curve must describe the windows the band counts."""
    s = _series_s()
    carved = _carve(s)
    with pytest.raises(ValueError, match=message):
        _build(
            carved,
            thresholds=_display(_LADDER_S),
            series=s * _SCALE,
            km=_km(other()),
        )


def test_a_window_table_from_another_record_is_refused() -> None:
    """Oracle: construction. Reads of record A beside windows relabelled as record C."""
    s = _series_s()
    carved = _carve(s)
    carved.windows_in_spec = carved.windows_in_spec.assign(dataset_id="C")
    with pytest.raises(ValueError, match=r"mix records \['A', 'C'\]"):
        _build(carved, thresholds=_display(_LADDER_S), series=s * _SCALE)


@pytest.mark.parametrize("phase", [4.0, 6.0])
def test_a_window_table_from_another_carve_of_the_same_record_is_refused(phase) -> None:
    """Oracle: construction. Reads of one carve beside windows and Kaplan-Meier set of
    another carve under the same record id, which every record-level check accepts."""
    s = _series_s()
    carved = _carve(s)
    other = _carve(_series_s(phase=phase))
    carved.windows_in_spec = other.windows_in_spec
    with pytest.raises(ValueError, match="are not one carve"):
        _build(carved, thresholds=_display(_LADDER_S), series=s * _SCALE, km=_km(other))


def test_a_window_table_with_the_same_windows_at_other_reads_is_refused() -> None:
    """Oracle: construction. Every key and the window count agree; one window starts one
    read late, so only the per-window comparison of first read can see it."""
    s = _series_s()
    carved = _carve(s)
    shifted = carved.windows_in_spec.copy()
    shifted.loc[shifted.index[0], "t_birth_s"] += 900.0  # one read spacing
    carved.windows_in_spec = shifted
    with pytest.raises(ValueError, match="are not one carve"):
        _build(carved, thresholds=_display(_LADDER_S), series=s * _SCALE)


def test_a_window_table_carved_under_another_gap_policy_is_refused() -> None:
    """Oracle: construction. The same reads carved twice: with a gap (gap_mult 10) and
    without one (gap_mult 1e6). Both tile the reads alike, so first read, last read and
    read count agree; the gap ends a window as censored in one carve and not the other."""
    t_h = _T_H.copy()
    # A 10 h hole (40 read spacings) right after the last read of an in-spec window at
    # rung x: one carve censors that window at the gap, the other sees it die at the next
    # read, and no window is split, so the tiling is the same.
    t_h[199:] += 10.0
    s = _series_s()

    def carve(gap_mult):
        return windows.run(
            windows.WindowsInputs(
                t_rel_s=t_h * 3600.0,
                values=s,
                thresholds=_LADDER_S,
                dataset_id="A",
                gap_mult=gap_mult,
            )
        )

    with_gap, without_gap = carve(10.0), carve(1e6)
    assert with_gap.diagnostics["gap_spans_s"], "the fixture has no gap"
    tiling = ["threshold_label", "window_index", "t_birth_s", "t_last_s", "n_reads"]
    assert with_gap.windows_in_spec[tiling].equals(
        without_gap.windows_in_spec[tiling]
    ), "the fixture's gap splits a window, so it no longer isolates the death"
    with pytest.raises(ValueError, match="are not one carve"):
        build_within_calibration_panel_data(
            t_h=t_h,
            primary_series=s * _SCALE,
            primary_label="T2* (µs)",
            thresholds=_display(_LADDER_S),
            meta={"dataset": "A"},
            windows=without_gap.windows_in_spec,
            reads=with_gap.reads,
            gap_spans_s=with_gap.diagnostics["gap_spans_s"],
            kaplan_meier=_km(without_gap),
        )


def test_both_bands_refuse_the_out_of_spec_table() -> None:
    """Oracle: construction. The out-of-spec table, which shares the in-spec schema.

    Each band refuses it itself, so the test calls each band directly.
    """
    s = _series_s()
    carved = _carve(s)
    out = carved.windows_out_of_spec
    assert len(out), "the fixture has no out-of-spec window"
    with pytest.raises(ValueError, match="expected the in-spec window table"):
        distinguish_band.run(
            distinguish_band.make_inputs_from_windows(
                reads=carved.reads,
                windows=out,
                thresholds=_display(_LADDER_S),
                median_read_spacing_s=900.0,
                gap_spans_h=[],
                sigma_display=None,
            )
        )
    with pytest.raises(ValueError, match="expected the in-spec window table"):
        reliability_band.run(
            reliability_band.make_inputs_from_windows(
                t_h=_T_H,
                values=s * _SCALE,
                reads=carved.reads,
                windows=out,
                thresholds=_display(_LADDER_S),
                gap_spans_h=[],
                kaplan_meier=_km(carved),
            )
        )


def test_the_in_spec_fraction_is_nan_where_no_time_was_observed() -> None:
    """Oracle: the definition. In-spec over observed time has no value with none observed.

    Two ways to have none: every interval inside a read gap, and a record with one finite
    read, through the whole builder. A 0.0 here would read as "never in spec".
    """
    t_h = np.arange(3.0)
    series = np.array([5.0, 3.0, 5.0])
    ladder = [("4", 4.0, True)]
    reads = windows.run(
        windows.WindowsInputs(t_rel_s=t_h * 3600.0, values=series, thresholds=ladder)
    ).reads
    observed = _threshold_in_spec_frac(t_h, series, reads, ladder, [])
    assert observed["4"] == pytest.approx(0.5, abs=1e-12)
    all_gap = _threshold_in_spec_frac(
        t_h, series, reads, ladder, [(0.0, 1.0), (1.0, 2.0)]
    )
    assert np.isnan(all_gap["4"])

    one_read = np.array([np.nan, 5.0, np.nan])
    carved = windows.run(
        windows.WindowsInputs(t_rel_s=t_h * 3600.0, values=one_read, thresholds=ladder)
    )
    data = build_within_calibration_panel_data(
        t_h=t_h,
        primary_series=one_read,
        primary_label="metric",
        thresholds=ladder,
        meta={"dataset": "A"},
        windows=carved.windows_in_spec,
        reads=carved.reads,
        gap_spans_s=carved.diagnostics["gap_spans_s"],
        kaplan_meier=_km(carved),
    )
    assert np.isnan(data.reliability.occupancy["4"])


@pytest.mark.parametrize(
    ("values_s", "mean_us", "std_us"),
    [
        ([np.nan, np.nan], np.nan, np.nan),
        ([2e-6, np.nan], 2.0, np.nan),
        ([1e-6, 3e-6], 2.0, np.sqrt(2.0)),
    ],
    ids=["no_valid_read", "one_valid_read", "two_valid_reads"],
)
def test_the_t2star_summary_is_nan_where_undefined(values_s, mean_us, std_us) -> None:
    """Oracle: analytic. Sample std of 1 and 3 µs is sqrt(2) µs; undefined below two reads.

    The mean of no read is undefined too, and is stated without numpy's empty-slice warning.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        result = t2star.run(
            t2star.T2StarInputs(
                t_rel_s=np.array([0.0, 60.0]),
                t2star_s=np.asarray(values_s, dtype=float),
            )
        )
    for key, expected in (("mean_us", mean_us), ("std_us", std_us)):
        got = result.diagnostics[key]
        if np.isnan(expected):
            assert np.isnan(got), (key, got)
        else:
            assert got == pytest.approx(expected, rel=1e-12), (key, got)


@pytest.mark.parametrize(
    "meta",
    [{}, {"dataset": ""}, {"dataset": "  "}, {"dataset": None}],
    ids=["missing", "empty", "blank", "none"],
)
def test_panel_data_without_a_dataset_name_is_refused(meta) -> None:
    """Oracle: construction. The caption names the record it draws, so an artifact that
    cannot name one is refused before anything is computed."""
    s = _series_s()
    carved = _carve(s)
    with pytest.raises(ValueError, match="must name the record"):
        build_within_calibration_panel_data(
            t_h=_T_H,
            primary_series=s * _SCALE,
            primary_label="T2* (µs)",
            thresholds=_display(_LADDER_S),
            meta=meta,
            windows=carved.windows_in_spec,
            reads=carved.reads,
            gap_spans_s=carved.diagnostics["gap_spans_s"],
            kaplan_meier=_km(carved),
        )


_N_NORM = 20
_T2STAR_NORM = {
    "t_rel_s": np.arange(_N_NORM, dtype=float),
    "T2star_s": np.full(_N_NORM, 3e-6),
}
_FIDELITY_NORM = {
    "t_rel_s": np.arange(_N_NORM, dtype=float) * 60.0,
    "delta_hz": np.full(_N_NORM, 1e3),
    "rabi_hz": np.full(_N_NORM, 1e6),
}
_ADAPTERS = {
    "t2star": lambda norm: t2star.run(t2star.make_inputs_from_norm(norm)),
    "fidelity": lambda norm: fidelity.run(fidelity.make_inputs_from_norm(norm, {})),
}


@pytest.mark.parametrize("adapter", ["t2star", "fidelity"])
@pytest.mark.parametrize(
    "meta", [None, {}, {"dataset_id": ""}], ids=["no_meta", "no_key", "empty"]
)
def test_both_adapters_refuse_a_record_without_a_dataset_id(adapter, meta) -> None:
    """Oracle: specification. A missing id is refused, never drawn as an unnamed caption."""
    norm = dict(_T2STAR_NORM if adapter == "t2star" else _FIDELITY_NORM)
    if meta is not None:
        norm["meta"] = meta
    with pytest.raises(KeyError, match="non-empty meta\\['dataset_id'\\]"):
        _ADAPTERS[adapter](norm)


def _panel_through_the_recipe(adapter: str, norm: dict) -> object:
    """The adapter's panel data built by the steps the recipe wires, in order."""
    if adapter == "t2star":
        result = recipes._t2star_run(norm)
        ladder = [("x", 3e-6, True)]
        carved = recipes._windows_run(result, 10.0, 1.0, False, ladder)
        km = recipes._kaplan_meier_set(recipes._event_tables(carved), 0.95)
        return recipes._t2star_panel_data(result, carved, km, 5, False, 0, ladder)
    profile = "longrun" if adapter == "fidelity_longrun" else ""
    config = {"dataset_profile": profile}
    result = fidelity.run(fidelity.make_inputs_from_norm(norm, config))
    carved = recipes._fidelity_windows(result, 10.0)
    km = recipes._kaplan_meier_set(recipes._event_tables(carved), 0.95)
    return recipes._fidelity_panel_data(result, carved, km, 0)


_WAVE = np.sin(np.arange(_N_NORM) / 3.0)


@pytest.mark.parametrize("adapter", ["t2star", "fidelity", "fidelity_longrun"])
def test_the_dataset_id_reaches_the_caption_through_every_adapter(adapter) -> None:
    """Oracle: the norm. The id the norm names is the one the panel's caption shows, on
    both fidelity return paths (`longrun` and the default) and on T2*."""
    if adapter == "t2star":
        norm = {"t_rel_s": _T2STAR_NORM["t_rel_s"], "T2star_s": 3e-6 + 1e-6 * _WAVE}
    else:
        norm = {**_FIDELITY_NORM, "delta_hz": 2e4 * _WAVE}
        norm["raw_frequency_hz"] = 5e9 + 2e4 * _WAVE
    norm["meta"] = {"dataset_id": "Q7"}
    assert _panel_through_the_recipe(adapter, norm).meta["dataset"] == "Q7"
