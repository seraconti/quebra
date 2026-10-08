"""Regression test for the WithinCalibration panel split (builder vs renderer).

The output-builder must produce a COMPLETE typed artifact (all derived fields
populated) and the renderer must be a pure function of it. These tests pin the
numeric invariants and the damage-seam behavior; the array-equality-vs-pre-split
parity was verified separately against golden baselines captured on the known-good
commit.
"""

from __future__ import annotations

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from quebra.analyzers import event_table, kaplan_meier, windows
from quebra.analyzers.within_calibration_compute import (
    _threshold_in_spec_frac,
    build_within_calibration_panel_data,
)
from quebra.panels.within_calibration import (
    WithinCalibrationPanel,
    WithinCalibrationPanelData,
)

pytestmark = pytest.mark.unit


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
    return result.windows_in_spec, result.reads


def _km_of(carved):
    """The job's one Kaplan-Meier node, built the way the recipe builds it."""
    return kaplan_meier.kaplan_meier_set(event_table.event_tables_from_carve(carved))


def _km_for(t_h, series, thresholds):
    return _km_of(
        windows.run(
            windows.WindowsInputs(
                t_rel_s=np.asarray(t_h, dtype=float) * 3600.0,
                values=np.asarray(series, dtype=float),
                thresholds=thresholds,
                dataset_id="unit",
            )
        )
    )


_THRESHOLDS = [("2 µs", 2.0, True), ("3 µs", 3.0, True), ("4 µs", 4.0, True)]


def _t2star_like() -> WithinCalibrationPanelData:
    t_h = np.linspace(0.0, 12.0, 500)
    s = 3.0 + 0.6 * np.sin(t_h) + 0.3 * np.cos(3.0 * t_h)
    return build_within_calibration_panel_data(
        t_h=t_h,
        primary_series=s,
        primary_label="T2* (µs)",
        thresholds=_THRESHOLDS,
        meta={"dataset": "unit-test"},
        windows=_carved(t_h, s, _THRESHOLDS)[0],
        reads=_carved(t_h, s, _THRESHOLDS)[1],
        gap_spans_s=[],  # uniform spacing: the carve records no gap
        kaplan_meier=_km_for(t_h, s, _THRESHOLDS),
    )


def test_builder_populates_all_derived_fields() -> None:
    d = _t2star_like()
    labels = {lbl for lbl, _, _ in d.thresholds}
    assert set(d.reliability.cumulative_time_per_threshold) == labels
    assert set(d.reliability.cumulative_damage_per_threshold) == labels
    assert set(d.reliability.ttf_per_threshold) == labels
    assert {key for key in d.reliability.kaplan_meier.curves} == {
        (label, side) for label in labels for side in event_table.SIDES
    }
    assert set(d.reliability.occupancy) == labels
    assert all(0.0 <= d.reliability.occupancy[label] <= 1.0 for label in labels)
    assert np.isfinite(d.signal.cv)
    # band 2 and band 3 each carry a per-threshold entry for every rung
    assert set(d.distinguish.state_series_per_threshold) == labels
    assert set(d.distinguish.state_counts_per_threshold) == labels
    assert set(d.reliability.compliance_state_series) == labels


def test_cumulative_time_is_monotonic_and_bounded() -> None:
    d = _t2star_like()
    total_h = float(d.signal.t_h[-1] - d.signal.t_h[0])
    for arr in d.reliability.cumulative_time_per_threshold.values():
        assert len(arr) == len(d.signal.t_h)
        assert np.all(np.diff(arr) >= -1e-12)  # non-decreasing
        assert arr[0] == 0.0
        assert arr[-1] <= total_h + 1e-9


def test_occupancy_equals_a_hand_computed_fraction_of_observed_time() -> None:
    """Oracle: a four-interval record whose occupancy is 0.5 by inspection.

    Four one-hour intervals on a uniform grid. `_threshold_in_spec_frac` charges an
    interval to the reading at its LEFT endpoint (`oos[:-1]`), so with readings
    [10, 5, 1, 1, 10] against a threshold of 5 the out-of-spec intervals are the third
    and fourth, giving 2 of 4 observed hours and an occupancy of exactly 0.5.

    The reading of exactly 5.0 is the point of the case. `AGENTS.md` section 5 fixes
    in spec as `T2* >= threshold`, and `windows.in_spec_mask` applies it as
    `margin >= 0`, so a value AT the threshold is in spec and must not be charged.
    Flipping `windows.in_spec_mask`'s `>=` to `>` moves this to 0.25.
    """
    t_h = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    series = np.array([10.0, 5.0, 1.0, 1.0, 10.0])
    thresholds = [("5", 5.0, True)]
    _, reads = _carved(t_h, series, thresholds)
    frac = _threshold_in_spec_frac(t_h, series, reads, thresholds, [])
    assert frac["5"] == pytest.approx(0.5, abs=1e-12)


def test_default_damage_is_excess_integral_not_noop() -> None:
    # A series that spends time out of spec must accumulate non-zero default damage,
    # and a non-identity damage_fn must change the curve (the seam is live).
    t_h = np.linspace(0.0, 10.0, 400)
    s = np.linspace(1.0, 6.0, 400)  # rises through the T2* thresholds
    thr = [("3 µs", 3.0, True)]
    win, rd = _carved(t_h, s, thr)
    default = build_within_calibration_panel_data(
        t_h=t_h,
        primary_series=s,
        primary_label="T2* (µs)",
        thresholds=thr,
        meta={"dataset": "unit-test"},
        windows=win,
        reads=rd,
        gap_spans_s=[],
        kaplan_meier=_km_for(t_h, s, thr),
    )
    squared = build_within_calibration_panel_data(
        t_h=t_h,
        primary_series=s,
        primary_label="T2* (µs)",
        thresholds=thr,
        meta={"dataset": "unit-test"},
        windows=win,
        reads=rd,
        gap_spans_s=[],
        kaplan_meier=_km_for(t_h, s, thr),
        damage_fn=lambda x: x**2,
    )
    d_default = default.reliability.cumulative_damage_per_threshold["3 µs"]
    d_squared = squared.reliability.cumulative_damage_per_threshold["3 µs"]
    assert d_default[-1] > 0.0
    assert not np.allclose(d_default, d_squared)


def test_renderer_produces_figure_without_data_arithmetic() -> None:
    d = _t2star_like()
    fig = WithinCalibrationPanel(name="unit").build_matplotlib(d)
    assert len(fig.axes) >= 5
    matplotlib.pyplot.close(fig)


def test_the_caption_names_the_dataset_and_refuses_an_artifact_without_one() -> None:
    """Oracle: AGENTS.md section 3, errors are raised, not swallowed; the name set by hand.

    Both adapters put the record's name in `meta["dataset"]`, so an artifact without it is
    broken, and a caption that drops the name without a word is a wrong-but-plausible figure.
    """
    import copy

    import matplotlib.pyplot as plt

    pd_ = _t2star_like()
    broken = copy.deepcopy(pd_)
    del broken.meta["dataset"]
    # Built around the builder, which refuses an empty name: the renderer refuses it too.
    blank = copy.deepcopy(pd_)
    blank.meta["dataset"] = "  "
    fig = plt.figure()
    try:
        WithinCalibrationPanel._draw_caption(fig, pd_)
        caption = fig.get_suptitle()
        assert caption.startswith("unit-test  ·  T2* (µs): mean "), caption
        with pytest.raises(KeyError, match="dataset"):
            WithinCalibrationPanel._draw_caption(fig, broken)
        with pytest.raises(ValueError, match="caption must name the record"):
            WithinCalibrationPanel._draw_caption(fig, blank)
    finally:
        plt.close(fig)


def test_renderer_rejects_wrong_type() -> None:
    import pytest

    with pytest.raises(TypeError):
        WithinCalibrationPanel(name="unit").build_matplotlib(object())


def _with_read_gaps() -> WithinCalibrationPanelData:
    """Three observed stretches separated by two unobserved holes, 5-8 h and 12-15 h.

    The median read spacing is about 0.025 h, so a 3 h hole is far past any gap multiple the
    carve uses: it records both holes in `signal.gap_spans_h`, through the real analyzer, the
    same way a job does. Two holes, because a split that honoured only the first would pass a
    one-gap record.
    """
    t_h = np.concatenate(
        [
            np.linspace(0.0, 5.0, 200),
            np.linspace(8.0, 12.0, 160),
            np.linspace(15.0, 20.0, 200),
        ]
    )
    s = 3.0 + 0.6 * np.sin(t_h) + 0.3 * np.cos(3.0 * t_h)
    carved = windows.run(
        windows.WindowsInputs(
            t_rel_s=t_h * 3600.0, values=s, thresholds=_THRESHOLDS, dataset_id="unit"
        )
    )
    # From the carve's diagnostics, exactly as `recipes.py` passes it for every production
    # job: the window and read tables do not hold every gap.
    return build_within_calibration_panel_data(
        t_h=t_h,
        primary_series=s,
        primary_label="T2* (µs)",
        thresholds=_THRESHOLDS,
        meta={"dataset": "unit-test-gap"},
        windows=carved.windows_in_spec,
        reads=carved.reads,
        gap_spans_s=carved.diagnostics["gap_spans_s"],
        kaplan_meier=_km_of(carved),
    )


def test_the_read_gaps_are_excluded_from_observed_time_and_timelines() -> None:
    """Oracle: observed time summed by hand, 5 + 4 + 5 = 14 h, on the two-gap fixture.

    A missing gap list changes NUMBERS, not only ink. The fixture's series never reaches 4 µs,
    so it is out of spec at 4 µs for every observed hour: cumulative time out of spec must end
    at the observed 14 h, not the 20 h of wall clock. And every state timeline must mark the
    two holes, 5-8 h and 12-15 h, as unobserved rather than as a state.
    """
    from quebra.analyzers.windows import STATE_UNOBSERVED

    pd_ = _with_read_gaps()
    assert np.max(pd_.signal.values) < 4.0, (
        "the fixture reaches 4 µs; the oracle is void"
    )
    cum_h = pd_.reliability.cumulative_time_per_threshold["4 µs"]
    assert cum_h[-1] == pytest.approx(14.0, abs=1e-9), (
        f"cumulative time out of spec ends at {cum_h[-1]} h, not the observed 14 h"
    )
    for label, segments in pd_.distinguish.state_series_per_threshold.items():
        holes = [(lo, hi) for lo, hi, state in segments if state == STATE_UNOBSERVED]
        np.testing.assert_allclose(holes, [(5.0, 8.0), (12.0, 15.0)], err_msg=label)


@pytest.mark.parametrize(
    "method", ["_draw_primary", "_draw_cumulative_time", "_draw_cumulative_damage"]
)
def test_no_curve_on_the_scan_clock_is_drawn_across_a_read_gap(method):
    """Oracle: the gap spans the carve recorded, checked against every line actually drawn.

    docs/FIGURE_STANDARD.md:138: nothing is drawn across a gap; the trace breaks. Every curve
    on the scan clock is drawn one line per observed stretch. `_draw_primary` is the reference
    case: a failure there points at `observed_slices` or the fixture rather than at the
    cumulative curves.

    Asserted: no line crosses a gap; the gaps lie inside the drawn span, so gaps and lines in
    different units cannot make the crossing check vacuous; and every read is drawn. For the
    two cumulative drawers, also per threshold: its lines together cover every read, at that
    threshold's own values, and carry its label exactly once.

    Asserted against the drawn lines, not the artifact: `_observed_dt_h` zeroes gap intervals,
    so the cumulative arrays are already flat across a gap, and only the axes show ink over
    unobserved time.
    """
    import matplotlib.pyplot as plt

    from quebra.panels import _within_calibration_render as render

    pd_ = _with_read_gaps()
    gaps_h = pd_.signal.gap_spans_h
    assert len(gaps_h) == 2, f"the fixture should carve two read gaps, got {gaps_h}"
    slices = render.observed_slices(pd_)
    assert len(slices) == 3, (
        f"expected the gaps to split the record in three, got {slices}"
    )

    fig, ax = plt.subplots()
    try:
        panel = WithinCalibrationPanel(name="gap_probe")
        if method == "_draw_primary":
            getattr(panel, method)(ax, pd_, color="black")
        else:
            getattr(panel, method)(ax, pd_)
        # Data lines only: a threshold `axhline` lives in axes coordinates on x.
        lines = [ln for ln in ax.get_lines() if ln.get_transform() == ax.transData]
        assert lines, (
            f"{method} drew no line in data coordinates, so this asserts nothing"
        )
        xs = [np.asarray(ln.get_xdata(), dtype=float) for ln in lines]
        for x in xs:
            assert np.all(np.isfinite(x)), f"{method} drew a line with non-finite x"

        drawn = np.concatenate(xs)
        for before_h, after_h in gaps_h:
            assert drawn.min() <= before_h and drawn.max() >= after_h, (
                f"gap {before_h}-{after_h} h lies outside the drawn span "
                f"{drawn.min()}-{drawn.max()}: gaps and lines are not in the same unit"
            )
            for x in xs:
                crosses = x.min() <= before_h and x.max() >= after_h
                assert not crosses, (
                    f"{method} drew one line from {x.min():.2f} h to {x.max():.2f} h, "
                    f"straight across the read gap {before_h:.2f}-{after_h:.2f} h"
                )

        t_h = np.asarray(pd_.signal.t_h, dtype=float)
        undrawn = t_h[~np.isin(t_h, drawn)]
        assert undrawn.size == 0, (
            f"{method} left {undrawn.size} observed reads undrawn, first at {undrawn[0]:.2f} h"
        )

        labels = [ln.get_label() for ln in lines if not ln.get_label().startswith("_")]
        assert len(labels) == len(set(labels)), (
            f"{method} repeats a legend label: {labels}"
        )

        if method != "_draw_primary":
            _assert_each_threshold_drawn_whole(method, pd_, lines, t_h)
    finally:
        plt.close(fig)


def _assert_each_threshold_drawn_whole(method, pd_, lines, t_h) -> None:
    """Each threshold's lines, found by its colour, cover every read at its own values."""
    from matplotlib.colors import to_rgba

    from quebra.plots import theme

    per_threshold = (
        pd_.reliability.cumulative_time_per_threshold
        if method == "_draw_cumulative_time"
        else pd_.reliability.cumulative_damage_per_threshold
    )
    n = len(pd_.thresholds)
    colours = [to_rgba(theme.threshold_color(i, n)) for i in range(n)]
    assert len(set(colours)) == n, (
        "thresholds share a colour, so lines cannot be told apart"
    )
    for i, (label, _, _) in enumerate(pd_.thresholds):
        mine = [ln for ln in lines if to_rgba(ln.get_color()) == colours[i]]
        x = np.concatenate([np.asarray(ln.get_xdata(), dtype=float) for ln in mine])
        y = np.concatenate([np.asarray(ln.get_ydata(), dtype=float) for ln in mine])
        np.testing.assert_array_equal(
            np.sort(x), t_h, err_msg=f"{method} {label}: reads missing or drawn twice"
        )
        expected = np.asarray(per_threshold[label], dtype=float)[
            np.searchsorted(t_h, x)
        ]
        np.testing.assert_array_equal(
            y, expected, err_msg=f"{method} {label}: wrong values"
        )
        shown = [ln.get_label() for ln in mine if not ln.get_label().startswith("_")]
        assert shown == [label], f"{method} {label}: legend entries {shown}"


def test_a_missing_gap_list_is_refused_not_read_as_no_gaps() -> None:
    """Oracle: AGENTS.md section 3, errors are raised, checked at every layer the list crosses.

    A record read as gap-free draws through its holes and counts unobserved hours as observed,
    and nothing downstream can tell. So each function that hands the carve's gap list on takes
    it as a required argument, and the builder refuses None, which is what a soft
    `diagnostics.get("gap_spans_s")` would hand it.
    """
    import dataclasses
    import inspect
    import types

    import pandas as pd

    from quebra.analyzers import fidelity, t2star
    from quebra.analyzers import within_calibration_compute as compute
    from quebra.analyzers.signal_band import SignalBand

    for fn, param in (
        (build_within_calibration_panel_data, "gap_spans_s"),
        (t2star.make_panel_data, "gap_spans_s"),
        (fidelity.make_panel_data, "gap_spans_s"),
        (compute._cumulative_time_out_of_spec, "gap_spans_h"),
        (compute._cumulative_damage, "gap_spans_h"),
        (compute._threshold_in_spec_frac, "gap_spans_h"),
    ):
        default = inspect.signature(fn).parameters[param].default
        assert default is inspect.Parameter.empty, (
            f"{fn.__module__}.{fn.__name__} defaults {param} to {default!r}"
        )
    (gap_field,) = [
        f for f in dataclasses.fields(SignalBand) if f.name == "gap_spans_h"
    ]
    assert gap_field.default is dataclasses.MISSING, "SignalBand defaults its gap list"
    assert gap_field.default_factory is dataclasses.MISSING, (
        "SignalBand defaults its gap list"
    )

    # Each adapter must hand None through for the builder to refuse, not soften it to [].
    t_rel_s = np.arange(10, dtype=float) * 60.0
    frame = pd.DataFrame(
        {
            "t_rel_s": t_rel_s,
            "t2star_s": np.full(10, 3e-6),
            "infidelity": np.full(10, 1e-3),
        }
    )
    stand_in = types.SimpleNamespace(frame=frame, meta={"dataset_id": "unit-test"})
    for adapter in (t2star.make_panel_data, fidelity.make_panel_data):
        with pytest.raises(ValueError, match="gap_spans_s is required"):
            adapter(
                stand_in,  # type: ignore[arg-type]
                windows=pd.DataFrame(),
                reads=pd.DataFrame(),
                gap_spans_s=None,  # type: ignore[arg-type]
                kaplan_meier=None,  # type: ignore[arg-type]
            )

    t_h = np.concatenate([np.linspace(0.0, 5.0, 200), np.linspace(8.0, 12.0, 160)])
    s = 3.0 + 0.6 * np.sin(t_h)
    win, rd = _carved(t_h, s, _THRESHOLDS)
    with pytest.raises(ValueError, match="gap_spans_s is required"):
        build_within_calibration_panel_data(
            t_h=t_h,
            primary_series=s,
            primary_label="T2* (µs)",
            thresholds=_THRESHOLDS,
            meta={"dataset": "unit-test"},
            windows=win,
            reads=rd,
            gap_spans_s=None,  # type: ignore[arg-type]
            kaplan_meier=_km_for(t_h, s, _THRESHOLDS),
        )


def test_a_cumulative_curve_longer_than_the_scan_clock_is_refused() -> None:
    """Oracle: artifacts whose last cumulative array is one point too long, or empty, against the reads.

    Drawing one line per observed stretch slices each array by read index, which would drop
    an extra point without a word, and an empty array would be skipped as "no data". So both
    drawers check every threshold's length before anything else. Only the LAST threshold is
    corrupted, so a check that stops after the first catches nothing.
    """
    import copy

    import matplotlib.pyplot as plt

    panel = WithinCalibrationPanel(name="length_probe")
    for corrupt in (lambda a: np.append(a, 0.0), lambda a: a[:0]):
        pd_ = copy.deepcopy(_t2star_like())
        for per_threshold in (
            pd_.reliability.cumulative_time_per_threshold,
            pd_.reliability.cumulative_damage_per_threshold,
        ):
            last = list(per_threshold)[-1]
            per_threshold[last] = corrupt(np.asarray(per_threshold[last]))
        for method in ("_draw_cumulative_time", "_draw_cumulative_damage"):
            fig, ax = plt.subplots()
            try:
                with pytest.raises(ValueError, match="points against"):
                    getattr(panel, method)(ax, pd_)
            finally:
                plt.close(fig)


def test_the_timeline_cull_counts_an_unmeasured_threshold_apart_from_one_below_5pct() -> (
    None
):
    """Oracle: occupancy values set by hand on the artifact, one per threshold.

    `_threshold_in_spec_frac` gives NaN where there is no observed time to measure on, so
    that threshold was never measured below 5% in spec. The fixture's 3 µs is drawn and its
    4 µs is measured below 5%; 2 µs is set to NaN. The cull drops both, and the title must
    name them apart. With every threshold NaN the empty axes must not claim a measurement.
    A missing occupancy key is a broken artifact (`check_thresholds` requires every label),
    so it raises rather than being read as 0.
    """
    import copy

    import matplotlib.pyplot as plt

    base = _t2star_like()
    occupancy = base.reliability.occupancy
    assert occupancy["3 µs"] >= 0.05 > occupancy["4 µs"], f"premise: {occupancy}"
    panel = WithinCalibrationPanel(name="cull_probe")

    def draw(pd_):
        fig, ax = plt.subplots()
        try:
            panel._draw_threshold_timeline(
                ax, pd_, pd_.reliability.compliance_state_series, title="T"
            )
            ticks = [tick.get_text() for tick in ax.get_yticklabels()]
            return ax.get_title(), ticks, [text.get_text() for text in ax.texts]
        finally:
            plt.close(fig)

    one_unmeasured = copy.deepcopy(base)
    one_unmeasured.reliability.occupancy["2 µs"] = float("nan")
    title, ticks, _ = draw(one_unmeasured)
    assert title == "T (1 thresholds below 5% in spec, 1 not measured)", title
    assert ticks == ["3 µs"], ticks

    none_measured = copy.deepcopy(base)
    for label in none_measured.reliability.occupancy:
        none_measured.reliability.occupancy[label] = float("nan")
    title, _, notes = draw(none_measured)
    assert title == "T (3 not measured)", title
    assert notes == ["In-spec time not measured: no observed time"], notes

    missing = copy.deepcopy(base)
    del missing.reliability.occupancy["3 µs"]
    with pytest.raises(KeyError, match="3 µs"):
        draw(missing)


def test_the_survival_panel_states_what_kaplan_meier_counted_and_left_out(
    capsys,
) -> None:
    """Oracle: counts per threshold from the carve's own window table.

    docs/FIGURE_STANDARD.md: a panel that drops data says how much it dropped, in the panel,
    not in the log. Kaplan-Meier keeps censored windows and leaves out windows whose birth
    was not observed, so every threshold that carved a window has a legend entry with its
    windows, deaths, censorings and left-out count, including one whose every window had an
    unobserved birth and which therefore has no curve (its entry draws no line sample); and
    the reliability band step prints nothing.
    """
    import matplotlib.pyplot as plt

    t_h = np.concatenate(
        [
            np.linspace(0.0, 5.0, 200),
            np.linspace(8.0, 12.0, 160),
            np.linspace(15.0, 20.0, 200),
        ]
    )
    s = 3.0 + 0.6 * np.sin(t_h) + 0.3 * np.cos(3.0 * t_h)
    carved = windows.run(
        windows.WindowsInputs(
            t_rel_s=t_h * 3600.0, values=s, thresholds=_THRESHOLDS, dataset_id="unit"
        )
    )
    pd_ = build_within_calibration_panel_data(
        t_h=t_h,
        primary_series=s,
        primary_label="T2* (µs)",
        thresholds=_THRESHOLDS,
        meta={"dataset": "unit-test"},
        windows=carved.windows_in_spec,
        reads=carved.reads,
        gap_spans_s=carved.diagnostics["gap_spans_s"],
        kaplan_meier=_km_of(carved),
    )
    from quebra.analyzers import reliability_band

    capsys.readouterr()
    reliability_band.run(
        reliability_band.make_inputs_from_windows(
            t_h=t_h,
            values=s,
            reads=carved.reads,
            windows=carved.windows_in_spec,
            thresholds=_THRESHOLDS,
            gap_spans_h=pd_.signal.gap_spans_h,
            kaplan_meier=_km_of(carved),
            damage_fn=None,
        )
    )
    printed = capsys.readouterr()
    assert printed.out == "" and printed.err == "", f"the band step printed: {printed}"

    fig, ax = plt.subplots()
    try:
        WithinCalibrationPanel(name="survival_probe")._draw_survival(ax, pd_)
        entries = ax.get_legend()
        legend = [text.get_text() for text in entries.get_texts()]
        samples = {
            text.get_text(): handle.get_linestyle()
            for text, handle in zip(entries.get_texts(), entries.legend_handles)
        }
        notes = [text.get_text() for text in ax.texts]
    finally:
        plt.close(fig)
    assert "independence checks: NOT ASSESSED for this band" in notes, notes
    no_curve = [text for text in samples if text.endswith(", no curve")]
    assert no_curve, "no no-curve entry to check"
    for text in no_curve:
        assert samples[text] == "None", f"{text!r} draws a line sample"

    table = carved.windows_in_spec
    with_curve = without_curve = 0
    for label, _, _ in _THRESHOLDS:
        mine = table[table["threshold_label"] == label]
        if len(mine) == 0:
            assert not any(t.startswith(label) for t in legend), (label, legend)
            continue
        seen = mine[mine["birth_observed"]]
        n_censored = int(seen["censored"].sum())
        expected = (
            f"{label}: n={len(seen)}, d={len(seen) - n_censored}, c={n_censored}, "
            f"u={len(mine) - len(seen)}"
        )
        if len(seen) == 0:
            expected += ", no curve"
            without_curve += 1
        else:
            with_curve += 1
        assert expected in legend, (expected, legend)
    assert with_curve and without_curve, (
        "the fixture no longer has both kinds of threshold"
    )

    # Every window with an unobserved birth: no curve at all, and the counts still show.
    only = [thr for thr in _THRESHOLDS if thr[0] == "2 µs"]
    carved_only = windows.run(
        windows.WindowsInputs(
            t_rel_s=t_h * 3600.0, values=s, thresholds=only, dataset_id="unit"
        )
    )
    pd_only = build_within_calibration_panel_data(
        t_h=t_h,
        primary_series=s,
        primary_label="T2* (µs)",
        thresholds=only,
        meta={"dataset": "unit-test"},
        windows=carved_only.windows_in_spec,
        reads=carved_only.reads,
        gap_spans_s=carved_only.diagnostics["gap_spans_s"],
        kaplan_meier=_km_of(carved_only),
    )
    fig, ax = plt.subplots()
    try:
        WithinCalibrationPanel(name="survival_probe")._draw_survival(ax, pd_only)
        assert ax.get_legend() is not None, "no legend when no curve is drawn"
        legend = [text.get_text() for text in ax.get_legend().get_texts()]
        notes = [text.get_text() for text in ax.texts]
        assert notes == [
            "No in-spec windows with an observed birth for defined thresholds"
        ], notes
    finally:
        plt.close(fig)
    n_only = len(carved_only.windows_in_spec)
    assert n_only, "the single-threshold carve produced no window"
    assert legend == [f"2 µs: n=0, d=0, c=0, u={n_only}, no curve"], legend


def test_the_survival_line_and_band_cover_the_curve_out_to_its_longest_window() -> None:
    """Oracle: durations set by hand, checked against the curve the artifact carries.

    `KaplanMeierCurve` holds S(t) as the left ends of its segments, defined out to
    `max_observed_min`, with a log-log band wherever both bounds are finite. So the drawn
    line must end at `max_observed_min` at the curve's last S, and the band must fill the
    middle of every segment whose bounds are finite and of no other. One curve per threshold:
    a death at 1 then windows censored at 5 and 10 (the tail [1, 10) at S = 2/3); no death,
    censored at 4 and 8 (a flat line, no band); deaths at 1, 2 and 3 (the segment [2, 3) at
    S = 1/3 ends where S reaches 0, so its right end is NaN and its left value is not).
    """
    import copy

    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba

    from quebra.plots import theme

    by_threshold = {
        "2 µs": ([1.0, 5.0, 10.0], [True, False, False]),
        "3 µs": ([4.0, 8.0], [False, False]),
        "4 µs": ([1.0, 2.0, 3.0], [True, True, True]),
    }
    pd_ = copy.deepcopy(_t2star_like())
    km = pd_.reliability.kaplan_meier
    for label, (duration_min, death_observed) in by_threshold.items():
        km.curves[(label, windows.SIDE_IN_SPEC)] = kaplan_meier.run(
            kaplan_meier.KaplanMeierInputs(
                duration_min=np.array(duration_min),
                death_observed=np.array(death_observed),
                label=label,
                threshold_label=label,
                side=windows.SIDE_IN_SPEC,
                n_windows_carved=len(duration_min),
            )
        )
    last_before_zero = km.curve("4 µs", windows.SIDE_IN_SPEC)
    assert np.isfinite(last_before_zero.band_lower[-2]), "premise: [2, 3) has a band"
    assert np.isnan(last_before_zero.band_lower[-1]), "premise: S = 0 has no band"

    fig, ax = plt.subplots()
    try:
        WithinCalibrationPanel(name="survival_probe")._draw_survival(ax, pd_)
        n = len(pd_.thresholds)
        for i, (label, _, _) in enumerate(pd_.thresholds):
            curve = km.curve(label, windows.SIDE_IN_SPEC)
            colour = to_rgba(theme.threshold_color(i, n))
            (line,) = [ln for ln in ax.get_lines() if to_rgba(ln.get_color()) == colour]
            x = np.asarray(line.get_xdata(), dtype=float)
            y = np.asarray(line.get_ydata(), dtype=float)
            assert x[0] == 0.0 and x[-1] == curve.max_observed_min, (label, x)
            assert y[-1] == curve.survival[-1], (label, y)

            paths = [
                path
                for band in ax.collections
                if np.allclose(band.get_facecolor()[0][:3], colour[:3])
                for path in band.get_paths()
            ]
            ends = np.append(curve.time_min[1:], curve.max_observed_min)
            for left, right, lo, hi in zip(
                curve.time_min, ends, curve.band_lower, curve.band_upper
            ):
                if right <= left:
                    continue
                mid = 0.5 * (left + right)
                if np.isfinite(lo) and np.isfinite(hi):
                    point = (mid, 0.5 * (lo + hi))
                    assert any(path.contains_point(point) for path in paths), (
                        f"{label}: no band over [{left}, {right})"
                    )
                else:
                    filled = [
                        y_probe
                        for y_probe in np.linspace(0.0, 1.0, 51)
                        if any(path.contains_point((mid, y_probe)) for path in paths)
                    ]
                    assert not filled, f"{label}: band over [{left}, {right}), {filled}"
    finally:
        plt.close(fig)
