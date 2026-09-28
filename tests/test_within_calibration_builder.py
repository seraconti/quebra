"""Regression gate for the WithinCalibration panel split (builder vs renderer).

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

from quebra.analyzers import windows
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
    return result.windows, result.reads


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
    )


def test_builder_populates_all_derived_fields() -> None:
    d = _t2star_like()
    labels = {lbl for lbl, _, _ in d.thresholds}
    assert set(d.reliability.cumulative_time_per_threshold) == labels
    assert set(d.reliability.cumulative_damage_per_threshold) == labels
    assert set(d.reliability.ttf_per_threshold) == labels
    assert set(d.reliability.threshold_window_stats) == labels
    assert set(d.reliability.survival_curve_min) == labels
    assert set(d.reliability.occupancy) == labels
    assert set(d.reliability.threshold_summary) == labels
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
    in-spec as `T2* >= threshold`, so a value AT the threshold is in spec and must not be
    charged. Flipping `_out_of_spec_mask`'s `<` to `<=` moves this to 0.25, which is the
    boundary convention that no other test in the suite pins.
    """
    t_h = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    series = np.array([10.0, 5.0, 1.0, 1.0, 10.0])
    frac = _threshold_in_spec_frac(t_h, series, [("5", 5.0, True)], [])
    assert frac["5"] == pytest.approx(0.5, abs=1e-12)


def test_in_spec_frac_matches_summary() -> None:
    """Consistency between two consumers of one convention, NOT an oracle for it.

    Both sides descend from `_out_of_spec_mask`, so a mutation to that helper moves them
    together and this identity still holds. The oracle for the convention itself is the
    hand-computed case above; this test guards only that the summary and the band do not
    drift apart.
    """
    d = _t2star_like()
    for label, summ in d.reliability.threshold_summary.items():
        assert summ is not None  # dense synthetic series
        expected_frac_oos = 100.0 * (1.0 - d.reliability.occupancy[label])
        assert abs(summ["frac_oos_pct"] - expected_frac_oos) < 1e-9
        assert 0.0 <= d.reliability.occupancy[label] <= 1.0


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
        meta={},
        windows=win,
        reads=rd,
        gap_spans_s=[],
    )
    squared = build_within_calibration_panel_data(
        t_h=t_h,
        primary_series=s,
        primary_label="T2* (µs)",
        thresholds=thr,
        meta={},
        windows=win,
        reads=rd,
        gap_spans_s=[],
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
        windows=carved.windows,
        reads=carved.reads,
        gap_spans_s=carved.diagnostics["gap_spans_s"],
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
        (compute._threshold_summary, "gap_spans_h"),
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
    stand_in = types.SimpleNamespace(frame=frame, meta={})
    for adapter in (t2star.make_panel_data, fidelity.make_panel_data):
        with pytest.raises(ValueError, match="gap_spans_s is required"):
            adapter(
                stand_in,  # type: ignore[arg-type]
                windows=pd.DataFrame(),
                reads=pd.DataFrame(),
                gap_spans_s=None,  # type: ignore[arg-type]
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
            meta={},
            windows=win,
            reads=rd,
            gap_spans_s=None,  # type: ignore[arg-type]
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
