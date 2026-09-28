"""Regression gate for the AcrossCalibration panel split (builder vs renderer).

The output-builder must produce a COMPLETE typed artifact and the renderer must be a
pure function of it. Elapsed days, hours and the histogram are the pre-split arithmetic;
the binned statistics keep every bin, an empty one as NaN.
"""

from __future__ import annotations

import matplotlib
import numpy as np

matplotlib.use("Agg")

from quebra.panels._across_calibration_compute import (
    build_across_calibration_panel_data,
)
from quebra.panels.across_calibration import (
    AcrossCalibrationPanel,
    AcrossCalibrationPanelData,
)
import pytest

pytestmark = pytest.mark.unit


def _case(
    n: int = 320, span_days: float = 180.0, seed: int = 1
) -> AcrossCalibrationPanelData:
    rng = np.random.default_rng(seed)
    gaps = rng.gamma(2.0, 1.0, n)
    times = np.cumsum(gaps)
    times = times / times[-1] * (span_days * 86400.0)
    ev = 1.6e9 + times
    iv = np.diff(ev)
    ev = ev[1:]
    stats = {
        "count": int(len(iv)),
        "mean_s": float(np.mean(iv)),
        "std_s": float(np.std(iv)),
        "min_s": float(np.min(iv)),
        "max_s": float(np.max(iv)),
    }
    return build_across_calibration_panel_data(
        intervals_s=iv,
        event_times_unix_s=ev,
        stats=stats,
        meta={"qubit": "2", "device": "6D2S", "dataset_id": "unit"},
    )


def test_builder_populates_derived_fields() -> None:
    d = _case()
    assert len(d.elapsed_days) == len(d.event_times_unix_s)
    assert len(d.intervals_h) == len(d.intervals_s)
    assert d.elapsed_days[0] == 0.0
    assert len(d.binned_interval_stats) == 5
    centers = d.binned_interval_stats[0]
    assert len(centers) > 1  # 180-day span → multiple 14-day bins
    assert len(d.histogram_counts) > 0
    assert len(d.histogram_edges) == len(d.histogram_counts) + 1


def test_unit_conversions_are_exact() -> None:
    d = _case()
    np.testing.assert_array_equal(d.intervals_h, d.intervals_s / 3600.0)
    t0 = float(d.event_times_unix_s[0])
    np.testing.assert_array_equal(d.elapsed_days, (d.event_times_unix_s - t0) / 86400.0)


def test_histogram_counts_conserve_positive_intervals() -> None:
    d = _case()
    assert int(d.histogram_counts.sum()) == int(np.sum(d.intervals_s > 0))


def test_short_span_single_bin() -> None:
    d = _case(n=40, span_days=9.0, seed=2)
    centers = d.binned_interval_stats[0]
    assert len(centers) == 1  # span < 14 days → single bin
    assert centers[0] == 7.0, (
        f"the single bin is centred at {centers[0]}, not t0 + 7 days"
    )


def test_renderer_produces_figure() -> None:
    d = _case()
    fig = AcrossCalibrationPanel(name="unit").build_matplotlib(d)
    assert len(fig.axes) == 4
    matplotlib.pyplot.close(fig)


def test_renderer_rejects_wrong_type() -> None:
    import pytest

    with pytest.raises(TypeError):
        AcrossCalibrationPanel(name="unit").build_matplotlib(object())


def _from_event_days(days, intervals_h) -> AcrossCalibrationPanelData:
    """An artifact built by the real builder from event times given in elapsed days."""
    iv = np.asarray(intervals_h, dtype=float) * 3600.0
    return build_across_calibration_panel_data(
        intervals_s=iv,
        event_times_unix_s=1.6e9 + np.asarray(days, dtype=float) * 86400.0,
        stats={
            "count": int(len(iv)),
            "mean_s": float(np.mean(iv)),
            "std_s": float(np.std(iv)),
            "min_s": float(np.min(iv)),
            "max_s": float(np.max(iv)),
        },
        meta={"qubit": "2", "device": "6D2S", "dataset_id": "unit"},
    )


def test_an_empty_bin_is_nan_in_the_artifact_not_absent():
    """Oracle: bins constructed by hand, with one deliberately left empty.

    A line joining two points asserts something about the interval between them, so the
    median, the IQR band and the p90 line must not join the bins either side of a stretch with
    no calibration events. The empty bin is in the artifact as NaN, where every consumer
    inherits it; matplotlib breaks a line and a fill at NaN.

    Not justified by the read-gap rule. FIGURE_STANDARD's `gap` is a within-calibration
    inter-read interval, and this is the across-calibration tier.
    """
    # Edges 0, 14, 28, 42: events in the first and third bins, none in days 14-28.
    d = _from_event_days([0.0, 1.0, 30.0, 31.0], [10.0, 12.0, 20.0, 22.0])
    centers, medians, q1s, q3s, p90s = d.binned_interval_stats

    np.testing.assert_allclose(centers, [7.0, 21.0, 35.0])
    for name, arr in (("median", medians), ("q1", q1s), ("q3", q3s), ("p90", p90s)):
        assert np.isnan(arr[1]), f"the empty bin's {name} is {arr[1]}, not NaN"
        assert np.isfinite(arr[0]) and np.isfinite(arr[2]), (
            f"{name} lost a populated bin"
        )


def test_every_event_lands_in_exactly_one_bin():
    """Oracle: bin membership counted by hand, and each bin's statistics computed by hand.

    Bins are closed on the left and the final bin is also closed on the right, as in
    `np.histogram`. So an event on the first edge and one on the last edge are both binned,
    whether or not the span is an exact multiple of the bin width, and a record whose events
    all share one time still gets its bin. The second bin is skewed, so its mean differs from
    its median and every percentile sits at a distinct value: the statistics are pinned, not
    just the membership. Percentiles are numpy's default linear interpolation.
    """
    # Span 28 d = 2 x 14: edges 0, 14, 28. Day 14 opens bin 1; day 28 closes it.
    d = _from_event_days([0.0, 5.0, 14.0, 20.0, 28.0], [1.0, 2.0, 3.0, 4.0, 50.0])
    centers, medians, q1s, q3s, p90s = d.binned_interval_stats
    np.testing.assert_allclose(centers, [7.0, 21.0])
    # Bin 0 holds {1, 2}; bin 1 holds {3, 4, 50}, whose mean is 19.
    np.testing.assert_allclose(medians, [1.5, 4.0])
    np.testing.assert_allclose(q1s, [1.25, 3.5])
    np.testing.assert_allclose(q3s, [1.75, 27.0])
    np.testing.assert_allclose(p90s, [1.9, 40.8])

    one = _from_event_days([0.0], [7.0])
    centers, medians, _q1, _q3, _p90 = one.binned_interval_stats
    np.testing.assert_allclose(centers, [7.0])
    np.testing.assert_allclose(medians, [7.0])


def test_a_populated_bin_with_no_populated_neighbour_is_marked():
    """Oracle: bins built by hand so that the first, one interior and the last bin stand alone.

    A line needs two consecutive finite points, so a populated bin between two empty ones, or
    at either end beside one, draws nothing on the median or p90 line. The renderer marks
    exactly those bins on both, above the IQR fill, with no legend entry. A record of two or
    more bins with no empty bin gets no marker.
    """
    import matplotlib.pyplot as plt

    def rolling(d):
        fig = AcrossCalibrationPanel(name="unit").build_matplotlib(d)
        (ax,) = [a for a in fig.axes if a.get_title().startswith("14-day")]
        return fig, ax

    def markers(ax):
        return [
            ln
            for ln in ax.get_lines()
            if ln.get_marker() == "o" and ln.get_linestyle() == "None"
        ]

    # Edges 0, 14, ..., 112. Populated: 0-14, 28-42, 42-56, 70-84, 98-112. Empty: the rest.
    days = [0.0, 2.0, 30.0, 31.0, 45.0, 46.0, 72.0, 73.0, 100.0, 112.0]
    d = _from_event_days(days, np.arange(1.0, 11.0))
    centers, medians, _q1, _q3, p90s = d.binned_interval_stats
    assert np.isnan(medians).sum() == 3, "fixture lost its empty bins"
    lone = np.isin(centers, [7.0, 77.0, 105.0])

    fig, ax = rolling(d)
    try:
        drawn = markers(ax)
        assert len(drawn) == 2, (
            f"expected a median and a p90 marker set, got {len(drawn)}"
        )
        by_y = {}
        for ln in drawn:
            np.testing.assert_allclose(np.asarray(ln.get_xdata(), float), centers[lone])
            by_y[tuple(np.round(np.asarray(ln.get_ydata(), float), 9))] = ln
            assert ln.get_markersize() > 0, "an invisible marker marks nothing"
            assert ln.get_label().startswith("_"), (
                f"marker in the legend: {ln.get_label()}"
            )
            fill_z = max(c.get_zorder() for c in ax.collections)
            assert ln.get_zorder() > fill_z, "marker drawn under the IQR fill"
        assert set(by_y) == {
            tuple(np.round(medians[lone], 9)),
            tuple(np.round(p90s[lone], 9)),
        }, "markers are not at the lone bins' median and p90"
        legend = [t.get_text() for t in ax.get_legend().get_texts()]
        assert legend == ["Median", "IQR", "p90"], legend
    finally:
        plt.close(fig)

    full = _case()
    assert np.isfinite(full.binned_interval_stats[1]).all(), "fixture has an empty bin"
    fig, ax = rolling(full)
    try:
        assert markers(ax) == [], "a record with no empty bin gained a marker"
    finally:
        plt.close(fig)
