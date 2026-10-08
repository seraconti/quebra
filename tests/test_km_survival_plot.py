"""The Kaplan-Meier figure's legend: what the image alone says about each curve."""

from __future__ import annotations

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from quebra.analyzers import kaplan_meier as km
from quebra.plots.km_survival_plot import KMSurvivalPlot

pytestmark = pytest.mark.unit


def test_the_legend_states_how_many_windows_each_curve_rests_on() -> None:
    """Oracle: the number of windows handed to the estimator, set by hand per curve.

    The module docstring records a deviation from FIGURE_STANDARD: the figure does not state
    what it excluded, and the legend's `n` is the stated mitigation, the only trace in the
    image of how much data backs each curve. So each entry carries its own curve's `n`.
    """
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(3)
    sizes = {"early": 12, "late": 7}
    # Every other count on the curve differs from n, so a legend quoting the wrong one fails:
    # two of "late"'s windows are censored (deaths 5), and both curves carve more windows
    # than they estimate from, the difference being windows whose birth was not observed.
    observed = {
        "early": np.ones(12, dtype=bool),
        "late": np.array([1, 1, 0, 1, 0, 1, 1], dtype=bool),
    }
    carved = {"early": 15, "late": 9}
    curves = [
        km.run(
            km.KaplanMeierInputs(
                duration_min=rng.exponential(30.0, n) + 1.0,
                death_observed=observed[label],
                label=label,
                threshold_label="3 µs",
                n_windows_carved=carved[label],
                n_unobserved_birth_dropped=carved[label] - n,
            )
        )
        for label, n in sizes.items()
    ]
    for curve, n in zip(curves, sizes.values()):
        assert curve.n_windows == n, (
            "the estimator did not use every window handed to it"
        )
        assert curve.n_windows_carved > n, (
            "the fixture no longer separates n from carved"
        )
    assert curves[1].n_deaths == 5, (
        "the fixture no longer separates n from the death count"
    )
    fig = KMSurvivalPlot(name="unit").build_matplotlib(km.compare(curves, "3 µs"))
    try:
        texts = [
            text.get_text()
            for ax in fig.axes
            if ax.get_legend() is not None
            for text in ax.get_legend().get_texts()
        ]
    finally:
        plt.close(fig)
    for label, n in sizes.items():
        assert any(t.startswith(f"{label},") and f"n = {n}," in t for t in texts), texts


def test_the_band_fills_every_segment_whose_bounds_are_defined() -> None:
    """Oracle: hand-set durations, checked against the band each curve carries.

    `KaplanMeierCurve` gives the log-log band on every segment whose two bounds are finite,
    out to `max_observed_min`. Deaths at 1, 2 and 3 put S = 1/3 on [2, 3) with finite
    bounds and S = 0 at 3 with none, so the last segment before S reaches 0 must be filled.
    One death at 2 then windows censored at 20 and 40 leaves a tail [2, 40) that must be
    filled too. Each curve's band is one `fill_between` collection, drawn in pair order.
    """
    import matplotlib.pyplot as plt

    hand_set = {
        "to_zero": ([1.0, 2.0, 3.0], [True, True, True]),
        "censored_tail": ([2.0, 20.0, 40.0], [True, False, False]),
    }
    curves = [
        km.run(
            km.KaplanMeierInputs(
                duration_min=np.array(duration_min),
                death_observed=np.array(death_observed),
                label=label,
                threshold_label="3 µs",
                n_windows_carved=len(duration_min),
            )
        )
        for label, (duration_min, death_observed) in hand_set.items()
    ]
    to_zero = curves[0]
    assert np.isfinite(to_zero.band_lower[-2]), "premise: [2, 3) has a band"
    assert np.isnan(to_zero.band_lower[-1]), "premise: S = 0 has no band"

    comparison = km.compare(curves, "3 µs")
    fig = KMSurvivalPlot(name="unit").build_matplotlib(comparison)
    try:
        (ax,) = fig.axes
        bands = ax.collections
        assert len(bands) == 2, f"expected one band per curve, got {len(bands)}"
        for curve, band in zip(comparison.pair_curves(), bands, strict=True):
            ends = np.append(curve.time_min[1:], curve.max_observed_min)
            checked = 0
            for left, right, lo, hi in zip(
                curve.time_min, ends, curve.band_lower, curve.band_upper
            ):
                if right <= left or not (np.isfinite(lo) and np.isfinite(hi)):
                    continue
                point = (0.5 * (left + right), 0.5 * (lo + hi))
                assert any(path.contains_point(point) for path in band.get_paths()), (
                    f"{curve.label}: no band over [{left}, {right})"
                )
                checked += 1
            assert checked, f"{curve.label}: no defined segment was checked"
    finally:
        plt.close(fig)
