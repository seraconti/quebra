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
