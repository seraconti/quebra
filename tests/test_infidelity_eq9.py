"""The infidelity of Capannelli et al. Supp. V Eq. (9) as a function of x = Delta / Omega.

Subject: the helpers in `analyzers/fidelity.py` that turn a rung on the infidelity into a
rung on |x| (`x_threshold`, `detuning_ladder`, `DETUNING_LADDER`), the refusals that keep
"I <= T" an interval in |x|, the peak, the plotting bounds and the log-axis floor. Each test
names its oracle: the defining equation, the issue's table, the Taylor series, a
cancellation-free rearrangement of Eq. (9), the published Fig. 7b, or a brute-force grid.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize import brentq

from quebra.analyzers.fidelity import (
    ABS_X_FLOOR,
    DETUNING_LADDER,
    FIDELITY_LADDER,
    MAX_RUNG_INFIDELITY,
    MIN_RUNG_INFIDELITY,
    PEAK_X,
    detuning_ladder,
    gate_fidelity,
    infidelity_at,
    infidelity_bounds,
    x_threshold,
)

RUNGS = [value for _, value, _ in FIDELITY_LADDER]


def _stable_infidelity(x: np.ndarray) -> np.ndarray:
    """1 - Eq. (9) with no subtraction near 1, from u = sqrt(1 + x^2), e = u - 1:
    3 I = (x^2 / u^2) sin^2(u pi/4) + (e + 2 sin^2(e pi/4)) / u, every term >= 0."""
    x = np.asarray(x, dtype=float)
    u = np.sqrt(1.0 + x * x)
    e = x * x / (u + 1.0)
    return (
        (x * x / (u * u)) * np.sin(u * np.pi / 4.0) ** 2
        + (e + 2.0 * np.sin(e * np.pi / 4.0) ** 2) / u
    ) / 3.0


def _root(f, hi: float) -> float:
    return float(
        brentq(f, 0.0, hi, xtol=np.finfo(float).tiny, rtol=4.0 * np.finfo(float).eps)
    )


# --- the rung on |x| -----------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("rung", RUNGS)
def test_x_threshold_solves_the_defining_equation(rung: float) -> None:
    """Oracle: the defining equation I(x_T) = T, to 1e-12. A round trip through
    `gate_fidelity`, so a transcription error in Eq. (9) cannot turn it red; the series,
    the cancellation-free form and Fig. 7b below can."""
    x_t = x_threshold(rung)
    assert abs(1.0 - float(gate_fidelity(np.array(x_t), np.array(1.0))) - rung) <= 1e-12


@pytest.mark.unit
def test_x_threshold_reproduces_the_issue_table() -> None:
    """Oracle: issue #18's table, x_T(1e-2, 1e-3, 1e-4) = 0.173660, 0.054787, 0.017321 to
    6 decimals."""
    for rung, expected in ((1e-2, 0.173660), (1e-3, 0.054787), (1e-4, 0.017321)):
        assert abs(x_threshold(rung) - expected) < 5e-7, rung


@pytest.mark.unit
def test_infidelity_follows_its_taylor_series() -> None:
    """Oracle: the series I(x) = x^2/3 + c4 x^4 + O(x^6), c4 = (pi^2 + 4 pi - 28)/96, and
    I(2.2769) = 2/3 to 1e-4. Held on [0.02, 0.2]: below that 1 - F cancels."""
    c4 = (np.pi**2 + 4.0 * np.pi - 28.0) / 96.0
    x = np.geomspace(0.02, 0.2, 25)
    residual = infidelity_at(x) - x**2 / 3.0 - c4 * x**4
    assert np.all(np.abs(residual) <= 0.01 * x**6)
    assert abs(float(infidelity_at(2.2769)) - 2.0 / 3.0) < 1e-4


@pytest.mark.unit
def test_x_threshold_matches_a_cancellation_free_eq9() -> None:
    """Oracle: an independent evaluation of Eq. (9), rearranged so nothing cancels.
    x_T agrees to 1e-12 relative on the shipped rungs and to 1e-9 at MIN_RUNG_INFIDELITY,
    the precision `fidelity.MIN_RUNG_INFIDELITY`'s comment states."""
    for rung, tolerance in [(r, 1e-12) for r in RUNGS] + [(MIN_RUNG_INFIDELITY, 1e-9)]:
        stable = _root(lambda x: float(_stable_infidelity(x)) - rung, PEAK_X)
        assert abs(x_threshold(rung) - stable) <= tolerance * stable, rung


# Capannelli et al., npj Quantum Information 11:192 (2025), doi:10.1038/s41534-025-01134-6,
# Fig. 7b, on page 7 of the published PDF: the minimum
# Rabi frequency (MHz) that keeps each infidelity at each detuning (kHz), AS PRINTED.
_FIG_7B_DETUNING_KHZ = (50.0, 100.0, 150.0, 200.0)
_FIG_7B_MHZ = {
    1e-2: (0.29, 0.58, 0.86, 1.15),
    1e-3: (0.91, 1.83, 2.73, 3.65),
    1e-4: (2.89, 5.77, 8.66, 11.55),
}
# The one printed cell that is truncated, not rounded: Delta / x_T = 2.738 MHz.
_FIG_7B_TRUNCATED = (1e-3, 150.0)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("rung", "detuning_khz"),
    [(r, d) for r in _FIG_7B_MHZ for d in _FIG_7B_DETUNING_KHZ],
)
def test_minimum_rabi_reproduces_capannelli_fig_7b(rung, detuning_khz) -> None:
    """Oracle: Capannelli et al. Fig. 7b, the twelve printed cells. Omega_min =
    Delta / x_T rounds to each cell; the 0.1% / 150 kHz cell is printed truncated (2.73
    for 2.738), so it is held to within one unit of its last digit."""
    printed = _FIG_7B_MHZ[rung][_FIG_7B_DETUNING_KHZ.index(detuning_khz)]
    omega_mhz = detuning_khz / x_threshold(rung) / 1e3
    if (rung, detuning_khz) == _FIG_7B_TRUNCATED:
        assert 0.0 <= omega_mhz - printed < 0.01
    else:
        assert round(omega_mhz, 2) == printed


# --- where the rung stops being an interval in |x| -----------------------------------


@pytest.mark.unit
@pytest.mark.parametrize(
    "rung", [MAX_RUNG_INFIDELITY, 0.5, 2.0 / 3.0, 1e-8, 0.0, -1e-3, np.nan, np.inf]
)
def test_x_threshold_refuses_rungs_outside_its_range(rung: float) -> None:
    """Oracle: the stated range [MIN_RUNG_INFIDELITY, MAX_RUNG_INFIDELITY)."""
    with pytest.raises(ValueError, match="x_threshold needs"):
        x_threshold(rung)


@pytest.mark.unit
def test_x_threshold_accepts_both_ends_of_its_range() -> None:
    """Oracle: the stated range; its lower end is in, a value just under its upper end is in."""
    for rung in (MIN_RUNG_INFIDELITY, 0.3149):
        assert 0.0 < x_threshold(rung) < PEAK_X


@pytest.mark.unit
def test_the_infidelity_stays_above_the_largest_rung_past_the_peak() -> None:
    """Oracle: Cauchy-Schwarz, (u cos + sin)^2 <= u^2 + 1, gives I >= 1/3 - 1/(3 u^2) for
    every x, above MAX_RUNG_INFIDELITY from x = 4.2 on; a grid covers [PEAK_X, 4.2]."""
    x = np.linspace(0.0, 50.0, 200_001)
    u2 = 1.0 + x * x
    assert np.all(infidelity_at(x) >= 1.0 / 3.0 - 1.0 / (3.0 * u2) - 1e-15)
    assert 1.0 / 3.0 - 1.0 / (3.0 * (1.0 + 4.2**2)) > MAX_RUNG_INFIDELITY
    assert float(np.min(infidelity_at(np.linspace(PEAK_X, 4.2, 200_001)))) > (
        MAX_RUNG_INFIDELITY
    )


@pytest.mark.unit
def test_peak_x_is_the_first_root_of_the_rearranged_eq9() -> None:
    """Oracle: Eq. (9) as F = 1/3 + (u cos(u pi/4) + sin(u pi/4))^2 / (3 u^2), whose first
    zero of u cos(u pi/4) + sin(u pi/4) on (2, 4) is the peak, where I = 2/3."""
    x = np.linspace(0.0, 10.0, 10_001)
    u = np.sqrt(1.0 + x * x)
    rearranged = 1.0 / 3.0 + (
        u * np.cos(u * np.pi / 4) + np.sin(u * np.pi / 4)
    ) ** 2 / (3.0 * u * u)
    assert np.max(np.abs(gate_fidelity(x, np.ones_like(x)) - rearranged)) < 1e-14
    u_star = float(
        brentq(
            lambda v: v * np.cos(v * np.pi / 4) + np.sin(v * np.pi / 4),
            2.0,
            4.0,
            xtol=np.finfo(float).tiny,
            rtol=4.0 * np.finfo(float).eps,
        )
    )
    assert abs(PEAK_X - np.sqrt(u_star**2 - 1.0)) < 1e-12
    assert abs(float(infidelity_at(PEAK_X)) - 2.0 / 3.0) < 1e-15


# --- the ladder on |x| and the log-axis floor ----------------------------------------


@pytest.mark.unit
def test_the_detuning_ladder_literal_is_the_helpers_ladder() -> None:
    """Oracle: `detuning_ladder(FIDELITY_LADDER)`. Same labels in the same order, every
    rung small-good, each value equal to the helper's to 1e-12 relative."""
    computed = detuning_ladder(FIDELITY_LADDER)
    assert [label for label, _, _ in DETUNING_LADDER] == [
        label for label, _, _ in FIDELITY_LADDER
    ]
    assert not any(big for _, _, big in DETUNING_LADDER)
    for (_, literal, _), (_, value, _) in zip(DETUNING_LADDER, computed, strict=True):
        assert abs(literal - value) <= 1e-12 * value


@pytest.mark.unit
def test_detuning_ladder_refuses_a_big_good_rung() -> None:
    """Oracle: infidelity is small-good, so a big-good rung is a mistake, not a rung."""
    with pytest.raises(ValueError, match="big-values-good"):
        detuning_ladder([("99%", 1e-2, True)])


@pytest.mark.unit
def test_the_floor_is_an_infidelity_of_1e_16_far_below_every_rung() -> None:
    """Oracle: I = x^2/3 at small x, so I(ABS_X_FLOOR) = 1e-16; and the floor is three
    decades below the smallest rung, so clipping at it moves no read across a rung."""
    assert abs(ABS_X_FLOOR**2 / 3.0 - 1e-16) <= 1e-12 * 1e-16
    assert ABS_X_FLOOR < 1e-3 * min(value for _, value, _ in DETUNING_LADDER)


# --- the plotting bounds -------------------------------------------------------------


@pytest.mark.unit
def test_infidelity_bounds_are_the_extremes_over_the_interval() -> None:
    """Oracle: brute force, the min and max of I on a dense grid of each read's interval
    [max(x - k s, 0), x + k s], including reads whose interval reaches x = 0."""
    rng = np.random.default_rng(20261030)
    for _ in range(200):
        k = float(rng.uniform(0.5, 3.0))
        x = float(rng.uniform(0.0, 0.6))
        s = float(rng.uniform(0.0, 0.3))
        lower, upper = infidelity_bounds(np.array([x]), np.array([s]), k)
        grid = infidelity_at(np.linspace(max(x - k * s, 0.0), x + k * s, 2001))
        assert lower[0] == pytest.approx(float(np.min(grid)), abs=1e-15)
        assert upper[0] == pytest.approx(float(np.max(grid)), abs=1e-15)


@pytest.mark.unit
def test_infidelity_bounds_refuse_what_they_cannot_bound() -> None:
    """Oracle: the stated domain. Past PEAK_X the endpoints are not the extremes; a
    negative x or sigma is no |Delta| / Omega or error bar; a non-finite sigma, NaN or
    infinite, gives NaN bounds on that read and leaves the others alone."""
    with pytest.raises(ValueError, match="PEAK_X"):
        infidelity_bounds(np.array([2.0]), np.array([0.2]), 2.0)
    with pytest.raises(ValueError, match="negative"):
        infidelity_bounds(np.array([-0.1]), np.array([0.01]), 1.0)
    with pytest.raises(ValueError, match="negative"):
        infidelity_bounds(np.array([0.1]), np.array([-0.01]), 1.0)
    for bad in (np.nan, np.inf):
        lower, upper = infidelity_bounds(
            np.array([0.1, 0.1]), np.array([bad, 0.01]), 1.0
        )
        assert np.isnan(lower[0]) and np.isnan(upper[0]), bad
        assert lower[1] == pytest.approx(float(infidelity_at(0.09)), rel=1e-12), bad
        assert upper[1] == pytest.approx(float(infidelity_at(0.11)), rel=1e-12), bad
