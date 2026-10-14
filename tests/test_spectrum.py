"""The spectrum, spectral window and G10 of an irregular record, against known answers.

Oracles:
- the spectral window: exactly 1 at the alias frequencies of a perfect lattice, and 1 / N in
  expectation for scattered read times (N unit phasors with independent phases);
- the white level: white noise of variance sigma^2 has the one-sided density 2 sigma^2
  Delta_bar, on a gapped irregular record too (a SCALE test, named as such);
- a sinusoid peaks in its own frequency bin;
- a power law S ∝ f^-alpha gives G10 = 10^(3 - alpha) exactly.
"""

from __future__ import annotations

import numpy as np
import pytest

from quebra.analyzers import spectrum as sp


def _gapped_irregular_times(rng: np.random.Generator, n: int = 4000) -> np.ndarray:
    """Spacings uniform on 10-40 s, with three gaps of 20,000 s added."""
    steps = rng.uniform(10.0, 40.0, size=n - 1)
    steps[[n // 4, n // 2, 3 * n // 4]] += 20000.0
    return np.r_[0.0, np.cumsum(steps)]


def _band(t: np.ndarray) -> tuple[np.ndarray, float, float]:
    spacing = float(np.median(np.diff(t)))
    edge = 1.0 / (2.0 * spacing)
    return sp.frequency_grid(float(t[-1] - t[0]), edge), spacing, edge


@pytest.mark.unit
def test_the_window_is_one_at_every_alias_of_a_perfect_lattice() -> None:
    """Oracle: analytic. On t_n = n Delta every phasor at f = k / Delta is 1."""
    t = np.arange(500, dtype=float) * 24.0
    np.testing.assert_allclose(
        sp.spectral_window(t, [1 / 24.0, 2 / 24.0, 3 / 24.0]), 1.0
    )


@pytest.mark.statistical
def test_the_window_is_one_over_n_in_expectation_for_scattered_times() -> None:
    """Oracle: analytic. N phasors with independent uniform phases have E|sum|^2 = N, so
    E[N W] = 1. N W is exponential with sd 1, so 400 draws put the mean within 0.2 at
    four standard errors."""
    rng = np.random.default_rng(3)
    n = 200
    draws = [
        n * sp.spectral_window(rng.uniform(0.0, 1e5, size=n), 0.01)[0]
        for _ in range(400)
    ]
    assert abs(np.mean(draws) - 1.0) < 0.2


@pytest.mark.unit
def test_a_sinusoid_peaks_in_its_own_bin_on_a_gapped_record() -> None:
    """Oracle: the sinusoid's frequency, placed exactly on the grid frequency k / T nearest
    1 / 600 s, so the peak must be that bin. The gaps throw sidelobes a few bins away; a
    frequency halfway between two bins could lose its peak to one of them, so it is not
    used. A spectrum evaluated at f instead of 2 pi f puts the peak a factor 2 pi away."""
    rng = np.random.default_rng(5)
    t = _gapped_irregular_times(rng, 2000)
    f, spacing, _ = _band(t)
    k = int(round((t[-1] - t[0]) / 600.0))
    f0 = f[k - 1]
    s = sp.lomb_scargle_spectrum(t, 2.0 * np.sin(2.0 * np.pi * f0 * t), f, spacing)
    assert int(np.argmax(s)) == k - 1


@pytest.mark.statistical
def test_white_noise_sits_at_two_sigma_squared_delta_on_a_gapped_record() -> None:
    """Oracle: the one-sided density of white noise, 2 sigma^2 Delta_bar (a SCALE test).

    Gapped and irregular, so a scale taken from T / N (1.6 Delta_bar here) or one missing
    its factor 2 is far outside the tolerance. Over 200 seeds of this generator the ratio
    measured 0.9986 on average with sd 0.025 and a largest deviation of 0.076; 10% is about
    four sd. The integral of S over the band is sigma^2 as a CONSISTENCY check.
    """
    rng = np.random.default_rng(11)
    t = _gapped_irregular_times(rng)
    sigma = 1.7
    f, spacing, _ = _band(t)
    s = sp.lomb_scargle_spectrum(t, rng.normal(scale=sigma, size=t.size), f, spacing)
    assert abs(s.mean() / (2.0 * sigma**2 * spacing) - 1.0) < 0.10
    assert abs(np.trapezoid(s, f) / sigma**2 - 1.0) < 0.10


@pytest.mark.unit
@pytest.mark.parametrize("alpha", [0.0, 1.0, 2.0])
def test_g10_recovers_an_exact_power_law(alpha: float) -> None:
    """Oracle: analytic. For S = f^-alpha, L(f) ∝ f^(3 - alpha), so G10 = 10^(3 - alpha).
    The grid holds 10,000 frequencies; the trapezoid and the head piece below 1 / T leave a
    relative error under 1e-3."""
    span = 1e6
    edge = 1e-2
    f = sp.frequency_grid(span, edge)
    ratio, alpha_hat = sp.g10(f, f**-alpha, edge)
    assert ratio == pytest.approx(10.0 ** (3.0 - alpha), rel=1e-3)
    assert alpha_hat == pytest.approx(alpha, abs=1e-3)


@pytest.mark.unit
@pytest.mark.parametrize(
    "call, match",
    [
        (lambda: sp.frequency_grid(10.0, 0.01), "holds no frequency"),
        (lambda: sp.g10(np.array([0.01, 0.02]), np.ones(2), 0.05), "f_e / 10"),
        (
            lambda: sp.lomb_scargle_spectrum(
                np.arange(5.0), np.ones(4), np.ones(3), 1.0
            ),
            "aligned",
        ),
    ],
    ids=["no-frequency", "short-span-g10", "misaligned"],
)
def test_an_unanswerable_spectrum_raises(call, match) -> None:
    """Oracle: specification. None of these has a spectrum to report."""
    with pytest.raises(ValueError, match=match):
        call()


@pytest.mark.unit
def test_chunking_the_frequencies_changes_nothing_but_rounding(monkeypatch) -> None:
    """Oracle: one call over every frequency. Chunks of a few frequencies must reproduce it,
    so no chunk is dropped, repeated or misaligned."""
    rng = np.random.default_rng(2)
    t = _gapped_irregular_times(rng, 400)
    q = rng.normal(size=t.size)
    f, spacing, _ = _band(t)
    whole = sp.lomb_scargle_spectrum(t, q, f, spacing)
    monkeypatch.setattr(sp, "_CHUNK_ELEMENTS", 7 * t.size)
    np.testing.assert_allclose(
        sp.lomb_scargle_spectrum(t, q, f, spacing), whole, rtol=1e-10
    )
