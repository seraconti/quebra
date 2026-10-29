"""The spectrum, spectral window and G10 of an irregular record, against known answers.

Oracles:
- the spectral window: exactly 1 at the alias frequencies of a perfect lattice, and 1 / N in
  expectation for scattered read times (N unit phasors with independent phases); computed in
  blocks of frequencies, it equals the same window computed in one block;
- the white level: white noise of variance sigma^2 has the one-sided density 2 sigma^2
  Delta_bar, on a gapped irregular record too (a SCALE test, named as such);
- a sinusoid peaks in its own frequency bin;
- a power law S ∝ f^-alpha gives G10 = 10^(3 - alpha) for the continuous integral, and an
  exactly computable value on a ten-frequency grid;
- simulation truth for leakage: on a sparse, gapped schedule a red process (an
  Ornstein-Uhlenbeck process, a random walk) never reads a G10 below its analytic value, and a
  near-white process with a line near the band edge reads below it.
"""

from __future__ import annotations

import math

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
def test_the_window_in_blocks_equals_the_window_in_one(monkeypatch) -> None:
    """Oracle: the window's formula evaluated in one array, here in the test. Blocks of 7
    frequencies split 40 frequencies over 25 reads into five full blocks and a short one; a
    budget of 3 elements, under one read's worth, still evaluates one frequency per block.
    Every value must equal the formula's."""
    rng = np.random.default_rng(11)
    t = np.sort(rng.uniform(0.0, 1000.0, size=25))
    f = np.linspace(0.001, 0.05, 40)
    whole = np.abs(np.exp(-2j * np.pi * np.outer(f, t)).sum(axis=1)) ** 2 / t.size**2
    for budget in (7 * t.size + 3, 3):
        monkeypatch.setattr(sp, "_CHUNK_ELEMENTS", budget)
        np.testing.assert_array_equal(sp.spectral_window(t, f), whole)


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
    its factor 2 is far outside the tolerance. Over seeds 0-199, each drawing the times and
    then the noise from `default_rng(seed)` as this test does, the ratio measured 0.9986 on
    average with sd 0.025 and a largest deviation of 0.076; 10% is about four sd. The
    integral of S over the band is sigma^2 as a CONSISTENCY check.
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
    """Oracle: analytic. For S = f^-alpha, L(f) ∝ f^(3 - alpha), so G10 = 10^(3 - alpha)
    for the continuous integral. The grid holds 10,000 frequencies; the trapezoid and the
    head piece below 1 / T leave a relative error under 1e-3."""
    f = sp.frequency_grid(1e6, 1e-2)
    assert sp.g10(f, f**-alpha) == pytest.approx(10.0 ** (3.0 - alpha), rel=1e-3)


@pytest.mark.unit
def test_g10_on_the_shortest_grid_is_exactly_28_for_an_inverse_square_spectrum() -> (
    None
):
    """Oracle: hand arithmetic. Ten frequencies k / T below an edge at 10.9 / T; S = f^-2
    makes (2 pi f)^2 S constant, c = (2 pi)^2. The head piece, with S held flat below 1 / T,
    gives c / (3 T) where the true piece is c / T; nine trapezoid steps add 9 c / T. The
    decade is read below the top bin, f_top / 10 = 1 / T, so G10 = (28 / 3) / (1 / 3) = 28,
    where the continuous integral gives 10. Read below the edge instead (10.9 / T / 10) it
    would be 22.05; without the head piece L(1 / T) is 0 and g10 raises."""
    span = 1000.0
    f = sp.frequency_grid(span, 10.9 / span)
    assert f.size == 10
    assert sp.g10(f, f**-2.0) == pytest.approx(28.0, rel=1e-12)


@pytest.mark.unit
def test_g10_interpolates_the_lower_point_between_grid_frequencies() -> None:
    """Oracle: hand arithmetic. Twelve frequencies k / T (an edge at 12.5 / T); S = f^-2
    makes (2 pi f)^2 S constant, c. L(1 / T) = c / (3 T) and L(2 / T) = (1/3 + 1) c / T, so
    at f_top / 10 = 1.2 / T the linear reading is (1/3 + 0.2) c / T; L(f_top) = (12 - 2/3)
    c / T; G10 = (34 / 3) / (8 / 15) = 21.25. Reading the grid point below or above instead
    gives 34 or 8.5."""
    f = sp.frequency_grid(1000.0, 0.0125)
    assert f.size == 12
    assert sp.g10(f, f**-2.0) == pytest.approx(21.25, rel=1e-12)


@pytest.mark.unit
def test_ten_grid_frequencies_hold_a_decade_whatever_the_rounding() -> None:
    """Oracle: exact arithmetic. Ten frequencies k / T give f_top / 10 = 1 / T = f_0. In
    floating point f_top / 10 lands a bit below f_0 for some spans (78 of these 500 at
    f704c18), and a test that compares them exactly would refuse a grid that holds its
    decade."""
    spans = np.random.default_rng(2).uniform(100.0, 1e6, size=500)
    for span in spans:
        f = sp.frequency_grid(float(span), 10.5 / float(span))
        assert f.size == 10
        assert sp.has_decade(f), span
        assert sp.g10(f, f**-2.0) == pytest.approx(28.0, rel=1e-9)


@pytest.mark.unit
def test_the_grid_runs_from_one_over_t_to_the_last_frequency_below_the_edge() -> None:
    """Oracle: the definition. k / T from k = 1, the last within 1 / T below f_e, at a
    span that holds a non-integer number of frequencies."""
    span, edge = 1000.0, 0.01234
    f = sp.frequency_grid(span, edge)
    assert f[0] == 1.0 / span
    assert f[-1] <= edge < f[-1] + 1.0 / span
    np.testing.assert_allclose(np.diff(f), 1.0 / span)


@pytest.mark.unit
def test_a_large_constant_offset_leaves_the_spectrum_unchanged() -> None:
    """Oracle: invariance to a constant offset. The spectrum of q + 2e6 (the size of a
    fitted Ramsey frequency in hertz) is the spectrum of q; with the offset kept, its power
    leaks across the band. This pins the removal of the mean: centring the values once
    would pass it too, so it does not tell the floating mean from a single centring."""
    rng = np.random.default_rng(6)
    t = _gapped_irregular_times(rng, 2000)
    q = rng.normal(size=t.size) + np.sin(2.0 * np.pi * t / 5e4)
    f, spacing, _ = _band(t)
    np.testing.assert_allclose(
        sp.lomb_scargle_spectrum(t, q + 2e6, f, spacing),
        sp.lomb_scargle_spectrum(t, q, f, spacing),
        rtol=1e-6,
    )


def _sparse_gapped_times(rng: np.random.Generator, n: int = 700) -> np.ndarray:
    """A schedule like the sparser records: 60% of spacings 45-60 s, 40% 100-600 s, and
    twelve gaps of 3000 s."""
    short = rng.random(n - 1) < 0.6
    steps = np.where(
        short, rng.uniform(45.0, 60.0, n - 1), rng.uniform(100.0, 600.0, n - 1)
    )
    steps[rng.choice(n - 1, 12, replace=False)] += 3000.0
    return np.r_[0.0, np.cumsum(steps)]


def _ou(t_rel_s: np.ndarray, tau_s: float, rng: np.random.Generator) -> np.ndarray:
    """An Ornstein-Uhlenbeck process of variance 1, exact at any read times."""
    x = np.empty(t_rel_s.size)
    x[0] = rng.normal()
    a = np.exp(-np.diff(t_rel_s) / tau_s)
    z = rng.normal(size=t_rel_s.size - 1)
    for i in range(1, t_rel_s.size):
        x[i] = a[i - 1] * x[i - 1] + math.sqrt(1.0 - a[i - 1] ** 2) * z[i - 1]
    return x


# The line of the near-white case, as a fraction of the top grid frequency: near the edge.
_LINE_AT = 0.9


def _g10_truth(t_rel_s: np.ndarray, f_top_hz: float, process: str) -> float:
    """G10 of the process's own spectrum, by the analytic L(F) = integral of (2 pi f)^2 S.

    OU of variance 1 and time tau: S = 4 tau / (1 + (2 pi f tau)^2), L = (2 / (pi tau^2))
    (X - atan X), X = 2 pi F tau. Random walk: S ∝ f^-2, L ∝ F. White noise of variance 1
    on the median spacing d: S = 2 d, L = (2 pi)^2 2 d F^3 / 3. A line of variance 1 at f0
    adds (2 pi f0)^2 to L for F >= f0."""
    span_s = float(t_rel_s[-1] - t_rel_s[0])
    spacing_s = float(np.median(np.diff(t_rel_s)))
    tau_s = span_s / 10.0
    f0_hz = _LINE_AT * f_top_hz

    def L(f_hz: float) -> float:
        if process == "ou":
            x = 2.0 * math.pi * f_hz * tau_s
            return 2.0 / (math.pi * tau_s**2) * (x - math.atan(x))
        if process == "random walk":
            return f_hz
        white = (2.0 * math.pi) ** 2 * 2.0 * spacing_s * f_hz**3 / 3.0
        return white + ((2.0 * math.pi * f0_hz) ** 2 if f_hz >= f0_hz else 0.0)

    return L(f_top_hz) / L(f_top_hz / 10.0)


@pytest.mark.statistical
@pytest.mark.parametrize("process", ["ou", "random walk", "white and edge line"])
def test_leakage_moves_g10_up_for_red_spectra_on_a_sparse_gapped_schedule(
    process: str,
) -> None:
    """Oracle: simulation truth, the analytic G10 of each process (`_g10_truth`).

    On a sparse, gapped schedule every run of a red process (OU with tau = T / 10, a random
    walk) reads a G10 at or above its truth of about 10: leakage moves it up. A near-white
    process with a line at 0.9 of the top grid frequency is the exception: every run reads
    below its truth of about 3430, but stays far above the red truths. Ten seeds per
    process. Over seeds 0-29 of this recipe, each drawing from `default_rng(seed)` as the
    test does, the OU runs read 25.8 to 54.2 times their truth, the random walk 24.5 to
    46.5, and the line runs 0.43 to 0.70 of theirs (estimates 1460 to 2395), so neither
    assertion sits near its margin.
    """
    t_rel_s = _sparse_gapped_times(np.random.default_rng(0))
    spacing_s = float(np.median(np.diff(t_rel_s)))
    f_hz = sp.frequency_grid(float(t_rel_s[-1] - t_rel_s[0]), 1.0 / (2.0 * spacing_s))
    truth = _g10_truth(t_rel_s, float(f_hz[-1]), process)
    ratios = []
    for seed in range(10):
        rng = np.random.default_rng(seed)
        if process == "ou":
            q = _ou(t_rel_s, float(t_rel_s[-1] - t_rel_s[0]) / 10.0, rng)
        elif process == "random walk":
            steps = rng.normal(size=t_rel_s.size - 1) * np.sqrt(np.diff(t_rel_s))
            q = np.r_[0.0, np.cumsum(steps)]
        else:
            phase = rng.uniform(0.0, 2.0 * math.pi)
            f0_hz = _LINE_AT * f_hz[-1]
            line = math.sqrt(2.0) * np.sin(2.0 * math.pi * f0_hz * t_rel_s + phase)
            q = rng.normal(size=t_rel_s.size) + line
        estimate = sp.g10(f_hz, sp.lomb_scargle_spectrum(t_rel_s, q, f_hz, spacing_s))
        ratios.append(estimate / truth)
    if process == "white and edge line":
        assert max(ratios) < 1.0
        assert min(ratios) * truth > 100.0
    else:
        assert min(ratios) >= 1.0


@pytest.mark.unit
@pytest.mark.parametrize(
    "call, match",
    [
        (lambda: sp.frequency_grid(10.0, 0.01), "holds no frequency"),
        (lambda: sp.g10(np.array([0.01, 0.02]), np.ones(2)), "f_top / 10"),
        (lambda: sp.g10(np.arange(1.0, 21.0), np.zeros(20)), "not positive"),
        (
            lambda: sp.lomb_scargle_spectrum(
                np.arange(5.0), np.ones(4), np.ones(3), 1.0
            ),
            "aligned",
        ),
        (
            lambda: sp.lomb_scargle_spectrum(
                np.arange(5.0), np.full(5, 0.3), np.ones(3), 1.0
            ),
            "zero variance",
        ),
        (
            lambda: sp.lomb_scargle_spectrum(
                np.arange(5.0), np.arange(5.0), np.array([0.1, -0.1]), 1.0
            ),
            "finite and positive",
        ),
        (
            lambda: sp.lomb_scargle_spectrum(
                np.arange(5.0), np.arange(5.0), np.array([]), 1.0
            ),
            "non-empty",
        ),
        (lambda: sp.spectral_window(np.array([0.0, np.nan]), 0.1), "finite read times"),
        (lambda: sp.spectral_window(np.arange(3.0), [[0.1]]), "1-D"),
        (lambda: sp.spectral_window(np.arange(3.0), np.inf), "finite and positive"),
        (lambda: sp.spectral_window(np.ones((2, 2)), 0.1), "finite read times"),
    ],
    ids=[
        "no-frequency",
        "short-span-g10",
        "zero-derivative-variance",
        "misaligned",
        "constant-values",
        "negative-frequency",
        "no-frequencies",
        "nan-time",
        "2d-frequencies",
        "infinite-frequency",
        "2d-times",
    ],
)
def test_an_unanswerable_spectrum_raises(call, match) -> None:
    """Oracle: specification. None of these has a spectrum to report."""
    with pytest.raises(ValueError, match=match):
        call()


@pytest.mark.unit
@pytest.mark.parametrize("last_chunk", ["seven", "one"])
def test_chunking_the_frequencies_changes_nothing_but_rounding(
    monkeypatch, last_chunk: str
) -> None:
    """Oracle: one call over every frequency, pinned by a chunk size above the whole grid.
    Chunks of seven frequencies, or chunks leaving a single frequency last, must reproduce
    it, so no chunk is dropped, repeated or misaligned."""
    rng = np.random.default_rng(2)
    t = _gapped_irregular_times(rng, 400)
    q = rng.normal(size=t.size)
    f, spacing, _ = _band(t)
    monkeypatch.setattr(sp, "_CHUNK_ELEMENTS", (f.size + 1) * t.size)
    whole = sp.lomb_scargle_spectrum(t, q, f, spacing)
    per_chunk = 7 if last_chunk == "seven" else f.size - 1
    monkeypatch.setattr(sp, "_CHUNK_ELEMENTS", per_chunk * t.size)
    np.testing.assert_allclose(
        sp.lomb_scargle_spectrum(t, q, f, spacing), whole, rtol=1e-10
    )
