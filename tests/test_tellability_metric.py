"""Tellability's per-metric rows against known answers: spacing, Q1.1, Q1.2, Q2.1, Q2.3 and
the clipped fraction.

Oracles:
- arithmetic on hand-built read times: the effective spacing sum(delta^2) / sum(delta) is
  the spacing on a regular grid and 550 / 19 s on seven 10 s and three 40 s spacings;
- analytic: the band maximum of the spectral window is 1, at f_s, on a perfect lattice, and f_s
  is a grid point whatever the span; on a jittered lattice W at the lattice rate averages the
  jitter's characteristic function squared, and the band maximum 0.983 of it, the mean of
  sinc^2 over an offset within half a grid step; an alias at 2 f_s and the main lobe at 0 lie
  outside the cell, and an alias exactly on its lower end, f_s / 2, is read there;
- constructed truth: on reads near a lattice at a rate off f_s, the band maximum finds the peak
  beside f_s that W(f_s) misses;
- simulation truth: white noise whose per-read sigma is known gives a noise floor equal to
  the edge level (four times it when sigma is reported twice too large), a linear trend plus
  white noise gives a second-difference scale equal to sigma, and the second-difference
  ratio is 1 for white noise and small for a slow signal;
- hand arithmetic for the gap rule, repeated times, a three-read segment, the clipped
  fraction and the half records' band edge;
- the convention page for repeated times: the spectrum and the halves average them, the
  spectral window counts every read, and sigma collapses to the sigma of the mean.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from quebra.analyzers import spectrum
from quebra.analyzers import tellability as tl
from quebra.analyzers.windows import observed_reads


def _obs(t, q, sigma=None, gap_mult=10.0):
    return observed_reads(
        np.asarray(t, dtype=float),
        np.asarray(q, dtype=float),
        sigma=sigma,
        gap_mult=gap_mult,
    )


@pytest.mark.unit
def test_the_effective_spacing_of_a_regular_grid_with_a_gap_is_the_spacing() -> None:
    """Oracle: arithmetic. Every within-segment spacing is 10 s; the 1000 s gap is none."""
    t = np.r_[np.arange(8) * 10.0, 1070.0 + np.arange(5) * 10.0]
    row = tl.spacing_row(_obs(t, np.ones(t.size)))
    assert (row.effective_spacing_s, row.mean_spacing_s, row.spacing_cv) == (
        10.0,
        10.0,
        0.0,
    )


@pytest.mark.unit
def test_the_effective_spacing_is_length_biased() -> None:
    """Oracle: arithmetic. Seven 10 s and three 40 s spacings, and a 1000 s gap: the mean is
    19 s and Delta_eff = (7 x 100 + 3 x 1600) / 190 = 550 / 19 s. Delta_eff = mean x (1 + CV^2)
    holds with the 1/n variance, as a CONSISTENCY check."""
    steps = [10.0] * 4 + [40.0] * 2 + [1000.0] + [10.0] * 3 + [40.0]
    t = np.r_[0.0, np.cumsum(steps)]
    row = tl.spacing_row(_obs(t, np.ones(t.size)))
    assert row.mean_spacing_s == pytest.approx(19.0, rel=1e-15)
    assert row.effective_spacing_s == pytest.approx(550.0 / 19.0, rel=1e-15)
    assert row.effective_spacing_s == pytest.approx(
        row.mean_spacing_s * (1.0 + row.spacing_cv**2), rel=1e-12
    )


@pytest.mark.unit
def test_the_window_band_maximum_is_one_at_f_s_on_a_lattice() -> None:
    """Oracle: analytic. f_s = 1 / Delta_bar is an alias of a perfect lattice, where W is 1.
    At 0.1 s, 4 T f_s rounds to 1196.0000000000018, which an odd step count would put half a
    step off f_s."""
    for spacing_s in (24.0, 0.1):
        t = np.arange(300) * spacing_s
        band = tl.window_band_maximum(_obs(t, np.ones(t.size)))
        assert band.window_band_max == pytest.approx(1.0), spacing_s
        assert band.window_band_peak_over_fs == 1.0, spacing_s


@pytest.mark.unit
def test_f_s_is_a_grid_point_when_the_span_is_not_a_multiple_of_the_spacing() -> None:
    """Oracle: analytic. 501 reads 24 s apart and one more 6 s after the last: the median
    spacing stays 24 s and the span is 12006 s, 500.25 spacings. The step count is even, so
    f_s is the middle grid point and the maximum is W(f_s) itself."""
    t = np.r_[np.arange(501) * 24.0, 500 * 24.0 + 6.0]
    obs = _obs(t, np.ones(t.size))
    first_alias, _ = tl.band_edges(obs)
    band = tl.window_band_maximum(obs)
    assert band.window_band_peak_over_fs == 1.0
    assert band.window_band_max == pytest.approx(
        spectrum.spectral_window(t, first_alias)[0], rel=1e-12
    )


@pytest.mark.statistical
def test_the_window_band_maximum_of_a_jittered_lattice_reads_its_coherence() -> None:
    """Oracle: analytic. On t_n = n Delta + N(0, sigma^2), W at the lattice rate 1 / Delta has
    mean 1/N + (1 - 1/N) exp(-(2 pi sigma / Delta)^2), the jitter's characteristic function
    squared. The lattice rate is not f_s = 1 / median spacing, so it falls at an offset from
    the grid that varies with the seed, uniform within half a step, 1/(8T). A 1/T-wide peak
    read there keeps on average the mean of sinc^2(u) over u in [0, 1/8], 0.983. Over 200
    seeds, W at the rate agrees with the analytic mean, and the band maximum with 0.983 of
    it, each within 4 standard errors.
    """
    n, delta, sigma, seeds = 500, 24.0, 0.05 * 24.0, 200
    expected = 1.0 / n + (1.0 - 1.0 / n) * math.exp(
        -((2.0 * math.pi * sigma / delta) ** 2)
    )
    at_rate, band = [], []
    for seed in range(seeds):
        rng = np.random.default_rng(seed)
        t = np.sort(np.arange(n) * delta + rng.normal(0.0, sigma, n))
        at_rate.append(spectrum.spectral_window(t, 1.0 / delta)[0])
        band.append(tl.window_band_maximum(_obs(t, np.ones(n))).window_band_max)
    se_rate = float(np.std(at_rate, ddof=1)) / math.sqrt(seeds)
    se_band = float(np.std(band, ddof=1)) / math.sqrt(seeds)
    assert abs(float(np.mean(at_rate)) - expected) <= 4.0 * se_rate
    u = np.linspace(0.0, 0.125, 100_001)[1:]
    grid_loss = float(np.mean((np.sin(np.pi * u) / (np.pi * u)) ** 2))
    assert abs(float(np.mean(band)) - grid_loss * expected) <= 4.0 * se_band


@pytest.mark.unit
def test_the_window_band_maximum_finds_the_peak_beside_f_s() -> None:
    """Oracle: constructed truth. Nine spacings of 24 s, then one of 28.8 s, repeated: the
    median spacing is 24 s, but the reads sit close to a lattice at the mean spacing, 24.48
    s. W(f_s) reads the tail of that lattice's peak; the band maximum finds the peak, at
    24 / 24.48 of f_s, within one grid step, which is at most 1/(4T)."""
    steps = np.tile(np.r_[np.full(9, 24.0), 28.8], 50)[:499]
    t = np.r_[0.0, np.cumsum(steps)]
    obs = _obs(t, np.ones(t.size))
    first_alias, _ = tl.band_edges(obs)
    band = tl.window_band_maximum(obs)
    max_step_over_fs = 1.0 / (4.0 * (t[-1] - t[0]) * first_alias)
    assert spectrum.spectral_window(t, first_alias)[0] < 0.01
    assert band.window_band_max >= 0.5
    assert abs(band.window_band_peak_over_fs - 24.0 / 24.48) <= max_step_over_fs


@pytest.mark.unit
def test_the_alias_at_twice_f_s_lies_outside_the_cell() -> None:
    """Oracle: analytic. Reads on a 10 s lattice, two thirds of the spacings 20 s: the median
    spacing is 20 s, so f_s = 0.05 Hz, and every read is a multiple of 10 s, so W(0.1 Hz) =
    W(2 f_s) = 1. That alias is the next harmonic's, outside the cell, and the band maximum
    stays near the noise: the limit the second validity assumption names."""
    rng = np.random.default_rng(21)
    steps = rng.choice([10.0, 20.0], size=599, p=[1.0 / 3.0, 2.0 / 3.0])
    t = np.r_[0.0, np.cumsum(steps)]
    obs = _obs(t, np.ones(t.size))
    first_alias, _ = tl.band_edges(obs)
    assert first_alias == 1.0 / 20.0
    assert spectrum.spectral_window(t, 2.0 * first_alias)[0] == pytest.approx(1.0)
    assert tl.window_band_maximum(obs).window_band_max < 0.2


@pytest.mark.unit
def test_the_cell_keeps_the_main_lobe_out() -> None:
    """Oracle: analytic. W tends to 1 as f tends to 0 for any read times. On 1000 uniform
    random times nothing else is coherent, so a maximum near 1, or one below f_s / 2, could
    only come from the main lobe the cell's lower edge excludes."""
    rng = np.random.default_rng(22)
    t = np.sort(rng.uniform(0.0, 50_000.0, size=1000))
    band = tl.window_band_maximum(_obs(t, np.ones(t.size)))
    assert band.window_band_max < 0.1
    assert 0.5 < band.window_band_peak_over_fs < 1.5


@pytest.mark.unit
def test_an_alias_on_the_cells_lower_end_is_read_there() -> None:
    """Oracle: analytic. Alternating 10 s and 30 s spacings: the median spacing is 20 s, so
    f_s = 0.05 Hz, and the reads repeat every 40 s, so W peaks at 0.025 Hz = f_s / 2, where
    each pair's two phasors sit a quarter turn apart and W = |1 + e^(-i pi/2)|^2 / 4 = 0.5.
    The pairs give the same 0.5 at 3 f_s / 2; one last read 35 s after the final pair, in
    phase with them at f_s / 2 and against them at 3 f_s / 2, makes the lower end the
    maximum, (300 sqrt(2) + 1)^2 / 601^2 = 0.5007. Both ends belong to the cell."""
    steps = np.r_[np.tile([10.0, 30.0], 300)[:-1], 35.0]
    t = np.r_[0.0, np.cumsum(steps)]
    obs = _obs(t, np.ones(t.size))
    first_alias, _ = tl.band_edges(obs)
    assert first_alias == 0.05
    band = tl.window_band_maximum(obs)
    assert band.window_band_peak_over_fs == 0.5
    assert band.window_band_max == pytest.approx(
        (300.0 * math.sqrt(2.0) + 1.0) ** 2 / 601.0**2, rel=1e-9
    )


@pytest.mark.statistical
def test_the_noise_floor_matches_the_edge_for_white_noise_with_known_sigma() -> None:
    """Oracle: simulation truth. White noise whose per-read sigma is reported exactly has
    S = 2 mean(sigma^2) Delta_bar everywhere, so floor / S(f_e) is 1.

    Half the reads have sigma 0.2 and half 2.0, so a floor built from mean(sigma)^2 instead
    of mean(sigma^2) reads 1.21 / 2.02, about 0.6. Over seeds 0-59, each drawing the times,
    sigma and the noise from `default_rng(seed)` as this test does, the ratio measured
    0.9935 on average with sd 0.064 and a largest deviation of 0.21, and that mutant at most
    0.71; 0.25 is about four sd. The mean over about 300 edge ordinates also
    carries the +1/K bias of a reciprocal, near 0.3% here and inside the tolerance.
    """
    rng = np.random.default_rng(4)
    n = 6000
    t = np.r_[0.0, np.cumsum(rng.uniform(10.0, 40.0, size=n - 1))]
    sigma = np.where(rng.random(n) < 0.5, 0.2, 2.0)
    obs = _obs(t, rng.normal(scale=sigma), sigma=sigma)
    _, ratio = tl.noise_floor_ratio(obs, tl.spectrum_rows(obs), sigma)
    assert abs(ratio - 1.0) < 0.25


@pytest.mark.statistical
def test_the_second_difference_scale_recovers_sigma_under_a_steep_trend() -> None:
    """Oracle: simulation truth. Second differences weighted for the spacing are zero on a
    linear trend, so a trend plus white noise of sd 0.5 gives sigma_hat = 0.5. Equal-spacing
    weights leave part of the trend in on irregular reads: over seeds 0-99, each drawing the
    times and the noise from `default_rng(seed)` as this test does, the right scale measured
    0.9989 sigma on average (sd 0.011) and the equal-weight one 1.118 sigma; 5% is about
    four sd."""
    rng = np.random.default_rng(6)
    n = 6000
    t = np.r_[0.0, np.cumsum(rng.uniform(10.0, 40.0, size=n - 1))]
    obs = _obs(t, 0.05 * t + rng.normal(scale=0.5, size=n))
    assert tl.second_difference_scale(obs) == pytest.approx(0.5, rel=0.05)


def _skewed_times(rng: np.random.Generator, n: int) -> np.ndarray:
    """Spacings of 10 s (70%) and 40 s (30%): median 10 s, Delta_eff about 28.9 s."""
    return np.r_[0.0, np.cumsum(np.where(rng.random(n - 1) < 0.7, 10.0, 40.0))]


@pytest.mark.statistical
def test_a_sigma_reported_twice_too_large_reads_a_noise_floor_ratio_of_four() -> None:
    """Oracle: arithmetic on the white level. White noise of sd 1 reported as 2 has a floor
    of 2 x 4 x Delta_bar against an edge level of 2 x 1 x Delta_bar: a ratio of 4, read
    the same whatever Delta_bar is, since floor and spectrum share it. Inverted it reads
    1/4; on these skewed spacings a floor taken on Delta_eff instead of the median spacing
    reads about 11.6. The tolerance is the white-noise test's 25%."""
    rng = np.random.default_rng(31)
    t = _skewed_times(rng, 6000)
    obs = _obs(t, rng.normal(size=t.size), sigma=np.full(t.size, 2.0))
    _, ratio = tl.noise_floor_ratio(obs, tl.spectrum_rows(obs), obs.sigma)
    assert ratio == pytest.approx(4.0, rel=0.25)


@pytest.mark.statistical
def test_the_edge_level_is_the_top_tenth_of_the_band_not_the_whole_band() -> None:
    """Oracle: arithmetic on a known spectrum. White noise of sd 1, reported as 1, plus two
    sines of variance 8 on grid frequencies, one at the fifth and one near 0.7 f_e: the top
    tenth of the band holds the white level alone, so the ratio is near 1. The band's mean
    level also holds both sines, (1 + 16) / f_e against 1 / f_e, a ratio near 1/17, and a
    top half of the band holds the second, (1 + 16) / f_e over half the band, near 1/17
    too. The reads are evenly spaced: on irregular reads a sine would leak into the top
    tenth (`spectrum.g10`)."""
    rng = np.random.default_rng(32)
    t = np.arange(6000) * 24.0
    span_s = t[-1] - t[0]
    f_low_hz = 5.0 / span_s
    f_high_hz = np.floor(0.7 * span_s / (2.0 * 24.0)) / span_s
    q = (
        rng.normal(size=t.size)
        + 4.0 * np.sin(2.0 * np.pi * f_low_hz * t)
        + 4.0 * np.sin(2.0 * np.pi * f_high_hz * t)
    )
    obs = _obs(t, q, sigma=np.ones(t.size))
    _, ratio = tl.noise_floor_ratio(obs, tl.spectrum_rows(obs), obs.sigma)
    assert ratio == pytest.approx(1.0, rel=0.25)


@pytest.mark.statistical
@pytest.mark.parametrize("signal", ["white", "slow"])
def test_the_second_difference_ratio_reads_where_the_variance_lives(
    signal: str,
) -> None:
    """Oracle: simulation truth. For white noise sigma_hat_D2 estimates sd(q), so the ratio
    is 1 (5%, the scale test's tolerance). For a slow sine of amplitude 3 plus noise of sd
    0.1, sigma_hat_D2 is 0.1 and sd(q) about 2.1, a ratio near 0.05; inverted it would read
    about 21."""
    rng = np.random.default_rng(33)
    t = np.r_[0.0, np.cumsum(rng.uniform(10.0, 40.0, size=5999))]
    if signal == "white":
        q = rng.normal(size=t.size)
    else:
        q = 3.0 * np.sin(2.0 * np.pi * t / 5e4) + rng.normal(scale=0.1, size=t.size)
    ratio = tl.second_difference_row(_obs(t, q)).second_difference_ratio
    if signal == "white":
        assert ratio == pytest.approx(1.0, rel=0.05)
    else:
        assert ratio < 0.1


@pytest.mark.unit
def test_a_three_read_segment_gives_its_one_residual() -> None:
    """Oracle: hand arithmetic. Reads at 0, 10 and 30 s of 0, 5 and 1: a = 20 / 30,
    b = 10 / 30, residual a 0 + b 1 - 5 = -14/3, weight 1 / (a^2 + b^2 + 1) = 9 / 14, so
    sigma_hat^2 = (196 / 9) (9 / 14) = 14."""
    obs = _obs([0.0, 10.0, 30.0], [0.0, 5.0, 1.0])
    assert tl.second_difference_scale(obs) == pytest.approx(math.sqrt(14.0), rel=1e-12)


@pytest.mark.unit
def test_repeated_times_collapse_for_the_spectrum_and_halves_but_not_the_window() -> (
    None
):
    """Oracle: the convention page. The spectrum and the halves see repeated reads averaged
    to one; the spectral window counts every read."""
    rng = np.random.default_rng(15)
    times = np.r_[0.0, np.cumsum(rng.uniform(10.0, 40.0, size=399))]
    t = np.sort(np.r_[times, times[::7]])
    q = rng.normal(size=t.size)
    obs = _obs(t, q)
    collapsed = _obs(*tl.collapse_repeated_times(t, q))
    np.testing.assert_allclose(
        tl.spectrum_rows(obs).spectrum_per_hz,
        tl.spectrum_rows(collapsed).spectrum_per_hz,
        rtol=1e-12,
    )
    assert tl.g10_halves(obs) == tl.g10_halves(collapsed)
    first_alias, _ = tl.band_edges(obs)
    band = tl.window_band_maximum(obs)
    peak_hz = band.window_band_peak_over_fs * first_alias
    every_read = spectrum.spectral_window(t, peak_hz)[0]
    assert band.window_band_max == pytest.approx(every_read, rel=1e-9)
    collapsed_read = spectrum.spectral_window(collapsed.t_rel_s, peak_hz)[0]
    assert collapsed_read != pytest.approx(every_read, rel=1e-6)


@pytest.mark.unit
def test_sigma_collapses_to_the_sigma_of_the_mean() -> None:
    """Oracle: hand arithmetic. Reads of sigma 3 and 4 at one time make a read of sigma
    sqrt(9 + 16) / 2 = 2.5; a group holding a read with no sigma has none."""
    t = np.array([0.0, 10.0, 10.0, 20.0, 30.0, 30.0])
    sigma = np.array([1.0, 3.0, 4.0, 2.0, 1.0, np.nan])
    out = tl.collapse_repeated_sigma(t, sigma)
    np.testing.assert_array_equal(out[:3], [1.0, 2.5, 2.0])
    assert np.isnan(out[3])


@pytest.mark.statistical
def test_doubled_reads_keep_the_noise_floor_and_the_sigma_bound_at_one() -> None:
    """Oracle: simulation truth. Every time read twice, white noise of sd 1 reported as 1:
    the averaged reads have sd 1 / sqrt(2), and so does the sigma of their mean, so Q1.2
    and Q3.3 both read 1 (the white tolerance, 25%). Taking sigma read by read would put
    Q1.2 near 2 and Q3.3 near 1.4."""
    rng = np.random.default_rng(34)
    times = np.r_[0.0, np.cumsum(rng.uniform(10.0, 40.0, size=2999))]
    t = np.repeat(times, 2)
    sigma = np.ones(t.size)
    obs = _obs(t, rng.normal(size=t.size), sigma=sigma)
    _, ratio = tl.noise_floor_ratio(obs, tl.spectrum_rows(obs), sigma)
    assert ratio == pytest.approx(1.0, rel=0.25)
    assert tl.sigma_bound_ratio(obs, sigma) == pytest.approx(1.0, rel=0.1)


@pytest.mark.unit
def test_second_differences_follow_the_gap_rule_they_are_given() -> None:
    """Oracle: hand arithmetic. Two noiseless lines of slope 0.1 per s, the second starting
    at 50 after a gap, off the first line's extension (100). Under gap_mult 10 the gap
    splits them and every residual is 0. Under gap_mult 200 the two triples across the gap
    have residuals -500/970 and +500/970 and weight 1 / (1 + (960^2 + 10^2) / 970^2), and
    the four other triples 0, so sigma_hat^2 is their weighted sum over 6."""
    t = np.r_[np.arange(5) * 10.0, 1000.0 + np.arange(3) * 10.0]
    q = np.r_[np.arange(5.0), 50.0 + np.arange(3.0)]
    assert tl.second_difference_scale(_obs(t, q, gap_mult=10.0)) == 0.0
    weight = 1.0 / (1.0 + (960.0**2 + 10.0**2) / 970.0**2)
    expected = math.sqrt(2.0 * (500.0 / 970.0) ** 2 * weight / 6.0)
    assert tl.second_difference_scale(_obs(t, q, gap_mult=200.0)) == pytest.approx(
        expected, rel=1e-12
    )


@pytest.mark.unit
def test_repeated_times_collapse_to_their_mean_value() -> None:
    """Oracle: hand arithmetic. Reads at 10 s of 1 and 3 collapse to one read of 2, which
    puts every read on the line q = t / 5: the scale is 0."""
    obs = _obs([0.0, 10.0, 10.0, 20.0, 30.0], [0.0, 1.0, 3.0, 4.0, 6.0])
    assert tl.second_difference_scale(obs) == 0.0


@pytest.mark.unit
def test_the_half_records_read_the_whole_records_band_edge() -> None:
    """Oracle: the definition. Each half's G10 is taken up to the WHOLE record's f_e. The
    second half is read every 30 s, so its own edge would be a third of the whole one's.
    A strong line at 0.004 Hz sits inside the top decade of that own edge (from 0.0017 Hz)
    but below the whole edge's (from 0.005 Hz), so the two G10 differ; the test checks
    both that the halves use the whole edge and that the choice matters."""
    rng = np.random.default_rng(8)
    t = np.r_[np.arange(3000) * 10.0, 30000.0 + np.arange(1, 1001) * 30.0]
    q = rng.normal(size=t.size) + 5.0 * np.sin(2.0 * np.pi * 0.004 * t)
    obs = _obs(t, q)
    first, second = tl.g10_halves(obs)
    _, edge = tl.band_edges(obs)
    half = (t.size + 1) // 2
    for row, ts, qs in ((first, t[:half], q[:half]), (second, t[half:], q[half:])):
        f = spectrum.frequency_grid(ts[-1] - ts[0], edge)
        expected = spectrum.g10(
            f, spectrum.lomb_scargle_spectrum(ts, qs, f, obs.median_spacing_s)
        )
        assert row.g10 == pytest.approx(expected, rel=1e-12)
    own_edge = 1.0 / (2.0 * 30.0)
    ts, qs = t[half:], q[half:]
    f_own = spectrum.frequency_grid(ts[-1] - ts[0], own_edge)
    own = spectrum.g10(f_own, spectrum.lomb_scargle_spectrum(ts, qs, f_own, 30.0))
    assert not math.isclose(second.g10, own, rel_tol=0.05)


@pytest.mark.unit
@pytest.mark.parametrize("value", [5.0, 0.1])
def test_a_constant_record_has_no_spectrum(value: float) -> None:
    """Oracle: specification. Equal values have a zero spectrum, so the rows read off it
    are not computed rather than computed from rounding. 0.1 at repeated times also
    checks the range is read before the repeated reads are averaged: three reads of 0.1
    average to 0.10000000000000002."""
    t = np.r_[np.arange(400) * 10.0, np.full(2, 3990.0)]
    obs = _obs(np.sort(t), np.full(t.size, value))
    with pytest.raises(tl.NotComputed) as excinfo:
        tl.spectrum_rows(obs)
    assert excinfo.value.reason == "zero variance"


@pytest.mark.unit
def test_a_half_with_constant_values_has_no_g10() -> None:
    """Oracle: specification. The first half of the reads is exactly 0: the halves group
    is not computed, and nothing escapes as a plain error."""
    rng = np.random.default_rng(9)
    t = np.arange(400) * 10.0
    q = np.r_[np.zeros(200), rng.normal(size=200)]
    with pytest.raises(tl.NotComputed) as excinfo:
        tl.g10_halves(_obs(t, q))
    assert excinfo.value.reason == "zero variance"


@pytest.mark.unit
def test_g10_is_not_computed_only_for_a_short_span(monkeypatch) -> None:
    """Oracle: specification. Fewer than ten grid frequencies is the record's own limit;
    any other error from `spectrum.g10` is a bug and is raised as itself."""
    rng = np.random.default_rng(10)
    short = _obs(np.arange(15) * 10.0, rng.normal(size=15))
    with pytest.raises(tl.NotComputed) as excinfo:
        tl.g10_row(short, tl.spectrum_rows(short))
    assert excinfo.value.reason == "span shorter than 20 median spacings"

    obs = _obs(np.arange(400) * 10.0, rng.normal(size=400))
    rows = tl.spectrum_rows(obs)

    def broken(*args, **kwargs):
        raise ValueError("not a property of the record")

    monkeypatch.setattr(tl.spectrum, "g10", broken)
    with pytest.raises(ValueError, match="not a property") as excinfo:
        tl.g10_row(obs, rows)
    assert not isinstance(excinfo.value, tl.NotComputed)


@pytest.mark.unit
def test_the_clipped_fraction_counts_reads_at_the_floor() -> None:
    """Oracle: count. Two of four reads sit at the 1e-16 floor; a metric with no floor has
    no clipped fraction."""
    assert tl.clipped_fraction(np.array([1e-16, 1e-16, 1e-3, 0.5]), 1e-16) == 0.5
    with pytest.raises(tl.NotComputed, match="no floor for this metric"):
        tl.clipped_fraction(np.array([1.0, 2.0]), None)


@pytest.mark.unit
@pytest.mark.parametrize(
    "call, reason",
    [
        (
            lambda: tl.spacing_row(_obs([0.0], [1.0])),
            "no segment with at least 2 reads",
        ),
        (
            lambda: tl.second_difference_scale(
                _obs([0.0, 10.0, 1000.0, 1010.0], np.ones(4))
            ),
            "no segment with at least 3 reads",
        ),
        (
            lambda: tl.second_difference_row(_obs(np.arange(9.0), np.ones(9))),
            "zero variance",
        ),
        (
            lambda: tl.mean_sigma_squared(np.array([np.nan, np.inf])),
            "no read with a finite sigma",
        ),
    ],
    ids=["one-read", "short-segments", "constant", "no-sigma"],
)
def test_an_unanswerable_row_names_its_reason(call, reason) -> None:
    """Oracle: specification. Each degenerate record carries its reason, from the set."""
    with pytest.raises(tl.NotComputed) as caught:
        call()
    assert caught.value.reason == reason
    assert caught.value.status == f"not computed: {reason}"
