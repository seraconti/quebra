"""The spectrum of an irregularly read metric record, its spectral window and its G10.

Everything is computed on the reads as observed, never on a resampled grid: the periodogram
is Lomb-Scargle, which takes the read times as they are. With the median read spacing
Delta_bar, the first alias frequency is f_s = 1 / Delta_bar and the band edge f_e =
1 / (2 Delta_bar), and the spectrum is evaluated on the independent frequencies k / T of the
record's span T, from 1 / T up to f_e.

- `spectral_window(t, f)` is W(f) = |sum_n exp(-2 pi i f t_n)|^2 / N^2: 1 at every alias of a
  perfect lattice, about 1 / N for scattered read times.
- `lomb_scargle_spectrum(t, q, f, spacing_s)` is S(f) = 2 Delta_bar P(f), with P the
  floating-mean Lomb-Scargle power of the values. For white noise of variance sigma^2, P
  averages sigma^2, so S is the one-sided spectral density 2 sigma^2 Delta_bar, which
  integrates to sigma^2 over [0, f_e]. For a coloured process on irregular reads, S is the
  density seen through the spectral window, which leaks power between frequencies; it is not
  the density itself. Unweighted: the reported error bars do not enter.
- `derivative_variance(f, S)` is L(f) = integral from 0 to f of (2 pi f')^2 S(f') df', and
  `g10` is G10 = L(f_top) / L(f_top / 10), f_top the last grid frequency. For the continuous
  integral, a flat spectrum gives 1000 and S ∝ f^-2 gives 10. No exponent is derived from
  G10: on irregular reads leakage biases it (see `g10`).

These are the spectral rows of Sera's L0 specification (Q1.1, Q1.2, Q2.3); the specification
names the estimator of S as a free choice, and Lomb-Scargle without interpolation is a
deliberate one.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.signal import lombscargle

# Frequencies x reads evaluated per call to `lombscargle`, which holds arrays of that size.
_CHUNK_ELEMENTS = 4_000_000


def spectral_window(
    t_rel_s: np.ndarray, frequency_hz: float | np.ndarray
) -> np.ndarray:
    """W(f) = |sum_n exp(-2 pi i f t_n)|^2 / N^2 at each frequency.

    Sampling convolves the true spectrum with W. W(f_s) near 0 means nothing above the band
    edge can fold coherently onto it; W(f_s) near 1 means it could.

    Validity assumptions:

    - Assumption: none about the signal; W depends on the read times only.
      Diagnostic: not applicable.
      Consequence of violation: not applicable.
      Reference: Deeming (1975), Astrophysics and Space Science 36, 137-158; Shapiro and
      Silverman (1960), Journal of the Society for Industrial and Applied Mathematics 8(2),
      225-248; Sera's L0 specification, section 3.1. No equation locator.
    """
    t = np.asarray(t_rel_s, dtype=float)
    if t.ndim != 1 or t.size == 0 or not np.all(np.isfinite(t)):
        raise ValueError("the spectral window needs a 1-D array of finite read times")
    # A scalar frequency is accepted: `window_at_first_alias` asks for one.
    f = _check_frequencies(np.atleast_1d(np.asarray(frequency_hz, dtype=float)))
    phase = np.exp(-2j * np.pi * np.outer(f, t))
    return np.abs(phase.sum(axis=1)) ** 2 / t.size**2


def _check_frequencies(f: np.ndarray) -> np.ndarray:
    if f.ndim != 1 or f.size == 0:
        raise ValueError(
            f"frequencies must be a non-empty 1-D array, got shape {f.shape}"
        )
    if not (np.all(np.isfinite(f)) and np.all(f > 0.0)):
        raise ValueError("frequencies must be finite and positive")
    return f


def frequency_grid(span_s: float, band_edge_hz: float) -> np.ndarray:
    """The independent frequencies k / T, k = 1, 2, ..., up to the band edge.

    Raises ValueError when the span holds no frequency below the edge.
    """
    if not (math.isfinite(span_s) and span_s > 0.0):
        raise ValueError(f"the record span must be finite and positive, got {span_s!r}")
    if not (math.isfinite(band_edge_hz) and band_edge_hz > 0.0):
        raise ValueError(
            f"the band edge must be finite and positive, got {band_edge_hz!r}"
        )
    n = int(math.floor(band_edge_hz * span_s))
    if n < 1:
        raise ValueError(
            f"a span of {span_s:.6g} s holds no frequency k / T below the band edge "
            f"{band_edge_hz:.6g} Hz"
        )
    return np.arange(1, n + 1, dtype=float) / span_s


def lomb_scargle_spectrum(
    t_rel_s: np.ndarray,
    values: np.ndarray,
    frequency_hz: np.ndarray,
    spacing_s: float,
) -> np.ndarray:
    """S(f) = 2 Delta_bar P(f): for white noise, the one-sided spectral density of the values.

    P is `scipy.signal.lombscargle` with a floating mean, evaluated at the angular
    frequencies 2 pi f. Units: the values' unit squared per hertz. Raises ValueError on
    values of zero range: their spectrum is zero, and any power computed is rounding.

    Validity assumptions:

    - Assumption: the process is stationary over the record.
      Diagnostic: the split-half G10 of the tellability rows.
      Consequence of violation: S averages a changing spectrum over the record.
      Reference: Lomb (1976), Astrophysics and Space Science 39, 447-462; Scargle (1982),
      The Astrophysical Journal 263, 835-853; the floating mean is Zechmeister and Kurster
      (2009), Astronomy and Astrophysics 496, 577-584, as `scipy.signal.lombscargle` cites
      them; no equation locator.
    - Assumption: P averages sigma^2 for white noise of variance sigma^2, so 2 Delta_bar P
      is its one-sided density.
      Diagnostic: the white-level test of this module, on a gapped irregular record.
      Consequence of violation: the noise floor of the Q1.2 row is off by the same factor.
      Reference: no source located for the scaling on irregular reads; measured in the test.
    """
    t = np.asarray(t_rel_s, dtype=float)
    q = np.asarray(values, dtype=float)
    f = _check_frequencies(np.asarray(frequency_hz, dtype=float))
    if t.shape != q.shape or t.ndim != 1:
        raise ValueError(
            f"times and values must be 1-D and aligned, got {t.shape}, {q.shape}"
        )
    if not (np.all(np.isfinite(t)) and np.all(np.isfinite(q))):
        raise ValueError("the spectrum needs finite times and values")
    if q.size == 0 or np.ptp(q) == 0.0:
        raise ValueError("zero variance: the values have zero range")
    if not (math.isfinite(spacing_s) and spacing_s > 0.0):
        raise ValueError(f"the spacing must be finite and positive, got {spacing_s!r}")
    # In chunks of frequencies, so the working arrays stay near _CHUNK_ELEMENTS elements
    # whatever the record length; chunking changes the result only by rounding.
    chunk = max(1, _CHUNK_ELEMENTS // max(t.size, 1))
    power = np.concatenate(
        [
            np.atleast_1d(
                np.asarray(
                    lombscargle(
                        t, q, 2.0 * np.pi * f[i : i + chunk], floating_mean=True
                    ),
                    dtype=float,
                )
            )
            for i in range(0, f.size, chunk)
        ]
    )
    return 2.0 * spacing_s * power


def derivative_variance(frequency_hz: np.ndarray, spectrum: np.ndarray) -> np.ndarray:
    """L(f) = integral from 0 to f of (2 pi f')^2 S(f') df', at every grid frequency.

    The trapezoid rule over the grid, plus the piece from 0 to the first frequency with S
    held at its first value.

    Validity assumptions:

    - Assumption: S is flat below the first grid frequency 1 / T.
      Diagnostic: none shipped.
      Consequence of violation: the head piece misstates L near the bottom of the grid. For
      a spectrum falling as f^-2 it gives a third of the true piece, so at the shortest span
      `g10` accepts (ten grid frequencies) G10 reads 28 where the continuous integral gives
      10; with thousands of frequencies the error falls below 1e-3 (both pinned in
      `test_spectrum.py`).
      Reference: no source located.
    """
    f = np.asarray(frequency_hz, dtype=float)
    s = np.asarray(spectrum, dtype=float)
    if f.shape != s.shape or f.ndim != 1 or f.size == 0:
        raise ValueError("frequencies and spectrum must be 1-D, aligned and non-empty")
    integrand = (2.0 * np.pi * f) ** 2 * s
    head = (2.0 * np.pi) ** 2 * s[0] * f[0] ** 3 / 3.0
    steps = 0.5 * (integrand[1:] + integrand[:-1]) * np.diff(f)
    return head + np.concatenate(([0.0], np.cumsum(steps)))


def g10(frequency_hz: np.ndarray, spectrum: np.ndarray) -> float:
    """G10 = L(f_top) / L(f_top / 10), f_top the last grid frequency.

    f_top lies within 1 / T below the band edge, and the decade below it is exact. L at
    f_top / 10 is read by linear interpolation on the grid. Raises ValueError when f_top / 10
    lies below the grid's first frequency (fewer than ten grid frequencies), or when
    L(f_top / 10) is not positive.

    G10 is reported as an indicator, never turned into a spectral exponent.

    Validity assumptions:

    - Assumption: the periodogram does not move power between frequencies, as on evenly
      spaced reads.
      Diagnostic: none shipped; W(f_s) measures folding from above the band edge, not
      leakage inside the band.
      Consequence of violation: on irregular reads, leakage moves G10 up for a spectrum
      that falls across the band, so a small G10 means a red spectrum. A large G10 can be a
      red spectrum read through leakage; on sparse, gappy schedules a red spectrum can read
      close to white. A near-white spectrum with power near the band edge can read below
      its truth, but stays large. `test_spectrum.py` pins both directions on a sparse,
      gapped schedule.
      Reference: no source located. Sera's L0 specification, section 4.3, reports its own
      calibration on synthetic series (saturating for spectra steeper than about f^-2.5);
      the read times of that calibration are not stated, and it is not reproduced here.
    """
    f = _check_frequencies(np.asarray(frequency_hz, dtype=float))
    lower = float(f[-1]) / 10.0
    if lower < f[0]:
        raise ValueError(
            f"f_top / 10 = {lower:.6g} Hz lies below the first frequency {f[0]:.6g} Hz; "
            f"the record spans too few read spacings for G10"
        )
    cumulative = derivative_variance(f, spectrum)
    at_top = float(cumulative[-1])
    at_decade = float(np.interp(lower, f, cumulative))
    if not at_decade > 0.0:
        raise ValueError(
            f"L(f_top / 10) = {at_decade!r} is not positive; G10 is undefined"
        )
    return at_top / at_decade
