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
  integrates to sigma^2 over [0, f_e]. Unweighted: the reported error bars do not enter.
- `derivative_variance(f, S)` is L(f) = integral from 0 to f of (2 pi f')^2 S(f') df', and
  `g10` is G10 = L(f_e) / L(f_e / 10) with alpha = 3 - log10 G10, which a power law
  S ∝ f^-alpha (alpha < 3) gives exactly.

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
      Reference: Deeming (1975), Fourier analysis with unequally-spaced data, and Shapiro
      and Silverman (1960), as cited in Sera's L0 specification; no equation locator.
    """
    t = np.asarray(t_rel_s, dtype=float)
    if t.size == 0:
        raise ValueError("the spectral window needs at least one read time")
    f = np.atleast_1d(np.asarray(frequency_hz, dtype=float))
    phase = np.exp(-2j * np.pi * np.outer(f, t))
    return np.abs(phase.sum(axis=1)) ** 2 / t.size**2


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
    """S(f) = 2 Delta_bar P(f), the one-sided spectral density of the values.

    P is `scipy.signal.lombscargle` with a floating mean, evaluated at the angular
    frequencies 2 pi f. Units: the values' unit squared per hertz.

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
    f = np.asarray(frequency_hz, dtype=float)
    if t.shape != q.shape or t.ndim != 1:
        raise ValueError(
            f"times and values must be 1-D and aligned, got {t.shape}, {q.shape}"
        )
    if not (np.all(np.isfinite(t)) and np.all(np.isfinite(q))):
        raise ValueError("the spectrum needs finite times and values")
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
    held at its first value. That piece is exact for a spectrum flat below 1 / T.
    """
    f = np.asarray(frequency_hz, dtype=float)
    s = np.asarray(spectrum, dtype=float)
    if f.shape != s.shape or f.ndim != 1 or f.size == 0:
        raise ValueError("frequencies and spectrum must be 1-D, aligned and non-empty")
    integrand = (2.0 * np.pi * f) ** 2 * s
    head = (2.0 * np.pi) ** 2 * s[0] * f[0] ** 3 / 3.0
    steps = 0.5 * (integrand[1:] + integrand[:-1]) * np.diff(f)
    return head + np.concatenate(([0.0], np.cumsum(steps)))


def g10(
    frequency_hz: np.ndarray, spectrum: np.ndarray, band_edge_hz: float
) -> tuple[float, float]:
    """G10 = L(f_e) / L(f_e / 10) and alpha = 3 - log10 G10.

    L is read at the two frequencies by linear interpolation on the grid; the grid ends at
    the last k / T below f_e, within 1 / T of it, and L(f_e) is read there. Raises
    ValueError when f_e / 10 lies below the grid's first frequency, which needs a span of
    at least 10 / f_e, or when L(f_e / 10) is not positive.

    Validity assumptions:

    - Assumption: S follows a single power law f^-alpha across the top decade of the band,
      with alpha < 3.
      Diagnostic: the split-half G10 values beside it.
      Consequence of violation: alpha is an effective exponent of that decade only; above
      alpha near 2.5 the estimate saturates (Sera's L0 specification, Q2.3).
      Reference: Sera's L0 specification, Q2.3; no equation locator.
    """
    f = np.asarray(frequency_hz, dtype=float)
    lower = band_edge_hz / 10.0
    if lower < f[0]:
        raise ValueError(
            f"f_e / 10 = {lower:.6g} Hz lies below the first frequency {f[0]:.6g} Hz; "
            f"the record spans too few read spacings for G10"
        )
    cumulative = derivative_variance(f, spectrum)
    at_edge = float(cumulative[-1])
    at_decade = float(np.interp(lower, f, cumulative))
    if not at_decade > 0.0:
        raise ValueError(
            f"L(f_e / 10) = {at_decade!r} is not positive; G10 is undefined"
        )
    ratio = at_edge / at_decade
    return ratio, 3.0 - math.log10(ratio)
