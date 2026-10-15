"""Tellability: can a record, read on its schedule with its error bars, tell when a crossing
happened?

The rows follow Sera's L0 specification, by question: Q1 (is the band edge trustworthy), Q2
(do crossings belong to the signal or to the read grid), Q3 (can a crossing be placed, given
the error bars) and Q4 (do the durations carry more than the in-spec fraction). The
definitions, and where they depart from L0, are in `docs/TELLABILITY_SEMANTICS.md`.

No row is turned into a label and nothing downstream reads a row to change what runs: every
number ships with a status, `computed`, `not computed: <reason>` or `not declared`, and the
reader reads the panel. A row that cannot be computed raises `NotComputed` with a reason from
the fixed set below; the step records it as that row's status.

Every row is computed on the reads as observed, under the carve's gap rule
(`windows.observed_reads`). A spacing that crosses a gap is never treated as observed time,
and no pair, difference or crossing spans a gap.
"""

from __future__ import annotations

import math
import re
import dataclasses
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from typing import TypeVar

import dcor
import numpy as np
import pandas as pd
from scipy import stats

from quebra.analyzers import spectrum
from quebra.analyzers.min_p import min_p
from quebra.analyzers.windows import (
    ObservedReads,
    WindowsResult,
    confident_mask,
    margin,
    observed_reads,
)
from quebra.core._artifact_guard import StaleArtifactGuard
from quebra.core.types import LOAD_COUNT_KEYS

# ---------------------------------------------------------------------------
# Status vocabulary
# ---------------------------------------------------------------------------

COMPUTED = "computed"
NOT_DECLARED = "not declared"

# The reasons a row can carry, as patterns: `<n>` is an integer. A reason outside this set is
# refused, so the vocabulary cannot drift one call site at a time.
REASON_TEMPLATES = (
    "sigma propagation not derived",
    "fewer than <n> reads",
    "no segment with at least <n> reads",
    "span shorter than <n> median spacings",
    "no down-crossing at this threshold",
    "no in-to-out transition within a segment",
    "fewer than two rungs with a down-crossing",
    "zero variance",
    "no floor for this metric",
    "no read with a finite sigma",
)
_REASON_PATTERNS = tuple(
    re.compile("^" + re.escape(template).replace(re.escape("<n>"), r"\d+") + "$")
    for template in REASON_TEMPLATES
)


def check_reason(reason: str) -> str:
    """Return `reason` if it matches one of `REASON_TEMPLATES`, or raise ValueError."""
    if not any(pattern.fullmatch(reason) for pattern in _REASON_PATTERNS):
        raise ValueError(
            f"{reason!r} is not a tellability reason; the set is {REASON_TEMPLATES}"
        )
    return reason


class NotComputed(ValueError):
    """A row the record cannot answer. `reason` is one of `REASON_TEMPLATES`."""

    def __init__(self, reason: str) -> None:
        self.reason = check_reason(reason)
        super().__init__(reason)

    @property
    def status(self) -> str:
        return f"not computed: {self.reason}"


# ---------------------------------------------------------------------------
# Shared preparation
# ---------------------------------------------------------------------------


def collapse_repeated_times(
    t_rel_s: np.ndarray, values: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Reads at one time collapsed to one read at their mean value.

    The spectrum and the second differences need distinct times. Repeated times are legal
    (`windows.spacing`), and averaging them is the one rule both rows use.
    """
    t = np.asarray(t_rel_s, dtype=float)
    q = np.asarray(values, dtype=float)
    unique, inverse = np.unique(t, return_inverse=True)
    if unique.size == t.size:
        return t, q
    sums = np.bincount(inverse, weights=q)
    counts = np.bincount(inverse)
    return unique, sums / counts


def collapse_repeated_sigma(t_rel_s: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    """The sigma of each read `collapse_repeated_times` makes: that of a mean of m reads,
    sqrt(sum sigma_i^2) / m, and unknown (NaN) when any read of the group has no sigma, as
    a read with no sigma is never sure (`windows.confident_mask`)."""
    t = np.asarray(t_rel_s, dtype=float)
    s = np.asarray(sigma, dtype=float)
    unique, inverse = np.unique(t, return_inverse=True)
    if unique.size == t.size:
        return s
    known = np.isfinite(s)
    counts = np.bincount(inverse)
    sum_sq = np.bincount(inverse, weights=np.where(known, s * s, 0.0))
    out = np.sqrt(sum_sq) / counts
    out[np.bincount(inverse, weights=~known) > 0] = np.nan
    return out


def mean_sigma_squared(sigma: np.ndarray) -> float:
    """sigma_bar^2 = mean(sigma_i^2) over the reads whose sigma is finite."""
    s = np.asarray(sigma, dtype=float)
    s = s[np.isfinite(s)]
    if s.size == 0:
        raise NotComputed("no read with a finite sigma")
    return float(np.mean(s * s))


# ---------------------------------------------------------------------------
# Per metric: spacing, Q1.1, Q1.2, Q2.1, Q2.3, clipped
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class SpacingRow:
    effective_spacing_s: float
    mean_spacing_s: float
    spacing_cv: float


def spacing_row(obs: ObservedReads) -> SpacingRow:
    """Delta_eff = sum(delta^2) / sum(delta) over the within-segment spacings, the mean
    spacing and the CV beside it.

    The CV takes the 1/n variance, so Delta_eff = mean x (1 + CV^2) holds exactly. Delta_eff
    is the mean length of the interval a crossing lands in when a stationary process crosses
    at a constant rate (Sera's L0 specification, section 8.3).
    """
    within = np.diff(obs.t_rel_s)[obs.within_segment()]
    if within.size == 0:
        raise NotComputed("no segment with at least 2 reads")
    mean = float(np.mean(within))
    if not mean > 0.0:
        raise NotComputed("zero variance")
    return SpacingRow(
        effective_spacing_s=obs.effective_spacing_s(),
        mean_spacing_s=mean,
        spacing_cv=float(np.std(within) / mean),
    )


def band_edges(obs: ObservedReads) -> tuple[float, float]:
    """(f_s, f_e) = (1 / Delta_bar, 1 / (2 Delta_bar)), Delta_bar the median positive
    spacing of the whole record, the carve's own."""
    if obs.t_rel_s.size < 2 or not math.isfinite(obs.median_spacing_s):
        raise NotComputed("fewer than 2 reads")
    return 1.0 / obs.median_spacing_s, 1.0 / (2.0 * obs.median_spacing_s)


def window_at_first_alias(obs: ObservedReads) -> float:
    """Q1.1: W(f_s), the spectral window at the first alias frequency."""
    first_alias_hz, _ = band_edges(obs)
    return float(spectrum.spectral_window(obs.t_rel_s, first_alias_hz)[0])


@dataclass(frozen=True, slots=True)
class SpectrumRows:
    frequency_hz: np.ndarray
    spectrum_per_hz: np.ndarray
    # mean of S over the top 10% of the band, [0.9 f_e, f_e]
    edge_level_per_hz: float


def spectrum_rows(obs: ObservedReads) -> SpectrumRows:
    """The spectrum on [1 / T, f_e], and its level at the edge (mean over the top 10%).

    Values of zero range have no spectrum: any power computed from them is rounding. The
    range is read on the reads themselves, before repeated times are averaged, because an
    average of equal values can differ from them in the last bit.
    """
    _, edge = band_edges(obs)
    if np.ptp(obs.values) == 0.0:
        raise NotComputed("zero variance")
    t, q = collapse_repeated_times(obs.t_rel_s, obs.values)
    span = float(t[-1] - t[0])
    if span * edge < 1.0:
        raise NotComputed("span shorter than 2 median spacings")
    f = spectrum.frequency_grid(span, edge)
    s = spectrum.lomb_scargle_spectrum(t, q, f, obs.median_spacing_s)
    top = f >= 0.9 * edge
    if not np.any(top):
        raise NotComputed("span shorter than 20 median spacings")
    return SpectrumRows(
        frequency_hz=f, spectrum_per_hz=s, edge_level_per_hz=float(np.mean(s[top]))
    )


def noise_floor_ratio(
    obs: ObservedReads, rows: SpectrumRows, sigma: np.ndarray
) -> tuple[float, float]:
    """Q1.2: (floor, floor / S(f_e)), floor = 2 sigma_bar^2 Delta_bar.

    Independent per-read errors of variance sigma_i^2 contribute the flat one-sided density
    2 sigma_bar^2 Delta_bar; a ratio near 1 says the edge is instrument noise, one above 1
    that the reported sigma exceeds what the spectrum shows there. sigma_bar^2 is taken on
    the reads the spectrum sees: repeated times are collapsed, each to the sigma of its mean.

    Validity assumptions:

    - Assumption: the per-read errors are independent, with the reported variances.
      Diagnostic: the self-consistency row (Q3.2) and the dependence panel (Q3.0).
      Consequence of violation: correlated errors are not flat, and a misreported sigma
      moves the floor by the same factor as the ratio.
      Reference: Sera's L0 specification, section 3.2; no equation locator.
    - Assumption: the edge level is the spectrum at the edge.
      Diagnostic: none shipped.
      Consequence of violation: on irregular reads the edge level holds power leaked from
      lower frequencies (`spectrum.g10`), so a ratio below 1 can be leakage, not signal.
      Reference: no source located.
    """
    sigma_c = collapse_repeated_sigma(obs.t_rel_s, sigma)
    floor = 2.0 * mean_sigma_squared(sigma_c) * obs.median_spacing_s
    if not rows.edge_level_per_hz > 0.0:
        raise NotComputed("zero variance")
    return floor, floor / rows.edge_level_per_hz


def second_difference_scale(obs: ObservedReads) -> float:
    """sigma_hat_D2, the Gasser-Sroka-Jennen-Steinmetz noise scale, within segments.

    For each read with a neighbour on both sides in its segment, the pseudo-residual
    e_i = a_i q_(i-1) + b_i q_(i+1) - q_i, with a_i = (t_(i+1) - t_i) / (t_(i+1) - t_(i-1))
    and b_i = 1 - a_i, is zero for a locally linear signal; sigma_hat^2 is the mean of
    e_i^2 / (a_i^2 + b_i^2 + 1). On equal spacing it is (1/6) E[(q_(i+1) - 2 q_i +
    q_(i-1))^2], L0's form. Repeated times are collapsed to their mean first.

    Validity assumptions:

    - Assumption: the signal is locally linear over three consecutive reads, and the
      errors are independent with a common variance.
      Diagnostic: the ratio to sd(q) beside it (Q2.1) and sigma_bar / sigma_hat (Q3.3).
      Consequence of violation: curvature between reads is counted as noise, so the scale
      is an upper bound on the error sd.
      Reference: Gasser, Sroka and Jennen-Steinmetz (1986), Biometrika 73, 625-633, as cited
      in Sera's L0 specification; no equation locator.
    """
    terms: list[np.ndarray] = []
    for s, e in obs.segments():
        t, q = collapse_repeated_times(obs.t_rel_s[s:e], obs.values[s:e])
        if t.size < 3:
            continue
        left, right = t[1:-1] - t[:-2], t[2:] - t[1:-1]
        a = right / (left + right)
        b = 1.0 - a
        resid = a * q[:-2] + b * q[2:] - q[1:-1]
        terms.append(resid * resid / (a * a + b * b + 1.0))
    if not terms:
        raise NotComputed("no segment with at least 3 reads")
    return float(math.sqrt(np.mean(np.concatenate(terms))))


@dataclass(frozen=True, slots=True)
class SecondDifferenceRow:
    second_difference_scale: float
    value_sd: float
    second_difference_ratio: float


def second_difference_row(obs: ObservedReads) -> SecondDifferenceRow:
    """Q2.1: sigma_hat_D2 / sd(q), the share of the variance between adjacent reads."""
    scale = second_difference_scale(obs)
    # Zero range, not sd > 0: equal values can give an sd of rounding, 1e-17 and not 0.
    if np.ptp(obs.values) == 0.0:
        raise NotComputed("zero variance")
    sd = float(np.std(obs.values, ddof=1))
    return SecondDifferenceRow(
        second_difference_scale=scale,
        value_sd=sd,
        second_difference_ratio=scale / sd,
    )


@dataclass(frozen=True, slots=True)
class G10Row:
    g10: float


def g10_row(obs: ObservedReads, rows: SpectrumRows) -> G10Row:
    """Q2.3: G10 on the whole record. Fewer than ten grid frequencies (a span under 20
    median spacings) is the one condition the record cannot answer; any other error from
    `spectrum.g10` is raised."""
    f = rows.frequency_hz
    if not spectrum.has_decade(f):
        raise NotComputed("span shorter than 20 median spacings")
    return G10Row(g10=spectrum.g10(f, rows.spectrum_per_hz))


def g10_halves(obs: ObservedReads) -> tuple[G10Row, G10Row]:
    """Q2.3 split in time: G10 of the first and second half of the reads, each on its own
    span but up to the WHOLE record's band edge, so the halves read the same band. A half
    whose reads have zero range is not computed, as the whole record's spectrum is not."""
    _, edge = band_edges(obs)
    t, q = collapse_repeated_times(obs.t_rel_s, obs.values)
    half = (t.size + 1) // 2
    out = []
    for ts, qs in ((t[:half], q[:half]), (t[half:], q[half:])):
        span = float(ts[-1] - ts[0]) if ts.size > 1 else 0.0
        if span * edge < 1.0:
            raise NotComputed("span shorter than 20 median spacings")
        f = spectrum.frequency_grid(span, edge)
        # The decade test is g10's own, so a half `g10` would refuse never reaches it.
        if not spectrum.has_decade(f):
            raise NotComputed("span shorter than 20 median spacings")
        in_half = (obs.t_rel_s >= ts[0]) & (obs.t_rel_s <= ts[-1])
        if np.ptp(obs.values[in_half]) == 0.0:
            raise NotComputed("zero variance")
        s = spectrum.lomb_scargle_spectrum(ts, qs, f, obs.median_spacing_s)
        out.append(G10Row(g10=spectrum.g10(f, s)))
    return out[0], out[1]


def clipped_fraction(values: np.ndarray, floor: float | None) -> float:
    """The fraction of reads at the metric's floor, the value it is clipped at."""
    if floor is None:
        raise NotComputed("no floor for this metric")
    q = np.asarray(values, dtype=float)
    if q.size == 0:
        raise NotComputed("fewer than 1 reads")
    return float(np.mean(q <= floor))


# ---------------------------------------------------------------------------
# Per metric: Q3.0 the dependence panel, Q3.2, Q3.3
# ---------------------------------------------------------------------------

# Every member is oriented so that a larger value is stronger dependence.
DEPENDENCE_MEMBERS = (
    "pearson_abs",
    "spearman_abs",
    "chatterjee_xi",
    "distance_correlation",
    "decile_spread",
)


@dataclass(frozen=True, slots=True)
class DependencePanel:
    # Each member's value on the observed pairs.
    statistics: dict[str, float]
    member_p_values: dict[str, float]
    p_value: float
    argmin_labels: tuple[str, ...]
    n_permutations: int


def paired_values(
    values: np.ndarray, sigma: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """The (q, sigma) pairs whose sigma is finite. A negative sigma is refused: it is no
    error bar, and the decile spread divides by the mean sigma."""
    q = np.asarray(values, dtype=float)
    s = np.asarray(sigma, dtype=float)
    keep = np.isfinite(q) & np.isfinite(s)
    if np.any(s[keep] < 0.0):
        raise ValueError("a reported sigma is negative")
    if not np.any(keep):
        raise NotComputed("no read with a finite sigma")
    return q[keep], s[keep]


def _decile_of(q: np.ndarray) -> np.ndarray:
    """The decile of each read in the ranking of q, 0 to 9; ties are split by position."""
    rank = np.empty(q.size, dtype=int)
    rank[np.argsort(q, kind="stable")] = np.arange(q.size)
    return (10 * rank) // q.size


def dependence_statistics(
    q: np.ndarray, s: np.ndarray, permutations: np.ndarray
) -> dict[str, np.ndarray]:
    """Every member on the observed pairs and on each row of `permutations`.

    Row 0 of the result is the observed pairing; row b + 1 pairs q with
    s[permutations[b]]. Each member's array has B + 1 values, ready for `min_p`.

    - `pearson_abs`, `spearman_abs`: |r| of q and sigma, of their values and of their
      ranks.
    - `chatterjee_xi`: xi with q as the explanatory variable and sigma as the response, so
      it measures how far sigma is a function of q (`scipy.stats.chatterjeexi(q, sigma)`).
    - `distance_correlation`: `dcor.distance_correlation`, its O(n log n) merge-sort method.
    - `decile_spread`: (largest - smallest mean sigma over the ten deciles of q) / mean sigma.
      Its decile means are sums of floats, so two rows that tie as sets of reads can differ
      in the last bit and a tie is counted against only up to rounding; a continuous sigma,
      as the records report, does not tie.
    """
    n = q.size
    rows = np.vstack([np.arange(n), np.asarray(permutations, dtype=int)])
    s_rows = s[rows]

    def abs_corr(x: np.ndarray, y_rows: np.ndarray) -> np.ndarray:
        xc = x - x.mean()
        yc = y_rows - y_rows.mean(axis=1, keepdims=True)
        return np.abs(yc @ xc) / (np.linalg.norm(xc) * np.linalg.norm(yc, axis=1))

    s_ranks = stats.rankdata(s)
    decile = _decile_of(q)
    counts = np.bincount(decile, minlength=10)
    decile_means = np.vstack(
        [np.bincount(decile, weights=row, minlength=10) / counts for row in s_rows]
    )
    return {
        "pearson_abs": abs_corr(q, s_rows),
        "spearman_abs": abs_corr(stats.rankdata(q), s_ranks[rows]),
        "chatterjee_xi": np.array(
            [float(stats.chatterjeexi(q, row).statistic) for row in s_rows]
        ),
        "distance_correlation": np.array(
            [
                float(
                    dcor.distance_correlation(
                        q, row, method=dcor.DistanceCovarianceMethod.MERGESORT
                    )
                )
                for row in s_rows
            ]
        ),
        "decile_spread": (decile_means.max(axis=1) - decile_means.min(axis=1))
        / s.mean(),
    }


def dependence_panel(
    values: np.ndarray, sigma: np.ndarray, *, n_permutations: int, seed: int
) -> DependencePanel:
    """Q3.0: is the reported sigma independent of the value it accompanies?

    Each null row shuffles sigma against q, which destroys any dependence and keeps both
    marginals, and every member is recomputed on the same rows; `min_p` calibrates the
    smallest member p over the panel. Nothing here decides anything: every member's value
    and p ship, with the panel's p and the members that attain it. `seed` is required: a
    permutation p that changes from run to run while the run identity stays the same is
    not reproducible.

    Validity assumptions:

    - Assumption: under the null the (q, sigma) pairs are exchangeable, so shuffling sigma
      against q leaves their joint law unchanged.
      Diagnostic: none here. Both q and sigma are time series, and a shared slow drift
      makes them dependent without either being a function of the other.
      Consequence of violation: the panel's p is too small, reporting dependence that is
      only shared drift.
      Reference: Westfall and Young (1993), Resampling-Based Multiple Testing, as Sera's L0
      specification cites it for the panel; no section locator located.
    - Assumption: each member measures what its name says.
      Diagnostic: the known-answer tests of the members.
      Consequence of violation: the panel's arg-min names the wrong kind of dependence.
      Reference: Chatterjee (2021), Journal of the American Statistical Association 116,
      2009-2022, for xi; Szekely, Rizzo and Bakirov (2007), The Annals of Statistics 35,
      2769-2794, for the distance correlation; no equation locator.
    """
    if n_permutations < 1:
        raise ValueError(f"n_permutations must be >= 1, got {n_permutations}")
    q, s = paired_values(values, sigma)
    if q.size < 10:
        raise NotComputed("fewer than 10 reads")
    if np.ptp(q) == 0.0 or np.ptp(s) == 0.0:
        raise NotComputed("zero variance")
    rng = np.random.default_rng(seed)
    permutations = np.vstack([rng.permutation(q.size) for _ in range(n_permutations)])
    table = dependence_statistics(q, s, permutations)
    result = min_p(table)
    return DependencePanel(
        statistics={name: float(table[name][0]) for name in DEPENDENCE_MEMBERS},
        member_p_values=result.member_p_values,
        p_value=result.p_value,
        argmin_labels=result.argmin_labels,
        n_permutations=result.n_permutations,
    )


def self_consistency(values: np.ndarray, sigma: np.ndarray) -> float:
    """Q3.2: (v - mean(sigma^2)) / v, with v the observed variance of q.

    v is the variance with ddof = 1 over every read with a finite value; mean(sigma^2) is
    over the reads whose sigma is finite. Independent additive noise can only add
    variance, so a negative value says the reported error bars are larger than the scatter
    they describe. It is a statement about the reporting, not a measurement of the noise.

    Validity assumptions:

    - Assumption: sigma is independent of the value it accompanies, and the noise is
      additive.
      Diagnostic: the dependence panel (Q3.0).
      Consequence of violation: the row is computed but not interpretable.
      Reference: Sera's L0 specification, section 5.3, which makes the row conditional on
      the independence of sigma (section 5.2); no equation locator.
    """
    q = np.asarray(values, dtype=float)
    q = q[np.isfinite(q)]
    if q.size < 2:
        raise NotComputed("fewer than 2 reads")
    if np.ptp(q) == 0.0:
        raise NotComputed("zero variance")
    v = float(np.var(q, ddof=1))
    return (v - mean_sigma_squared(sigma)) / v


def sigma_bound_ratio(obs: ObservedReads, sigma: np.ndarray) -> float:
    """Q3.3: sigma_bar / sigma_hat_D2, with sigma_bar^2 = mean(sigma_i^2).

    sigma_hat_D2 counts curvature between reads as noise, so it bounds the error sd from
    above; a ratio above 1 says the reported sigma exceeds that bound. Both sides are taken
    on the reads with repeated times collapsed, each to the sigma of its mean.

    Validity assumptions:

    - Assumption: sigma is independent of the value it accompanies.
      Diagnostic: the dependence panel (Q3.0).
      Consequence of violation: the row is computed but not interpretable.
      Reference: Sera's L0 specification, section 5.4, which makes the row conditional on
      the independence of sigma (section 5.2); no equation locator.
    """
    bound = second_difference_scale(obs)
    if np.ptp(obs.values) == 0.0 or not bound > 0.0:
        raise NotComputed("zero variance")
    sigma_c = collapse_repeated_sigma(obs.t_rel_s, sigma)
    return math.sqrt(mean_sigma_squared(sigma_c)) / bound


# ---------------------------------------------------------------------------
# Per threshold: Q2.2 the count exponent
# ---------------------------------------------------------------------------

# L0's ladder of decimation factors.
DECIMATION_LADDER = (1, 2, 3, 5, 8, 13, 20)


@dataclass(frozen=True, slots=True)
class CountExponentRow:
    decimation_ladder: tuple[int, ...]
    # n_avg(m) and P_avg(m), averaged over the m offsets, one per rung.
    down_crossings: tuple[float, ...]
    kept_pairs: tuple[float, ...]
    # [n(m) / P(m)] / [n(1) / P(1)]; None on a rung that keeps no pair.
    count_ratio: tuple[float | None, ...]
    count_exponent: float
    # The rungs left out of the fit: no kept pair, or no down-crossing.
    rungs_excluded: tuple[int, ...]


def count_exponent_row(
    obs: ObservedReads, in_spec: np.ndarray, ladder: tuple[int, ...] = DECIMATION_LADDER
) -> CountExponentRow:
    """Q2.2: how the down-crossing count falls as the reads are thinned.

    Thinning by m keeps every m-th read of each segment, starting again in every segment,
    from each of the m offsets in turn. Over the m offsets, the kept consecutive pairs are
    exactly the within-segment pairs at lag m, so n_avg(m) = (lag-m pairs that go from in
    spec to out of spec) / m and P_avg(m) = sum over segments of (L_s - min(m, L_s)) / m.
    gamma = 1 - the least-squares slope of log(n_avg / P_avg) on log m.

    Under independent identically distributed reads E[n] / P is the same at every m, on any
    layout, so gamma is 1 to first order: the fit is on the logs of the counts, not of
    their expectations. A smooth signal whose crossings land at random points of the
    within-segment intervals, more than the largest m reads apart, gives 0 in the same
    sense. When every crossing lies at least the largest m reads from both ends of its
    segment, n_avg(m) is the number of crossings and gamma is the slope of
    log sum_s (L_s - min(m, L_s)) on log m, which is negative; a crossing nearer a segment
    end is kept by fewer offsets, which pulls gamma toward 1. On a gapless record gamma
    equals L0's count exponent up to O(m / N), in expectation.

    Validity assumptions:

    - Assumption: log(n / P) is linear in log m over the ladder.
      Diagnostic: `count_ratio` per rung, beside the fitted exponent.
      Consequence of violation: gamma is an average slope over the ladder.
      Reference: Sera's L0 specification, Q2.2, for the count exponent; the averaging over
      offsets follows subsampling-and-averaging, Ait-Sahalia, Mykland and Zhang (2005), as
      listed in L0's references; no equation locator.
    """
    rungs = tuple(int(m) for m in ladder)
    if not rungs or rungs[0] != 1 or any(b <= a for a, b in zip(rungs, rungs[1:])):
        raise ValueError(
            f"the ladder must start at 1 and increase strictly, got {tuple(ladder)}"
        )
    inside = np.asarray(in_spec, dtype=bool)
    lengths = [e - s for s, e in obs.segments()]
    n_avg, p_avg = [], []
    for m in ladder:
        down = 0
        for s, e in obs.segments():
            if e - s > m:
                seg = inside[s:e]
                down += int(np.count_nonzero(seg[:-m] & ~seg[m:]))
        n_avg.append(down / m)
        p_avg.append(sum(length - min(m, length) for length in lengths) / m)
    if not (p_avg[0] > 0.0 and n_avg[0] > 0.0):
        raise NotComputed("no down-crossing at this threshold")
    base = n_avg[0] / p_avg[0]
    ratio = tuple(n / p / base if p > 0.0 else None for n, p in zip(n_avg, p_avg))
    fit = [i for i, (n, p) in enumerate(zip(n_avg, p_avg)) if n > 0.0 and p > 0.0]
    if len(fit) < 2:
        raise NotComputed("fewer than two rungs with a down-crossing")
    x = np.log([ladder[i] for i in fit])
    y = np.log([n_avg[i] / p_avg[i] for i in fit])
    slope = float(np.polyfit(x, y, 1)[0])
    return CountExponentRow(
        decimation_ladder=tuple(ladder),
        down_crossings=tuple(n_avg),
        kept_pairs=tuple(p_avg),
        count_ratio=ratio,
        count_exponent=1.0 - slope,
        rungs_excluded=tuple(ladder[i] for i in range(len(ladder)) if i not in fit),
    )


# ---------------------------------------------------------------------------
# Per threshold: Q3.4, Q4.1, the in-spec fraction bounds, rho_k
# ---------------------------------------------------------------------------

CONFIDENCE_LEVELS = (1, 2)


@dataclass(frozen=True, slots=True)
class CrossingCounts:
    k: float
    clean: int
    merged: int
    unbracketed: int


def crossing_counts(
    obs: ObservedReads, margin_v: np.ndarray, sigma: np.ndarray, k: float
) -> CrossingCounts:
    """Q3.4 at one level k: how many down-crossings the error bars can place.

    A raw down-crossing is a step, inside a segment, from a read in spec to one out of
    spec. Walking the segment on confident reads only (`windows.confident_mask`: sigma
    finite and |margin| >= k sigma), its bracket is the nearest confident read at or before
    it and the nearest one at or after it, in the same segment. The crossing is CLEAN when
    that bracket runs from a confident in-spec read to a confident out-of-spec one and holds
    exactly one raw crossing of either direction; MERGED when such a bracket holds more;
    UNBRACKETED otherwise: a segment edge or a gap comes first, or the bracket stays on one
    side or runs from out of spec to in spec, and the data does not assert the crossing. It is a reading of resolvability, not an estimator:
    the bracket's ends depend on the process, so it is not Turnbull's interval.

    Validity assumptions:

    - Assumption: the reported sigma is a Gaussian 1-sigma, so k = 1 and 2 read as 68% and
      95% confidence.
      Diagnostic: the dependence panel (Q3.0) and the self-consistency row (Q3.2).
      Consequence of violation: the levels are nominal.
      Reference: Sera's L0 specification, Q3.4; no equation locator.
    """
    m = np.asarray(margin_v, dtype=float)
    if not np.any(np.isfinite(np.asarray(sigma, dtype=float))):
        raise NotComputed("no read with a finite sigma")
    inside = m >= 0.0
    confident = confident_mask(m, sigma, k)
    clean = merged = unbracketed = 0
    for s, e in obs.segments():
        seg_in = inside[s:e]
        seg_conf = np.flatnonzero(confident[s:e])
        changes = np.flatnonzero(seg_in[:-1] != seg_in[1:])
        for i in np.flatnonzero(seg_in[:-1] & ~seg_in[1:]):
            before = seg_conf[seg_conf <= i]
            after = seg_conf[seg_conf >= i + 1]
            if before.size == 0 or after.size == 0:
                unbracketed += 1
                continue
            p, q = int(before[-1]), int(after[0])
            if not (seg_in[p] and not seg_in[q]):
                unbracketed += 1
            elif np.count_nonzero((changes >= p) & (changes < q)) == 1:
                clean += 1
            else:
                merged += 1
    return CrossingCounts(k=k, clean=clean, merged=merged, unbracketed=unbracketed)


@dataclass(frozen=True, slots=True)
class RunLengthRow:
    read_fraction_in_spec: float
    mean_in_spec_run_reads: float
    run_length_ratio: float


def run_length_row(windows_in_spec: pd.DataFrame, in_spec: np.ndarray) -> RunLengthRow:
    """Q4.1: do the in-spec durations say more than the in-spec fraction?

    The mean in-spec run in reads is read from the carve's in-spec window table for this
    threshold: sum(n_reads - censored) / sum(1 - censored), the reads that have an observed
    successor in their segment over the observed in-to-out transitions. Runs cut by a gap
    or by an end of the record count their reads and no transition. A is the per-read
    fraction in spec, and the ratio mean run x (1 - A) is 1 when the reads are independent:
    then run lengths are geometric with mean 1 / (1 - A).

    Validity assumptions:

    - Assumption: the in/out sequence is stationary.
      Diagnostic: rho_k beside it.
      Consequence of violation: the ratio averages runs of different regimes.
      Reference: Feller, An Introduction to Probability Theory and Its Applications, Vol. 1,
      on runs, as cited in Sera's L0 specification; no section locator located.
    """
    inside = np.asarray(in_spec, dtype=bool)
    if inside.size == 0:
        raise NotComputed("fewer than 1 reads")
    a = float(inside.mean())
    censored = windows_in_spec["censored"].to_numpy(dtype=bool)
    n_reads = windows_in_spec["n_reads"].to_numpy(dtype=int)
    exits = int(np.count_nonzero(~censored))
    if exits == 0:
        raise NotComputed("no in-to-out transition within a segment")
    mean_run = float(np.sum(n_reads - censored)) / exits
    return RunLengthRow(
        read_fraction_in_spec=a,
        mean_in_spec_run_reads=mean_run,
        run_length_ratio=mean_run * (1.0 - a),
    )


@dataclass(frozen=True, slots=True)
class InSpecBounds:
    in_spec_fraction_lower: float
    in_spec_fraction_upper: float
    unsure_k: float


def in_spec_fraction_bounds(
    obs: ObservedReads, margin_v: np.ndarray, sigma: np.ndarray, k: float
) -> InSpecBounds:
    """Lower and upper bounds on the time-weighted in-spec fraction.

    Each within-segment interval [t_i, t_(i+1)) is held by its left read, the convention of
    the carve's durations. An interval is UNSETTLED when it holds a crossing (the margin
    changes sign across it) or its left read is not sure at level k: its sigma is unknown,
    or |margin| < k sigma (`windows.confident_mask`). The lower bound counts every unsettled
    interval out of spec, the upper bound counts it in spec; observed time is the sum of
    the within-segment spacings, so no gap enters. The carve's own in-spec fraction lies
    between the two. A record with no finite sigma has no bounds.

    Validity assumptions:

    - Assumption: the reported sigma is a Gaussian 1-sigma, so k reads as a confidence
      level.
      Diagnostic: the dependence panel (Q3.0) and the self-consistency row (Q3.2).
      Consequence of violation: the level is nominal.
      Reference: no source located; the bounds are not an L0 row.
    - Assumption: no excursion to the other side happens between two reads on the same
      side.
      Diagnostic: rho_k and the count exponent beside it.
      Consequence of violation: the bounds bound the fraction of time held by reads in
      spec, not the true time-weighted fraction; an excursion between reads lies outside
      both.
      Reference: no source located.
    """
    m = np.asarray(margin_v, dtype=float)
    within = obs.within_segment()
    delta = np.diff(obs.t_rel_s)[within]
    observed = float(np.sum(delta))
    if not observed > 0.0:
        raise NotComputed("no segment with at least 2 reads")
    s = np.asarray(sigma, dtype=float)
    if not np.any(np.isfinite(s)):
        raise NotComputed("no read with a finite sigma")
    inside = m >= 0.0
    unsure = ~confident_mask(m, s, k)
    left_in, right_in, left_unsure = inside[:-1], inside[1:], unsure[:-1]
    unsettled = ((left_in != right_in) | left_unsure)[within]
    settled_in = (left_in[within]) & ~unsettled
    lower = float(np.sum(delta[settled_in])) / observed
    upper = lower + float(np.sum(delta[unsettled])) / observed
    return InSpecBounds(
        in_spec_fraction_lower=lower, in_spec_fraction_upper=upper, unsure_k=float(k)
    )


MAX_AUTOCORR_LAG = 10


def in_spec_autocorr(obs: ObservedReads, in_spec: np.ndarray) -> tuple[float, ...]:
    """rho_k of the in-spec indicator, k = 1, 2, ..., up to 10, within segments.

    Centred on the record's per-read fraction A and scaled by A (1 - A); the numerator is
    the mean over the within-segment pairs at lag k. The array stops at the last lag that
    has such a pair, so it never carries an undefined value.

    Validity assumptions:

    - Assumption: the in-spec indicator is stationary across the record, so one A centres
      every segment.
      Diagnostic: the split-half G10 and the run-length row beside it.
      Consequence of violation: segments at different in-spec fractions add a shared
      offset, which raises rho_k at every lag.
      Reference: no source located.
    """
    x = np.asarray(in_spec, dtype=float)
    a = float(x.mean()) if x.size else 0.0
    variance = a * (1.0 - a)
    if not variance > 0.0:
        raise NotComputed("zero variance")
    centred = x - a
    out = []
    for lag in range(1, MAX_AUTOCORR_LAG + 1):
        products = [
            centred[s : e - lag] * centred[s + lag : e]
            for s, e in obs.segments()
            if e - s > lag
        ]
        if not products:
            break
        out.append(float(np.mean(np.concatenate(products))) / variance)
    if not out:
        raise NotComputed("no segment with at least 2 reads")
    return tuple(out)


# ---------------------------------------------------------------------------
# The step: inputs, result, run
# ---------------------------------------------------------------------------

# How the dataset declares sigma: given and to be used, absent, or known but not usable.
SIGMA_DECLARED = "declared"
# The status a declared sigma takes when no read carries a finite one: decided once, by
# `run`, so every row that reads sigma says the same.
NO_FINITE_SIGMA = "not computed: no read with a finite sigma"
SIGMA_STATUSES = (
    SIGMA_DECLARED,
    NOT_DECLARED,
    "not computed: sigma propagation not derived",
    NO_FINITE_SIGMA,
)

# Each field group ships with one status, and its fields are all None unless it is computed.
METRIC_GROUPS: dict[str, tuple[str, ...]] = {
    "spacing": (
        "median_spacing_s",
        "gap_threshold_s",
        "effective_spacing_s",
        "mean_spacing_s",
        "spacing_cv",
    ),
    "window": ("first_alias_hz", "band_edge_hz", "spectral_window_at_first_alias"),
    "spectrum": ("spectrum_frequency_hz", "spectrum_per_hz", "edge_level_per_hz"),
    "noise_floor": ("noise_floor_per_hz", "edge_noise_ratio"),
    "second_difference": (
        "second_difference_scale",
        "value_sd",
        "second_difference_ratio",
    ),
    "g10": ("g10",),
    "g10_halves": ("g10_first_half", "g10_second_half"),
    "dependence": ("dependence",),
    "self_consistency": ("self_consistency",),
    "sigma_bound": ("sigma_bound_ratio",),
    "clipped": ("clipped_fraction",),
}
THRESHOLD_GROUPS: dict[str, tuple[str, ...]] = {
    "count_exponent": (
        "decimation_ladder",
        "down_crossings",
        "kept_pairs",
        "count_ratio",
        "count_exponent",
        "rungs_excluded",
    ),
    "crossing_counts": ("crossing_counts", "clean_fraction_k1"),
    "run_length": (
        "read_fraction_in_spec",
        "mean_in_spec_run_reads",
        "run_length_ratio",
    ),
    "bounds": ("in_spec_fraction_lower", "in_spec_fraction_upper", "unsure_k"),
    "autocorr": ("in_spec_autocorr",),
}
# The groups that read sigma; they carry the sigma status when sigma is not declared.
SIGMA_GROUPS = frozenset(
    {
        "noise_floor",
        "dependence",
        "self_consistency",
        "sigma_bound",
        "crossing_counts",
        "bounds",
    }
)


def check_sigma_status(status: str) -> str:
    """One of `SIGMA_STATUSES`; any other status, a reason from the set included, raises."""
    if status not in SIGMA_STATUSES:
        raise ValueError(
            f"unknown sigma status {status!r}; the set is {SIGMA_STATUSES}"
        )
    return status


def _finite(value: object) -> bool:
    """Every number inside `value` is finite. None inside a tuple marks an entry its row
    defines as absent, as a rung of `count_ratio` that keeps no pair."""
    if isinstance(value, (str, bool, np.bool_)):
        return True
    if isinstance(value, (int, float, np.integer, np.floating)):
        return bool(np.isfinite(value))
    if isinstance(value, np.ndarray):
        return bool(np.all(np.isfinite(value)))
    if isinstance(value, (tuple, list)):
        return all(item is None or _finite(item) for item in value)
    if isinstance(value, dict):
        return all(_finite(item) for item in value.values())
    if dataclasses.is_dataclass(value):
        return all(_finite(getattr(value, f.name)) for f in dataclasses.fields(value))
    raise TypeError(f"no finiteness rule for {type(value).__name__}")


def _check_statuses(row: object, groups: Mapping[str, tuple[str, ...]]) -> None:
    """A row's statuses name exactly its groups, each from the vocabulary; a computed group
    has every field set and finite, any other has every field None."""
    status: dict[str, str] = getattr(row, "status")
    if set(status) != set(groups):
        raise ValueError(
            f"statuses {sorted(status)} do not name the groups {sorted(groups)}"
        )
    for group, names in groups.items():
        word = status[group]
        if word not in (COMPUTED, NOT_DECLARED):
            if not word.startswith("not computed: "):
                raise ValueError(f"group {group!r} has an unknown status {word!r}")
            check_reason(word[len("not computed: ") :])
        values = [getattr(row, name) for name in names]
        if word == COMPUTED:
            bad = [n for n, v in zip(names, values) if v is None or not _finite(v)]
            if bad:
                raise ValueError(
                    f"group {group!r} is computed, but {bad} are missing or not finite"
                )
        elif any(v is not None for v in values):
            raise ValueError(f"group {group!r} is {word!r} but carries a value")


@dataclass(frozen=True, slots=True)
class TellabilityInputs:
    dataset_id: str
    metric_key: str
    value_key: str
    t_rel_s: np.ndarray
    values: np.ndarray
    sigma: np.ndarray | None
    sigma_status: str
    # (label, value, big_values_good), carved by `windows` into the `windows` table.
    thresholds: tuple[tuple[str, float, bool], ...]
    windows: WindowsResult | None
    unsure_k: float
    gap_mult: float
    clip_floor: float | None
    n_permutations: int
    seed: int
    # Rows of the file the load dropped before any step ran, keyed by LOAD_COUNT_KEYS.
    load_counts: Mapping[str, int] = field(default_factory=dict)


@dataclass
class MetricTellability(StaleArtifactGuard):
    metric_key: str
    value_key: str
    sigma_status: str
    n_reads: int
    n_segments: int
    # Reads dropped for a non-finite time or value, and reads kept with no finite sigma
    # (None when no sigma is declared).
    n_reads_dropped_nonfinite: int
    n_reads_sigma_unknown: int | None
    # Rows the load dropped before this step (core.types.LOAD_COUNT_KEYS); None when the
    # load did not record the count.
    n_rows_dropped_by_schema: int | None
    n_rows_dropped_invalid: int | None
    median_spacing_s: float | None
    gap_threshold_s: float | None
    effective_spacing_s: float | None
    mean_spacing_s: float | None
    spacing_cv: float | None
    first_alias_hz: float | None
    band_edge_hz: float | None
    spectral_window_at_first_alias: float | None
    spectrum_frequency_hz: np.ndarray | None
    spectrum_per_hz: np.ndarray | None
    edge_level_per_hz: float | None
    noise_floor_per_hz: float | None
    edge_noise_ratio: float | None
    second_difference_scale: float | None
    value_sd: float | None
    second_difference_ratio: float | None
    g10: float | None
    g10_first_half: float | None
    g10_second_half: float | None
    dependence: DependencePanel | None
    self_consistency: float | None
    sigma_bound_ratio: float | None
    clipped_fraction: float | None
    # One status per field group, keyed as in METRIC_GROUPS.
    status: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        _check_statuses(self, METRIC_GROUPS)


@dataclass
class ThresholdTellability(StaleArtifactGuard):
    threshold_label: str
    threshold_value: float
    big_values_good: bool
    decimation_ladder: tuple[int, ...] | None
    down_crossings: tuple[float, ...] | None
    kept_pairs: tuple[float, ...] | None
    count_ratio: tuple[float | None, ...] | None
    count_exponent: float | None
    rungs_excluded: tuple[int, ...] | None
    crossing_counts: tuple[CrossingCounts, ...] | None
    clean_fraction_k1: float | None
    read_fraction_in_spec: float | None
    mean_in_spec_run_reads: float | None
    run_length_ratio: float | None
    in_spec_fraction_lower: float | None
    in_spec_fraction_upper: float | None
    unsure_k: float | None
    in_spec_autocorr: tuple[float, ...] | None
    # One status per field group, keyed as in THRESHOLD_GROUPS.
    status: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        _check_statuses(self, THRESHOLD_GROUPS)


@dataclass
class TellabilityResult(StaleArtifactGuard):
    dataset_id: str
    metric: MetricTellability
    thresholds: tuple[ThresholdTellability, ...]
    # The settings that made it, so the artifact reads without its job file.
    settings: dict[str, object] = field(default_factory=dict)


T = TypeVar("T")


def _attempt(status: dict[str, str], group: str, compute: Callable[[], T]) -> T | None:
    """Run one group's computation; record `computed` or the reason it was not.

    The contract with the row functions: `NotComputed` means the record cannot answer the
    row, and becomes its status. Any other exception is a contradiction in the inputs or a
    bug, and is raised; turning it into a status would ship a wrong but plausible row.
    """
    try:
        value = compute()
    except NotComputed as exc:
        status[group] = exc.status
        return None
    status[group] = COMPUTED
    return value


def _check_inputs(inputs: TellabilityInputs) -> None:
    check_sigma_status(inputs.sigma_status)
    if not (math.isfinite(inputs.unsure_k) and inputs.unsure_k > 0.0):
        raise ValueError(
            f"unsure_k must be finite and positive, got {inputs.unsure_k!r}"
        )
    if int(inputs.n_permutations) < 1:
        raise ValueError(f"n_permutations must be >= 1, got {inputs.n_permutations}")
    if int(inputs.seed) < 0:
        raise ValueError(f"seed must be >= 0, got {inputs.seed}")
    unknown = sorted(set(inputs.load_counts) - set(LOAD_COUNT_KEYS))
    if unknown:
        raise ValueError(
            f"unknown load counts {unknown}; the keys are {LOAD_COUNT_KEYS}"
        )
    if inputs.sigma is not None:
        s = np.asarray(inputs.sigma, dtype=float)
        if np.any(s[np.isfinite(s)] < 0.0):
            raise ValueError("a reported sigma is negative")
    labels = [label for label, _, _ in inputs.thresholds]
    if len(set(labels)) != len(labels):
        raise ValueError(f"a threshold label is repeated: {labels}")
    if (inputs.sigma is not None) != (inputs.sigma_status == SIGMA_DECLARED):
        raise ValueError(
            f"sigma is {'given' if inputs.sigma is not None else 'absent'} but its status "
            f"is {inputs.sigma_status!r}; a sigma array goes with 'declared' and only it"
        )
    if not inputs.thresholds:
        return
    if inputs.windows is None:
        raise ValueError("thresholds were given without the carve's windows for them")
    carved_gap_mult = float(inputs.windows.meta["gap_mult"])
    if carved_gap_mult != float(inputs.gap_mult):
        raise ValueError(
            f"the windows were carved at gap_mult {carved_gap_mult} but tellability runs at "
            f"{inputs.gap_mult}; one record, one gap rule"
        )
    carved = set(inputs.windows.windows_in_spec["threshold_label"]) | set(
        inputs.windows.diagnostics.get("per_threshold", {})
    )
    missing = sorted(label for label, _, _ in inputs.thresholds if label not in carved)
    if missing:
        raise ValueError(f"the windows were not carved at threshold(s) {missing}")
    # A label is only a name: the value and the side it was carved at must match too.
    tables = (inputs.windows.windows_in_spec, inputs.windows.windows_out_of_spec)
    for label, value, big_values_good in inputs.thresholds:
        for table in tables:
            rows = table[table["threshold_label"] == label]
            if len(rows) and (
                float(rows["threshold_value"].iloc[0]) != float(value)
                or bool(rows["big_values_good"].iloc[0]) != bool(big_values_good)
            ):
                raise ValueError(
                    f"threshold {label!r} was carved at "
                    f"({rows['threshold_value'].iloc[0]}, "
                    f"{rows['big_values_good'].iloc[0]}), not at ({value}, "
                    f"{big_values_good})"
                )
    carved_k = float(inputs.windows.meta["k"])
    if carved_k != float(inputs.unsure_k):
        raise ValueError(
            f"the windows were carved at k {carved_k} but tellability reads unsure reads "
            f"at {inputs.unsure_k}; one record, one confidence level"
        )


def _check_windows_reads(inputs: TellabilityInputs, obs: ObservedReads) -> None:
    """The windows must have been carved on these reads: as many, at the same median
    spacing. Windows of another record, or of part of this one, are refused."""
    if inputs.windows is None or not inputs.thresholds:
        return
    diagnostics = inputs.windows.diagnostics
    n_carved = int(diagnostics["n_reads_finite"])
    median = float(diagnostics["median_spacing_s"])
    same = median == obs.median_spacing_s or (
        math.isnan(median) and math.isnan(obs.median_spacing_s)
    )
    if n_carved != obs.t_rel_s.size or not same:
        raise ValueError(
            f"the windows were carved on {n_carved} reads at a median spacing of "
            f"{median} s, but tellability reads {obs.t_rel_s.size} at "
            f"{obs.median_spacing_s} s; one record, one set of reads"
        )


def _load_count(inputs: TellabilityInputs, key: str) -> int | None:
    return int(inputs.load_counts[key]) if key in inputs.load_counts else None


def _metric_rows(
    inputs: TellabilityInputs,
    obs: ObservedReads,
    sigma: np.ndarray | None,
    sigma_status: str,
) -> MetricTellability:
    st: dict[str, str] = {}
    declared = sigma_status == SIGMA_DECLARED

    def sigma_group(group: str, compute: Callable[[], T]) -> T | None:
        if not declared:
            st[group] = sigma_status
            return None
        return _attempt(st, group, compute)

    spacing = _attempt(st, "spacing", lambda: spacing_row(obs))
    edges = _attempt(st, "window", lambda: band_edges(obs))
    window = window_at_first_alias(obs) if edges is not None else None
    rows = _attempt(st, "spectrum", lambda: spectrum_rows(obs))
    if rows is None:
        # A sigma row says first that it has no sigma, and only then that it has no spectrum.
        st["noise_floor"] = st["spectrum"] if declared else sigma_status
        st["g10"] = st["spectrum"]
        floor = g10 = None
    else:
        floor = sigma_group("noise_floor", lambda: noise_floor_ratio(obs, rows, sigma))
        g10 = _attempt(st, "g10", lambda: g10_row(obs, rows))
    second = _attempt(st, "second_difference", lambda: second_difference_row(obs))
    halves = _attempt(st, "g10_halves", lambda: g10_halves(obs))
    panel = sigma_group(
        "dependence",
        lambda: dependence_panel(
            obs.values, sigma, n_permutations=inputs.n_permutations, seed=inputs.seed
        ),
    )
    consistency = sigma_group(
        "self_consistency", lambda: self_consistency(obs.values, sigma)
    )
    bound = sigma_group("sigma_bound", lambda: sigma_bound_ratio(obs, sigma))
    clipped = _attempt(
        st, "clipped", lambda: clipped_fraction(obs.values, inputs.clip_floor)
    )
    return MetricTellability(
        metric_key=inputs.metric_key,
        value_key=inputs.value_key,
        sigma_status=sigma_status,
        n_reads=int(obs.t_rel_s.size),
        n_segments=len(obs.segments()),
        n_reads_dropped_nonfinite=int(obs.n_dropped_nonfinite),
        n_reads_sigma_unknown=(
            int(np.count_nonzero(~np.isfinite(sigma))) if sigma is not None else None
        ),
        n_rows_dropped_by_schema=_load_count(inputs, "n_rows_dropped_by_schema"),
        n_rows_dropped_invalid=_load_count(inputs, "n_rows_dropped_invalid"),
        median_spacing_s=obs.median_spacing_s if spacing else None,
        gap_threshold_s=obs.gap_threshold_s if spacing else None,
        effective_spacing_s=spacing.effective_spacing_s if spacing else None,
        mean_spacing_s=spacing.mean_spacing_s if spacing else None,
        spacing_cv=spacing.spacing_cv if spacing else None,
        first_alias_hz=edges[0] if edges else None,
        band_edge_hz=edges[1] if edges else None,
        spectral_window_at_first_alias=window,
        spectrum_frequency_hz=rows.frequency_hz if rows else None,
        spectrum_per_hz=rows.spectrum_per_hz if rows else None,
        edge_level_per_hz=rows.edge_level_per_hz if rows else None,
        noise_floor_per_hz=floor[0] if floor else None,
        edge_noise_ratio=floor[1] if floor else None,
        second_difference_scale=second.second_difference_scale if second else None,
        value_sd=second.value_sd if second else None,
        second_difference_ratio=second.second_difference_ratio if second else None,
        g10=g10.g10 if g10 else None,
        g10_first_half=halves[0].g10 if halves else None,
        g10_second_half=halves[1].g10 if halves else None,
        dependence=panel,
        self_consistency=consistency,
        sigma_bound_ratio=bound,
        clipped_fraction=clipped,
        status=st,
    )


def _clean_fraction(
    counts: tuple[CrossingCounts, ...],
) -> tuple[tuple[CrossingCounts, ...], float]:
    first = next((c for c in counts if c.k == 1), None)
    if first is None:
        raise ValueError("the crossing counts hold no k = 1 level")
    total = first.clean + first.merged + first.unbracketed
    if total == 0:
        raise NotComputed("no down-crossing at this threshold")
    return counts, first.clean / total


def _threshold_rows(
    inputs: TellabilityInputs,
    obs: ObservedReads,
    sigma: np.ndarray | None,
    windows_in_spec: pd.DataFrame,
    label: str,
    value: float,
    big_values_good: bool,
    sigma_status: str,
) -> ThresholdTellability:
    st: dict[str, str] = {}
    declared = sigma_status == SIGMA_DECLARED
    margin_v = margin(obs.values, value, big_values_good)
    inside = margin_v >= 0.0
    table = windows_in_spec[windows_in_spec["threshold_label"] == label]
    exponent = _attempt(st, "count_exponent", lambda: count_exponent_row(obs, inside))
    if declared:
        crossings = _attempt(
            st,
            "crossing_counts",
            lambda: _clean_fraction(
                tuple(
                    crossing_counts(obs, margin_v, sigma, k) for k in CONFIDENCE_LEVELS
                )
            ),
        )
        bounds = _attempt(
            st,
            "bounds",
            lambda: in_spec_fraction_bounds(obs, margin_v, sigma, inputs.unsure_k),
        )
    else:
        st["crossing_counts"] = st["bounds"] = sigma_status
        crossings = bounds = None
    runs = _attempt(st, "run_length", lambda: run_length_row(table, inside))
    autocorr = _attempt(st, "autocorr", lambda: in_spec_autocorr(obs, inside))
    return ThresholdTellability(
        threshold_label=label,
        threshold_value=float(value),
        big_values_good=bool(big_values_good),
        decimation_ladder=exponent.decimation_ladder if exponent else None,
        down_crossings=exponent.down_crossings if exponent else None,
        kept_pairs=exponent.kept_pairs if exponent else None,
        count_ratio=exponent.count_ratio if exponent else None,
        count_exponent=exponent.count_exponent if exponent else None,
        rungs_excluded=exponent.rungs_excluded if exponent else None,
        crossing_counts=crossings[0] if crossings else None,
        clean_fraction_k1=crossings[1] if crossings else None,
        read_fraction_in_spec=runs.read_fraction_in_spec if runs else None,
        mean_in_spec_run_reads=runs.mean_in_spec_run_reads if runs else None,
        run_length_ratio=runs.run_length_ratio if runs else None,
        in_spec_fraction_lower=bounds.in_spec_fraction_lower if bounds else None,
        in_spec_fraction_upper=bounds.in_spec_fraction_upper if bounds else None,
        unsure_k=bounds.unsure_k if bounds else None,
        in_spec_autocorr=autocorr,
        status=st,
    )


def run(inputs: TellabilityInputs) -> TellabilityResult:
    """Every tellability row of one metric record, each with its status.

    The reads are taken as observed under the carve's gap rule (`windows.observed_reads`,
    at `gap_mult`), and the windows passed in must have been carved under the same rule and
    at every threshold asked for: one record, one gap rule. Raises ValueError on
    contradictory inputs; a row the record cannot answer is NOT an error and ships as `not
    computed: <reason>`.
    """
    _check_inputs(inputs)
    obs = observed_reads(
        inputs.t_rel_s, inputs.values, sigma=inputs.sigma, gap_mult=inputs.gap_mult
    )
    _check_windows_reads(inputs, obs)
    sigma = obs.sigma
    sigma_status = inputs.sigma_status
    if sigma_status == SIGMA_DECLARED and not np.any(np.isfinite(sigma)):
        sigma_status = NO_FINITE_SIGMA
    metric = _metric_rows(inputs, obs, sigma, sigma_status)
    table = inputs.windows.windows_in_spec if inputs.windows is not None else None
    rows = tuple(
        _threshold_rows(inputs, obs, sigma, table, label, value, big, sigma_status)
        for label, value, big in inputs.thresholds
        if table is not None
    )
    return TellabilityResult(
        dataset_id=inputs.dataset_id,
        metric=metric,
        thresholds=rows,
        settings={
            "gap_mult": float(inputs.gap_mult),
            "unsure_k": float(inputs.unsure_k),
            "clip_floor": inputs.clip_floor,
            "n_permutations": int(inputs.n_permutations),
            "seed": int(inputs.seed),
            "decimation_ladder": DECIMATION_LADDER,
            "confidence_levels": CONFIDENCE_LEVELS,
            "max_autocorr_lag": MAX_AUTOCORR_LAG,
        },
    )
