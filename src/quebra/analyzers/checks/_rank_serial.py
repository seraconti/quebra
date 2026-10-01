"""Rank serial autocorrelation: the shared core of C5 and C6.

C5 and C6 ask the same question - does the ORDER of the durations carry information? -
and differ only in how they aggregate across lags, so the lag-wise autocorrelations are
computed once here.

Three design points that the results depend on:

**Ranks, not values.** A duration vector with a heavy right tail lets one long window
dominate a Pearson autocorrelation. Ranks bound each observation's leverage, which matters
at N = 20 where a single outlier is 5% of the sample.

**Pairs never cross a segment boundary.** A lag-1 pair spanning a read gap would relate two
windows separated by unobserved hours, which is exactly the relation the carve refuses to
assert when it ends a window with `gap_start`. Ranks are still computed GLOBALLY (over the
concatenated vector) so that segments are on one comparable scale; only the PAIRING is
within-segment. This keeps the permutation test exact: the multiset of values in each block
is fixed by a within-block permutation, so the multiset of global ranks in each block is
fixed too, and permuting the values is identical to permuting the ranks.

**Studentization uses the permutation ensemble's own moments, and how much it matters was
measured rather than assumed.** Every lag shares one denominator (the full centred sum of
squares) while its numerator sums only the within-segment pairs at that lag, so the null
spread of `r_h` GROWS with the pair count - the opposite of the usual intuition that fewer
observations means more noise. Measured null SD per lag:

    one segment of 60      lags 1-5, pairs 59..55   SD 0.1269..0.1224   ratio 1.04
    10 segments of 5       lags 1-4, pairs 40..10   SD 0.0985..0.0600   ratio 1.64

So on an ungapped record studentization is very nearly a no-op, and the two C5 variants
should agree. It bites when segments are short and numerous - the censored/gapped case -
where an unstudentized `max_h |r_h|` picks lag 1 in 43% of replicates against 21-29% once
each lag is put on its own scale. Shipping both variants is what turns that from a claim
into a number the bench reports.

The moments are computed from the observed value TOGETHER WITH the permutations, not from
the permutations alone: that makes the scaling a symmetric function of the augmented
sample, which is what keeps the resulting permutation test exactly valid at finite B
rather than valid only to O(1/B).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats

# Cap on the number of lags. Beyond a handful of lags the pair count collapses and each
# extra coordinate costs power in the max/portmanteau without plausibly carrying signal.
MAX_LAG_CAP = 5

# A lag needs at least this many within-segment pairs to be included at all.
MIN_PAIRS_PER_LAG = 2


@dataclass(frozen=True)
class LagLayout:
    """Which positions pair with which, per lag, over the concatenated vector."""

    lags: tuple[int, ...]
    left: tuple[np.ndarray, ...]
    right: tuple[np.ndarray, ...]
    n_pairs: tuple[int, ...]


def lag_layout(
    sizes: list[int] | tuple[int, ...], max_lag: int = MAX_LAG_CAP
) -> LagLayout:
    """Within-segment `(i, i+h)` index pairs for each usable lag `h`."""
    sizes = [int(s) for s in sizes]
    lags: list[int] = []
    left: list[np.ndarray] = []
    right: list[np.ndarray] = []
    counts: list[int] = []
    for lag in range(1, int(max_lag) + 1):
        lo_idx: list[np.ndarray] = []
        offset = 0
        for size in sizes:
            if size > lag:
                base = offset + np.arange(size - lag)
                lo_idx.append(base)
            offset += size
        if not lo_idx:
            break
        flat = np.concatenate(lo_idx)
        if len(flat) < MIN_PAIRS_PER_LAG:
            break
        lags.append(lag)
        left.append(flat)
        right.append(flat + lag)
        counts.append(len(flat))
    if not lags:
        raise ValueError(
            f"no usable lag: segment sizes {sizes} give fewer than "
            f"{MIN_PAIRS_PER_LAG} within-segment pairs even at lag 1. A record of "
            "very short segments carries no ordering information to test."
        )
    return LagLayout(
        lags=tuple(lags),
        left=tuple(left),
        right=tuple(right),
        n_pairs=tuple(counts),
    )


def global_ranks(x: np.ndarray) -> np.ndarray:
    """Midrank-averaged ranks. Ties get their average rank, which is what keeps a
    quantised duration vector (many windows exactly one read long) from being ordered
    arbitrarily by input order."""
    x = np.asarray(x, dtype=float)
    ranks = stats.rankdata(x, method="average")
    if float(np.var(ranks)) == 0.0:
        raise ValueError(
            "every duration is tied, so the rank vector is constant and no serial "
            "statistic is defined. This is a degenerate cell, not a result."
        )
    return ranks


def autocorrelations(rank_matrix: np.ndarray, layout: LagLayout) -> np.ndarray:
    """Lag-wise rank autocorrelations for a `(rows, total)` matrix of rank vectors.

    Returns `(rows, n_lags)`. The denominator is the FULL centred sum of squares, shared
    by every lag, so the values are comparable across lags in the usual autocorrelation
    sense; the pair-count differences between lags are then handled by studentization
    (C5) or by explicit weighting (C6).
    """
    rank_matrix = np.atleast_2d(np.asarray(rank_matrix, dtype=float))
    centred = rank_matrix - rank_matrix.mean(axis=1, keepdims=True)
    denominator = (centred**2).sum(axis=1)
    if np.any(denominator <= 0.0):
        raise ValueError("a rank vector is constant; no serial statistic is defined")
    out = np.empty((rank_matrix.shape[0], len(layout.lags)), dtype=float)
    for column, (left, right) in enumerate(zip(layout.left, layout.right)):
        out[:, column] = (centred[:, left] * centred[:, right]).sum(
            axis=1
        ) / denominator
    return out


def studentize(observed: np.ndarray, null: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Centre and scale each lag by the AUGMENTED sample's own moments.

    `observed` is `(n_lags,)`, `null` is `(B, n_lags)`. Returns the studentized observed
    row and the studentized null matrix. Pooling the observed row into the moments is
    what makes the transformation symmetric across the augmented sample - see the module
    docstring.
    """
    observed = np.asarray(observed, dtype=float).reshape(1, -1)
    null = np.asarray(null, dtype=float)
    if null.shape[1] != observed.shape[1]:
        raise ValueError(
            f"lag count mismatch: observed has {observed.shape[1]}, null has "
            f"{null.shape[1]}"
        )
    augmented = np.vstack([observed, null])
    centre = augmented.mean(axis=0, keepdims=True)
    scale = augmented.std(axis=0, keepdims=True)
    if np.any(scale <= 0.0):
        raise ValueError(
            "a lag has zero spread across the permutation ensemble, so it cannot be "
            "studentized. Every permutation gave the same autocorrelation at that lag, "
            "which means the duration vector is effectively constant."
        )
    return ((observed - centre) / scale)[0], (null - centre) / scale
