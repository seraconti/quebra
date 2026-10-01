"""Excursion shape statistics for one in-spec window.

Ported from `monoliths/v13fig/v13_shape.py:119-178`. The question each statistic asks is
the same: inside a window, does the margin above the threshold predict how much longer
the window has left to live? If it does, an excursion has a shape worth modelling; if it
does not, the window is a coin flip and only its duration means anything.

Primary statistic is Chatterjee's xi (Chatterjee 2021, JASA 116(536):2009-2022), which
tends to 1 iff forward time is a function of margin. Spearman `rho2` and `limb_rho` are
kept alongside deliberately: xi is weaker than Spearman against smooth monotone
alternatives at small n, and the windows here are short: over the eligible rows of the
shipped artifacts the median window is 11 reads on 0704 and 14 on 1004. (An earlier draft
said "~7", which is the 3 us cell rather than the median over the eligible set - a pooled
number labelled as though it were the whole.)

SIX DIFFERENCES FROM THE REFERENCE, all deliberate:
1. `rho2_excess` in the monolith guards only `n < 3` and then evaluates `rho ** 2`, while
   `rho2` also guards `isnan(rho)`. Both are guarded here.
2. `chatterjee_xi` on a constant `x` returns a value decided by stable-argsort read
   order - an artefact, not a statistic. Returns NaN here.
3. `dcor` silently returned NaN above 400 reads (its cost is O(n^2) in memory), with no
   way to tell that from "too few reads" or "constant input". Still skipped, but the
   reason is now counted: `dcor_size_skipped`, `dcor_too_short`, `dcor_constant`.
4. Forward time is taken from the read table's `forward_time_s` (time to the recorded
   death) rather than the monolith's time-to-last-in-spec-read. The two differ by a
   constant per window, and every statistic here is rank- or distance-based, so all are
   invariant to that shift. Using the table keeps one definition of death in the repo.
5. The xi p-value is an ADDITION, not a port - the reference has no null model anywhere.
   The calibration is chosen by measuring the ties: tie-free data gets the closed form
   `sqrt(n) * xi -> N(0, 2/5)`, tied data gets a permutation against the independence
   null, because the 2/5 variance is derived for continuous observations only. Which one
   produced a given p-value is recorded in `xi_p_method` beside it. Both are asymptotic or
   exact in READS PER WINDOW, and a window here is short - median 11 reads on 0704, 14 on
   1004 - so the p is reported beside `n` and never alone.

6. The reference `profile` requires 6 reads per window; `shape_curve` requires only 2,
   because eligibility is already decided upstream by `shape_min_reads`.

Known selection effect, reported rather than hidden: `limb_rho` is undefined whenever the
peak sits in the last three reads, and those are exactly the right-skewed windows. Its
median is therefore a median over a left-peaked subset, which is why every median here
ships with its defined-count.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata, spearmanr

DCOR_MAX_READS = 400
XI_NULL_VARIANCE = 2.0 / 5.0

STAT_COLUMNS = [
    "window_index",
    "n_reads",
    "duration_s",
    "peak_frac",
    "rho",
    "rho2",
    "rho2_excess",
    "limb_rho",
    "xi_m_to_f",
    "xi_f_to_m",
    "xi_p_value",
    "xi_p_method",
    "xi_tie_frac",
    "dcor",
]

# Columns of the shape table that are NOT statistics: an identifier and a categorical
# label. Aggregating over them is a category error, and `to_numpy(dtype=float)` on the
# label raises outright - which is how this list came to exist.
NON_NUMERIC_STAT_COLUMNS = frozenset({"window_index", "xi_p_method"})


def chatterjee_xi(x: np.ndarray, y: np.ndarray) -> float:
    """Chatterjee's rank correlation, TIE-CORRECTED form; -> 1 iff y is a function of x.

    Chatterjee (2021), JASA 116(536):2009-2022, in the general form that admits ties,
    called the tie-corrected form throughout this repository. Its page and equation
    locators are recorded in `docs/GOLD_STANDARD.md` section 6, NOT RE-CHECKED:

        xi = 1 - n * sum_i |r_{i+1} - r_i| / (2 * sum_i l_i (n - l_i))

    with the pairs sorted by x, `r_i = #{j : y_j <= y_(i)}` and `l_i = #{j : y_j >= y_(i)}`.
    With no ties in y this reduces algebraically to `1 - 3*sum|dr|/(n^2 - 1)`, which is
    the tie-free reduction - correct on tie-free data and wrong the moment a
    quantised metric produces tied responses.

    Measured on the shipped artifacts, which is where the number has to come from: the
    T2* ladder carries 234 window-rows on 0704 and 75 on 1004, 309 in total, and NOT ONE
    has a tie on either axis - so the correction moved no shipped number. (An earlier draft
    of this docstring said "0 of 279"; no artifact holds 279, and the count is recorded here
    per dataset because pooling two datasets into one figure is what produced the wrong
    one.) The correction is not decoration: a fidelity ladder, where most windows are a
    single read long, ties heavily.

    NOTE ON THE FINITE-n CEILING. `xi -> 1` for a functional relation only asymptotically.
    On TIE-FREE data the maximum attainable is `1 - 3/(n+1)` - 0.857 at n = 20. With ties
    in y the ceiling is HIGHER, because the denominator shrinks: measured, n = 12 with two
    distinct y values reaches 0.833 where the tie-free ceiling would be 0.769. Either way
    xi must never be read against 1 at the window sizes here, which run from 3 reads up.

    NaN for n < 3, for constant x (the stable argsort would otherwise return a number
    describing acquisition order) and for constant y (the denominator vanishes).

    DEVIATION FROM THE SOURCE, recorded rather than hidden: Chatterjee breaks ties in x
    "uniformly at random". This uses a deterministic stable sort, because a seeded shuffle
    inside a reproducibility tool is its own problem - the same decision `tie_fraction`
    documents. The cost is measured in the tie experiment rather than assumed to be nil.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    if n < 3:
        return float("nan")
    finite = np.isfinite(x) & np.isfinite(y)
    if not finite.all():
        x, y, n = x[finite], y[finite], int(finite.sum())
        if n < 3:
            return float("nan")
    if float(np.max(x)) == float(np.min(x)):
        return float("nan")
    if float(np.max(y)) == float(np.min(y)):
        # Constant y makes l_i = n for every i, so the tie-corrected denominator is 0.
        # It is also the degenerate case: a constant response is trivially a function of
        # anything, and reporting xi = 1 for it would be an artefact, not a finding.
        return float("nan")
    order = np.argsort(x, kind="stable")
    # r_i = #{j : y_j <= y_(i)} and l_i = #{j : y_j >= y_(i)}, the max-rank convention.
    r = rankdata(y, method="max")[order]
    ell = rankdata(-y, method="max")
    denominator = 2.0 * float(np.sum(ell * (n - ell)))
    if denominator <= 0.0:
        return float("nan")
    return float(1.0 - n * float(np.abs(np.diff(r)).sum()) / denominator)


def tie_fraction(x: np.ndarray) -> float:
    """Fraction of values in `x` that share their value with another read.

    chatterjee_xi sorts by `x`; tied values are ordered by acquisition sequence, so
    part of xi comes from measurement order rather than from the data. Deterministic
    and reported rather than randomised away - a seeded shuffle inside a
    reproducibility tool is its own problem.
    """
    v = np.asarray(x, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) < 2:
        return float("nan")
    _, counts = np.unique(v, return_counts=True)
    return float(np.sum(counts[counts > 1]) / len(v))


XI_METHOD_ASYMPTOTIC = "asymptotic"
XI_METHOD_PERMUTATION = "permutation"
XI_METHOD_NONE = "none"

# Above this fraction of tied observations the closed form is refused outright. The 2/5
# null is Chatterjee (2021) Theorem 2.1, p. 2011, which assumes Y continuous; with ties in
# Y the limiting variance is the tau^2 of Theorem 2.2 (formula (3)), not 2/5. Any tie at
# all makes the closed form approximate, so the cutoff is zero and tied data goes to a
# permutation null instead.
XI_TIE_CUTOFF = 0.0


def xi_p_value_asymptotic(xi: float, n: int) -> float:
    """One-sided p under independence via `sqrt(n) * xi -> N(0, 2/5)`.

    VALID ONLY FOR TIE-FREE DATA. The 2/5 variance is derived for continuous observations;
    with ties the null variance is different and this over- or under-states significance
    with no warning. Use `xi_p_value`, which picks the calibration by measuring the ties.

    Asymptotic in n = READS IN THE WINDOW. Chatterjee calls the normal approximation
    "roughly valid even for n as small as 20"; at the n this repo sees (~7) it is
    indicative at best, which is why n travels beside it.
    """
    if not np.isfinite(xi) or n < 3:
        return float("nan")
    return float(norm.sf(np.sqrt(n) * xi / np.sqrt(XI_NULL_VARIANCE)))


def xi_p_value(
    x: np.ndarray,
    y: np.ndarray,
    *,
    seed: int,
    n_resamples: int = 999,
    tie_cutoff: float = XI_TIE_CUTOFF,
) -> tuple[float, str]:
    """One-sided p for `xi(x, y)`, with the calibration chosen by measuring the ties.

    Returns `(p_value, method)`. The method is part of the answer, not metadata: a
    closed-form p and a permutation p on the same window are different objects, and a
    reader comparing two rows has to know which is which.

    Tie-free -> the closed form, which is what the T2* ladder gets (measured: 0 of 309
    windows tied on either axis) and costs nothing. Tied -> a permutation against the
    independence null, which is exact whatever the ties. The choice is this repository's,
    recorded in `.claude/qre_checks_reference.tex`. Switching on a measurement rather than picking one globally is what
    keeps the cheap path cheap without leaving the tied path invalid.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    n = int(len(x))
    value = chatterjee_xi(x, y)
    if not np.isfinite(value) or n < 3:
        return float("nan"), XI_METHOD_NONE

    ties = max(tie_fraction(x), tie_fraction(y))
    if np.isfinite(ties) and ties > tie_cutoff:
        from quebra.analyzers.permutation import paired_permutation_test

        result = paired_permutation_test(
            x, y, chatterjee_xi, seed=seed, n_resamples=n_resamples
        )
        return result.p_value, XI_METHOD_PERMUTATION
    return xi_p_value_asymptotic(value, n), XI_METHOD_ASYMPTOTIC


DCOR_TOO_SHORT = "too_short"
DCOR_CONSTANT = "constant"
DCOR_SIZE_SKIPPED = "size_skipped"
DCOR_OK = "ok"


def dcor_with_reason(x: np.ndarray, y: np.ndarray) -> tuple[float, str]:
    """Distance correlation, plus WHY it is NaN when it is.

    Three different conditions produce NaN - too few usable reads, a constant input, and
    the O(n^2) memory skip above DCOR_MAX_READS - and a bare defined-count cannot tell a
    data limitation from a compute limitation.

    Non-finite pairs are masked FIRST so the count that decides "too short" is the
    usable count. Attributing a NaN caused by non-finite input to "constant" would put
    a wrong reason in the contract, which is worse than no reason at all.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    n = len(x)
    if n > DCOR_MAX_READS:
        return float("nan"), DCOR_SIZE_SKIPPED
    if n < 4:
        return float("nan"), DCOR_TOO_SHORT
    value = dcor(x, y)
    return (value, DCOR_OK) if np.isfinite(value) else (float("nan"), DCOR_CONSTANT)


def dcor(x: np.ndarray, y: np.ndarray) -> float:
    """Distance correlation (Szekely, Rizzo & Bakirov 2007); 0 iff independent."""
    x = np.asarray(x, dtype=float)[:, None]
    y = np.asarray(y, dtype=float)[:, None]
    if len(x) < 4:
        return float("nan")
    a, b = np.abs(x - x.T), np.abs(y - y.T)
    A = a - a.mean(0) - a.mean(1)[:, None] + a.mean()
    B = b - b.mean(0) - b.mean(1)[:, None] + b.mean()
    vx, vy = (A * A).mean(), (B * B).mean()
    if vx <= 0 or vy <= 0:
        return float("nan")
    # The clamp fixes a negative sample dCov^2 before the sqrt; it also means an exact
    # 0.0 is ambiguous between "independent" and "clamped".
    return float(np.sqrt(max((A * B).mean(), 0.0) / np.sqrt(vx * vy)))


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or float(np.max(x)) == float(np.min(x)):
        return float("nan")
    if float(np.max(y)) == float(np.min(y)):
        return float("nan")
    return float(spearmanr(x, y).statistic)


def window_stats(
    margin: np.ndarray,
    forward_time_s: np.ndarray,
    duration_s: float,
    window_index: int,
    *,
    xi_seed: int = 0,
) -> tuple[dict[str, float], str]:
    """Every shape statistic for one window, plus why dcor is NaN if it is."""
    m = np.asarray(margin, dtype=float)
    f = np.asarray(forward_time_s, dtype=float)
    n = int(len(m))
    peak = int(np.argmax(m)) if n else 0

    _xi_p, _xi_method = xi_p_value(m, f, seed=xi_seed)

    rho = _spearman(m, f)
    rho2 = float("nan") if not np.isfinite(rho) else rho**2
    xi_mf = chatterjee_xi(m, f)
    dcor_value, dcor_reason = dcor_with_reason(m, f)
    return {
        "window_index": int(window_index),
        "n_reads": n,
        "duration_s": float(duration_s),
        "peak_frac": float(peak / max(n - 1, 1)),
        "rho": rho,
        "rho2": rho2,
        # raw rho2 carries a 1/(n-1) null inflation, so short windows score higher for
        # free; the excess is the length-comparable form.
        "rho2_excess": float("nan") if not np.isfinite(rho2) else rho2 - 1.0 / (n - 1),
        "limb_rho": _spearman(m[peak:], f[peak:]) if n - peak >= 4 else float("nan"),
        "xi_m_to_f": xi_mf,
        "xi_f_to_m": chatterjee_xi(f, m),
        "xi_p_value": _xi_p,
        "xi_p_method": _xi_method,
        # How much of xi_m_to_f rests on tie-breaking rather than on the data.
        "xi_tie_frac": tie_fraction(m),
        "dcor": dcor_value,
    }, dcor_reason


def shape_curve(
    ages: list[np.ndarray], margins: list[np.ndarray], n_bins: int = 25
) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None, int]:
    """Mean margin against normalised window age, each window scaled by its own peak.

    Returns (grid, mean, standard error, n_windows). Mean and stderr are None when no
    window qualifies, so the renderer can say so rather than draw an empty axis.
    """
    grid = np.linspace(0.0, 1.0, n_bins)
    rows: list[np.ndarray] = []
    for age, margin in zip(ages, margins):
        if len(margin) < 2 or float(np.max(margin)) <= 0.0:
            continue
        span = float(age[-1] - age[0])
        if span <= 0.0:
            continue
        normalised_age = (age - age[0]) / span
        # np.interp requires increasing xp; duplicate timestamps make age non-strict.
        order = np.argsort(normalised_age, kind="stable")
        rows.append(
            np.interp(grid, normalised_age[order], (margin / np.max(margin))[order])
        )
    if not rows:
        return grid, None, None, 0
    stacked = np.asarray(rows)
    return (
        grid,
        stacked.mean(0),
        stacked.std(0) / np.sqrt(len(stacked)),
        int(len(stacked)),
    )


def pooled_xi(
    reads: pd.DataFrame, windows: pd.DataFrame, *, xi_seed: int = 0
) -> tuple[float, float, int, int, str]:
    """Chatterjee's xi over the reads of ALL windows pooled together.

    A different question from the median of per-window xi, and worth both. The median
    asks "inside a typical excursion, does margin predict remaining life"; the pooled
    statistic asks "across every excursion at this threshold, taken as one sample". Pooling
    mixes windows of different lengths, so a pooled effect can come from between-window
    structure rather than within - which is exactly why the two are shown side by side
    and never collapsed into one number.

    Deliberately NOT filtered by `shape_min_reads`: pooling is what gives coverage at
    thresholds where no single window is long enough, and dropping short windows here would
    throw away the reason to pool. Complete windows only, as everywhere else.

    Returns (xi, p_value, n_reads_pooled, n_windows_pooled).
    """
    if not len(windows) or not len(reads):
        return float("nan"), float("nan"), 0, 0, XI_METHOD_NONE
    in_spec = reads[reads["in_spec"]]
    keep = in_spec[in_spec["window_index"].isin(windows["window_index"])]
    if len(keep) < 3:
        return (
            float("nan"),
            float("nan"),
            int(len(keep)),
            int(len(windows)),
            XI_METHOD_NONE,
        )
    margin = keep["margin"].to_numpy(dtype=float)
    forward = keep["forward_time_s"].to_numpy(dtype=float)
    xi = chatterjee_xi(margin, forward)
    # The method is part of the answer here too: a pooled sample runs to thousands of
    # reads and is far more tie-prone than any single window, so which calibration
    # produced this p is not a detail.
    p_pooled, method = xi_p_value(margin, forward, seed=xi_seed)
    return xi, p_pooled, int(len(keep)), int(len(windows)), method


def pooled_rho(reads: pd.DataFrame, windows: pd.DataFrame) -> tuple[float, int, int]:
    """Spearman rho over the reads of ALL windows pooled, per the same rule as pooled_xi.

    Kept beside xi rather than folded into it: rho carries a SIGN (positive = a larger
    margin predicts more time left) and is stronger than xi against smooth monotone
    dependence, while xi catches non-monotone structure rho is blind to. Reading either
    alone is the mistake the claim-discipline rule exists to prevent.

    Returns (rho, n_reads_pooled, n_windows_pooled).
    """
    if not len(windows) or not len(reads):
        return float("nan"), 0, 0
    in_spec = reads[reads["in_spec"]]
    keep = in_spec[in_spec["window_index"].isin(windows["window_index"])]
    if len(keep) < 3:
        return float("nan"), int(len(keep)), int(len(windows))
    return (
        _spearman(
            keep["margin"].to_numpy(dtype=float),
            keep["forward_time_s"].to_numpy(dtype=float),
        ),
        int(len(keep)),
        int(len(windows)),
    )


def for_windows(
    *,
    reads: pd.DataFrame,
    windows: pd.DataFrame,
    threshold_value: float,
    big_values_good: bool,
    xi_seed: int = 0,
) -> tuple[pd.DataFrame, dict[str, float], dict[str, int], tuple]:
    """Shape statistics for every window in `windows`, plus medians and defined-counts.

    `windows` must already be filtered to eligible windows (complete, and at or above
    `shape_min_reads`) - this function does not decide eligibility.
    """
    del threshold_value, big_values_good  # margin already carries the direction
    if len(windows) and not len(reads):
        # Same disagreement as a window with no reads, one level up.
        raise ValueError(
            f"{len(windows)} eligible window(s) but an empty read table - the window "
            f"and read tables disagree."
        )
    if not len(windows):
        empty = pd.DataFrame(columns=STAT_COLUMNS)
        return empty, {}, {}, (np.linspace(0.0, 1.0, 25), None, None, 0)

    in_spec = reads[reads["in_spec"]]
    rows: list[dict[str, float]] = []
    ages: list[np.ndarray] = []
    margins: list[np.ndarray] = []
    dcor_reasons: dict[str, int] = {}
    for _, w in windows.iterrows():
        index = int(w["window_index"])
        sub = in_spec[in_spec["window_index"] == index].sort_values("t_read_s")
        if not len(sub):
            # The window and read tables come from one carve; a window with no reads
            # means they disagree, which is a bug upstream, not a data condition.
            raise ValueError(
                f"window {index} of threshold "
                f"{w.get('threshold_label', '?')!r} has no in-spec reads in the read "
                f"table - the window and read tables disagree."
            )
        margin = sub["margin"].to_numpy(dtype=float)
        stats_row, dcor_reason = window_stats(
            margin=margin,
            forward_time_s=sub["forward_time_s"].to_numpy(dtype=float),
            duration_s=float(w["duration_s"]),
            window_index=index,
            xi_seed=xi_seed,
        )
        rows.append(stats_row)
        dcor_reasons[dcor_reason] = dcor_reasons.get(dcor_reason, 0) + 1
        ages.append(sub["window_age_s"].to_numpy(dtype=float))
        margins.append(margin)

    stats = pd.DataFrame(rows, columns=STAT_COLUMNS)
    medians: dict[str, float] = {}
    defined: dict[str, int] = {}
    for column in STAT_COLUMNS:
        # `window_index` is an identifier and `xi_p_method` is a categorical label naming
        # which calibration produced the p-value beside it; a median over either is
        # meaningless, and `to_numpy(dtype=float)` on the label raises.
        if column in NON_NUMERIC_STAT_COLUMNS:
            continue
        values = stats[column].to_numpy(dtype=float) if len(stats) else np.array([])
        finite = values[np.isfinite(values)]
        defined[column] = int(len(finite))
        medians[column] = float(np.median(finite)) if len(finite) else float("nan")
    # Why dcor is missing where it is missing, so a low defined-count is readable.
    for reason in (DCOR_TOO_SHORT, DCOR_CONSTANT, DCOR_SIZE_SKIPPED):
        defined[f"dcor_{reason}"] = int(dcor_reasons.get(reason, 0))
    return stats, medians, defined, shape_curve(ages, margins)
