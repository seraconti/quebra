"""Survival estimators on an event table: Nelson-Aalen, curve summaries, and more.

Every estimator here is a step: it takes an `event_table.EventTable` and nothing else,
works in seconds, and returns a typed result, with a typed empty result for an empty
table. The product-limit arithmetic, the quantile rule and the log-log band are
`kaplan_meier.product_limit`, `survival_quantile` and `loglog_band`, so there is one
implementation of each.

Each result names its standing (`check_attachment`): what a check rejection would do to
it. Bands, intervals and standard errors rest on exchangeability and carry the check
outcome; point values of the Kaplan-Meier curve and its functionals are descriptive, with
`n_censored` beside them.

The fixtures in `jobs/reference/` pin the numbers to R `survival`: the Kaplan-Meier curve
and its band, H and its standard error, the quantiles with their intervals, and the
restricted mean with its standard error. R reports no band on H, so a hand computation in
`tests/test_survival_r_reference.py` pins the Nelson-Aalen band.

Each estimator's docstring states its formula and the assumptions that are its own. The
one below is shared by Nelson-Aalen, the curve summaries and the Kaplan-Meier curve.

Validity assumptions

- Assumption: the multiplicative intensity model. Every window at risk at age t shares one
  hazard, and censoring is independent of the window's future.
  Diagnostic: the attached check outcome (`a1_renewal_durations`); censoring here happens
  only at read gaps and at the end of the scan.
  Consequence of violation: Kaplan-Meier and Nelson-Aalen estimate a mixture whose weights
  shift as windows leave the risk set; with censored windows the point values move too.
  Reference: no source located (the model is usually credited to Aalen, not opened here).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Generic, TypeVar

import numpy as np

from quebra.analyzers.check_attachment import (
    DESCRIPTIVE,
    RESTS_ON_EXCHANGEABILITY,
    SURVIVES,
    CheckAttachment,
)
from quebra.analyzers.event_table import EventTable, EventTables
from quebra.analyzers.kaplan_meier import (
    CONF_LEVEL,
    QUANTILE_TOL,
    loglog_band,
    product_limit,
    survival_quantile,
    z_two_sided,
)

DEFAULT_QUANTILES = (0.25, 0.5, 0.75)

T = TypeVar("T")


def _restriction(table: EventTable, rmst_tau_s: float | None) -> float:
    """The restriction time: `rmst_tau_s`, else the largest observed age (NaN if empty).

    Validated before any branch on the table, so an empty table refuses what a full one
    does.
    """
    if rmst_tau_s is None:
        return float(np.max(table.age_s)) if table.n_windows else float("nan")
    tau = float(rmst_tau_s)
    if not np.isfinite(tau) or tau < 0.0:
        raise ValueError(f"rmst_tau_s must be finite and non-negative; got {tau!r}")
    return tau


def _counts(table: EventTable) -> dict[str, object]:
    return {
        "dataset_id": table.dataset_id,
        "threshold_label": table.threshold_label,
        "side": table.side,
        "n_windows": table.n_windows,
        "n_deaths": table.n_deaths,
        "n_censored": table.n_windows - table.n_deaths,
        "n_unobserved_birth_dropped": table.n_unobserved_birth_dropped,
    }


# ---------------------------------------------------------------------------
# Nelson-Aalen
# ---------------------------------------------------------------------------


@dataclass
class NelsonAalen:
    """Cumulative hazard H(t) = sum d/n, its Aalen variance sum d/n^2, and a log band.

    `time_s` starts at 0 with H = 0 and holds one entry per age with a death. The band is
    H exp(+/- z sqrt(V) / H), NaN where H = 0. No step fills `checks`, so it is empty
    and reads NOT ASSESSED.
    """

    dataset_id: str
    threshold_label: str
    side: str
    conf_level: float
    time_s: np.ndarray
    cumulative_hazard: np.ndarray
    variance: np.ndarray
    band_lower: np.ndarray
    band_upper: np.ndarray
    n_windows: int
    n_deaths: int
    n_censored: int
    n_unobserved_birth_dropped: int
    standing: str = RESTS_ON_EXCHANGEABILITY
    checks: CheckAttachment = field(default_factory=CheckAttachment)


def nelson_aalen(table: EventTable, *, conf_level: float = CONF_LEVEL) -> NelsonAalen:
    """Nelson-Aalen cumulative hazard, its Aalen variance, and a band symmetric on log H.

    With d_j deaths at age t_j and n_j windows at risk there, H(t) is the sum of
    d_j / n_j over t_j <= t and V(t) the sum of d_j / n_j^2; tied deaths enter as one term.
    The band is H exp(+/- z sqrt(V) / H), so it stays positive, and is NaN where H = 0.
    R `survival` (`ctype = 1`, `std.chaz`) gives the same H and sqrt(V) on the fixtures.
    Reference: no source located for the estimator, the variance or the band (the
    estimator is usually credited to Nelson and to Aalen, not opened here).

    Validity assumptions

    - Assumption: a falling hazard is a property of each window, not of the mix.
      Diagnostic: none here. A mixture of laws each with a flat hazard gives a falling
      pooled hazard, so a concave H cannot be read as windows getting safer with age.
      Consequence of violation: a selection effect read as ageing.
      Reference: no source located (the frailty argument usually credited to Vaupel,
      Manton and Stallard, 1979, not opened here).
    - Assumption: log H is close to normal with standard error sqrt(V) / H.
      Diagnostic: none here; the approximation is poorest where few deaths have
      accumulated, at the youngest ages.
      Consequence of violation: pointwise coverage away from `conf_level`.
      Reference: no source located.
    """
    died = table.n_events > 0
    d = table.n_events[died].astype(float)
    n = table.n_at_risk[died].astype(float)
    time_s = np.concatenate([[0.0], table.event_age_s[died]])
    hazard = np.concatenate([[0.0], np.cumsum(d / n)])
    variance = np.concatenate([[0.0], np.cumsum(d / n**2)])
    lower = np.full_like(hazard, np.nan)
    upper = np.full_like(hazard, np.nan)
    ok = hazard > 0.0
    spread = z_two_sided(conf_level) * np.sqrt(variance[ok]) / hazard[ok]
    lower[ok] = hazard[ok] * np.exp(-spread)
    upper[ok] = hazard[ok] * np.exp(spread)
    return NelsonAalen(
        conf_level=conf_level,
        time_s=time_s,
        cumulative_hazard=hazard,
        variance=variance,
        band_lower=lower,
        band_upper=upper,
        **_counts(table),  # type: ignore[arg-type]
    )


# ---------------------------------------------------------------------------
# Curve summaries
# ---------------------------------------------------------------------------


@dataclass
class CurveSummaries:
    """Quantiles of the Kaplan-Meier curve with intervals, and the restricted mean.

    A quantile is None where the curve, or the band curve for an interval end, never
    reaches 1 - p. `rmst_tau_s` is recorded because the restricted mean means nothing
    without it. Point values are descriptive; intervals and the SE rest on exchangeability.

    `max_observed_age_s` is the largest age in the table, censored or not. A `rmst_tau_s`
    past it extends the last value of S flat to tau, as R's `summary(fit, rmean = tau)`
    does, so the restricted mean then counts time no window was watched for.

    No step fills `checks`, so it is empty and reads NOT ASSESSED.
    """

    dataset_id: str
    threshold_label: str
    side: str
    conf_level: float
    quantile_p: tuple[float, ...]
    quantile_s: tuple[float | None, ...]
    quantile_lower_s: tuple[float | None, ...]
    quantile_upper_s: tuple[float | None, ...]
    rmst_tau_s: float
    rmst_s: float
    rmst_se_s: float
    max_observed_age_s: float
    n_windows: int
    n_deaths: int
    n_censored: int
    n_unobserved_birth_dropped: int
    point_standing: str = DESCRIPTIVE
    interval_standing: str = RESTS_ON_EXCHANGEABILITY
    checks: CheckAttachment = field(default_factory=CheckAttachment)


def _restricted_mean(
    time: np.ndarray, survival: np.ndarray, d: np.ndarray, n: np.ndarray, tau: float
) -> tuple[float, float]:
    """Area under the step curve on [0, tau], and its Greenwood-type standard error.

    Var = sum over deaths before tau of A_j^2 d_j / (n_j (n_j - d_j)), with A_j the area
    from t_j to tau. A risk set fully consumed leaves S = 0 after it, so A_j = 0 and the
    term is 0 rather than 0 * inf.
    """
    edges = np.concatenate([np.minimum(time, tau), [tau]])
    pieces = survival * np.diff(edges)
    # area_after[j]: the area from time[j] to tau, the segment starting at time[j] included
    area_after = np.cumsum(pieces[::-1])[::-1]
    variance = 0.0
    for j in range(1, len(time)):
        if time[j] >= tau or area_after[j] == 0.0:
            continue
        variance += area_after[j] ** 2 * d[j - 1] / (n[j - 1] * (n[j - 1] - d[j - 1]))
    return float(np.sum(pieces)), float(np.sqrt(variance))


def curve_summaries(
    table: EventTable,
    *,
    quantiles: tuple[float, ...] = DEFAULT_QUANTILES,
    rmst_tau_s: float | None = None,
    conf_level: float = CONF_LEVEL,
) -> CurveSummaries:
    """Quantiles with Brookmeyer-Crowley intervals, and the restricted mean with its SE.

    The p quantile is the age where S first falls below 1 - p or, where S sits exactly
    at 1 - p over a stretch, that stretch's midpoint (R's rule,
    `kaplan_meier.survival_quantile`). Its interval applies the same rule to the two
    log-log band curves, which inverts the band. The restricted mean is the area under S
    on [0, tau]; its variance is the sum over death ages t_j < tau of
    A_j^2 d_j / (n_j (n_j - d_j)), with A_j the area under S from t_j to tau, and a term
    with A_j = 0 counts 0. R `survival` 3.8.6 (`quantile.survfit`,
    `summary(fit, rmean = tau)`) gives the same numbers on the fixtures. Reference: no
    source located for the interval (usually credited to Brookmeyer and Crowley, not
    opened here) or for the variance (Klein and Moeschberger's textbook, not opened here).

    `rmst_tau_s` defaults to the largest observed age. Every probability in `quantiles`
    must lie in (0, 1).

    Validity assumptions

    - Assumption: the quantile interval and the variance are the formulas above, which
      are R's.
      Diagnostic: agreement with R on the committed fixtures.
      Consequence of violation: an interval that is not the one its label names.
      Reference: R `survival` 3.8.6 `quantile.survfit` and `summary(fit, rmean = tau)`.
    - Assumption: the band the quantile interval inverts covers S pointwise at
      `conf_level`.
      Diagnostic: none here; the band is `kaplan_meier.loglog_band` on Greenwood's
      variance, an approximation poorest where few windows remain at risk.
      Consequence of violation: a quantile interval whose coverage is not `conf_level`.
      Reference: no source located.
    """
    if not all(0.0 < p < 1.0 for p in quantiles):
        raise ValueError(
            f"every quantile probability must lie in (0, 1); got {quantiles!r}"
        )
    tau = _restriction(table, rmst_tau_s)
    z = z_two_sided(conf_level)
    counts = _counts(table)
    if table.n_windows == 0:
        none = tuple(None for _ in quantiles)
        return CurveSummaries(
            conf_level=conf_level,
            quantile_p=tuple(quantiles),
            quantile_s=none,
            quantile_lower_s=none,
            quantile_upper_s=none,
            rmst_tau_s=tau,
            rmst_s=float("nan"),
            rmst_se_s=float("nan"),
            max_observed_age_s=float("nan"),
            **counts,  # type: ignore[arg-type]
        )
    time, survival, n_at_risk, greenwood = product_limit(
        table.event_age_s, table.n_events, table.n_at_risk
    )
    lower, upper = loglog_band(survival, greenwood, z)
    max_age = float(np.max(table.age_s))
    point = tuple(survival_quantile(time, survival, p, max_age) for p in quantiles)
    low = tuple(survival_quantile(time, lower, p, max_age) for p in quantiles)
    high = tuple(survival_quantile(time, upper, p, max_age) for p in quantiles)
    died = table.n_events > 0
    rmst, se = _restricted_mean(
        time,
        survival,
        table.n_events[died].astype(float),
        table.n_at_risk[died].astype(float),
        tau,
    )
    return CurveSummaries(
        conf_level=conf_level,
        quantile_p=tuple(quantiles),
        quantile_s=point,
        quantile_lower_s=low,
        quantile_upper_s=high,
        rmst_tau_s=tau,
        rmst_s=rmst,
        rmst_se_s=se,
        max_observed_age_s=max_age,
        **counts,  # type: ignore[arg-type]
    )


# ---------------------------------------------------------------------------
# Turnbull
# ---------------------------------------------------------------------------

TURNBULL_TOL = 1e-10
TURNBULL_MAX_ITER = 100_000


@dataclass
class Turnbull:
    """Turnbull's self-consistency estimate on the length intervals [age_lo_s, age_hi_s].

    The length interval is closed, and a single point when lo == hi: an exact length,
    which duplicate timestamps produce. The likelihood reads a non-degenerate interval as
    (age_lo_s, age_hi_s], as R's `Surv(type = "interval2")` does; `turnbull` states what
    that costs.

    The estimate puts mass on innermost intervals only: `innermost_lo_s` and
    `innermost_hi_s` bound each, `innermost_point` marks an exact length (lo == hi), and
    `mass` is the probability on it. Where inside an interval the mass sits is not
    identified, so no curve is drawn from this; the survival after each interval,
    `survival_after`, is the identified part.

    `converged` is True when an iteration raised the log-likelihood by `tol` or less:
    the EM stopped, which is not the same as the masses sitting at the maximum, since on
    a flat likelihood they can still be moving. It is False when `max_iter` ran out
    first, and True for an empty table, where there is nothing to iterate.

    Standalone: it plays no part in any other result. No step fills `checks`, so it is
    empty and reads NOT ASSESSED.
    """

    dataset_id: str
    threshold_label: str
    side: str
    innermost_lo_s: np.ndarray
    innermost_hi_s: np.ndarray
    innermost_point: np.ndarray
    mass: np.ndarray
    survival_after: np.ndarray
    log_likelihood: float
    converged: bool
    n_iter: int
    last_change: float
    tol: float
    max_iter: int
    n_windows: int
    n_deaths: int
    n_censored: int
    n_unobserved_birth_dropped: int
    standing: str = RESTS_ON_EXCHANGEABILITY
    checks: CheckAttachment = field(default_factory=CheckAttachment)


def _innermost(
    lo: np.ndarray, hi: np.ndarray, exact: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Turnbull's innermost intervals: a left end followed by a right end, nothing between.

    A non-degenerate interval is read half-open, (lo, hi]: at a tie a right end sorts
    before a left end, since x lies in (., x] but not in (x, .]. An exact length x is the
    point [x, x]: its left end sorts before every right end at x.
    """
    ends: list[tuple[float, int, bool]] = []
    for a, b, point in zip(lo, hi, exact):
        ends.append((float(a), 0 if point else 2, True))
        ends.append((float(b), 1, False))
    ends.sort(key=lambda e: (e[0], e[1]))
    out: set[tuple[float, float, bool]] = set()
    for (x1, k1, left1), (x2, _k2, left2) in zip(ends, ends[1:]):
        if left1 and not left2 and (x2 > x1 or k1 == 0):
            out.add((x1, x2, k1 == 0 and x1 == x2))
    rows = sorted(out)
    return (
        np.array([r[0] for r in rows], dtype=float),
        np.array([r[1] for r in rows], dtype=float),
        np.array([r[2] for r in rows], dtype=bool),
    )


def turnbull(
    table: EventTable,
    *,
    tol: float = TURNBULL_TOL,
    max_iter: int = TURNBULL_MAX_ITER,
) -> Turnbull:
    """Self-consistency EM (Turnbull) on the table's length intervals.

    The innermost intervals run from a left end to the next right end with no end
    between them (`_innermost`). With A_ij = 1 when innermost interval j lies in window
    i's interval, and N windows, the update from equal masses is
    m_j <- m_j (1/N) sum_i A_ij / (sum_k A_ik m_k). Reference: no source located for the
    construction or the update (Turnbull, 1976, not opened here).

    Stops when the log-likelihood rises by `tol` or less, or after `max_iter`
    iterations; then `converged` is False and the masses are the last iterate, reported
    as such rather than raised, because one slow cell must not stop a job. `tol` must be
    finite and non-negative (0 stops at the first iteration that does not raise it) and
    `max_iter` at least 1, checked before the table is read.

    Validity assumptions

    - Assumption: the interval on a window's length is non-informative, i.e. the chance of
      seeing this interval does not depend on where in it the true length lies.
      Diagnostic: none in the step. The read grid violates it: a length L shows n reads
      with a probability triangular in L across ((n-1)d, (n+1)d], so the estimate is
      biased, and the bias does not shrink as windows accumulate. A design-time
      simulation, not shipped, showed this; no tracked artifact measures its size.
      Consequence of violation: a curve that is not the length distribution of the windows
      seen; read it as an estimate under that assumption, never as a check on Kaplan-Meier.
      Reference: no source located (the condition is usually called the constant-sum
      condition, not opened here).
    - Assumption: birth and death are each censored to a read interval; one interval on the
      length is a simplification of that doubly censored structure.
      Diagnostic: none. Consequence of violation: as above. Reference: no source located.
    - Assumption: an end of a non-degenerate interval carries no probability, so reading
      [age_lo_s, age_hi_s] as (age_lo_s, age_hi_s] loses nothing.
      Diagnostic: none in the step. It holds when the crossing times have a continuous
      law. A duplicate timestamp at one crossing pins that crossing and leaves the length
      continuous across the interval; a duplicate at both makes a point, read as exact.
      Consequence of violation: windows whose intervals touch at an end keep separate
      innermost intervals where the closed reading would let them share a point, so mass
      can sit on the wrong side of that end.
      Reference: R `survival` 3.8.6 `?Surv`, argument `time2` (intervals "open on the
      left and closed on the right"); the R fixtures agree on where the mass sits.
    """
    if not (math.isfinite(tol) and tol >= 0.0):
        raise ValueError(f"tol must be finite and non-negative; got {tol!r}")
    if max_iter < 1:
        raise ValueError(f"max_iter must be at least 1; got {max_iter!r}")
    lo, hi = table.age_lo_s, table.age_hi_s
    exact = lo == hi
    q, p, point = _innermost(lo, hi, exact)
    counts = _counts(table)
    if table.n_windows == 0:
        empty = np.array([], dtype=float)
        return Turnbull(
            innermost_lo_s=empty,
            innermost_hi_s=empty.copy(),
            innermost_point=np.array([], dtype=bool),
            mass=empty.copy(),
            survival_after=empty.copy(),
            log_likelihood=float("nan"),
            converged=True,
            n_iter=0,
            last_change=0.0,
            tol=tol,
            max_iter=max_iter,
            **counts,  # type: ignore[arg-type]
        )
    # A[i, j]: innermost interval j lies inside window i's interval.
    inside_half_open = (q[None, :] >= lo[:, None]) & (p[None, :] <= hi[:, None])
    point_inside = (q[None, :] > lo[:, None]) & (q[None, :] <= hi[:, None])
    A = np.where(point[None, :], point_inside, inside_half_open)
    A = np.where(exact[:, None], point[None, :] & (q[None, :] == lo[:, None]), A)
    mass = np.full(len(q), 1.0 / len(q))
    old = -np.inf
    change = np.inf
    converged = False
    n_iter = 0
    for n_iter in range(1, max_iter + 1):
        denom = A @ mass
        mass = mass * (A / denom[:, None]).mean(axis=0)
        loglik = float(np.sum(np.log(A @ mass)))
        change = loglik - old
        old = loglik
        if change <= tol:
            converged = True
            break
    return Turnbull(
        innermost_lo_s=q,
        innermost_hi_s=p,
        innermost_point=point,
        mass=mass,
        survival_after=1.0 - np.cumsum(mass),
        log_likelihood=old,
        converged=converged,
        n_iter=n_iter,
        last_change=float(change),
        tol=tol,
        max_iter=max_iter,
        **counts,  # type: ignore[arg-type]
    )


# ---------------------------------------------------------------------------
# Placement bracket
# ---------------------------------------------------------------------------


@dataclass
class PlacementBracket:
    """How far the curve can move when each crossing moves inside its read interval.

    `lower` and `upper` are step curves on `time_s` (right-continuous from 0). Every
    Kaplan-Meier curve computed from crossing times consistent with the reads lies between
    them at every age, the shipped one included. `rmst_lower_s` and `rmst_upper_s` are
    their areas on [0, tau_s]; `width_s` is the difference. `width_censoring_only_s` is
    the same width with every window's length taken as its Kaplan-Meier age, so the part
    of `width_s` the grid accounts for is `width_s - width_censoring_only_s`.

    `median_lower_s` and `median_upper_s` are the medians of `lower` and `upper` by R's
    midpoint rule (`kaplan_meier.survival_quantile`), each on its own step times. The rule
    is monotone in the curve, so every consistent placement's median lies between them.
    `median_upper_s` is None when `upper` never falls below 1/2: at least half the windows
    are censored, and a censored placement can push the median past any finite value.

    Reported as numbers, with no threshold and no pass or fail. A statement about this
    record's estimate, needing no exchangeability, so it survives a check rejection.

    Why it bounds: redistribute-to-the-right moves a censored window's mass only onto
    later deaths, so any Kaplan-Meier curve lies between the curve counting censored
    windows as deaths at their censoring times and 1 - (deaths by t) / n. Widening each
    death to its length interval and each censoring to its lower end only widens those
    curves. No source located for redistribute-to-the-right (usually credited to Efron,
    not opened here).

    Validity assumptions

    - Assumption: each crossing lies inside its read interval, and the carve's windows are
      the true windows.
      Diagnostic: none yet; nothing measures excursions shorter than the read spacing or
      windows merged across one.
      Consequence of violation: the bracket says nothing about missed or merged windows.
      Reference: `docs/WINDOW_SEMANTICS.md`.
    - Assumption: a censored window's true censoring time is at least its last read.
      Diagnostic: holds by the carve's construction. Consequence of violation: the lower
      curve is no longer a bound. Reference: `analyzers/windows.py`, censoring convention.
    """

    dataset_id: str
    threshold_label: str
    side: str
    tau_s: float
    time_s: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    rmst_lower_s: float
    rmst_upper_s: float
    width_s: float
    width_censoring_only_s: float
    median_lower_s: float | None
    median_upper_s: float | None
    n_windows: int
    n_deaths: int
    n_censored: int
    n_unobserved_birth_dropped: int
    standing: str = SURVIVES


def _lower_at(lo: np.ndarray, at: np.ndarray) -> np.ndarray:
    """lower(t) = #{lo > t} / n over all windows."""
    return np.array([np.count_nonzero(lo > t) for t in at], dtype=float) / len(lo)


def _upper_at(hi: np.ndarray, event: np.ndarray, at: np.ndarray) -> np.ndarray:
    """upper(t) = 1 - #{deaths, hi <= t} / n."""
    dead = np.array([np.count_nonzero(event & (hi <= t)) for t in at], dtype=float)
    return 1.0 - dead / len(hi)


def _bracket_curves(
    lo: np.ndarray, hi: np.ndarray, event: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Both curves on the joint grid of their step times."""
    time_s = np.unique(np.concatenate([[0.0], lo, hi[event]]))
    return time_s, _lower_at(lo, time_s), _upper_at(hi, event, time_s)


def _bracket_medians(
    lo: np.ndarray, hi: np.ndarray, event: np.ndarray
) -> tuple[float | None, float | None]:
    """Each curve's median by R's convention, on that curve's own step times.

    On the joint grid the midpoint rule would end a flat stretch at the other curve's
    step. The upper median is None when `upper` never falls below 1/2: at least half the
    windows are censored, and a censored window placed later moves a consistent
    placement's median without bound.
    """
    t_lower = np.unique(np.concatenate([[0.0], lo]))
    t_upper = np.unique(np.concatenate([[0.0], hi[event]]))
    lower = _lower_at(lo, t_lower)
    upper = _upper_at(hi, event, t_upper)
    median_lower = survival_quantile(t_lower, lower, 0.5, float(t_lower[-1]))
    if upper[-1] >= 0.5 - QUANTILE_TOL:
        return median_lower, None
    return median_lower, survival_quantile(t_upper, upper, 0.5, float(t_upper[-1]))


def _bracket_areas(
    lo: np.ndarray, hi: np.ndarray, event: np.ndarray, tau: float
) -> tuple[float, float]:
    """Exact areas of the two curves on [0, tau], term by term."""
    n = len(lo)
    rmst_lower = math.fsum(np.minimum(lo, tau)) / n
    rmst_upper = (
        math.fsum(np.minimum(hi[event], tau)) + tau * int(np.count_nonzero(~event))
    ) / n
    return rmst_lower, rmst_upper


def placement_bracket(
    table: EventTable, *, rmst_tau_s: float | None = None
) -> PlacementBracket:
    """The bracket on the Kaplan-Meier curve set by placing crossings within the reads.

    `rmst_tau_s` defaults to the largest observed age, as in `curve_summaries`.
    """
    tau = _restriction(table, rmst_tau_s)
    counts = _counts(table)
    if table.n_windows == 0:
        empty = np.array([], dtype=float)
        return PlacementBracket(
            tau_s=tau,
            time_s=empty,
            lower=empty.copy(),
            upper=empty.copy(),
            rmst_lower_s=float("nan"),
            rmst_upper_s=float("nan"),
            width_s=float("nan"),
            width_censoring_only_s=float("nan"),
            median_lower_s=None,
            median_upper_s=None,
            **counts,  # type: ignore[arg-type]
        )
    lo, hi, event = table.age_lo_s, table.age_hi_s, table.event
    time_s, lower, upper = _bracket_curves(lo, hi, event)
    rmst_lower, rmst_upper = _bracket_areas(lo, hi, event, tau)
    km_lo, km_hi = _bracket_areas(table.age_s, table.age_s, event, tau)
    median_lower, median_upper = _bracket_medians(lo, hi, event)
    return PlacementBracket(
        tau_s=tau,
        time_s=time_s,
        lower=lower,
        upper=upper,
        rmst_lower_s=rmst_lower,
        rmst_upper_s=rmst_upper,
        width_s=rmst_upper - rmst_lower,
        width_censoring_only_s=km_hi - km_lo,
        median_lower_s=median_lower,
        median_upper_s=median_upper,
        **counts,  # type: ignore[arg-type]
    )


# ---------------------------------------------------------------------------
# Every threshold and side at once: the steps `recipes.wire_survival` wires
# ---------------------------------------------------------------------------


@dataclass
class ByThresholdSide(Generic[T]):
    """One result per (threshold label, side), on the carve's whole ladder."""

    ladder: tuple[str, ...]
    items: dict[tuple[str, str], T]

    def get(self, threshold_label: str, side: str) -> T:
        key = (threshold_label, side)
        if key not in self.items:
            raise KeyError(f"no result for {key}; ladder {list(self.ladder)}")
        return self.items[key]


def nelson_aalen_set(
    tables: EventTables, *, conf_level: float = CONF_LEVEL
) -> ByThresholdSide[NelsonAalen]:
    return ByThresholdSide(
        ladder=tables.ladder,
        items={
            k: nelson_aalen(t, conf_level=conf_level) for k, t in tables.tables.items()
        },
    )


def curve_summaries_set(
    tables: EventTables,
    *,
    quantiles: tuple[float, ...] = DEFAULT_QUANTILES,
    rmst_tau_s: float | None = None,
    conf_level: float = CONF_LEVEL,
) -> ByThresholdSide[CurveSummaries]:
    return ByThresholdSide(
        ladder=tables.ladder,
        items={
            k: curve_summaries(
                t, quantiles=quantiles, rmst_tau_s=rmst_tau_s, conf_level=conf_level
            )
            for k, t in tables.tables.items()
        },
    )


def turnbull_set(
    tables: EventTables,
    *,
    tol: float = TURNBULL_TOL,
    max_iter: int = TURNBULL_MAX_ITER,
) -> ByThresholdSide[Turnbull]:
    return ByThresholdSide(
        ladder=tables.ladder,
        items={
            k: turnbull(t, tol=tol, max_iter=max_iter) for k, t in tables.tables.items()
        },
    )


def placement_bracket_set(
    tables: EventTables, *, rmst_tau_s: float | None = None
) -> ByThresholdSide[PlacementBracket]:
    return ByThresholdSide(
        ladder=tables.ladder,
        items={
            k: placement_bracket(t, rmst_tau_s=rmst_tau_s)
            for k, t in tables.tables.items()
        },
    )
