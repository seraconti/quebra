"""Comparisons across records: a pair's log-rank test and restricted-mean difference, and a
random-effects summary of many records' restricted means.

Every step here takes the event tables of several records at one threshold and side, where
`analyzers/survival.py` takes one. The record is named by its table's `dataset_id`; two
tables of one comparison must agree on threshold label and side and name different records,
or the step raises.

The restricted means are compared up to one restriction time, the comparison tau, fixed by
a rule passed as a step argument (`tau_rule`), so the rule reaches the run identity and the
provenance label. The one rule, `min_largest_age`, takes the smallest of the records'
largest observed ages, so no record's curve is extended past its own data. Each result
records the tau it used and the record that set it.

A result that cannot be computed says so in a `*_not_computed_reason` field and carries None
in its numbers, never NaN: a missing number must not read as a measured one.

The fixtures in `jobs/reference/` pin the numbers to R: `survival::survdiff` for the
log-rank statistic, `summary(survfit, rmean = tau)` for each record's restricted mean, and
metafor 5.2.1 `rma(method = "DL")`, `confint` and `predict(predtype = "Riley")` for the
summary. All are written by `jobs/rscripts/comparison_reference.R`.

Validity assumptions

- Assumption: the windows of each record are exchangeable draws from one law, with
  independent censoring.
  Diagnostic: the check outcome attached to each record (`checks_a` and `checks_b` on a
  pair, each row's summary in a pooled result); empty reads NOT ASSESSED, and no step here
  fills it.
  Consequence of violation: the log-rank p, the SE of the difference and every interval of
  the summary are miscalibrated; the point difference is still the difference of the two
  areas.
  Reference: `analyzers/assumptions.py`, `a1_renewal_durations`.
- Assumption: the records are independent of each other.
  Diagnostic: none here. Two qubits read in the same night share an environment, and two
  halves of one scan share a qubit and a night.
  Consequence of violation: the log-rank variance, the SE of the difference and the
  summary's intervals are wrong, in a direction no step here measures.
  Reference: no source located.
- Assumption: the records were read on comparable schedules.
  Diagnostic: none here. A window's age lives on its record's read grid
  (`docs/WINDOW_SEMANTICS.md`), so two read spacings shift the ages apart by construction.
  Consequence of violation: a difference in read spacing reads as a difference in hazard,
  and as between-record variance in a summary.
  Reference: `docs/WINDOW_SEMANTICS.md`.
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
from scipy import optimize, stats

from quebra.analyzers.check_attachment import (
    DESCRIPTIVE,
    RESTS_ON_EXCHANGEABILITY,
    CheckAttachment,
)
from quebra.analyzers.event_table import EventTable
from quebra.analyzers.kaplan_meier import CONF_LEVEL, z_two_sided
from quebra.analyzers.survival import CurveSummaries, curve_summaries

TAU_RULE_MIN_LARGEST_AGE = "min_largest_age"
TAU_RULES = (TAU_RULE_MIN_LARGEST_AGE,)

NO_WINDOWS = "no windows at this threshold and side"
ZERO_VARIANCE = "no death before tau: zero variance, not pooled"
MIN_POOLED = 3
NO_DEATHS = "no deaths in either record"
ZERO_LOGRANK_VARIANCE = (
    "zero log-rank variance: at every death age one record had no window at risk, "
    "or every window at risk died"
)


# ---------------------------------------------------------------------------
# The comparison tau
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ComparisonTau:
    """The restriction time of one comparison, the record that set it, and who took part.

    `excluded` holds (record, reason) for every record left out before tau was taken;
    `included` names the rest, in the order given. `rmst_tau_s` and `tau_set_by` are None
    when no record is left.
    """

    tau_rule: str
    rmst_tau_s: float | None
    tau_set_by: str | None
    included: tuple[str, ...]
    excluded: tuple[tuple[str, str], ...]


def _require_comparable(tables: Sequence[EventTable]) -> None:
    """One threshold, one side, and every record named once; raise otherwise."""
    if not tables:
        raise ValueError("a comparison needs at least one event table")
    labels = {t.threshold_label for t in tables}
    sides = {t.side for t in tables}
    if len(labels) > 1 or len(sides) > 1:
        raise ValueError(
            f"records compared at different thresholds or sides: {sorted(labels)}, "
            f"{sorted(sides)}"
        )
    names = [t.dataset_id for t in tables]
    if len(set(names)) != len(names):
        raise ValueError(f"a record appears twice in one comparison: {names}")


def comparison_tau(tables: Sequence[EventTable], *, tau_rule: str) -> ComparisonTau:
    """The comparison tau of `tables` under `tau_rule`.

    `min_largest_age`: a record with no windows is excluded and named; tau is the smallest,
    over the rest, of each record's largest observed age, censored or not, and `tau_set_by`
    is the first record in the given order that attains it. Every record's curve is then
    observed up to tau, so no restricted mean counts time no window was watched for
    (`survival.CurveSummaries`).

    Raises on an unknown rule, and on tables that `_require_comparable` refuses.
    """
    if tau_rule not in TAU_RULES:
        raise ValueError(f"tau_rule must be one of {TAU_RULES}; got {tau_rule!r}")
    _require_comparable(tables)
    excluded = tuple((t.dataset_id, NO_WINDOWS) for t in tables if t.n_windows == 0)
    kept = [t for t in tables if t.n_windows > 0]
    if not kept:
        return ComparisonTau(tau_rule, None, None, (), excluded)
    largest = [float(np.max(t.age_s)) for t in kept]
    first = int(np.argmin(largest))
    return ComparisonTau(
        tau_rule=tau_rule,
        rmst_tau_s=largest[first],
        tau_set_by=kept[first].dataset_id,
        included=tuple(t.dataset_id for t in kept),
        excluded=excluded,
    )


# ---------------------------------------------------------------------------
# A pair: log-rank and the restricted-mean difference
# ---------------------------------------------------------------------------


@dataclass
class PairComparison:
    """Two records at one threshold and side: log-rank and restricted-mean difference.

    `logrank_z` is signed: positive when record a has more deaths than expected under one
    shared hazard. `logrank_chi2` is its square, R's `survdiff` statistic, and `logrank_p`
    the two-sided p, unadjusted for any other comparison. The difference is a - b, with
    SE sqrt(se_a^2 + se_b^2) from the two records' `curve_summaries` at `tau.rmst_tau_s`.

    `not_computed_reason` is set when a record has no windows: then nothing is computed, and
    `logrank_not_computed_reason` carries the same reason. That field alone is set when the
    log-rank statistic is undefined while the restricted means are not. The SE of the
    difference is 0 when neither record has a death before tau: the Greenwood-type variance
    has no term then, which is no claim of certainty.

    The point difference is descriptive; the p and the SE rest on exchangeability.
    `checks_a` / `checks_b` are each record's check-outcome field; no step here fills them,
    so they read NOT ASSESSED.
    """

    record_a: str
    record_b: str
    threshold_label: str
    side: str
    conf_level: float
    tau: ComparisonTau
    summary_a: CurveSummaries | None
    summary_b: CurveSummaries | None
    logrank_z: float | None
    logrank_chi2: float | None
    logrank_p: float | None
    rmst_difference_s: float | None
    rmst_difference_se_s: float | None
    not_computed_reason: str | None
    logrank_not_computed_reason: str | None
    point_standing: str = DESCRIPTIVE
    interval_standing: str = RESTS_ON_EXCHANGEABILITY
    checks_a: CheckAttachment = field(default_factory=CheckAttachment)
    checks_b: CheckAttachment = field(default_factory=CheckAttachment)

    def check_summary(self) -> str:
        return (
            f"{self.record_a}: {self.checks_a.summary()}; "
            f"{self.record_b}: {self.checks_b.summary()}"
        )


def _censored_data(table: EventTable) -> stats.CensoredData:
    return stats.CensoredData(
        uncensored=table.age_s[table.event], right=table.age_s[~table.event]
    )


def _logrank(
    a: EventTable, b: EventTable
) -> tuple[float | None, float | None, str | None]:
    """(z, p, reason): `scipy.stats.logrank`, or None with the reason it is undefined."""
    with warnings.catch_warnings(), np.errstate(divide="ignore", invalid="ignore"):
        warnings.simplefilter("ignore", RuntimeWarning)
        result = stats.logrank(_censored_data(a), _censored_data(b))
    z = float(result.statistic)
    if not math.isfinite(z):
        return (
            None,
            None,
            NO_DEATHS if a.n_deaths + b.n_deaths == 0 else ZERO_LOGRANK_VARIANCE,
        )
    return z, float(result.pvalue), None


def compare_pair(
    table_a: EventTable,
    table_b: EventTable,
    *,
    tau_rule: str,
    conf_level: float = CONF_LEVEL,
) -> PairComparison:
    """Log-rank test and restricted-mean difference of two records, up to the comparison tau.

    The log-rank statistic is (O_a - E_a) / sqrt(V), summed over the distinct death ages of
    the two records together, with the hypergeometric variance
    V = sum n_a n_b d (n - d) / (n^2 (n - 1)); scipy computes it and drops the terms where
    one window is at risk, whose variance is 0. R `survival` 3.8.6 `survdiff` gives the same
    chi-square on the fixtures. Each restricted mean and its SE are `curve_summaries` at the
    comparison tau; the difference's variance is the sum of the two, which holds for
    independent records. Reference: no source located for the test (usually credited to
    Mantel, not opened here).

    Validity assumptions

    - Assumption: the hazards of the two records are proportional.
      Diagnostic: none here. Consequence of violation: the test loses power against curves
      that cross, so a large p is not evidence that two curves agree.
      Reference: no source located.
    - The assumptions of the module docstring.
    """
    z_two_sided(conf_level)
    tau = comparison_tau((table_a, table_b), tau_rule=tau_rule)
    common = dict(
        record_a=table_a.dataset_id,
        record_b=table_b.dataset_id,
        threshold_label=table_a.threshold_label,
        side=table_a.side,
        conf_level=conf_level,
        tau=tau,
    )
    if tau.excluded:
        reason = "; ".join(f"{r}: {why}" for r, why in tau.excluded)
        return PairComparison(
            **common,  # type: ignore[arg-type]
            summary_a=None,
            summary_b=None,
            logrank_z=None,
            logrank_chi2=None,
            logrank_p=None,
            rmst_difference_s=None,
            rmst_difference_se_s=None,
            not_computed_reason=reason,
            logrank_not_computed_reason=reason,
        )
    summary_a = curve_summaries(
        table_a, rmst_tau_s=tau.rmst_tau_s, conf_level=conf_level
    )
    summary_b = curve_summaries(
        table_b, rmst_tau_s=tau.rmst_tau_s, conf_level=conf_level
    )
    z, p, reason = _logrank(table_a, table_b)
    return PairComparison(
        **common,  # type: ignore[arg-type]
        summary_a=summary_a,
        summary_b=summary_b,
        logrank_z=z,
        logrank_chi2=None if z is None else z * z,
        logrank_p=p,
        rmst_difference_s=summary_a.rmst_s - summary_b.rmst_s,
        rmst_difference_se_s=math.hypot(summary_a.rmst_se_s, summary_b.rmst_se_s),
        not_computed_reason=None,
        logrank_not_computed_reason=reason,
        checks_a=summary_a.checks,
        checks_b=summary_b.checks,
    )


# ---------------------------------------------------------------------------
# Many records: a DerSimonian-Laird summary of the restricted means
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RecordRow:
    """One record's restricted mean at the comparison tau, as a forest row.

    `lower_s` and `upper_s` are the normal interval rmst +/- z se, None when the SE is 0.
    `not_pooled_reason` names why a row with a restricted mean stays out of the pool.
    """

    record: str
    summary: CurveSummaries
    lower_s: float | None
    upper_s: float | None
    pooled: bool
    not_pooled_reason: str | None


@dataclass
class RandomEffectsSummary:
    """A DerSimonian-Laird summary of the records' restricted means at one tau.

    `rows` holds every record with windows, pooled or not; `tau.excluded` names the records
    with none. When fewer than `MIN_POOLED` rows have a variance, every pooled number is
    None and `not_pooled_reason` says so: the rows are listed, not pooled.

    `between_var_s2` is tau^2, the between-record variance of the true restricted means, in
    s^2; its interval is the Q-profile. `i2_fraction` is I^2 as a fraction, not a percent.
    The prediction interval is where the true restricted mean of one more record from the
    same population should fall; a measured one scatters wider, by its own variance. Every
    pooled number rests on exchangeability within each record and independence between them.
    Each row's summary carries that record's check-outcome field; no step here fills it, so
    it reads NOT ASSESSED.
    """

    threshold_label: str
    side: str
    conf_level: float
    tau: ComparisonTau
    rows: tuple[RecordRow, ...]
    k_pooled: int
    pooled_s: float | None
    pooled_se_s: float | None
    pooled_lower_s: float | None
    pooled_upper_s: float | None
    between_var_s2: float | None
    between_var_lower_s2: float | None
    between_var_upper_s2: float | None
    q_statistic: float | None
    q_df: int | None
    i2_fraction: float | None
    prediction_lower_s: float | None
    prediction_upper_s: float | None
    not_pooled_reason: str | None
    standing: str = RESTS_ON_EXCHANGEABILITY

    def check_summary(self) -> str:
        if not self.rows:
            return f"no record with windows; {CheckAttachment().summary()}"
        lines = [row.summary.checks.summary() for row in self.rows]
        if len(set(lines)) == 1:
            return f"every record: {lines[0]}"
        return "; ".join(f"{row.record}: {line}" for row, line in zip(self.rows, lines))


def _generalised_q(y: np.ndarray, v: np.ndarray, between_var: float) -> float:
    """sum (y_i - mu)^2 / (v_i + tau^2), mu the mean weighted by 1 / (v_i + tau^2)."""
    w = 1.0 / (v + between_var)
    mu = float(np.sum(w * y) / np.sum(w))
    return float(np.sum(w * (y - mu) ** 2))


def _q_root(y: np.ndarray, v: np.ndarray, target: float, start: float) -> float:
    """The tau^2 where the generalised Q, decreasing in tau^2, falls to `target`.

    The bracket's upper end doubles from `start` until Q is below the target: Q tends to 0
    as tau^2 grows, so it always crosses, and no fixed cap clips a bound measured in s^2.
    """
    hi = start
    for _ in range(2000):
        if _generalised_q(y, v, hi) < target:
            break
        hi *= 2.0
    else:
        raise ValueError(f"the Q-profile bracket never crossed {target!r}")
    return float(
        optimize.brentq(
            lambda t2: _generalised_q(y, v, t2) - target,
            0.0,
            hi,
            xtol=np.finfo(float).tiny,
            rtol=4.0 * np.finfo(float).eps,
            maxiter=500,
        )
    )


def _q_profile(
    y: np.ndarray, v: np.ndarray, conf_level: float, start: float
) -> tuple[float, float]:
    """The Q-profile interval for tau^2: where the generalised Q meets two chi-square quantiles.

    A bound whose quantile lies above Q at tau^2 = 0 is 0, so both are 0 when Q is below the
    lower-tail quantile, as metafor's `confint` reports.
    """
    df = len(y) - 1
    alpha = 1.0 - conf_level
    q_upper_tail = float(stats.chi2.ppf(1.0 - alpha / 2.0, df))
    q_lower_tail = float(stats.chi2.ppf(alpha / 2.0, df))
    q_zero = _generalised_q(y, v, 0.0)
    lower = 0.0 if q_zero < q_upper_tail else _q_root(y, v, q_upper_tail, start)
    upper = 0.0 if q_zero < q_lower_tail else _q_root(y, v, q_lower_tail, start)
    return lower, upper


def _row(summary: CurveSummaries, record: str, z: float) -> RecordRow:
    se = summary.rmst_se_s
    if se > 0.0:
        lower, upper = summary.rmst_s - z * se, summary.rmst_s + z * se
        return RecordRow(
            record, summary, lower, upper, pooled=True, not_pooled_reason=None
        )
    return RecordRow(
        record, summary, None, None, pooled=False, not_pooled_reason=ZERO_VARIANCE
    )


def random_effects_summary(
    tables: Sequence[EventTable],
    *,
    tau_rule: str,
    conf_level: float = CONF_LEVEL,
) -> RandomEffectsSummary:
    """DerSimonian-Laird summary of the records' restricted means at the comparison tau.

    Order: records with no windows are excluded (`comparison_tau`); tau is taken over the
    rest; each record's restricted mean y_i and variance v_i = se_i^2 are `curve_summaries`
    at tau. A record with v_i = 0 has no death before tau, so its restricted mean is tau and
    its inverse-variance weight infinite: it is listed with no interval and not pooled.
    With k pooled records, w_i = 1 / v_i, Q = sum w_i (y_i - y_w)^2 about the w-weighted
    mean, C = sum w_i - sum w_i^2 / sum w_i, and tau^2 = max(0, (Q - (k - 1)) / C). C is
    summed as 2 sum_{i<j} w_i w_j / sum w_i, which equals it and takes no difference of
    near-equal terms when one weight dominates. The pooled mean weights by 1 / (v_i + tau^2),
    with SE 1 / sqrt(sum of those weights) and a normal interval. I^2 = tau^2 / (tau^2 +
    (k - 1) / C), which is max(0, (Q - (k - 1)) / Q) here. The prediction interval is
    mu +/- t_{k-2} sqrt(SE^2 + tau^2), for a new record's true restricted mean. metafor
    5.2.1 (`rma(method = "DL")`, `confint`, `predict(predtype = "Riley")`) gives the same
    numbers on the fixtures. Reference: no source located for the estimator (usually
    credited to DerSimonian and Laird, not opened here); metafor's `?confint.rma` credits
    the Q-profile to Viechtbauer (2007) and `?predict.rma` the k - 2 degrees of freedom to
    Riley et al. (2011), neither opened here.

    Validity assumptions

    - Assumption: each restricted mean is normal about its record's true mean, with
      variance known and equal to v_i.
      Diagnostic: none here; the SE is Greenwood-type and poorest with few windows.
      Consequence of violation: the CI, the Q-profile interval and the prediction interval
      miss their stated coverage.
      Reference: metafor 5.2.1 `?confint.rma`, Note (the Q-profile is exact under
      normal effects with known sampling variances).
    - Assumption: the records' true restricted means are normal about one pooled mean.
      Diagnostic: none here; with three to ten records nothing here could reject it.
      Consequence of violation: the prediction interval is not where a new record's true
      restricted mean falls.
      Reference: metafor 5.2.1 `?predict.rma`, Note.
    - Assumption: tau^2 is known when the pooled CI is formed.
      Diagnostic: none here. The interval is normal, metafor's default; with few records
      it ignores how uncertain tau^2 is, and how far its coverage falls short here is not
      measured (OPEN).
      Consequence of violation: a pooled CI narrower than its stated level.
      Reference: no source located.
    - Assumption: the records left out for zero variance are not the ones that matter.
      Diagnostic: they are listed and counted. They are the records with no death before
      tau, the longest-lived, so the pooled mean leans short whenever any are left out.
      Consequence of violation: a pooled restricted mean biased down, by an amount not
      corrected here.
      Reference: no source located.
    - The assumptions of the module docstring.
    """
    z = z_two_sided(conf_level)
    tau = comparison_tau(tables, tau_rule=tau_rule)
    common = dict(
        threshold_label=tables[0].threshold_label,
        side=tables[0].side,
        conf_level=conf_level,
        tau=tau,
    )
    by_name = {t.dataset_id: t for t in tables}
    rows = tuple(
        _row(
            curve_summaries(
                by_name[name], rmst_tau_s=tau.rmst_tau_s, conf_level=conf_level
            ),
            name,
            z,
        )
        for name in tau.included
    )
    pooled = [row for row in rows if row.pooled]
    k = len(pooled)
    if k < MIN_POOLED:
        return RandomEffectsSummary(
            **common,  # type: ignore[arg-type]
            rows=rows,
            k_pooled=k,
            pooled_s=None,
            pooled_se_s=None,
            pooled_lower_s=None,
            pooled_upper_s=None,
            between_var_s2=None,
            between_var_lower_s2=None,
            between_var_upper_s2=None,
            q_statistic=None,
            q_df=None,
            i2_fraction=None,
            prediction_lower_s=None,
            prediction_upper_s=None,
            not_pooled_reason=(
                f"{k} records with a variance; at least {MIN_POOLED} are needed to pool"
            ),
        )
    y = np.array([row.summary.rmst_s for row in pooled], dtype=float)
    v = np.array([row.summary.rmst_se_s**2 for row in pooled], dtype=float)
    w = 1.0 / v
    q = _generalised_q(y, v, 0.0)
    df = k - 1
    c = 2.0 * float(np.sum(np.triu(np.outer(w, w), k=1))) / float(np.sum(w))
    between_var = max(0.0, (q - df) / c)
    w_star = 1.0 / (v + between_var)
    mu = float(np.sum(w_star * y) / np.sum(w_star))
    se = float(1.0 / np.sqrt(np.sum(w_star)))
    lower_var, upper_var = _q_profile(
        y, v, conf_level, start=max(between_var, float(np.mean(v)))
    )
    t_crit = float(stats.t.ppf(0.5 + conf_level / 2.0, k - 2))
    spread = t_crit * math.sqrt(se**2 + between_var)
    return RandomEffectsSummary(
        **common,  # type: ignore[arg-type]
        rows=rows,
        k_pooled=k,
        pooled_s=mu,
        pooled_se_s=se,
        pooled_lower_s=mu - z * se,
        pooled_upper_s=mu + z * se,
        between_var_s2=between_var,
        between_var_lower_s2=lower_var,
        between_var_upper_s2=upper_var,
        q_statistic=q,
        q_df=df,
        i2_fraction=between_var / (between_var + df / c),
        prediction_lower_s=mu - spread,
        prediction_upper_s=mu + spread,
        not_pooled_reason=None,
    )
