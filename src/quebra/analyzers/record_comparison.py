"""Comparisons across records: a pair's log-rank test and restricted-mean difference, and a
random-effects summary of many records' restricted means.

Every step here takes the event tables of several records at one threshold and side, where
`analyzers/survival.py` takes one. The record is named by its table's `dataset_id`; two
tables of one comparison must agree on threshold label and side and name different records,
or the step raises. A record never sits beside one of its halves (`HALF_SUFFIXES`): the
halves are row splits of one scan, so the same reads would be compared with themselves.

The per-ladder steps at the end take the records' carves, refuse carves made on different
ladders or under different gap rules (a label alone does not say which threshold it names),
and run the comparison at every threshold and side. `pair_table` and
`random_effects_table` flatten their results into one long table, one number per row, with
the job's evidence line for each pair or group: what the records' metadata shows was held
fixed between them. The figure-data builders at the end hand those results to
`plots/record_comparison_plot.py`; the pair figure's builder adds each record's Kaplan-Meier curve,
taken on the same event tables.

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
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace

import numpy as np
import pandas as pd
from scipy import optimize, stats

from quebra.analyzers.check_attachment import (
    DESCRIPTIVE,
    RESTS_ON_EXCHANGEABILITY,
    CheckAttachment,
)
from quebra.analyzers.event_table import (
    SIDES,
    EventTable,
    EventTables,
    event_tables_from_carve,
)
from quebra.analyzers.kaplan_meier import (
    CONF_LEVEL,
    KaplanMeierCurve,
    curve_from_event_table,
    z_two_sided,
)
from quebra.analyzers.survival import ByThresholdSide, CurveSummaries, curve_summaries
from quebra.analyzers.windows import WindowsResult
from quebra.core._artifact_guard import StaleArtifactGuard

TAU_RULE_MIN_LARGEST_AGE = "min_largest_age"
TAU_RULES = (TAU_RULE_MIN_LARGEST_AGE,)
# Where a forest's tau is taken: over each group's own records, or over every record of
# every group, so that groups drawn on one axis share one tau.
TAU_PER_GROUP = "per_group"
TAU_ACROSS_GROUPS = "across_groups"
TAU_SCOPES = (TAU_PER_GROUP, TAU_ACROSS_GROUPS)
# A record named `<record>_before` or `<record>_after` is a half of `<record>`.
HALF_SUFFIXES = ("_before", "_after")

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


def require_no_record_beside_its_halves(names: Sequence[str]) -> None:
    """Raise when a record and one of its halves are both in `names`."""
    present = set(names)
    for name in names:
        for suffix in HALF_SUFFIXES:
            if name.endswith(suffix) and name[: -len(suffix)] in present:
                raise ValueError(
                    f"{name} is a half of {name[: -len(suffix)]}: a record and its "
                    f"halves never sit in one comparison"
                )


def require_disjoint_groups(groups: Mapping[str, Sequence[str]]) -> list[str]:
    """Every record of `groups` in order; raise when one is named twice or beside its half.

    The halves rule spans all the groups, since every group may share one tau.
    """
    everyone = [name for names in groups.values() for name in names]
    if len(set(everyone)) != len(everyone):
        raise ValueError(f"a record appears twice across the groups: {everyone}")
    require_no_record_beside_its_halves(everyone)
    return everyone


def _require_comparable(tables: Sequence[EventTable]) -> None:
    """One threshold, one side, every record named once and none beside its halves."""
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
    require_no_record_beside_its_halves(names)


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

    `n_windows_*` and `n_deaths_*` are each record's counts at this threshold and side, read
    off its event table, so they are set even when nothing else is computed.

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
    n_windows_a: int
    n_deaths_a: int
    n_windows_b: int
    n_deaths_b: int
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
        n_windows_a=table_a.n_windows,
        n_deaths_a=table_a.n_deaths,
        n_windows_b=table_b.n_windows,
        n_deaths_b=table_b.n_deaths,
    )
    if tau.excluded:
        reason: str | None = "; ".join(f"{r}: {why}" for r, why in tau.excluded)
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


def _tau_over(
    own: ComparisonTau,
    tables: Sequence[EventTable],
    tau_from: Sequence[EventTable],
    tau_rule: str,
) -> ComparisonTau:
    """`own` with its tau taken over `tau_from`, which must hold every record of `tables`."""
    wider = comparison_tau(tau_from, tau_rule=tau_rule)
    if (tau_from[0].threshold_label, tau_from[0].side) != (
        tables[0].threshold_label,
        tables[0].side,
    ):
        raise ValueError(
            "tau taken at another threshold or side than the records pooled"
        )
    missing = sorted({t.dataset_id for t in tables} - {t.dataset_id for t in tau_from})
    if missing:
        raise ValueError(f"tau taken over records that leave out {missing}")
    if not own.included:
        return own
    return replace(own, rmst_tau_s=wider.rmst_tau_s, tau_set_by=wider.tau_set_by)


def random_effects_summary(
    tables: Sequence[EventTable],
    *,
    tau_rule: str,
    conf_level: float = CONF_LEVEL,
    tau_from: Sequence[EventTable] | None = None,
) -> RandomEffectsSummary:
    """DerSimonian-Laird summary of the records' restricted means at the comparison tau.

    `tau_from`, when given, holds every record of `tables` and more: tau is then the
    smallest largest age over all of them, and `tau.tau_set_by` may name a record that is
    not pooled here. `tau.included` and `tau.excluded` still name this summary's records.

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
    if tau_from is not None:
        tau = _tau_over(tau, tables, tau_from, tau_rule)
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


# ---------------------------------------------------------------------------
# Every threshold and side at once: the steps `comparison_recipe` wires
# ---------------------------------------------------------------------------


def _tables_on_one_ladder(
    carves: Sequence[WindowsResult],
) -> tuple[tuple[str, ...], list[EventTables]]:
    """The ladder every carve shares and each carve's event tables, or raise.

    A table carries a threshold's label and not its value, so records carved on two
    ladders could meet under one label at two thresholds. Every carve must have the same
    labels in the same order, each label one value and one direction wherever a carve has a
    window at it, and one gap rule, since the gap rule moves window boundaries.
    """
    tables = [event_tables_from_carve(carve) for carve in carves]
    ladders = {t.ladder for t in tables}
    if len(ladders) != 1:
        raise ValueError(f"records carved on different ladders: {sorted(ladders)}")
    gap_mults = {float(carve.meta["gap_mult"]) for carve in carves}  # type: ignore[arg-type]
    if len(gap_mults) != 1:
        raise ValueError(
            f"records carved under different gap rules: {sorted(gap_mults)}"
        )
    meanings: dict[str, set[tuple[float, bool]]] = {}
    for carve in carves:
        for frame in (carve.windows_in_spec, carve.windows_out_of_spec):
            columns = frame[["threshold_label", "threshold_value", "big_values_good"]]
            for label, value, big in columns.drop_duplicates().itertuples(index=False):
                meanings.setdefault(str(label), set()).add((float(value), bool(big)))
    clashes = {label: sorted(m) for label, m in meanings.items() if len(m) > 1}
    if clashes:
        raise ValueError(f"one label names two thresholds across records: {clashes}")
    (ladder,) = ladders
    return ladder, tables


def compare_pair_set(
    carve_a: WindowsResult,
    carve_b: WindowsResult,
    *,
    tau_rule: str,
    conf_level: float = CONF_LEVEL,
) -> ByThresholdSide[PairComparison]:
    """`compare_pair` of two records at every threshold and side of their one ladder.

    Validity assumptions

    - Assumption: each p is read at its own threshold and side.
      Diagnostic: none here; the comparison table notes every p as unadjusted.
      Consequence of violation: the smallest p over the ladder and both sides, read as one
      test's p, overstates the evidence, since each cell had its own chance to fall low.
      Reference: no source located.
    - The assumptions of `compare_pair`, at each threshold and side.
    """
    ladder, (a, b) = _tables_on_one_ladder((carve_a, carve_b))
    return ByThresholdSide(
        ladder=ladder,
        items={
            (label, side): compare_pair(
                a.table(label, side),
                b.table(label, side),
                tau_rule=tau_rule,
                conf_level=conf_level,
            )
            for label in ladder
            for side in SIDES
        },
    )


@dataclass
class RandomEffectsSet:
    """One random-effects summary per group, threshold and side.

    `groups` maps each group to its records, as the job declared them. `tau_scope` says
    whether a group's tau is taken over its own records (`TAU_PER_GROUP`) or over every
    record of every group (`TAU_ACROSS_GROUPS`), at each threshold and side.
    """

    tau_rule: str
    tau_scope: str
    conf_level: float
    groups: dict[str, tuple[str, ...]]
    ladder: tuple[str, ...]
    items: dict[tuple[str, str, str], RandomEffectsSummary]

    def get(self, group: str, threshold_label: str, side: str) -> RandomEffectsSummary:
        key = (group, threshold_label, side)
        if key not in self.items:
            raise KeyError(
                f"no summary for {key}; groups {list(self.groups)}, ladder "
                f"{list(self.ladder)}"
            )
        return self.items[key]


def random_effects_set(
    *carves: WindowsResult,
    groups: Mapping[str, Sequence[str]],
    tau_scope: str,
    tau_rule: str,
    conf_level: float = CONF_LEVEL,
) -> RandomEffectsSet:
    """`random_effects_summary` of each group at every threshold and side of one ladder.

    A carve is matched to its record by `meta["dataset_id"]`, so two carves naming one record
    raise. The groups pass `require_disjoint_groups`.

    Validity assumptions

    - Assumption: under `TAU_ACROSS_GROUPS`, the record that sets tau, the shortest-lived of
      any group, leaves each group's records time to die before it.
      Diagnostic: every group's `tau.tau_set_by` names that record, and each group's rows
      name the records left out for zero variance.
      Consequence of violation: more records have no death before tau and leave the pool,
      so every group's pooled restricted mean leans further short (the zero-variance
      assumption of `random_effects_summary`).
      Reference: no source located.
    - The assumptions of `random_effects_summary`, at each group, threshold and side.
    """
    if tau_scope not in TAU_SCOPES:
        raise ValueError(f"tau_scope must be one of {TAU_SCOPES}; got {tau_scope!r}")
    everyone = require_disjoint_groups(groups)
    ladder, tables = _tables_on_one_ladder(carves)
    carved = [str(c.meta["dataset_id"]) for c in carves]
    if len(set(carved)) != len(carved):
        raise ValueError(f"two carves name one record: {carved}")
    by_name = dict(zip(carved, tables))
    members = {group: tuple(names) for group, names in groups.items()}
    unknown = sorted(set(everyone) - set(by_name))
    if unknown:
        raise KeyError(f"groups name records with no carve: {unknown}")
    items: dict[tuple[str, str, str], RandomEffectsSummary] = {}
    for label in ladder:
        for side in SIDES:
            wider = (
                [by_name[name].table(label, side) for name in everyone]
                if tau_scope == TAU_ACROSS_GROUPS
                else None
            )
            for group, names in members.items():
                items[(group, label, side)] = random_effects_summary(
                    [by_name[name].table(label, side) for name in names],
                    tau_rule=tau_rule,
                    conf_level=conf_level,
                    tau_from=wider,
                )
    return RandomEffectsSet(
        tau_rule=tau_rule,
        tau_scope=tau_scope,
        conf_level=conf_level,
        groups=members,
        ladder=ladder,
        items=items,
    )


# ---------------------------------------------------------------------------
# The comparison table: one number per row
# ---------------------------------------------------------------------------

TABLE_COLUMNS = (
    "comparison",
    "group",
    "record",
    "threshold_label",
    "side",
    "quantity",
    "value",
    "note",
)
POOLED = "pooled"
UNADJUSTED = "unadjusted for any other comparison"
_RMST_QUANTITIES = ("rmst_s", "rmst_se_s")
_COUNT_QUANTITIES = ("n_windows", "n_deaths")
_RECORD_QUANTITIES = _RMST_QUANTITIES + _COUNT_QUANTITIES
_POOLED_QUANTITIES = (
    "k_pooled",
    "pooled_s",
    "pooled_se_s",
    "pooled_lower_s",
    "pooled_upper_s",
    "between_var_s2",
    "between_var_lower_s2",
    "between_var_upper_s2",
    "q_statistic",
    "q_df",
    "i2_fraction",
    "prediction_lower_s",
    "prediction_upper_s",
)


def _frame(rows: list[tuple[object, ...]]) -> pd.DataFrame:
    # Object columns, so a number that was not computed stays None and never reads as NaN.
    return pd.DataFrame(rows, columns=list(TABLE_COLUMNS), dtype=object)


def _tau_note(tau: ComparisonTau, scope: str | None = None) -> str:
    rule = tau.tau_rule if scope is None else f"{tau.tau_rule}, {scope}"
    if tau.rmst_tau_s is None:
        why = "; ".join(f"{record}: {reason}" for record, reason in tau.excluded)
        return f"{rule}; no record with windows ({why})"
    return f"{rule}, set by {tau.tau_set_by}"


def require_evidence(evidence: Mapping[str, str], groups: Sequence[str]) -> None:
    """Raise unless `evidence` holds one non-blank string per group and no other key."""
    missing = sorted(set(groups) - set(evidence))
    unknown = sorted(set(evidence) - set(groups))
    blank = sorted(
        g
        for g, line in evidence.items()
        if not isinstance(line, str) or not line.strip()
    )
    if missing or unknown or blank:
        raise ValueError(
            f"evidence needs one line per group: missing {missing}, unknown {unknown}, "
            f"blank {blank}"
        )


def require_flags(flags: Mapping[str, str], records: Sequence[str]) -> None:
    """Raise unless every flag names a record of `records` with a non-blank string."""
    unknown = sorted(set(flags) - set(records))
    blank = sorted(
        r for r, line in flags.items() if not isinstance(line, str) or not line.strip()
    )
    if unknown or blank:
        raise ValueError(
            f"a flag names a record of the comparison with a line of text: unknown "
            f"{unknown}, blank {blank}"
        )


def pair_table(
    pairs: ByThresholdSide[PairComparison],
    *,
    comparison: str,
    evidence: str,
    flags: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    """Every number of a pair at every threshold and side, one per row (`TABLE_COLUMNS`).

    `record` names the record a number belongs to, or "a vs b" for the pair's own numbers.
    A number not computed is None, and its row's note gives the reason. A pair not computed
    forms no record summary, so its restricted means are None, while each record's window
    and death counts, read off its own event table, are still given, with no note; its tau
    row keeps the tau the rule gave over the record with windows, and its note says nothing
    was computed at it.

    `evidence` says what the records' metadata shows was held fixed between them; it is
    the note of an `evidence` row at every threshold and side, and enters no number, which
    is why it is a kwarg of this step and not of the comparison. `flags` maps a record to
    a line said about it alone, such as how the paper's table marks it; a flagged record
    gets a `flag` row with that line as its note, at every threshold and side.
    """
    require_evidence({comparison: evidence}, [comparison])
    flags = dict(flags or {})
    first = next(iter(pairs.items.values()))
    require_flags(flags, [first.record_a, first.record_b])
    rows: list[tuple[object, ...]] = []
    for (label, side), pair in pairs.items.items():
        both = f"{pair.record_a} vs {pair.record_b}"
        tau_note = _tau_note(pair.tau)
        if pair.not_computed_reason is not None and pair.tau.rmst_tau_s is not None:
            tau_note += f"; nothing computed at it ({pair.not_computed_reason})"
        cells: list[tuple[str, str, object, str | None]] = [
            (both, "tau_s", pair.tau.rmst_tau_s, tau_note)
        ]
        for record, summary, counts in (
            (pair.record_a, pair.summary_a, (pair.n_windows_a, pair.n_deaths_a)),
            (pair.record_b, pair.summary_b, (pair.n_windows_b, pair.n_deaths_b)),
        ):
            cells += [
                (
                    record,
                    quantity,
                    None if summary is None else getattr(summary, quantity),
                    pair.not_computed_reason,
                )
                for quantity in _RMST_QUANTITIES
            ]
            cells += [
                (record, quantity, count, None)
                for quantity, count in zip(_COUNT_QUANTITIES, counts, strict=True)
            ]
            if record in flags:
                cells.append((record, "flag", None, flags[record]))
        cells += [
            (both, quantity, getattr(pair, quantity), pair.not_computed_reason)
            for quantity in ("rmst_difference_s", "rmst_difference_se_s")
        ]
        cells += [
            (
                both,
                quantity,
                getattr(pair, quantity),
                pair.logrank_not_computed_reason
                or (UNADJUSTED if quantity == "logrank_p" else None),
            )
            for quantity in ("logrank_z", "logrank_chi2", "logrank_p")
        ]
        cells.append((both, "check_outcome", None, pair.check_summary()))
        cells.append((both, "evidence", None, evidence))
        rows += [
            (comparison, "", record, label, side, quantity, value, note)
            for record, quantity, value, note in cells
        ]
    return _frame(rows)


def random_effects_table(
    summaries: RandomEffectsSet,
    *,
    comparison: str,
    evidence: Mapping[str, str],
    flags: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    """Every number of each group's summary at every threshold and side, one per row.

    `record` is "pooled" for the summary's own numbers: tau, with the rule, the scope and
    the record that set it in the note, and `conf_level`, the level of every interval in
    the group's rows. A record left out of the pool keeps its restricted mean with the
    reason in the note; a record with no windows has a None restricted mean, window and
    death counts of 0, and its reason. `evidence` holds one line per group and `flags` one
    line per flagged record, as in `pair_table`.
    """
    require_evidence(evidence, list(summaries.groups))
    flags = dict(flags or {})
    require_flags(flags, [r for records in summaries.groups.values() for r in records])
    rows: list[tuple[object, ...]] = []
    for (group, label, side), summary in summaries.items.items():
        cells: list[tuple[str, str, object, str | None]] = [
            (
                POOLED,
                "tau_s",
                summary.tau.rmst_tau_s,
                _tau_note(summary.tau, summaries.tau_scope),
            ),
            (POOLED, "conf_level", summary.conf_level, None),
        ]
        cells += [
            (POOLED, quantity, getattr(summary, quantity), summary.not_pooled_reason)
            for quantity in _POOLED_QUANTITIES
        ]
        for row in summary.rows:
            cells += [
                (
                    row.record,
                    quantity,
                    getattr(row.summary, quantity),
                    row.not_pooled_reason,
                )
                for quantity in _RECORD_QUANTITIES
            ]
            cells += [
                (row.record, "lower_s", row.lower_s, row.not_pooled_reason),
                (row.record, "upper_s", row.upper_s, row.not_pooled_reason),
            ]
            if row.record in flags:
                cells.append((row.record, "flag", None, flags[row.record]))
        for record, reason in summary.tau.excluded:
            if reason != NO_WINDOWS:
                raise ValueError(
                    f"{record} was excluded for {reason!r}: its counts are known here "
                    f"only when it has no windows"
                )
            cells += [
                (record, "rmst_s", None, reason),
                (record, "n_windows", 0, reason),
                (record, "n_deaths", 0, reason),
            ]
            if record in flags:
                cells.append((record, "flag", None, flags[record]))
        cells.append((POOLED, "check_outcome", None, summary.check_summary()))
        cells.append((POOLED, "evidence", None, evidence[group]))
        rows += [
            (comparison, group, record, label, side, quantity, value, note)
            for record, quantity, value, note in cells
        ]
    return _frame(rows)


# ---------------------------------------------------------------------------
# Figure data: what `plots/record_comparison_plot.py` draws, one artifact per figure
# ---------------------------------------------------------------------------


@dataclass
class PairPanel:
    """One threshold and side of a pair figure: both records' curves and the pair itself.

    The curves are `kaplan_meier.curve_from_event_table` on the event tables the comparison
    read, so a drawn curve and the numbers printed beside it rest on the same windows.
    """

    curve_a: KaplanMeierCurve
    curve_b: KaplanMeierCurve
    pair: PairComparison


@dataclass
class PairLadderData(StaleArtifactGuard):
    """Everything a pair figure draws: one `PairPanel` per threshold and side of the ladder."""

    comparison: str
    record_a: str
    record_b: str
    evidence: str
    ladder: tuple[str, ...]
    panels: dict[tuple[str, str], PairPanel]
    flags: dict[str, str]


def build_pair_ladder_data(
    pairs: ByThresholdSide[PairComparison],
    carve_a: WindowsResult,
    carve_b: WindowsResult,
    *,
    comparison: str,
    evidence: str,
    flags: Mapping[str, str] | None = None,
) -> PairLadderData:
    """The pair figure's data: each comparison of `pairs` beside both records' curves.

    The carves must be the two records `pairs` compared, a first, on the ladder `pairs`
    holds; anything else raises. `flags` is as in `pair_table`.
    """
    require_evidence({comparison: evidence}, [comparison])
    ladder, (a, b) = _tables_on_one_ladder((carve_a, carve_b))
    if ladder != pairs.ladder:
        raise ValueError(f"carves on the ladder {ladder}, pairs on {pairs.ladder}")
    panels: dict[tuple[str, str], PairPanel] = {}
    for (label, side), pair in pairs.items.items():
        table_a, table_b = a.table(label, side), b.table(label, side)
        carved = (table_a.dataset_id, table_b.dataset_id)
        if carved != (pair.record_a, pair.record_b):
            raise ValueError(
                f"carves of {carved} given for the pair {pair.record_a} vs "
                f"{pair.record_b}"
            )
        panels[(label, side)] = PairPanel(
            curve_a=curve_from_event_table(
                table_a, label=pair.record_a, conf_level=pair.conf_level
            ),
            curve_b=curve_from_event_table(
                table_b, label=pair.record_b, conf_level=pair.conf_level
            ),
            pair=pair,
        )
    if not panels:
        raise ValueError("the pairs hold no threshold and side to draw")
    first = next(iter(panels.values())).pair
    flags = dict(flags or {})
    require_flags(flags, [first.record_a, first.record_b])
    return PairLadderData(
        comparison=comparison,
        record_a=first.record_a,
        record_b=first.record_b,
        evidence=evidence,
        ladder=ladder,
        panels=panels,
        flags=flags,
    )


def _birth_unseen(summary: RandomEffectsSummary) -> tuple[int, int]:
    """(windows dropped for an unseen birth, windows carved) over the records with windows."""
    dropped = sum(row.summary.n_unobserved_birth_dropped for row in summary.rows)
    return dropped, dropped + sum(row.summary.n_windows for row in summary.rows)


@dataclass
class ForestData(StaleArtifactGuard):
    """Each group's summary at one threshold and side, in declared order, as a forest."""

    comparison: str
    threshold_label: str
    side: str
    tau_scope: str
    conf_level: float
    evidence: dict[str, str]
    summaries: dict[str, RandomEffectsSummary]
    # Per group: (windows dropped for an unseen birth, windows carved), records with
    # windows only; a record with none is named in its summary's `tau.excluded`.
    birth_unseen: dict[str, tuple[int, int]]
    flags: dict[str, str]


def build_forest_data(
    summaries: RandomEffectsSet,
    *,
    comparison: str,
    evidence: Mapping[str, str],
    threshold_label: str,
    side: str,
    flags: Mapping[str, str] | None = None,
) -> ForestData:
    """The forest at one threshold and side of `summaries`; an unknown one raises.

    Which threshold and side are kwargs of this step alone, so drawing another leaves
    the comparison step and its artifact unchanged. `flags` is as in `random_effects_table`.
    """
    require_evidence(evidence, list(summaries.groups))
    flags = dict(flags or {})
    require_flags(flags, [r for records in summaries.groups.values() for r in records])
    return ForestData(
        comparison=comparison,
        threshold_label=threshold_label,
        side=side,
        tau_scope=summaries.tau_scope,
        conf_level=summaries.conf_level,
        evidence={group: evidence[group] for group in summaries.groups},
        summaries={
            group: summaries.get(group, threshold_label, side)
            for group in summaries.groups
        },
        birth_unseen={
            group: _birth_unseen(summaries.get(group, threshold_label, side))
            for group in summaries.groups
        },
        flags=flags,
    )


@dataclass
class PooledLadderData(StaleArtifactGuard):
    """Each group's summary at every threshold and side, as the pooled ladder draws it."""

    comparison: str
    tau_scope: str
    conf_level: float
    ladder: tuple[str, ...]
    groups: tuple[str, ...]
    evidence: dict[str, str]
    summaries: dict[tuple[str, str, str], RandomEffectsSummary]
    # As `ForestData.birth_unseen`, per group, threshold and side.
    birth_unseen: dict[tuple[str, str, str], tuple[int, int]]


def build_pooled_ladder_data(
    summaries: RandomEffectsSet, *, comparison: str, evidence: Mapping[str, str]
) -> PooledLadderData:
    """`summaries` with each group's evidence line, for the pooled ladder figure."""
    require_evidence(evidence, list(summaries.groups))
    return PooledLadderData(
        comparison=comparison,
        tau_scope=summaries.tau_scope,
        conf_level=summaries.conf_level,
        ladder=summaries.ladder,
        groups=tuple(summaries.groups),
        evidence={group: evidence[group] for group in summaries.groups},
        summaries=dict(summaries.items),
        birth_unseen={
            key: _birth_unseen(summary) for key, summary in summaries.items.items()
        },
    )
