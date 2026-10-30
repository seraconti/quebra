"""Comparisons across records: a log-rank test and a restricted-mean difference for a pair.

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
log-rank statistic, `summary(survfit, rmean = tau)` for each record's restricted mean. Both
are written by `jobs/rscripts/comparison_reference.R`.

Validity assumptions

- Assumption: the windows of each record are exchangeable draws from one law, with
  independent censoring.
  Diagnostic: the check outcome attached to each record (`checks_a`, `checks_b`); empty
  reads NOT ASSESSED, and no step here fills it.
  Consequence of violation: the log-rank p and the SE of the difference are miscalibrated;
  the point difference is still the difference of the two areas.
  Reference: `analyzers/assumptions.py`, `a1_renewal_durations`.
- Assumption: the two records are independent of each other.
  Diagnostic: none here. Two qubits read in the same night share an environment, and two
  halves of one scan share a qubit and a night.
  Consequence of violation: the log-rank variance and the SE of the difference are wrong,
  in a direction no step here measures.
  Reference: no source located.
- Assumption: the two records were read on comparable schedules.
  Diagnostic: none here. A window's age lives on its record's read grid
  (`docs/WINDOW_SEMANTICS.md`), so two read spacings shift the ages apart by construction.
  Consequence of violation: a difference in read spacing reads as a difference in hazard.
  Reference: `docs/WINDOW_SEMANTICS.md`.
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
from scipy import stats

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
