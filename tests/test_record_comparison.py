"""What the comparisons across records refuse, and what they report as not computed.

Oracle, unless a test names another (a closed form, or exact rational arithmetic): the
specification in `analyzers/record_comparison.py`. The numbers themselves are pinned to R
in `tests/test_record_comparison_r_reference.py`.
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
import pytest

from quebra.analyzers import record_comparison
from quebra.analyzers.check_attachment import RESTS_ON_EXCHANGEABILITY
from quebra.analyzers.event_table import EventTable
from quebra.analyzers.windows import SIDE_IN_SPEC, SIDE_OUT_OF_SPEC

pytestmark = pytest.mark.unit

RULE = record_comparison.TAU_RULE_MIN_LARGEST_AGE


def _table(
    name: str,
    ages: list[float],
    deaths: list[bool],
    *,
    side: str = SIDE_IN_SPEC,
    threshold_label: str = "3 µs",
) -> EventTable:
    return EventTable(
        dataset_id=name,
        threshold_label=threshold_label,
        side=side,
        age_s=np.array(ages, dtype=float),
        event=np.array(deaths, dtype=bool),
    )


def test_records_that_cannot_be_compared_are_refused() -> None:
    """Oracle: specification. One threshold, one side, one rule, each record once."""
    a = _table("a", [1.0, 2.0], [True, False])
    with pytest.raises(ValueError, match="thresholds"):
        record_comparison.compare_pair(
            a, _table("b", [1.0], [True], threshold_label="4 µs"), tau_rule=RULE
        )
    with pytest.raises(ValueError, match="sides"):
        record_comparison.compare_pair(
            a, _table("b", [1.0], [True], side=SIDE_OUT_OF_SPEC), tau_rule=RULE
        )
    with pytest.raises(ValueError, match="appears twice"):
        record_comparison.compare_pair(a, a, tau_rule=RULE)
    with pytest.raises(ValueError, match="tau_rule"):
        record_comparison.compare_pair(
            a, _table("b", [1.0], [True]), tau_rule="max_largest_age"
        )


def test_an_empty_record_or_an_undefined_log_rank_is_not_computed() -> None:
    """Oracle: specification. None with a reason, never NaN.

    A record with no windows leaves the whole pair not computed. With windows on both sides,
    the log-rank statistic is undefined with no death at all, and when its variance is 0;
    the restricted means are computed in both.
    """
    empty = _table("empty", [], [])
    pair = record_comparison.compare_pair(
        _table("a", [1.0, 2.0], [True, False]), empty, tau_rule=RULE
    )
    assert pair.not_computed_reason == f"empty: {record_comparison.NO_WINDOWS}"
    assert pair.logrank_not_computed_reason == pair.not_computed_reason
    assert pair.tau.rmst_tau_s == 2.0 and pair.tau.tau_set_by == "a"
    assert pair.logrank_p is None and pair.rmst_difference_s is None

    no_deaths = record_comparison.compare_pair(
        _table("a", [1.0, 2.0], [False, False]),
        _table("b", [3.0], [False]),
        tau_rule=RULE,
    )
    assert no_deaths.logrank_not_computed_reason == record_comparison.NO_DEATHS
    assert no_deaths.logrank_z is None and no_deaths.logrank_p is None
    assert no_deaths.rmst_difference_s == 0.0

    # b's one window dies at age 5, when it alone is at risk: n - 1 = 0 in every term.
    lone = record_comparison.compare_pair(
        _table("a", [1.0, 2.0], [False, False]),
        _table("b", [5.0], [True]),
        tau_rule=RULE,
    )
    assert lone.logrank_not_computed_reason == record_comparison.ZERO_LOGRANK_VARIANCE
    assert lone.logrank_chi2 is None and lone.logrank_p is None
    assert lone.rmst_difference_se_s is not None

    # Both windows at risk at age 1 die there: n - d = 0, and every term is 0 again.
    all_die = record_comparison.compare_pair(
        _table("a", [1.0], [True]), _table("b", [1.0], [True]), tau_rule=RULE
    )
    assert (
        all_die.logrank_not_computed_reason == record_comparison.ZERO_LOGRANK_VARIANCE
    )
    assert all_die.rmst_difference_se_s == 0.0


def _records() -> list[EventTable]:
    """Three records with deaths before every tau the rule can pick."""
    return [
        _table("r1", [1.0, 2.0, 4.0, 6.0], [True, True, False, True]),
        _table("r2", [1.0, 3.0, 3.0, 5.0, 7.0], [True, False, True, True, False]),
        _table("r3", [2.0, 2.0, 4.0, 8.0], [True, True, True, False]),
    ]


def test_each_record_interval_is_its_restricted_mean_plus_or_minus_z_se() -> None:
    """Oracle: closed form, rmst +/- z se with z = 1.959963984540054 at 0.95."""
    summary = record_comparison.random_effects_summary(_records(), tau_rule=RULE)
    for row in summary.rows:
        half = 1.959963984540054 * row.summary.rmst_se_s
        assert row.lower_s == pytest.approx(row.summary.rmst_s - half, rel=1e-12)
        assert row.upper_s == pytest.approx(row.summary.rmst_s + half, rel=1e-12)


def test_fewer_than_three_records_with_a_variance_are_listed_not_pooled() -> None:
    """Oracle: specification. Below the floor every pooled number is None, with a reason."""
    two = _records()[:2] + [_table("empty", [], [])]
    summary = record_comparison.random_effects_summary(two, tau_rule=RULE)
    assert summary.tau.excluded == (("empty", record_comparison.NO_WINDOWS),)
    assert [row.record for row in summary.rows] == ["r1", "r2"]
    assert summary.k_pooled == 2
    assert summary.not_pooled_reason is not None
    assert summary.pooled_s is None and summary.between_var_s2 is None
    assert summary.prediction_lower_s is None and summary.i2_fraction is None


def test_a_record_with_no_death_before_tau_is_listed_at_tau_and_not_pooled() -> None:
    """Oracle: specification. Its weight would be infinite, so it is drawn with no interval.

    Its restricted mean is tau exactly, and the three other records are pooled.
    """
    lasting = _table("lasting", [6.5, 9.0], [True, False])
    summary = record_comparison.random_effects_summary(
        _records() + [lasting], tau_rule=RULE
    )
    assert summary.tau.rmst_tau_s == 6.0 and summary.tau.tau_set_by == "r1"
    row = next(r for r in summary.rows if r.record == "lasting")
    assert row.summary.rmst_s == 6.0 and row.summary.rmst_se_s == 0.0
    assert not row.pooled and row.not_pooled_reason == record_comparison.ZERO_VARIANCE
    assert row.lower_s is None and row.upper_s is None
    assert summary.k_pooled == 3 and summary.pooled_s is not None
    assert summary.standing == RESTS_ON_EXCHANGEABILITY


def test_a_summary_with_no_record_left_still_says_not_assessed() -> None:
    """Oracle: specification (AGENTS.md section 5). An empty attachment reads NOT ASSESSED."""
    summary = record_comparison.random_effects_summary(
        [_table("empty", [], [])], tau_rule=RULE
    )
    assert summary.rows == () and summary.k_pooled == 0
    assert "NOT ASSESSED" in summary.check_summary()


def _spread() -> list[EventTable]:
    """Three records whose restricted means sit far apart: tau^2 and its lower bound > 0."""
    return [
        _table("short", [1.0, 1.0, 2.0, 2.0, 3.0, 10.0], [True] * 5 + [False]),
        _table("middle", [3.0, 4.0, 5.0, 5.0, 6.0, 10.0], [True] * 5 + [False]),
        _table("long", [6.0, 7.0, 8.0, 9.0, 9.0, 10.0], [True] * 5 + [False]),
    ]


def test_every_interval_follows_conf_level_closed_form_at_90_percent() -> None:
    """Oracle: closed form at 0.90 with k = 3: z, t_1 = tan(0.45 pi), chi2_2 = -2 ln(p).

    The Q-profile bounds are checked through a generalised Q written here: at each bound it
    must equal its chi-square quantile with 2 degrees of freedom.
    """
    summary = record_comparison.random_effects_summary(
        _spread(), tau_rule=RULE, conf_level=0.90
    )
    z = 1.6448536269514722
    for row in summary.rows:
        half = z * row.summary.rmst_se_s
        assert row.upper_s - row.summary.rmst_s == pytest.approx(half, rel=1e-12)
    mu, se, between = summary.pooled_s, summary.pooled_se_s, summary.between_var_s2
    assert mu is not None and se is not None and between is not None and between > 0
    assert summary.pooled_upper_s - mu == pytest.approx(z * se, rel=1e-12)
    spread = math.tan(0.45 * math.pi) * math.sqrt(se**2 + between)
    assert summary.prediction_upper_s - mu == pytest.approx(spread, rel=1e-12)

    y = np.array([row.summary.rmst_s for row in summary.rows])
    v = np.array([row.summary.rmst_se_s**2 for row in summary.rows])

    def q(between_var: float) -> float:
        w = 1.0 / (v + between_var)
        return float(np.sum(w * (y - np.sum(w * y) / np.sum(w)) ** 2))

    lower, upper = summary.between_var_lower_s2, summary.between_var_upper_s2
    assert lower is not None and lower > 0 and upper is not None
    assert q(lower) == pytest.approx(-2.0 * math.log(0.05), rel=1e-9)
    assert q(upper) == pytest.approx(-2.0 * math.log(0.95), rel=1e-9)


def test_the_dl_variance_matches_exact_rationals_under_a_dominant_weight() -> None:
    """Oracle: DerSimonian-Laird in exact rationals (`fractions`) on the same y_i and v_i.

    Record a's one death falls 1e-6 s before tau, so its variance is about 1e-14 of the
    others' and sum w - sum w^2 / sum w would be a difference of near-equal floats.
    """
    tables = [
        _table("a", [10.0 - 1e-6] + [12.0] * 9, [True] + [False] * 9),
        _table("b", [1.0, 2.0, 3.0, 4.0, 10.0], [True] * 4 + [False]),
        _table("c", [5.0, 6.0, 7.0, 8.0, 11.0], [True] * 4 + [False]),
    ]
    summary = record_comparison.random_effects_summary(tables, tau_rule=RULE)
    y = [Fraction(row.summary.rmst_s) for row in summary.rows]
    v = [Fraction(row.summary.rmst_se_s**2) for row in summary.rows]
    assert max(v) / min(v) > 10**13
    w = [1 / vi for vi in v]
    total = sum(w)
    mean = sum(wi * yi for wi, yi in zip(w, y)) / total
    q = sum(wi * (yi - mean) ** 2 for wi, yi in zip(w, y))
    c = total - sum(wi * wi for wi in w) / total
    assert q > 2
    assert summary.between_var_s2 == pytest.approx(float((q - 2) / c), rel=1e-12)
