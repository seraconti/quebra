"""What the comparisons across records refuse, and what they report as not computed.

Oracle for every test: the specification in `analyzers/record_comparison.py`. The numbers
themselves are pinned to R in `tests/test_record_comparison_r_reference.py`.
"""

from __future__ import annotations

import numpy as np
import pytest

from quebra.analyzers import record_comparison
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
