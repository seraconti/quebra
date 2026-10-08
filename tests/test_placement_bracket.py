"""The placement bracket: every placement of the crossings inside the reads stays within it."""

from __future__ import annotations

import itertools
import math

import numpy as np
import pytest
from hypothesis import assume, example, given
from hypothesis import strategies as st

from quebra.analyzers import event_table, kaplan_meier, survival, windows
from quebra.analyzers.event_table import EventTable

SIDE = windows.SIDE_IN_SPEC


def _km_at(age: np.ndarray, event: np.ndarray, at: np.ndarray) -> np.ndarray:
    ages, n_events, _c, n_at_risk = event_table.risk_table(age, event)
    time, surv, _n, _g = kaplan_meier.product_limit(ages, n_events, n_at_risk)
    return kaplan_meier._step_eval(time, surv, at)


def _bracket_at(bracket: survival.PlacementBracket, at: np.ndarray):
    lower = kaplan_meier._step_eval(bracket.time_s, bracket.lower, at)
    upper = kaplan_meier._step_eval(bracket.time_s, bracket.upper, at)
    return lower, upper


def _inside(bracket, age, event, at) -> bool:
    lower, upper = _bracket_at(bracket, at)
    km = _km_at(age, event, at)
    return bool(np.all(lower <= km + 1e-12) and np.all(km <= upper + 1e-12))


@pytest.mark.unit
def test_with_no_censoring_the_width_is_both_spacings_on_a_regular_grid() -> None:
    """Oracle: analytic. Each death's interval is [age - 1, age + 1] on a one-second grid.

    Reads 9, 10 and 21 held in spec make the four windows unequal: lower ends 4, 2, 3, 2
    and upper ends 6, 4, 5, 4. With tau past every upper end the areas are the means of
    the ends, 11/4 and 19/4, and the width is the mean of hi - lo, so exactly 2.
    """
    t = np.arange(30.0)
    v = np.where((t // 3) % 2 == 0, 9.0, 1.0)
    v[[9, 10, 21]] = 9.0
    carved = windows.run(
        windows.WindowsInputs(t_rel_s=t, values=v, thresholds=[("5", 5.0, True)])
    )
    full = event_table.from_carve(carved, threshold_label="5", side=SIDE)
    deaths = full.event
    table = EventTable(
        dataset_id="",
        threshold_label="5",
        side=SIDE,
        age_s=full.age_s[deaths],
        event=np.ones(int(deaths.sum()), dtype=bool),
        age_lo_s=full.age_lo_s[deaths],
        age_hi_s=full.age_hi_s[deaths],
    )
    bracket = survival.placement_bracket(
        table, rmst_tau_s=float(table.age_hi_s.max()) + 1.0
    )
    np.testing.assert_array_equal(table.age_lo_s, [4.0, 2.0, 3.0, 2.0])
    np.testing.assert_array_equal(table.age_hi_s, [6.0, 4.0, 5.0, 4.0])
    assert bracket.rmst_lower_s == 11.0 / 4 and bracket.rmst_upper_s == 19.0 / 4
    assert bracket.width_s == 2.0
    assert bracket.width_censoring_only_s == 0.0


@pytest.mark.unit
def test_the_bracket_curves_and_areas_match_a_hand_computation() -> None:
    """Oracle: a hand computation of both curves and both areas, censoring included.

    Deaths in [1, 3] and [2, 4], a censoring at 2.5. The lower curve counts every window
    still short of its lower end: 1, 2/3 from 1, 1/3 from 2, 0 from 2.5. The upper curve
    loses only deaths past their upper ends: 2/3 from 3, 1/3 from 4. To tau = 5 the areas
    are (1 + 2 + 2.5) / 3 and (3 + 4 + 5) / 3. Containment alone cannot see a looser
    bound; these values can.
    """
    table = EventTable(
        dataset_id="",
        threshold_label="h",
        side=SIDE,
        age_s=np.array([2.0, 3.0, 2.5]),
        event=np.array([True, True, False]),
        age_lo_s=np.array([1.0, 2.0, 2.5]),
        age_hi_s=np.array([3.0, 4.0, np.inf]),
    )
    bracket = survival.placement_bracket(table, rmst_tau_s=5.0)
    at = np.array([0.0, 1.0, 2.0, 2.5, 3.0, 4.0, 4.5])
    lower, upper = _bracket_at(bracket, at)
    np.testing.assert_allclose(lower, [1, 2 / 3, 1 / 3, 0, 0, 0, 0])
    np.testing.assert_allclose(upper, [1, 1, 1, 1, 2 / 3, 1 / 3, 1 / 3])
    assert bracket.rmst_lower_s == pytest.approx(5.5 / 3)
    assert bracket.rmst_upper_s == pytest.approx(12.0 / 3)
    # Lengths at the Kaplan-Meier ages 2, 3 and 2.5: areas 7.5 / 3 and (2 + 3 + 5) / 3.
    assert bracket.width_censoring_only_s == pytest.approx(2.5 / 3)


def _two_windows(hi: list[float], event: list[bool]) -> EventTable:
    lo = np.array([1.0, 3.0])
    hi_s = np.array(hi)
    return EventTable(
        dataset_id="",
        threshold_label="m",
        side=SIDE,
        age_s=np.where(event, (lo + np.where(event, hi_s, 0.0)) / 2.0, lo),
        event=np.array(event),
        age_lo_s=lo,
        age_hi_s=hi_s,
    )


@pytest.mark.unit
def test_each_median_is_its_own_curves_median_by_r_s_midpoint_rule() -> None:
    """Oracle: a hand computation of R's quantile rule on each curve's own steps.

    Deaths in [1, 2] and [3, 4]. `lower` sits at 1/2 on [1, 3), so its median is 2;
    `upper` sits at 1/2 on [2, 4), so its median is 3. The placements at the two ends,
    (1, 3) and (2, 4), have Kaplan-Meier medians 2 and 3: both bounds are attained.
    """
    bracket = survival.placement_bracket(_two_windows([2.0, 4.0], [True, True]))
    assert bracket.median_lower_s == 2.0
    assert bracket.median_upper_s == 3.0
    for ages, median in (([1.0, 3.0], 2.0), ([2.0, 4.0], 3.0)):
        exact = EventTable(
            dataset_id="",
            threshold_label="m",
            side=SIDE,
            age_s=np.array(ages),
            event=np.array([True, True]),
        )
        assert survival.curve_summaries(exact, quantiles=(0.5,)).quantile_s == (median,)


@pytest.mark.unit
def test_the_upper_median_is_none_when_half_the_windows_are_censored() -> None:
    """Oracle: hand computation. A death in [1, 2] and a window censored at 3 or later.

    `upper` falls to 1/2 at 2 and stays there, so a censored placement at c gives a
    Kaplan-Meier median of (2 + c) / 2 by R's rule, unbounded in c. `lower` sits at 1/2
    on [1, 3), so its median is 2.
    """
    bracket = survival.placement_bracket(_two_windows([2.0, np.inf], [True, False]))
    assert bracket.median_lower_s == 2.0
    assert bracket.median_upper_s is None


@pytest.mark.unit
def test_every_placement_on_a_small_grid_stays_inside_the_bracket() -> None:
    """Oracle: exhaustive enumeration, on the case where Kaplan-Meier is not monotone.

    Deaths in [0.5, 2.5] and [3.0, 3.5], a censoring at 2.0 or later. Moving the first death
    from 1 to 2.5 LOWERS S(2.7), so Kaplan-Meier on the lower ends is not a lower bound.
    """
    table = EventTable(
        dataset_id="",
        threshold_label="f4",
        side=SIDE,
        age_s=np.array([1.5, 3.25, 2.0]),
        event=np.array([True, True, False]),
        age_lo_s=np.array([0.5, 3.0, 2.0]),
        age_hi_s=np.array([2.5, 3.5, np.inf]),
    )
    bracket = survival.placement_bracket(table)
    at = np.linspace(0.0, 4.0, 801)
    first = np.arange(0.75, 2.5001, 0.25)
    second = np.array([3.25, 3.5])
    censor = np.array([2.0, 2.5])
    for d1, d2, c in itertools.product(first, second, censor):
        age = np.array([d1, d2, c])
        assert _inside(bracket, age, table.event, at), (d1, d2, c)


def _known_crossings(rng: np.random.Generator):
    """A step-valued record with known crossing times, read on a grid with one 20 s gap.

    Every excursion lasts at least 2.5 read spacings, so outside the gap the grid misses
    none (the bracket assumes the carve's windows are the true ones). Excursions inside
    the gap go unseen; the gap policy censors the window open at the gap and starts the
    next one with an unobserved birth, which the event table drops.
    """
    durations = 2.5 + rng.exponential(3.0, 200)
    crossings = np.cumsum(durations) - 7.0
    t = rng.uniform(0.0, 1.0) + np.arange(300.0)
    t = t[(t < 120.0) | (t > 140.0)]
    index = np.searchsorted(crossings, t, side="right")
    v = np.where(index % 2 == 0, 9.0, 1.0)
    return t, v, crossings


def _exact_ages(frame, crossings):
    seen = frame[frame["birth_observed"]]
    k = np.searchsorted(crossings, seen["t_birth_s"].to_numpy(), side="right")
    birth = crossings[k - 1]
    death = crossings[k]
    censored = seen["censored"].to_numpy(dtype=bool)
    age = np.where(censored, seen["t_last_s"].to_numpy() - birth, death - birth)
    return age, ~censored


@pytest.mark.statistical
@pytest.mark.parametrize("side", event_table.SIDES)
def test_the_curve_on_the_true_crossing_times_lies_inside_the_bracket(side) -> None:
    """Oracle: simulation truth. The crossing times are known, the record is carved for real.

    Kaplan-Meier on the exact lengths (and exact censoring times) of the windows the carve
    saw must lie inside the bracket built from the reads, at every age, in every trial.
    """
    for seed in range(20):
        t, v, crossings = _known_crossings(np.random.default_rng(seed))
        carved = windows.run(
            windows.WindowsInputs(t_rel_s=t, values=v, thresholds=[("5", 5.0, True)])
        )
        frame = (
            carved.windows_in_spec
            if side == windows.SIDE_IN_SPEC
            else carved.windows_out_of_spec
        )
        table = event_table.from_carve(carved, threshold_label="5", side=side)
        bracket = survival.placement_bracket(table)
        age, event = _exact_ages(frame, crossings)
        at = np.linspace(0.0, float(np.max(table.age_s)) + 2.0, 600)
        assert _inside(bracket, age, event, at), seed


@pytest.mark.properties
@given(
    rows=st.lists(
        st.tuples(
            st.floats(0.0, 20.0),
            st.booleans(),
            st.floats(0.0, 2.0),
            st.floats(0.0, 2.0),
        ),
        min_size=1,
        max_size=25,
    )
)
@example(rows=[(3.95, True, 0.0, 0.0)])
@example(rows=[(0.9500000000000002, True, 0.0, 0.0)])
def test_the_shipped_kaplan_meier_lies_inside_the_bracket(rows) -> None:
    """Oracle: the bound in `placement_bracket`'s docstring, on arbitrary event tables.

    The curve is in minutes and the bracket in seconds. The shipped survival values are
    read at the table's own death ages in seconds, so no step moves by a division by 60
    (both explicit examples are ones that did). An example where that division merges two
    death ages has no seconds clock for its steps and is skipped.
    """
    age = np.array([r[0] for r in rows])
    event = np.array([r[1] for r in rows], dtype=bool)
    lo = np.where(event, np.maximum(age - [r[2] for r in rows], 0.0), age)
    hi = np.where(event, age + [r[3] for r in rows], np.inf)
    table = EventTable(
        dataset_id="",
        threshold_label="p",
        side=SIDE,
        age_s=age,
        event=event,
        age_lo_s=lo,
        age_hi_s=hi,
    )
    bracket = survival.placement_bracket(table)
    curve = kaplan_meier.run(
        kaplan_meier.make_inputs_from_event_table(table, label="p")
    )
    at = np.linspace(0.0, 25.0, 501)
    lower, upper = _bracket_at(bracket, at)
    step_s = np.r_[0.0, table.event_age_s[table.n_events > 0]]
    assume(len(step_s) == len(curve.time_min))
    km = kaplan_meier._step_eval(step_s, curve.survival, at)
    assert np.all(lower <= km + 1e-12) and np.all(km <= upper + 1e-12)


@pytest.mark.unit
def test_an_empty_table_gives_an_empty_bracket() -> None:
    """Oracle: specification. An empty side is data, not an error."""
    table = EventTable(
        dataset_id="",
        threshold_label="e",
        side=SIDE,
        age_s=np.array([]),
        event=np.array([], dtype=bool),
    )
    bracket = survival.placement_bracket(table)
    assert bracket.n_windows == 0 and math.isnan(bracket.width_s)
    assert bracket.median_lower_s is None


@pytest.mark.unit
@pytest.mark.parametrize("n_windows", [0, 1])
@pytest.mark.parametrize("tau", [-1.0, math.inf, math.nan])
def test_an_invalid_tau_raises_whether_or_not_the_table_is_empty(
    n_windows, tau
) -> None:
    """Oracle: specification. A restriction time must be finite and non-negative."""
    table = EventTable(
        dataset_id="",
        threshold_label="t",
        side=SIDE,
        age_s=np.ones(n_windows),
        event=np.ones(n_windows, dtype=bool),
    )
    with pytest.raises(ValueError, match="rmst_tau_s must be finite"):
        survival.placement_bracket(table, rmst_tau_s=tau)
