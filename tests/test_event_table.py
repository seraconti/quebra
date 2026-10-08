"""The event table every survival estimator reads: its risk table and its intervals."""

from __future__ import annotations

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from quebra.analyzers import event_table, windows


@pytest.mark.properties
@given(
    data=st.lists(
        st.tuples(st.sampled_from([0.0, 1.0, 2.0, 2.0, 3.0, 5.0]), st.booleans()),
        min_size=1,
        max_size=30,
    )
)
def test_the_risk_table_matches_a_brute_force_count(data) -> None:
    """Oracle: counting by definition, at every distinct age, ties included.

    Deaths at a are the observed windows of age a, censorings the others, and the number at
    risk counts every window of age a or more: a window censored at a is still at risk at a.
    """
    age = np.array([a for a, _ in data])
    event = np.array([e for _, e in data], dtype=bool)
    ages, n_events, n_censored, n_at_risk = event_table.risk_table(age, event)

    expected = sorted(set(age.tolist()))
    assert ages.tolist() == expected
    assert len(n_events) == len(n_censored) == len(n_at_risk) == len(expected)
    for a, d, c, n in zip(ages, n_events, n_censored, n_at_risk, strict=True):
        assert d == int(np.sum((age == a) & event))
        assert c == int(np.sum((age == a) & ~event))
        assert n == int(np.sum(age >= a))


@pytest.mark.unit
@pytest.mark.parametrize(
    "event", [np.array([1.0, 0.0]), [True, False]], ids=["floats", "a_list"]
)
def test_the_risk_table_refuses_an_event_flag_that_is_not_a_bool_array(event) -> None:
    """Oracle: specification. A bool cast reads any nonzero flag as a death."""
    with pytest.raises(TypeError, match="boolean"):
        event_table.risk_table(np.array([1.0, 2.0]), event)


def _regular_record():
    t = np.arange(20.0)
    v = np.array(
        [9, 9, 1, 1, 9, 9, 9, 1, 9, 1, 1, 1, 9, 9, 1, 9, 9, 9, 1, 9], dtype=float
    )
    return windows.run(
        windows.WindowsInputs(
            t_rel_s=t, values=v, thresholds=[("5", 5.0, True)], dataset_id="grid"
        )
    )


# Spacings that all differ, with one gap (index 9 to 10), so a birth-side spacing is never
# equal to a death-side one. In spec is a value of 9, out of spec 1.
_IRREGULAR_STEPS_S = [
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    1000,
    11,
    12,
    13,
    14,
    15,
    16,
    17,
    18,
    19,
]
_IRREGULAR_VALUES = [9, 1, 1, 9, 9, 9, 1, 9, 9, 9, 1, 1, 9, 9, 1, 9, 9, 9, 1, 1]
# Per side, each window whose birth is observed, read off the values by hand: (s, e, True)
# for a death, e the first read past the crossing; (s, last, False) for a censored window.
_IRREGULAR_WINDOWS = {
    windows.SIDE_IN_SPEC: [(3, 6, True), (7, 9, False), (12, 14, True), (15, 18, True)],
    windows.SIDE_OUT_OF_SPEC: [
        (1, 3, True),
        (6, 7, True),
        (14, 15, True),
        (18, 19, False),
    ],
}


@pytest.mark.unit
@pytest.mark.parametrize("side", event_table.SIDES)
def test_each_length_interval_is_read_off_the_window_s_read_indices(side) -> None:
    """Oracle: the definition, applied to read indices fixed by hand.

    A birth lies in [t[s-1], t[s]] and a death in [t[e-1], t[e]], so a death's length lies
    in [t[e-1] - t[s], t[e] - t[s-1]] and its KM age is t[e] - t[s]. A censored window's
    length lies in [t[last] - t[s], inf] and its age is the lower end. Every spacing
    differs, so a formula that swapped the birth-side and death-side spacing would move a
    bound; each side holds a censored window, at a gap on one and the scan end on the other.
    """
    t = np.concatenate([[0.0], np.cumsum(_IRREGULAR_STEPS_S, dtype=float)])
    carved = windows.run(
        windows.WindowsInputs(
            t_rel_s=t,
            values=np.asarray(_IRREGULAR_VALUES, dtype=float),
            thresholds=[("5", 5.0, True)],
        )
    )
    expected = _IRREGULAR_WINDOWS[side]
    death = np.array([d for _, _, d in expected])
    s = np.array([w[0] for w in expected])
    end = np.array([w[1] for w in expected])

    table = event_table.from_carve(carved, threshold_label="5", side=side)

    np.testing.assert_array_equal(table.event, death)
    np.testing.assert_array_equal(
        table.age_lo_s, np.where(death, t[end - 1], t[end]) - t[s]
    )
    np.testing.assert_array_equal(
        table.age_hi_s, np.where(death, t[end] - t[s - 1], np.inf)
    )
    np.testing.assert_array_equal(table.age_s, t[end] - t[s])


@pytest.mark.unit
def test_windows_whose_birth_was_not_seen_are_dropped_and_counted() -> None:
    """Oracle: the carve's own birth_observed column, side by side.

    The kept rows are asserted, not only counted: dropping a different window of the same
    count would otherwise pass.
    """
    carved = _regular_record()
    for side, frame in (
        (windows.SIDE_IN_SPEC, carved.windows_in_spec),
        (windows.SIDE_OUT_OF_SPEC, carved.windows_out_of_spec),
    ):
        table = event_table.from_carve(carved, threshold_label="5", side=side)
        kept = frame[frame["birth_observed"]]
        assert table.n_windows == len(kept)
        assert table.n_unobserved_birth_dropped == int((~frame["birth_observed"]).sum())
        np.testing.assert_array_equal(table.age_s, kept["duration_s"].to_numpy())
        np.testing.assert_array_equal(table.event, ~kept["censored"].to_numpy())
        np.testing.assert_array_equal(
            table.age_lo_s, (kept["t_last_s"] - kept["t_birth_s"]).to_numpy()
        )


@pytest.mark.unit
def test_an_empty_side_is_an_empty_table_but_a_label_off_the_ladder_raises() -> None:
    """Oracle: specification. A record always in spec has no out-of-spec windows."""
    t = np.arange(6.0)
    carved = windows.run(
        windows.WindowsInputs(
            t_rel_s=t, values=np.full(6, 9.0), thresholds=[("5", 5.0, True)]
        )
    )
    empty = event_table.from_carve(
        carved, threshold_label="5", side=windows.SIDE_OUT_OF_SPEC
    )
    assert empty.n_windows == 0 and len(empty.event_age_s) == 0
    assert empty.event.dtype == np.bool_ and empty.age_s.dtype == np.float64
    assert empty.age_lo_s.dtype == empty.age_hi_s.dtype == np.float64
    with pytest.raises(KeyError, match="not on the ladder"):
        event_table.from_carve(carved, threshold_label="6", side=windows.SIDE_IN_SPEC)
    # Without the carve's ladder the label must appear in the table, and on an empty
    # side it cannot, so a typo and an empty side are not confused.
    with pytest.raises(KeyError, match="not on the ladder"):
        event_table.from_windows(
            carved.windows_out_of_spec,
            threshold_label="5",
            side=windows.SIDE_OUT_OF_SPEC,
        )


@pytest.mark.unit
def test_a_record_too_short_to_carve_keeps_its_label_and_two_empty_tables() -> None:
    """Oracle: specification. The ladder is the carve's own, not the labels it carved.

    A one-read record carves nothing, and its threshold is still on the ladder with an
    empty table on each side.
    """
    carved = windows.run(
        windows.WindowsInputs(
            t_rel_s=np.array([0.0]),
            values=np.array([9.0]),
            thresholds=[("5", 5.0, True)],
        )
    )
    tables = event_table.event_tables_from_carve(carved)

    assert tables.ladder == ("5",)
    assert set(tables.tables) == {("5", side) for side in event_table.SIDES}
    for table in tables.tables.values():
        assert table.n_windows == 0 and table.n_unobserved_birth_dropped == 0


@pytest.mark.unit
def test_the_wrong_side_and_a_non_boolean_event_are_refused() -> None:
    """Oracle: specification. Either would put the wrong windows or deaths on a curve."""
    carved = _regular_record()
    with pytest.raises(ValueError, match="asked for the out_of_spec side"):
        event_table.from_windows(
            carved.windows_in_spec,
            threshold_label="5",
            side=windows.SIDE_OUT_OF_SPEC,
        )
    with pytest.raises(TypeError, match="boolean"):
        event_table.EventTable(
            dataset_id="",
            threshold_label="5",
            side=windows.SIDE_IN_SPEC,
            age_s=np.array([1.0, 2.0]),
            event=np.array([1.0, 0.0]),
        )


@pytest.mark.unit
@pytest.mark.parametrize("flag", ["birth_observed", "censored"])
@pytest.mark.parametrize("coded_as", ["strings", "objects_with_nan"])
def test_a_flag_column_that_is_not_bool_is_refused(flag, coded_as) -> None:
    """Oracle: specification. A bool cast reads NaN and the string "False" as True.

    Read that way, an unobserved birth would enter as a lifetime and a censored window as
    a death, so the column must already be bool.
    """
    frame = _regular_record().windows_in_spec.copy()
    if coded_as == "strings":
        frame[flag] = frame[flag].map({True: "True", False: "False"})
    else:
        frame[flag] = frame[flag].astype(object)
        frame.loc[0, flag] = np.nan
    with pytest.raises(TypeError, match=f"{flag!r} must be bool"):
        event_table.from_windows(frame, threshold_label="5", side=windows.SIDE_IN_SPEC)


def _table(**overrides):
    kwargs = dict(
        dataset_id="",
        threshold_label="5",
        side=windows.SIDE_IN_SPEC,
        age_s=np.array([1.0, 2.0]),
        event=np.array([True, False]),
        age_lo_s=np.array([0.5, 2.0]),
        age_hi_s=np.array([1.5, np.inf]),
    )
    kwargs.update(overrides)
    return event_table.EventTable(**kwargs)


_REFUSED = {
    "unknown_side": ({"side": "sideways"}, ValueError, "side must be one of"),
    "event_a_list": ({"event": [True, False]}, TypeError, "boolean"),
    "event_float": ({"event": np.array([1.0, 0.0])}, TypeError, "boolean"),
    "event_short": ({"event": np.array([True])}, ValueError, "differ in length"),
    "age_nan": ({"age_s": np.array([np.nan, 2.0])}, ValueError, "finite and non-neg"),
    "age_inf": ({"age_s": np.array([1.0, np.inf])}, ValueError, "finite and non-neg"),
    "age_negative": (
        {"age_s": np.array([-1.0, 2.0])},
        ValueError,
        "finite and non-neg",
    ),
    "lo_short": ({"age_lo_s": np.array([0.5])}, ValueError, "match age_s in length"),
    "hi_short": ({"age_hi_s": np.array([1.5])}, ValueError, "match age_s in length"),
    "lo_nan": (
        {"age_lo_s": np.array([np.nan, 2.0])},
        ValueError,
        "lo_s must be finite",
    ),
    "lo_minus_inf": (
        {"age_lo_s": np.array([-np.inf, 2.0])},
        ValueError,
        "lo_s must be fin",
    ),
    "hi_nan": ({"age_hi_s": np.array([np.nan, np.inf])}, ValueError, "must not be NaN"),
    "lo_above_age": ({"age_lo_s": np.array([1.2, 2.0])}, ValueError, "in its interval"),
    "hi_below_age": (
        {"age_hi_s": np.array([0.8, np.inf])},
        ValueError,
        "in its interval",
    ),
    "death_hi_inf": (
        {"age_hi_s": np.array([np.inf, np.inf])},
        ValueError,
        "finite age_hi",
    ),
    "censored_hi_finite": ({"age_hi_s": np.array([1.5, 3.0])}, ValueError, "= inf"),
}


@pytest.mark.unit
@pytest.mark.parametrize("case", list(_REFUSED))
def test_the_constructor_refuses_each_malformed_table(case) -> None:
    """Oracle: specification. Each check the constructor states, one bad field at a time.

    The base table is accepted, so each refusal is its own override's. Direct construction
    is an entry point, and these checks are all that stand between it and an estimator.
    """
    assert _table().n_windows == 2
    overrides, error, message = _REFUSED[case]
    with pytest.raises(error, match=message):
        _table(**overrides)


@pytest.mark.unit
def test_the_table_does_not_share_its_arrays_with_the_caller() -> None:
    """Oracle: specification. The risk table is computed once, at construction.

    So the arrays it was computed from must not change under it when the caller reuses
    its own arrays.
    """
    age, event = np.array([1.0, 2.0]), np.array([True, True])
    lo, hi = age - 0.5, age + 0.5
    table = _table(age_s=age, event=event, age_lo_s=lo, age_hi_s=hi)

    age[0], event[1], lo[0], hi[1] = 9.0, False, -1.0, np.inf
    assert table.age_s.tolist() == [1.0, 2.0]
    assert table.event.tolist() == [True, True]
    assert table.age_lo_s.tolist() == [0.5, 1.5]
    assert table.age_hi_s.tolist() == [1.5, 2.5]
    assert table.n_deaths == int(table.n_events.sum()) == 2
