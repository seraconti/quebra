"""Calibration cycles and their mean time between calibrations (MTBC), against hand tables
and R.

The point and the SE come from `survival.curve_summaries`, whose agreement with R is pinned
in `tests/test_survival_r_reference.py`. What these tests pin is what `mtbf.mtbc` adds or
chooses (tau as the longest cycle whether or not it is censored, and the interval) and what
`analyzers/cycles.py` builds from a record of spells: the cycles, the battery's segments and
the share of time per status.

The hand-built record: 1000 s, two outages, one pair of touching calibrating rows, a
`benchmarking` row that starts the instant a calibration ends and plays no part in the
cycles, and a last calibration that leaves a cycle censored at the record end. Its cycles
are the `qi_cycles_hand` case of `jobs/rscripts/survival_reference.R`.
"""

from __future__ import annotations

from statistics import NormalDist

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import quebra.analyzers.checks.c6_exchangeability as c6
from quebra.analyzers import cycles, mtbf
from quebra.analyzers.check_ledger import VERDICT_ABSENT, VERDICT_NOT_COMPUTED
from quebra.analyzers.check_outcome import PERMUTATION_KEYS
from quebra.analyzers.event_table import CYCLES, EventTable
from quebra.analyzers.kaplan_meier import curve_from_event_table
from quebra.core.types import Norm
from tests.fixtures import R_SURVIVAL_INPUTS, R_SURVIVAL_VALUES


def _cycles(age: list[float], event: list[bool]) -> EventTable:
    return EventTable(
        dataset_id="unit",
        threshold_label="",
        side=CYCLES,
        age_s=np.asarray(age, dtype=float),
        event=np.asarray(event, dtype=bool),
    )


@pytest.mark.statistical
def test_mtbc_hand_computed_censored_case() -> None:
    """Oracle: a hand table. Cycles 30 (death), 35 (censored), 48 (death), 60 (censored).

    tau = 60, the longest cycle though it is censored. S is 1 on [0, 30), 3/4 on [30, 48)
    (4 at risk, 1 death) and 3/8 on [48, 60] (2 at risk, 1 death), so
    RMST = 30 + 18 (3/4) + 12 (3/8) = 48. The areas after the two deaths are 18 and 4.5, so
    Var = 18^2 / (4 3) + 4.5^2 / (2 1) = 27 + 10.125 = 37.125. z comes from the standard
    library's normal, not from the code under test.
    """
    table = _cycles([30.0, 35.0, 48.0, 60.0], [True, False, True, False])
    m = mtbf.mtbc(table)
    se = float(np.sqrt(37.125))
    z = NormalDist().inv_cdf(0.975)
    assert m.tau_s == 60.0
    assert m.rmst_s == pytest.approx(48.0, rel=1e-12)
    assert m.rmst_se_s == pytest.approx(se, rel=1e-12)
    assert m.lower_s == pytest.approx(48.0 - z * se, rel=1e-12)
    assert m.upper_s == pytest.approx(48.0 + z * se, rel=1e-12)
    assert (m.n_cycles, m.n_deaths, m.n_censored) == (4, 2, 2)
    # Nothing was checked here: the interval carries the empty a2 attachment.
    assert m.checks.assumption_id == "a2_renewal_cycles"
    assert m.checks.summary(subject="interval") == (
        "independence checks: NOT ASSESSED for this interval"
    )
    z90 = NormalDist().inv_cdf(0.95)
    m90 = mtbf.mtbc(table, conf_level=0.90)
    assert m90.conf_level == 0.90
    assert m90.lower_s == pytest.approx(48.0 - z90 * se, rel=1e-12)
    assert m90.upper_s == pytest.approx(48.0 + z90 * se, rel=1e-12)


@pytest.mark.unit
def test_mtbc_refuses_a_table_of_windows() -> None:
    """A restricted mean of windows on a side of a threshold is a mean window length."""
    windows_table = EventTable(
        dataset_id="unit",
        threshold_label="5",
        side="in_spec",
        age_s=np.array([1.0, 2.0]),
        event=np.array([True, True]),
    )
    with pytest.raises(ValueError, match="event table of cycles"):
        mtbf.mtbc(windows_table)


@pytest.mark.statistical
def test_mtbc_without_censoring_is_the_arithmetic_mean() -> None:
    """Oracle: a closed form. With every cycle complete and tau the longest, the area under
    the empirical survival curve is the sample mean, ties and 0 s cycles included."""
    rng = np.random.default_rng(5)
    age = np.round(rng.gamma(2.0, 30.0, 57), 0)
    age[:3] = 0.0
    m = mtbf.mtbc(mtbf.intervals_event_table(age, dataset_id="unit"))
    assert m.n_censored == 0
    assert m.tau_s == float(np.max(age))
    assert m.rmst_s == pytest.approx(float(np.mean(age)), rel=1e-12)


# A unix-second base, so the hand values are exercised at the magnitude of a real record.
_T0 = 1_700_000_000.0
# (start, end, status), in seconds from the record start.
_HAND_ROWS = [
    (50, 60, "calibrating"),
    (60, 100, "benchmarking"),
    (130, 140, "calibrating"),
    (140, 150, "calibrating"),
    (200, 210, "calibrating"),
    (300, 400, "offline"),
    (420, 430, "calibrating"),
    (470, 475, "calibrating"),
    (560, 570, "calibrating"),
    (600, 610, "calibrating"),
    (690, 750, "offline"),
    (800, 810, "calibrating"),
]


def _spell_norm(rows: list[tuple[float, float, str]], length_s: float = 1000.0) -> Norm:
    return Norm(
        {
            "spell_start_unix_s": np.array([_T0 + r[0] for r in rows], dtype=float),
            "spell_end_unix_s": np.array([_T0 + r[1] for r in rows], dtype=float),
            "spell_status": np.array([r[2] for r in rows]),
            "meta": {
                "dataset_id": "hand",
                "record_start_unix_s": _T0,
                "record_end_unix_s": _T0 + length_s,
            },
        }
    )


@pytest.mark.unit
def test_cycles_hand_table_with_an_outage_and_a_merge() -> None:
    """Oracle: a hand table of the record in the module docstring.

    Segment 0 is [0, 300): 0-50 has no known start (dropped); 130-140 and 140-150 touch and
    merge; cycles 60-130 (70) and 150-200 (50) are complete, the first starting at 60, where
    calibrating ends, not at 100, where the benchmarking that follows ends; 210-300 (90) is
    censored by the outage. Segment 1 is [400, 690): 400-420 dropped; 40, 85 and 30
    complete; 610-690 (80) censored. Segment 2 is [750, 1000]: 750-800 dropped; 810-1000
    (190) censored at the end.
    """
    table = cycles.cycles(_spell_norm(_HAND_ROWS))
    assert table.length_s.tolist() == [70.0, 50.0, 90.0, 40.0, 85.0, 30.0, 80.0, 190.0]
    assert table.event.tolist() == [True, True, False, True, True, True, False, False]
    assert table.segment.tolist() == [0, 0, 0, 1, 1, 1, 1, 2]
    assert (table.start_unix_s - _T0).tolist() == [
        60,
        150,
        210,
        430,
        475,
        570,
        610,
        810,
    ]
    assert table.n_segments == 3
    assert table.n_unobserved_birth_dropped == 3
    assert table.n_calibrations_merged == 1
    event_table = cycles.cycle_event_table(table)
    assert event_table.side == CYCLES and event_table.n_unobserved_birth_dropped == 3


def _contract_breach(**changes: object) -> Norm:
    """The hand record with one key, array or meta field changed (None deletes it)."""
    norm = _spell_norm(_HAND_ROWS)
    meta = dict(norm["meta"])
    for key, value in changes.items():
        target = meta if key in meta or key == "dataset_id" else norm
        if value is None:
            del target[key]
        else:
            target[key] = value
    norm["meta"] = meta
    return norm


_STATUS = np.array([r[2] for r in _HAND_ROWS], dtype=object)
_START = np.array([_T0 + r[0] for r in _HAND_ROWS], dtype=float)
_END = np.array([_T0 + r[1] for r in _HAND_ROWS], dtype=float)


def _with(array: np.ndarray, index: int, value: object) -> np.ndarray:
    out = array.copy()
    out[index] = value
    return out


@pytest.mark.unit
def test_a_record_that_breaks_the_spell_contract_is_refused() -> None:
    """One case per breach the contract names: each raises, nothing is sorted, clipped,
    coerced or defaulted."""
    cases = [
        ({"spell_status": None}, KeyError, "needs"),
        ({"dataset_id": None}, KeyError, "needs"),
        ({"spell_status": _with(_STATUS, 5, None)}, ValueError, "non-empty string"),
        ({"spell_status": _with(_STATUS, 1, "")}, ValueError, "non-empty string"),
        ({"spell_status": _with(_STATUS, 1, cycles.NOT_LISTED)}, ValueError, "no row"),
        ({"spell_end_unix_s": _END[:-1]}, ValueError, "differ in length"),
        ({"spell_start_unix_s": _with(_START, 2, np.nan)}, ValueError, "finite"),
        ({"record_end_unix_s": np.inf}, ValueError, "record start and end"),
        ({"record_end_unix_s": _T0}, ValueError, "end after it starts"),
        ({"spell_end_unix_s": _with(_END, 0, _T0 + 40)}, ValueError, "end before"),
        ({"spell_start_unix_s": _with(_START, 3, _T0 + 135)}, ValueError, "overlap"),
        ({"spell_start_unix_s": _with(_START, 0, _T0 - 5)}, ValueError, "outside"),
        ({"record_end_unix_s": _T0 + 805}, ValueError, "outside"),
        ({"unobserved_start_unix_s": np.array([_T0 + 10])}, KeyError, "both"),
        (_span([10.0], [20.0, 30.0]), ValueError, "differ in length"),
        (_span([np.nan], [20.0]), ValueError, "finite"),
        (_span([20.0], [20.0]), ValueError, "positive length"),
        (_span([40.0, 10.0], [45.0, 20.0]), ValueError, "time order"),
        (_span([990.0], [1010.0]), ValueError, "outside"),
        (_span([55.0], [58.0]), ValueError, "overlap an unobserved span"),
    ]
    for changes, error, match in cases:
        with pytest.raises(error, match=match):
            cycles.cycles(_contract_breach(**changes))


def _span(lo: list[float], hi: list[float]) -> dict[str, np.ndarray]:
    """Unobserved spans, in seconds from `_T0`, as contract keys."""
    return {
        "unobserved_start_unix_s": _T0 + np.asarray(lo, dtype=float),
        "unobserved_end_unix_s": _T0 + np.asarray(hi, dtype=float),
    }


@pytest.mark.unit
def test_a_status_history_that_breaks_its_contract_is_refused() -> None:
    """One case per refusal of `spells_from_samples`: nothing is sorted or defaulted."""
    good = [(0, "idle"), (30, "calibrating")]
    cases = [
        (_samples(good), {"spacing_s": 0.0}, "spacing_s"),
        (_samples([]), {"spacing_s": 30.0}, "no samples"),
        (_samples([(30, "idle"), (0, "idle")]), {"spacing_s": 30.0}, "time order"),
    ]
    for samples, kwargs, match in cases:
        with pytest.raises(ValueError, match=match):
            cycles.spells_from_samples(samples, **kwargs)
    short = _samples(good)
    short["sample_status"] = np.array(["idle"])
    with pytest.raises(ValueError, match="differ in length"):
        cycles.spells_from_samples(short, spacing_s=30.0)
    nan = _samples(good)
    nan["sample_time_unix_s"] = np.array([_T0, np.nan])
    with pytest.raises(ValueError, match="finite"):
        cycles.spells_from_samples(nan, spacing_s=30.0)
    no_id = _samples(good)
    no_id["meta"] = {}
    with pytest.raises(KeyError, match="needs"):
        cycles.spells_from_samples(no_id, spacing_s=30.0)
    for blank in (None, "", np.nan):
        missing = _samples(good)
        missing["sample_status"] = np.array(["idle", blank], dtype=object)
        with pytest.raises(ValueError, match="non-empty string"):
            cycles.spells_from_samples(missing, spacing_s=30.0)


@pytest.mark.unit
def test_cycles_hand_table_at_the_record_and_outage_boundaries() -> None:
    """Oracle: a hand table of a second record, 500 s, whose every boundary is a corner.

    Offline 0-20 (at the record start: no stretch before it); calibrating 20-30 (at the
    segment start: no residual); calibrating 70-80 (cycle 30-70, 40 s); offline 80-120
    (touching the calibration: no censored tail) and 120-150 (touching outages: no segment
    between); idle 160-200 in segment 150-210, which holds no calibration (dropped);
    offline 210-240; calibrating 300-310 (residual 240-300 dropped) and 400-500 (cycle
    310-400, 90 s, and the record ends with it: no tail). Segments 20-80, 150-210, 240-500.
    """
    rows = [
        (0, 20, "offline"),
        (20, 30, "calibrating"),
        (70, 80, "calibrating"),
        (80, 120, "offline"),
        (120, 150, "offline"),
        (160, 200, "idle"),
        (210, 240, "offline"),
        (300, 310, "calibrating"),
        (400, 500, "calibrating"),
    ]
    table = cycles.cycles(_spell_norm(rows, length_s=500.0))
    assert table.length_s.tolist() == [40.0, 90.0]
    assert table.event.tolist() == [True, True]
    assert table.segment.tolist() == [0, 2]
    assert table.n_segments == 3
    assert table.n_unobserved_birth_dropped == 2
    assert table.n_calibrations_merged == 0
    segments, n_dropped = cycles.battery_segments(table, 2)
    assert (segments, n_dropped) == ([], 2)
    shares = cycles.status_shares(_spell_norm(rows, length_s=500.0))
    assert shares.statuses == ("offline", "calibrating", "idle", cycles.NOT_LISTED)
    assert shares.seconds == (120.0, 130.0, 40.0, 210.0)


@pytest.mark.unit
def test_status_shares_match_a_hand_table() -> None:
    """Oracle: a hand table. calibrating 10+10+10+10+10+5+10+10+10 = 85 s, benchmarking
    40 s, offline 100+60 = 160 s; the time no row covers, summed gap by gap, is
    50+30+50+90+20+40+85+30+80+50+190 = 715 s; together the 1000 s record."""
    shares = cycles.status_shares(_spell_norm(_HAND_ROWS))
    assert shares.statuses == (
        "calibrating",
        "benchmarking",
        "offline",
        cycles.NOT_LISTED,
    )
    assert shares.seconds == (85.0, 40.0, 160.0, 715.0)
    assert shares.record_length_s == 1000.0
    assert math.fsum(shares.seconds) == shares.record_length_s
    assert shares.fractions == pytest.approx((0.085, 0.040, 0.160, 0.715), rel=1e-12)


@pytest.mark.unit
def test_battery_segments_hand_table() -> None:
    """Oracle: a hand table. Segment 0: x = (70, 50), tau = 120 + 90 = 210. Segment 1:
    x = (40, 85, 30), tau = 155 + 80 = 235. Segment 2 holds only its censored cycle, no
    complete one, so it is dropped."""
    segments, n_dropped = cycles.battery_segments(
        cycles.cycles(_spell_norm(_HAND_ROWS)), 2
    )
    assert n_dropped == 1
    assert [s.x.tolist() for s in segments] == [[70.0, 50.0], [40.0, 85.0, 30.0]]
    assert [s.tau for s in segments] == [210.0, 235.0]
    assert [s.n_censored_dropped for s in segments] == [1, 1]


@pytest.mark.statistical
def test_cycle_mtbc_matches_r_survfit() -> None:
    """Oracle: R `survival` 3.8.6 (`survfit`, log-log band, `summary(fit, rmean = tau)`)
    on the `qi_cycles_hand` case: the cycles built from the spell record, through the event
    table, give R's curve, band, restricted mean and its SE."""
    inputs = pd.read_csv(R_SURVIVAL_INPUTS)
    case = inputs[inputs["case"] == "qi_cycles_hand"].sort_values("i")
    r = _r_hand

    table = cycles.cycle_event_table(cycles.cycles(_spell_norm(_HAND_ROWS)))
    assert table.age_s.tolist() == case["age"].tolist()
    assert table.event.tolist() == (case["event"] == 1).tolist()
    m = mtbf.mtbc(table)
    # The tau R was asked for; R's restricted mean below is at it.
    assert m.tau_s == float(case["tau"].iloc[0])
    assert m.rmst_s == pytest.approx(float(r("rmean")[0]), rel=1e-9)
    assert m.rmst_se_s == pytest.approx(float(r("rmean_se")[0]), rel=1e-9)
    curve = curve_from_event_table(table)
    np.testing.assert_allclose(curve.time_min[1:] * 60.0, r("km_time"), rtol=1e-9)
    np.testing.assert_allclose(curve.survival[1:], r("km_surv"), rtol=1e-9)
    for ours, quantity in (
        (curve.band_lower, "km_lower"),
        (curve.band_upper, "km_upper"),
    ):
        theirs = r(quantity)
        np.testing.assert_array_equal(np.isnan(ours[1:]), np.isnan(theirs))
        keep = ~np.isnan(theirs)
        np.testing.assert_allclose(ours[1:][keep], theirs[keep], rtol=1e-9)


def _r_hand(quantity: str) -> np.ndarray:
    """One quantity R wrote for the `qi_cycles_hand` case, in its own order."""
    values = pd.read_csv(R_SURVIVAL_VALUES)
    rows = values[
        (values["case"] == "qi_cycles_hand") & (values["quantity"] == quantity)
    ]
    return rows.sort_values("i")["value"].to_numpy(dtype=float)


def _bench_table() -> pd.DataFrame:
    return pd.read_csv(
        Path(__file__).resolve().parents[1]
        / "jobs"
        / "bench"
        / "results"
        / "size_table.csv"
    )


def _checks(table: cycles.CycleTable) -> cycles.CycleChecks:
    return cycles.cycle_checks(
        table,
        _bench_table(),
        run_set=PERMUTATION_KEYS,
        alpha=0.05,
        min_events_pass=35,
        tie_cutoff_distinct=5,
        lag_max=5,
        n_permutations=999,
        seed=20261030,
        min_events_per_segment=2,
    )


@pytest.mark.statistical
def test_c5_rejects_an_alternating_cycle_sequence() -> None:
    """Oracle: a known truth for a2_renewal_cycles: alternating cycles are dependent.

    Forty cycles alternate short and long with distinct values (10, 50, 11, 51, ...), one
    segment, 1 s calibrations, a 5 s censored tail. Their lag-1 rank correlation is strongly
    negative by construction, so C5 must reject and C6's lag-1 statistic must be negative;
    sorting the cycles would keep C5 rejecting (it is two-sided) but flip that sign.
    """
    lengths = [v for k in range(20) for v in (10.0 + k, 50.0 + k)]
    rows: list[tuple[float, float, str]] = [(0.0, 1.0, "calibrating")]
    for length in lengths:
        start = rows[-1][1] + length
        rows.append((start, start + 1.0, "calibrating"))
    table = cycles.cycles(_spell_norm(rows, length_s=rows[-1][1] + 5.0))
    assert table.length_s[table.event].tolist() == lengths
    checks = _checks(table)
    c5_row = checks.rows[checks.rows["check_id"] == "c5_rank_autocorr"].iloc[0]
    assert c5_row["verdict"] == "fail", c5_row["notes"]
    # No trend in an alternation: C1 does not reject it, as it would if the battery were
    # fed the cycles sorted.
    c1 = checks.rows[
        (checks.rows["check_id"] == "c1_lewis_robinson")
        & (checks.rows["calibration"] == "permutation")
    ].iloc[0]
    assert c1["verdict"] != "fail", c1["notes"]
    segments, _ = cycles.battery_segments(table, 2)
    z = c6.run(segments, rng=np.random.default_rng(0), max_lag=5).extra["z"]
    assert z[0] < 0.0


@pytest.mark.unit
def test_cycles_panel_data_carries_the_counts_and_the_check_summary() -> None:
    """Oracle: the hand table of the 1000 s record, the ledger cells, and R's restricted
    mean on the `qi_cycles_hand` case."""
    norm = _spell_norm(_HAND_ROWS)
    table = cycles.cycles(norm)
    checks = _checks(table)
    data = cycles.build_cycles_panel_data(table, checks, cycles.status_shares(norm))
    assert set(checks.rows["clock"]) == {"calendar"}
    assert data.battery_note == ""
    assert data.checks.checks_unanswered == ()
    for grid, (label, p_value, verdict, _note) in zip(
        checks.outcome.grids, data.verdict_rows, strict=True
    ):
        check_id, calibration, variant = grid.key
        cell = checks.rows[
            (checks.rows["check_id"] == check_id)
            & (checks.rows["calibration"] == calibration)
            & (checks.rows["variant"].fillna("") == variant)
        ]
        assert len(cell) == 1, grid.key
        assert label == grid.label
        assert verdict == cell["verdict"].iloc[0] != VERDICT_ABSENT
        assert p_value == cell["p_value"].iloc[0]
    assert data.mtbc.rmst_s == pytest.approx(_r_hand("rmean")[0], rel=1e-9)
    data90 = cycles.build_cycles_panel_data(
        table, checks, cycles.status_shares(norm), conf_level=0.90
    )
    se = _r_hand("rmean_se")[0]
    assert data90.mtbc.conf_level == 0.90
    assert data90.mtbc.lower_s == pytest.approx(
        _r_hand("rmean")[0] - NormalDist().inv_cdf(0.95) * se, rel=1e-9
    )
    assert (data.n_cycles, data.n_censored) == (8, 3)
    assert (data.n_unobserved_birth_dropped, data.n_calibrations_merged) == (3, 1)
    assert (data.n_segments, data.n_segments_dropped) == (3, 1)
    assert data.mtbc.checks.assumption_id == "a2_renewal_cycles"
    assert data.checks == data.mtbc.checks
    assert data.check_summary == data.checks.summary()
    assert data.check_summary.startswith("5 checks asked")
    assert [row[0] for row in data.verdict_rows] == [
        grid.label for grid in checks.outcome.grids
    ]
    assert len(data.verdict_rows) == 5
    assert data.bench_note == cycles.BENCH_TRANSFER_NOTE
    m = data.mtbc
    assert data.mtbc_min == (m.rmst_s / 60.0, m.lower_s / 60.0, m.upper_s / 60.0)
    assert data.share_h == tuple(s / 3600.0 for s in (85.0, 40.0, 160.0, 715.0))
    # The counts the samples step and the cut rule leave, each a different number: five
    # identical rows at 30 s (4 dropped), four conflicting at 90 s (3 dropped), a lone
    # offline sample just before a polling gap (a zero-length outage), and two gaps.
    sampled = cycles.spells_from_samples(
        _samples(
            [
                (0, "calibrating"),
                *[(30, "idle")] * 5,
                (60, "calibrating"),
                *[(90, "executing"), (90, "idle")] * 2,
                (120, "offline"),
                (500, "idle"),
                (530, "calibrating"),
                (560, "idle"),
                (900, "calibrating"),
                (930, "idle"),
            ]
        ),
        spacing_s=30.0,
    )
    sampled_table = cycles.cycles(sampled)
    sampled_data = cycles.build_cycles_panel_data(
        sampled_table, _checks(sampled_table), cycles.status_shares(sampled)
    )
    assert (
        sampled_data.n_duplicates_identical_dropped,
        sampled_data.n_duplicates_conflicting_dropped,
        sampled_data.n_zero_outages_dropped,
        sampled_data.n_unobserved_spans,
    ) == (4, 3, 1, 2)


@pytest.mark.unit
def test_a_record_with_nothing_to_score_still_gets_its_rows_and_says_why() -> None:
    """Oracle: the ledger's verdict vocabulary. One cycle, no complete one: no segment
    reaches the battery, nothing raises, every asked check reads unanswered, and the reason
    reaches the panel data. A run-set outside the permutation rows is refused."""
    rows = [(10.0, 20.0, "calibrating"), (60.0, 70.0, "calibrating")]
    norm = _spell_norm(rows, length_s=100.0)
    table = cycles.cycles(norm)
    assert table.event.tolist() == [True, False]
    checks = _checks(table)
    assert checks.rows["verdict"].tolist() == [VERDICT_NOT_COMPUTED]
    assert "no segment reached 2 events (1 dropped)" in checks.battery_note
    data = cycles.build_cycles_panel_data(table, checks, cycles.status_shares(norm))
    assert data.battery_note == checks.battery_note
    assert len(data.checks.checks_unanswered) == len(PERMUTATION_KEYS) == 5
    assert all(verdict == VERDICT_ABSENT for _l, _p, verdict, _n in data.verdict_rows)
    with pytest.raises(ValueError, match="permutation rows only"):
        cycles.cycle_checks(
            table,
            _bench_table(),
            run_set=(("c1_lewis_robinson", "asymptotic", ""),),
            alpha=0.05,
            min_events_pass=35,
            tie_cutoff_distinct=5,
            lag_max=5,
            n_permutations=99,
            seed=1,
            min_events_per_segment=2,
        )


@pytest.mark.unit
def test_cycle_checks_are_reproducible_from_their_seed() -> None:
    """Oracle: the permutation rng is a function of the seed. One seed twice gives the same
    rows; two seeds give different permutation p-values on the hand record."""
    table = cycles.cycles(_spell_norm(_HAND_ROWS))

    def rows(seed: int) -> pd.DataFrame:
        return cycles.cycle_checks(
            table,
            _bench_table(),
            run_set=PERMUTATION_KEYS,
            alpha=0.05,
            min_events_pass=35,
            tie_cutoff_distinct=5,
            lag_max=5,
            n_permutations=999,
            seed=seed,
            min_events_per_segment=2,
        ).rows

    pd.testing.assert_frame_equal(rows(11), rows(11))
    permutation = rows(11)["calibration"] == "permutation"
    assert not np.array_equal(
        rows(11)["p_value"][permutation].to_numpy(dtype=float),
        rows(12)["p_value"][permutation].to_numpy(dtype=float),
    )


def _samples(points: list[tuple[float, str]]) -> Norm:
    """A status history of point samples, times in seconds from `_T0`."""
    return Norm(
        {
            "sample_time_unix_s": np.array([_T0 + t for t, _ in points], dtype=float),
            "sample_status": np.array([s for _, s in points]),
            "meta": {"dataset_id": "samples"},
        }
    )


def _spell_rows(spells: Norm) -> list[tuple[float, float, str]]:
    return [
        (float(s - _T0), float(e - _T0), str(status))
        for s, e, status in zip(
            spells["spell_start_unix_s"],
            spells["spell_end_unix_s"],
            spells["spell_status"],
        )
    ]


@pytest.mark.unit
def test_a_spell_runs_from_its_first_sample_to_the_next_runs_first() -> None:
    """Oracle: a hand table. Samples every 30 s: idle, calibrating, calibrating,
    benchmarking, idle, idle. Spells: idle 0-30, calibrating 30-90, benchmarking 90-120,
    idle 120-150 (the last run ends at the record's last sample)."""
    spells = cycles.spells_from_samples(
        _samples(
            [
                (0, "idle"),
                (30, "calibrating"),
                (60, "calibrating"),
                (90, "benchmarking"),
                (120, "idle"),
                (150, "idle"),
            ]
        ),
        spacing_s=30.0,
    )
    assert _spell_rows(spells) == [
        (0.0, 30.0, "idle"),
        (30.0, 90.0, "calibrating"),
        (90.0, 120.0, "benchmarking"),
        (120.0, 150.0, "idle"),
    ]
    meta = spells["meta"]
    assert (meta["record_start_unix_s"] - _T0, meta["record_end_unix_s"] - _T0) == (
        0,
        150,
    )
    assert len(spells["unobserved_start_unix_s"]) == 0


@pytest.mark.unit
def test_a_polling_gap_is_unobserved_and_cuts_the_record_like_an_outage() -> None:
    """Oracle: a hand table. Samples: 0 calibrating, 30 and 60 idle, a 360 s step, 420
    idle, 450 calibrating, 480 executing, 510 calibrating, 540 idle.

    The step exceeds 2 x 30 s, so 60-420 is unobserved and the idle run is two spells,
    30-60 and 420-450. Segments 0-60 and 420-540: the cycle 30-60 runs into the span and is
    censored there; 420-450 has no known start and is dropped; 480-510 (30 s) is complete;
    the calibration 510-540 ends at the last sample, so no censored tail. Shares:
    calibrating 90 s, idle 60 s, executing 30 s, not listed 360 s.
    """
    spells = cycles.spells_from_samples(
        _samples(
            [
                (0, "calibrating"),
                (30, "idle"),
                (60, "idle"),
                (420, "idle"),
                (450, "calibrating"),
                (480, "executing"),
                (510, "calibrating"),
                (540, "idle"),
            ]
        ),
        spacing_s=30.0,
    )
    assert _spell_rows(spells) == [
        (0.0, 30.0, "calibrating"),
        (30.0, 60.0, "idle"),
        (420.0, 450.0, "idle"),
        (450.0, 480.0, "calibrating"),
        (480.0, 510.0, "executing"),
        (510.0, 540.0, "calibrating"),
        (540.0, 540.0, "idle"),
    ]
    assert (spells["unobserved_start_unix_s"] - _T0).tolist() == [60.0]
    assert (spells["unobserved_end_unix_s"] - _T0).tolist() == [420.0]
    table = cycles.cycles(spells)
    assert table.length_s.tolist() == [30.0, 30.0]
    assert table.event.tolist() == [False, True]
    assert (table.n_segments, table.n_unobserved_birth_dropped) == (2, 1)
    assert table.n_unobserved_spans == 1
    shares = cycles.status_shares(spells)
    assert dict(zip(shares.statuses, shares.seconds)) == {
        "calibrating": 90.0,
        "idle": 60.0,
        "executing": 30.0,
        cycles.NOT_LISTED: 360.0,
    }
    # At the threshold: a step of exactly 2 x 30 s is not a gap, one of 61 s is.
    edge = cycles.spells_from_samples(
        _samples(
            [(0, "idle"), (60, "calibrating"), (121, "calibrating"), (151, "idle")]
        ),
        spacing_s=30.0,
    )
    assert _spell_rows(edge) == [
        (0.0, 60.0, "idle"),
        (60.0, 60.0, "calibrating"),
        (121.0, 151.0, "calibrating"),
        (151.0, 151.0, "idle"),
    ]
    assert (edge["unobserved_start_unix_s"] - _T0).tolist() == [60.0]
    assert (edge["unobserved_end_unix_s"] - _T0).tolist() == [121.0]


@pytest.mark.unit
def test_a_lone_calibrating_sample_before_a_gap_ends_the_stretch_the_gap_closes() -> (
    None
):
    """Oracle: a hand table. Samples: 0 calibrating, 30 and 60 idle, 90 calibrating, a
    310 s step, 400 idle, 430 calibrating, 460 and 490 idle.

    The calibration at 90 is a spell of zero length at the span's start, so it lies in the
    stretch 0-90 that the span closes, not in the one after it. Segment 0-90: the cycle
    30-90 (60 s) is complete and the calibration ends at the segment's end, so no censored
    tail. Segment 400-490: 400-430 has no known start and is dropped; 460-490 (30 s) is
    censored at the record end. No cycle crosses the gap.
    """
    spells = cycles.spells_from_samples(
        _samples(
            [
                (0, "calibrating"),
                (30, "idle"),
                (60, "idle"),
                (90, "calibrating"),
                (400, "idle"),
                (430, "calibrating"),
                (460, "idle"),
                (490, "idle"),
            ]
        ),
        spacing_s=30.0,
    )
    table = cycles.cycles(spells)
    assert list(
        zip(
            (table.start_unix_s - _T0).tolist(),
            (table.end_unix_s - _T0).tolist(),
            table.event.tolist(),
        )
    ) == [(30.0, 90.0, True), (460.0, 490.0, False)]
    assert table.segment.tolist() == [0, 1]
    assert (table.n_segments, table.n_unobserved_birth_dropped) == (2, 1)


@pytest.mark.unit
def test_identical_duplicate_samples_are_one_sample_and_are_counted() -> None:
    """Oracle: a hand table. Two identical calibrating rows at 30 s are one sample."""
    spells = cycles.spells_from_samples(
        _samples([(0, "idle"), (30, "calibrating"), (30, "calibrating"), (60, "idle")]),
        spacing_s=30.0,
    )
    assert _spell_rows(spells) == [
        (0.0, 30.0, "idle"),
        (30.0, 60.0, "calibrating"),
        (60.0, 60.0, "idle"),
    ]
    meta = spells["meta"]
    assert meta["n_duplicates_identical_dropped"] == 1
    assert meta["n_duplicates_conflicting_dropped"] == 0
    assert cycles.cycles(spells).n_duplicates_identical_dropped == 1


@pytest.mark.unit
def test_conflicting_duplicate_samples_keep_the_last_in_file_order() -> None:
    """Oracle: a hand table. At 30 s, idle then calibrating: calibrating holds. At 90 s,
    executing, idle, executing: executing holds. Three rows dropped as conflicting."""
    spells = cycles.spells_from_samples(
        _samples(
            [
                (0, "idle"),
                (30, "idle"),
                (30, "calibrating"),
                (60, "idle"),
                (90, "executing"),
                (90, "idle"),
                (90, "executing"),
                (120, "idle"),
            ]
        ),
        spacing_s=30.0,
    )
    assert _spell_rows(spells) == [
        (0.0, 30.0, "idle"),
        (30.0, 60.0, "calibrating"),
        (60.0, 90.0, "idle"),
        (90.0, 120.0, "executing"),
        (120.0, 120.0, "idle"),
    ]
    meta = spells["meta"]
    assert meta["n_duplicates_conflicting_dropped"] == 3
    assert meta["n_duplicates_identical_dropped"] == 0
    assert cycles.cycles(spells).n_duplicates_conflicting_dropped == 3


@pytest.mark.unit
def test_a_zero_length_offline_row_is_dropped_and_counted_not_an_outage() -> None:
    """Oracle: a hand table. Calibrating 10-20 and 80-90 around an offline row at 50-50 in
    a 100 s record: the offline row cuts nothing, so 20-80 (60 s) is one complete cycle;
    0-10 is a residual and 90-100 (10 s) is censored at the record end."""
    rows = [(10, 20, "calibrating"), (50, 50, "offline"), (80, 90, "calibrating")]
    table = cycles.cycles(_spell_norm(rows, length_s=100.0))
    assert table.length_s.tolist() == [60.0, 10.0]
    assert table.event.tolist() == [True, False]
    assert table.n_zero_outages_dropped == 1
    assert (table.n_segments, table.n_unobserved_birth_dropped) == (1, 1)
