"""The ledger's group rows: which members answer, the draws they keep, and what each
group says.

A min-P over several checks on one permutation set keeps its size only when the set of
members that answer is a property of the record, never of the shuffle; a member's own
permutation draws are what the min-P ranks. Synthetic segments and windows only.
"""

from __future__ import annotations

import ast
import itertools
from collections import Counter
from dataclasses import replace
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

import quebra.analyzers.check_ledger as ledger  # noqa: E402
import quebra.analyzers.checks.c6_exchangeability as c6_module  # noqa: E402
from quebra.analyzers import windows  # noqa: E402
from quebra.analyzers.calibration_summary import bench_acceptance_at_n  # noqa: E402
from quebra.analyzers.independence_survey import (  # noqa: E402
    build_independence_survey,
    survey_summary,
)
from quebra.analyzers.min_p import min_p  # noqa: E402
from quebra.ledger_recipe import ledger_bag  # noqa: E402
from quebra.panels.check_ledger import CheckLedgerPanel  # noqa: E402

import quebra.analyzers.checks.c1_lewis_robinson as c1  # noqa: E402
import quebra.analyzers.checks.c2_anderson_darling as c2  # noqa: E402
import quebra.analyzers.checks.cvm_cramer_von_mises as cvm  # noqa: E402
from quebra.analyzers.checks._permutation import (  # noqa: E402
    PermutationSet,
    block_permutations,
    permutation_p_value,
    two_sided_p_value,
)
from quebra.analyzers.checks.battery import REFUSED, UNDEFINED, run_battery  # noqa: E402
from quebra.analyzers.checks._rank_serial import (  # noqa: E402
    MAX_LAG_CAP,
    autocorrelations,
    global_ranks,
    lag_layout,
    studentize,
)
from quebra.analyzers.checks.result import (  # noqa: E402
    CALIB_PERMUTATION,
    TAU_MARGIN,
    CheckRefused,
    Segment,
    last_event_time,
    validate_segment,
)

C1, C2, CVM = c1.CHECK_NAME, c2.CHECK_NAME, cvm.CHECK_NAME
C5, C6 = "c5_rank_autocorr", "c6_exchangeability"


def _seg(gaps: list[float], residual: float = 1.0) -> Segment:
    x = np.asarray(gaps, dtype=float)
    return Segment(x=x, tau=float(np.sum(x)) + residual)


# (site, segments, the checks it must leave refused, the checks that must still answer)
REFUSAL_SITES = [
    ("a zero gap", [_seg([1.0, 0.0, 2.0, 1.5, 0.8])], {C1, C2, CVM, C5, C6}, set()),
    ("no usable lag", [_seg([1.0, 2.0])], {C5, C6}, {C1, C2, CVM}),
    (
        "no spread at a lag",
        [_seg([40.0, 55.0]), _seg([47.0, 62.0]), _seg([51.0, 70.0])],
        {C5, C6},
        {C1, C2, CVM},
    ),
    (
        "every gap equal",
        [_seg([2.0, 2.0, 2.0, 2.0, 2.0])],
        {C1, C2, CVM, C5, C6},
        set(),
    ),
]


@pytest.mark.unit
@pytest.mark.parametrize(
    ("site", "segments", "refused", "answering"),
    REFUSAL_SITES,
    ids=[case[0] for case in REFUSAL_SITES],
)
def test_each_refusal_of_the_record_is_marked_as_one(
    site: str, segments: list[Segment], refused: set[str], answering: set[str]
) -> None:
    """Oracle: specification. Each admissible refusal reads only segment sizes, gap
    multisets and tau, or the augmented ensemble as a whole, and its undefined row says so
    (`extra[REFUSED]`), so a group may leave that member out and keep its size."""
    rows = {
        r.check: r
        for r in run_battery(
            segments,
            rng=np.random.default_rng(1),
            n_perm=99,
            undefined_as_rows=True,
        )
    }
    for check in refused:
        assert UNDEFINED in rows[check].extra, (site, check)
        assert rows[check].extra[REFUSED] is True, (site, check, rows[check].notes)
    for check in answering:
        assert rows[check].p_value is not None, (site, check, rows[check].notes)


@pytest.mark.unit
def test_a_raise_that_is_not_a_refusal_leaves_its_row_unmarked() -> None:
    """Oracle: specification. A raise other than `CheckRefused` leaves its row undefined
    and unmarked, so no group may leave that member out: here a `tau` below the last
    event, a producer's bug."""
    segment = Segment(x=np.array([1.0, 2.5, 0.7, 1.9, 3.1, 0.4]), tau=5.0)
    rows = run_battery(
        [segment], rng=np.random.default_rng(1), n_perm=99, undefined_as_rows=True
    )
    assert {row.check for row in rows} == {C1, C2, CVM, C5, C6}
    for row in rows:
        assert UNDEFINED in row.extra, row.check
        assert row.extra[REFUSED] is False, (row.check, row.notes)


# Every `raise CheckRefused` in the package, by module and function, each reviewed as
# reading only the record: a mark on a raise that reads the draws would let a group drop
# that member on some draws and keep it on others.
REVIEWED_REFUSALS = Counter(
    {
        ("analyzers/checks/result.py", "validate_segment"): 2,
        ("analyzers/checks/_multiprocess.py", "gamma_hat"): 1,
        ("analyzers/checks/_rank_serial.py", "lag_layout"): 1,
        ("analyzers/checks/_rank_serial.py", "global_ranks"): 1,
        ("analyzers/checks/_rank_serial.py", "autocorrelations"): 1,
        ("analyzers/checks/_rank_serial.py", "studentize"): 1,
    }
)


@pytest.mark.unit
def test_only_the_reviewed_raises_mark_a_refusal() -> None:
    """Oracle: the reviewed list above. A new mark, or a subclass of `CheckRefused`, must be
    reviewed as reading only the record before it joins the list."""
    package = Path(c1.__file__).resolve().parents[2]
    found: Counter[tuple[str, str]] = Counter()
    for path in sorted(package.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        parents = {
            child: node
            for node in ast.walk(tree)
            for child in ast.iter_child_nodes(node)
        }
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef):
                bases = {getattr(b, "id", getattr(b, "attr", "")) for b in node.bases}
                if "CheckRefused" in bases:
                    found[(path.relative_to(package).as_posix(), node.name)] += 1
            if not isinstance(node, ast.Raise) or node.exc is None:
                continue
            exc = node.exc.func if isinstance(node.exc, ast.Call) else node.exc
            if getattr(exc, "id", getattr(exc, "attr", "")) != "CheckRefused":
                continue
            owner: ast.AST = node
            while owner in parents and not isinstance(owner, ast.FunctionDef):
                owner = parents[owner]
            name = getattr(owner, "name", "<module>")
            found[(path.relative_to(package).as_posix(), name)] += 1
    assert found == REVIEWED_REFUSALS


@pytest.mark.unit
def test_the_strict_tau_refusal_reads_a_tau_on_the_margin_the_same_in_every_order() -> (
    None
):
    """Oracle: specification. `validate_segment`'s strict-tau refusal, which no check asks
    for since C1, C2 and CvM read n - 1 events, reads the gap multiset and `tau`, so a `tau`
    placed exactly on the margin is refused in every order of its gaps or in none."""
    rng = np.random.default_rng(11)
    for _ in range(20):
        gaps = rng.uniform(0.1, 3.0, size=40)
        tau = last_event_time(np.sort(gaps)) * (1.0 + TAU_MARGIN)
        outcomes = set()
        for _ in range(30):
            segment = Segment(x=rng.permutation(gaps), tau=tau)
            try:
                validate_segment(segment, require_strict_tau=True)
                outcomes.add("answers")
            except CheckRefused:
                outcomes.add("refused")
        assert len(outcomes) == 1, outcomes


@pytest.mark.unit
def test_each_permutation_row_keeps_the_draws_its_p_reads() -> None:
    """Oracle: each check's statistic recomputed row by row on the same permutation set
    (C1, C2 and CvM by their batch forms, C5 and C6 from the rank autocorrelations), and
    their own p rules. The draws are kept signed: C1's two-sided p takes |.| itself.
    Row order matters: a min-P pairs row j of every member."""
    segments = [
        _seg([1.0, 2.5, 0.7, 1.9, 3.1, 0.4, 2.2, 1.3]),
        _seg([0.9, 1.7, 2.8, 0.5, 1.1, 3.6, 0.8]),
    ]
    perm = block_permutations(
        [len(s.x) for s in segments], 199, np.random.default_rng(3)
    )
    rows = {
        r.check: r
        for r in run_battery(segments, perm=perm, include_asymptotic=True)
        if r.calibration == CALIB_PERMUTATION
    }
    for module, check in ((c1, C1), (c2, C2), (cvm, CVM)):
        np.testing.assert_array_equal(
            rows[check].null_statistics, module.statistic_batch(segments, perm)
        )
    layout = lag_layout([len(s.x) for s in segments], MAX_LAG_CAP)
    ranks = global_ranks(np.concatenate([s.x for s in segments]))
    _, z_null = studentize(
        autocorrelations(ranks, layout)[0], autocorrelations(perm.apply(ranks), layout)
    )
    np.testing.assert_array_equal(rows[C5].null_statistics, np.abs(z_null).max(axis=1))
    np.testing.assert_array_equal(rows[C6].null_statistics, (z_null**2).sum(axis=1))
    assert np.any(rows[C1].null_statistics < 0.0)
    assert rows[C1].p_value == two_sided_p_value(
        rows[C1].statistic, rows[C1].null_statistics
    )
    for check in (C2, CVM, C5, C6):
        null = rows[check].null_statistics
        assert null is not None and null.shape == (199,), check
        assert rows[check].p_value == permutation_p_value(rows[check].statistic, null)
    asymptotic = [
        r
        for r in run_battery(segments, perm=perm, include_asymptotic=True)
        if r.calibration != CALIB_PERMUTATION
    ]
    assert asymptotic and all(r.null_statistics is None for r in asymptotic)


# ----------------------------------------------------------------------------- group rows

ACCEPTANCE = bench_acceptance_at_n(
    pd.read_csv(
        Path(__file__).resolve().parents[1] / "jobs/bench/results/size_table.csv"
    ),
    alpha=0.05,
)
GRID = sorted({int(n) for n in ACCEPTANCE["n_target"]})
# Gaps that lengthen down the record: C1's statistic is negative there.
LENGTHENING = [0.6, 0.5, 0.9, 0.7, 1.1, 0.8, 1.4, 1.0, 1.2, 1.7, 1.3, 1.9, 1.5, 2.2]
# The same kind of gaps in no particular order.
SCATTERED = [0.7, 1.9, 0.4, 2.6, 1.3, 0.9, 2.2, 0.5, 1.6, 1.1, 2.9, 0.8, 1.4, 2.0]


def _context(**overrides: object) -> dict[str, object]:
    context: dict[str, object] = dict(
        dataset_id="t",
        label="3 µs",
        side="in_spec",
        clock="side",
        alpha=0.05,
        n_events=40,
        n_read=40,
        n_used=40,
        min_events_pass=0,
        n_distinct=99,
        tie_fraction=0.0,
        tie_cutoff=0,
        include_tau_checks=True,
        n_segments=1,
        n_censored_dropped=0,
    )
    context.update(overrides)
    return context


def _score(segments: list[Segment], **overrides: object) -> pd.DataFrame:
    kwargs: dict[str, object] = dict(
        n_dropped=0,
        n_degenerate=0,
        n_used=int(sum(len(s.x) for s in segments)),
        label="3 µs",
        side="in_spec",
        clock="side",
        dataset_id="t",
        acceptance=ACCEPTANCE,
        grid=GRID,
        alpha=0.05,
        min_events_pass=0,
        tie_cutoff_distinct=0,
        lag_max=5,
        n_permutations=199,
        seed=11,
        min_events_per_segment=2,
        include_tau_checks=True,
        include_c2_asymptotic=False,
        include_c3=False,
        c3_n_null_sim=0,
        groups=True,
    )
    kwargs.update(overrides)
    return pd.DataFrame(ledger.score_segments(segments, **kwargs)).set_index("check_id")


@pytest.mark.unit
def test_each_member_p_in_a_group_is_its_own_rows_p() -> None:
    """Oracle: each check's own permutation p. A member enters the min-P as its p reads it,
    C1 two-sided through |.|, so the group's row-0 p of each member is that check's p, as
    floats; C1's statistic is negative here, where a one-sided C1 would differ."""
    results = run_battery(
        [_seg(LENGTHENING, residual=0.6)],
        rng=np.random.default_rng(5),
        n_perm=199,
        undefined_as_rows=True,
    )
    by_check = {r.check: r for r in results}
    assert by_check[C1].statistic < 0.0
    for group, members in ledger.GROUPS.items():
        found = ledger.group_members(results, group, include_tau_checks=True)
        assert found.blocking == "" and found.absent == []
        member_p = min_p(found.columns).member_p_values
        for check, _calibration, _variant in members:
            assert member_p[ledger.MEMBER_SHORT[check]] == by_check[check].p_value, (
                check
            )


@pytest.mark.unit
def test_a_member_entered_twice_leaves_the_group_p_unchanged(monkeypatch) -> None:
    """Oracle: min-P is invariant to an increasing copy of a member (test_min_p.py on
    `min_p` itself); here through `group_rows`, so a builder that scaled the smallest
    member p by the member count would move."""
    results = run_battery(
        [_seg(LENGTHENING, residual=0.6)],
        rng=np.random.default_rng(5),
        n_perm=199,
        undefined_as_rows=True,
    )
    rows = ledger.group_rows(results, **_context())
    before = {r["check_id"]: r["p_value"] for r in rows}
    smallest = next(r["statistic"] for r in rows if r["check_id"] == ledger.GROUP_TREND)
    copy = replace(next(r for r in results if r.check == C1), check="c1_copy")
    monkeypatch.setitem(
        ledger.GROUPS,
        ledger.GROUP_TREND,
        (*ledger.GROUPS[ledger.GROUP_TREND], ("c1_copy", CALIB_PERMUTATION, "")),
    )
    monkeypatch.setitem(ledger.MEMBER_SHORT, "c1_copy", "C1 copy")
    monkeypatch.setattr(
        ledger, "TWO_SIDED_MEMBERS", ledger.TWO_SIDED_MEMBERS | {"c1_copy"}
    )
    after = {
        r["check_id"]: r["p_value"]
        for r in ledger.group_rows([*results, copy], **_context())
    }
    assert 3 * smallest < 1.0, "premise: the count would show"
    assert after == before


@pytest.mark.statistical
@pytest.mark.parametrize(
    "residual", [0.8, 0.0], ids=["time truncated", "ends at an event"]
)
def test_every_shuffle_of_a_small_record_keeps_each_group_within_its_level(
    residual: float,
) -> None:
    """Oracle: exhaustive enumeration. The 119 other orderings of 5 distinct gaps are the
    permutation set, so B + 1 is the whole orbit; over every ordering taken as the record,
    P(group p <= v) <= v at every attained v, for both groups, on a segment truncated in
    time and on one that ends at an event, whose first n - 1 events C1, C2 and CvM read."""
    x = np.array([0.7, 1.9, 0.4, 2.6, 1.3])
    tau = float(np.cumsum(x)[-1]) + residual
    assert Segment(x=x, tau=tau).ends_at_event is (residual == 0.0)
    orders = list(itertools.permutations(range(5)))
    perm = PermutationSet(sizes=(5,), indices=np.array(orders[1:]))
    p_values: dict[str, list[float]] = {group: [] for group in ledger.GROUPS}
    for order in orders:
        results = run_battery(
            [Segment(x=x[list(order)], tau=tau)],
            perm=perm,
            max_lag=3,
            undefined_as_rows=True,
        )
        context = _context(n_events=5, n_read=4 if residual == 0.0 else 5)
        for row in ledger.group_rows(results, **context):
            p_values[str(row["check_id"])].append(row["p_value"])
    for group, values in p_values.items():
        assert len(values) == 120 and None not in values, group
        ps = np.array(values, dtype=float)
        for v in np.unique(ps):
            assert np.mean(ps <= v) <= v + 1e-12, (group, v)


@pytest.mark.unit
def test_a_segment_that_ends_at_an_event_counts_the_trend_group_by_the_events_read() -> (
    None
):
    """Oracle: specification. At tau = T_N, C1, C2 and CvM read the first n - 1 events and
    all answer, so the trend group is their min-P; it counts and judges by n - 1 events, as
    their rows do, and says so, while the dependence group counts all n."""
    x = np.asarray(LENGTHENING, dtype=float)
    rows = _score([Segment(x=x, tau=float(np.cumsum(x)[-1]))])
    trend = rows.loc[ledger.GROUP_TREND]
    for check in (C1, C2, CVM):
        assert pd.notna(rows.loc[check, "p_value"]), check
        assert rows.loc[check, "n_events"] == len(x) - 1, check
    assert trend["notes"].startswith("members C1, C2, CvM; member p ")
    assert (
        "1 of 1 segment(s) end at an event: first n - 1 events read" in trend["notes"]
    )
    assert trend["n_events"] == len(x) - 1
    assert rows.loc[ledger.GROUP_DEPENDENCE, "n_events"] == len(x)

    s = np.asarray(SCATTERED, dtype=float)
    floor = _score([Segment(x=s, tau=float(np.cumsum(s)[-1]))], min_events_pass=len(s))
    assert (floor.loc[list(ledger.GROUPS), "p_value"] > 0.05).all()
    assert floor.loc[ledger.GROUP_TREND, "verdict"] == ledger.VERDICT_INCONCLUSIVE
    assert floor.loc[ledger.GROUP_DEPENDENCE, "verdict"] == ledger.VERDICT_PASS


@pytest.mark.unit
def test_a_member_failing_without_a_refusal_leaves_its_group_not_computed(
    monkeypatch,
) -> None:
    """Oracle: specification. A member undefined for a cause that could depend on the
    shuffle would make the member set do so: its group reads not computed with the reason,
    and the other group still answers."""

    def unforeseen(*_args: object, **_kwargs: object) -> None:
        raise ValueError("unforeseen")

    monkeypatch.setattr(c6_module, "run", unforeseen)
    rows = _score([_seg(LENGTHENING, residual=0.6)])
    dependence = rows.loc[ledger.GROUP_DEPENDENCE]
    assert dependence["verdict"] == ledger.VERDICT_NOT_COMPUTED
    assert pd.isna(dependence["p_value"])
    assert dependence["notes"] == (
        "C6 has no answer, and not by a refusal of the record: unforeseen"
    )
    assert pd.notna(rows.loc[ledger.GROUP_TREND, "p_value"])


@pytest.mark.unit
@pytest.mark.parametrize("refusal", [False, True], ids=["plain", "refusal"])
def test_c2_undefined_leaves_the_trend_group_to_the_cause(
    monkeypatch, refusal: bool
) -> None:
    """Oracle: specification. C2 undefined, C1 and CvM still answer. A refusal of the record
    leaves the trend group the min-P of C1 and CvM; any other cause leaves it not computed,
    naming C2. Either way every row counts the record's censored windows."""

    def undefined(*_args: object, **_kwargs: object) -> None:
        raise (CheckRefused if refusal else ValueError)("undefined here")

    monkeypatch.setattr(c2, "run", undefined)
    x = np.asarray(LENGTHENING, dtype=float)
    rows = _score([Segment(x=x, tau=float(x.sum()) + 0.6, n_censored_dropped=1)])
    trend = rows.loc[ledger.GROUP_TREND]
    assert pd.notna(rows.loc[C1, "p_value"]) and pd.isna(rows.loc[C2, "p_value"])
    if refusal:
        smallest = min(rows.loc[C1, "p_value"], rows.loc[CVM, "p_value"])
        assert trend["notes"].startswith(
            "members C1, CvM; no answer: C2: undefined here; "
        )
        assert trend["statistic"] == smallest <= trend["p_value"]
    else:
        assert trend["verdict"] == ledger.VERDICT_NOT_COMPUTED
        assert trend["notes"] == (
            "C2 has no answer, and not by a refusal of the record: undefined here"
        )
    assert pd.notna(rows.loc[ledger.GROUP_DEPENDENCE, "p_value"])
    counted = [C1, C2, ledger.GROUP_TREND, ledger.GROUP_DEPENDENCE]
    assert list(rows.loc[counted, "n_censored_dropped"]) == [1, 1, 1, 1]
    assert list(rows.loc[counted, "n_segments"]) == [1, 1, 1, 1]


@pytest.mark.unit
def test_a_run_set_without_c1_and_c2_leaves_the_trend_group_cvm_alone() -> None:
    """Oracle: specification. Members a run-set did not ask for read `not asked` and are
    left out; the trend group is then CvM's min-P, which is CvM's own p."""
    rows = _score([_seg(LENGTHENING, residual=0.6)], include_tau_checks=False)
    trend = rows.loc[ledger.GROUP_TREND]
    assert C1 not in rows.index and C2 not in rows.index
    assert trend["p_value"] == rows.loc[CVM, "p_value"]
    assert trend["notes"].startswith(
        "members CvM; no answer: C1 not asked; C2 not asked"
    )


@pytest.mark.unit
def test_a_ledger_reads_the_tie_rule_and_the_event_floor_into_its_group_rows() -> None:
    """Oracle: specification. Through `score_segments`: below the distinct-duration cutoff
    the dependence group reads ties and the trend group does not; below the event floor a
    group that does not reject reads inconclusive, and passes without the floor."""
    segments = [_seg(SCATTERED, residual=0.6)]
    ties = _score(segments, tie_cutoff_distinct=10**6)
    assert ties.loc[ledger.GROUP_DEPENDENCE, "verdict"] == ledger.VERDICT_TIES
    assert ties.loc[ledger.GROUP_TREND, "verdict"] != ledger.VERDICT_TIES
    groups = list(ledger.GROUPS)
    floor = _score(segments, min_events_pass=10**6).loc[groups]
    open_ = _score(segments).loc[groups]
    quiet = floor["p_value"] > 0.05
    assert quiet.any()
    assert set(floor.loc[quiet, "verdict"]) == {ledger.VERDICT_INCONCLUSIVE}
    assert set(open_.loc[quiet, "verdict"]) == {ledger.VERDICT_PASS}


@pytest.mark.unit
@pytest.mark.parametrize(
    ("p_value", "tie_sensitive", "n_distinct", "n_events", "floor", "verdict"),
    [
        (None, False, 99, 40, 35, ledger.VERDICT_NOT_COMPUTED),
        (None, True, 3, 40, 35, ledger.VERDICT_NOT_COMPUTED),
        (0.01, True, 3, 40, 35, ledger.VERDICT_TIES),
        (0.01, False, 3, 10, 35, ledger.VERDICT_FAIL),
        (0.3, False, 99, 10, 35, ledger.VERDICT_INCONCLUSIVE),
        (0.3, False, 99, 20, 10, ledger.VERDICT_PASS),
        (0.3, True, 99, 40, 35, ledger.VERDICT_PASS),
    ],
)
def test_a_group_verdict_keeps_the_tie_rule_and_the_event_floor_and_no_bench(
    p_value: float | None,
    tie_sensitive: bool,
    n_distinct: int,
    n_events: int,
    floor: int,
    verdict: str,
) -> None:
    """Oracle: specification. `_verdict`'s precedence without the bench: not computed, the
    tie rule for a tie-sensitive group, a rejection even below the floor, then the floor
    read from the knob."""
    got, _reason = ledger._group_verdict(
        p_value,
        "no member answered",
        alpha=0.05,
        n_events=n_events,
        min_events_pass=floor,
        tie_sensitive=tie_sensitive,
        n_distinct=n_distinct,
        tie_cutoff=5,
    )
    assert got == verdict


def _ledger(seed: int) -> ledger.CheckLedger:
    """A whole synthetic ledger: both bags at two rungs, both clocks, group rows included."""
    rng = np.random.default_rng(seed)
    values: list[float] = []
    for _ in range(70):
        values += [5.0] * int(rng.integers(2, 12))
        values += [1.0] * int(rng.integers(2, 8))
    series = np.asarray(values, dtype=float)
    ladder = [("3 µs", 3.0, True), ("4 µs", 4.0, True)]
    carve = windows.run(
        windows.WindowsInputs(
            t_rel_s=np.arange(len(series), dtype=float) * 60.0,
            values=series,
            thresholds=ladder,
            dataset_id=f"synthetic {seed}",
        )
    )
    return ledger.run(
        ledger.make_inputs_from_windows(
            carve,
            pd.read_csv(
                Path(__file__).resolve().parents[1]
                / "jobs/bench/results/size_table.csv"
            ),
            thresholds=ladder,
            alpha=0.05,
            min_events_pass=35,
            tie_cutoff_distinct=5,
            lag_max=5,
            n_permutations=99,
            seed=seed,
        )
    )


@pytest.mark.integration
def test_the_ledger_panel_draws_both_group_columns_unhatched_with_their_note() -> None:
    """Oracle: specification. Beside the checks, one column per group, labelled "trend
    group" and "dependence group", never hatched as miscalibrated, a legend entry marking
    them "no bench cell (permutation)", and the calibration sentence on the figure; on
    each bag's figure."""
    record = _ledger(3)
    for side in ("in_spec", "out_of_spec"):
        fig = CheckLedgerPanel(name="unit").build_matplotlib(
            ledger_bag(record, side=side)
        )
        try:
            table = fig.axes[0]
            columns = [t.get_text() for t in table.get_xticklabels()]
            assert "trend group" in columns and "dependence group" in columns
            group_x = {columns.index("trend group"), columns.index("dependence group")}
            group_cells = [
                patch
                for patch in table.patches
                if int(patch.get_x()) in group_x and patch.get_facecolor()[3] > 0
            ]
            assert group_cells and all(not patch.get_hatch() for patch in group_cells)
            legend = [t.get_text() for t in table.get_legend().get_texts()]
            assert (
                "trend group, dependence group: no bench cell (permutation)" in legend
            )
            assert ledger.GROUP_CALIBRATION_NOTE in [t.get_text() for t in fig.texts]
        finally:
            plt.close(fig)
    checks_only = replace(
        record, rows=record.rows[record.rows["calibration"] != ledger.CALIB_MIN_P]
    )
    fig = CheckLedgerPanel(name="unit").build_matplotlib(
        ledger_bag(checks_only, side="in_spec")
    )
    try:
        legend = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
        assert not any("group" in text for text in legend)
        assert ledger.GROUP_CALIBRATION_NOTE not in [t.get_text() for t in fig.texts]
    finally:
        plt.close(fig)


@pytest.mark.integration
def test_the_survey_carries_both_group_grids_and_their_table_rows() -> None:
    """Oracle: specification. One grid per group, bag and clock beside the checks' grids,
    and a summary-table row for each."""
    data = build_independence_survey(_ledger(3), _ledger(4), dataset_labels=("a", "b"))
    for key in ledger.GROUP_KEYS:
        found = {(g.side, g.clock) for g in data.grids if g.key == key}
        assert found == {(side, clock) for side in data.sides for clock in data.clocks}
    instruments = set(survey_summary(data)["instrument"])
    assert {
        "Trend group (min-P over those of C1, C2, CvM that answer)",
        "Dependence group (min-P over those of C5, C6 that answer)",
    } <= instruments
