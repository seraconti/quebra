"""The battery's rows as members of a group: which members answer, and the draws they keep.

A min-P over several checks on one permutation set keeps its size only when the set of
members that answer is a property of the record, never of the shuffle; a member's own
permutation draws are what the min-P ranks. Synthetic segments only.
"""

from __future__ import annotations

import ast
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

import quebra.analyzers.checks.c1_lewis_robinson as c1
import quebra.analyzers.checks.c2_anderson_darling as c2
import quebra.analyzers.checks.cvm_cramer_von_mises as cvm
from quebra.analyzers.checks._permutation import (
    block_permutations,
    permutation_p_value,
    two_sided_p_value,
)
from quebra.analyzers.checks.battery import REFUSED, UNDEFINED, run_battery
from quebra.analyzers.checks._rank_serial import (
    MAX_LAG_CAP,
    autocorrelations,
    global_ranks,
    lag_layout,
    studentize,
)
from quebra.analyzers.checks.result import (
    CALIB_PERMUTATION,
    TAU_MARGIN,
    CheckRefused,
    Segment,
    last_event_time,
    validate_segment,
)

pytestmark = pytest.mark.unit

C1, C2, CVM = c1.CHECK_NAME, c2.CHECK_NAME, cvm.CHECK_NAME
C5, C6 = "c5_rank_autocorr", "c6_exchangeability"


def _seg(gaps: list[float], residual: float = 1.0) -> Segment:
    x = np.asarray(gaps, dtype=float)
    return Segment(x=x, tau=float(np.sum(x)) + residual)


# (site, segments, the checks it must leave refused, the checks that must still answer)
REFUSAL_SITES = [
    (
        "C2's strict tau",
        [_seg([1.0, 2.5, 0.7, 1.9, 3.1, 0.4], residual=0.0)],
        {C2},
        {C1, CVM, C5, C6},
    ),
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


def test_c2_refuses_a_tau_on_the_margin_in_every_order_or_in_none() -> None:
    """Oracle: specification. C2's strict-tau refusal reads the gap multiset and `tau`, so
    a `tau` placed exactly on the margin is refused in every order of its gaps or in none."""
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
