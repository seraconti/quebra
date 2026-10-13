"""The single-step min-P panel p-value, against hand enumeration and analytic bounds.

Oracles:
- hand enumeration of a two-member panel at B = 2, row by row;
- a duplicated member: perfectly dependent members cost nothing, so the panel's p equals
  the member's own p exactly, however large B is;
- the checks' own permutation p (`checks._permutation.permutation_p_value`), which every
  member p must reproduce, ties included;
- the bounds that hold for every input: the smallest member p <= the panel's p <= the
  member count times it.
"""

from __future__ import annotations

import ast
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from quebra.analyzers import min_p as min_p_module
from quebra.analyzers.checks._permutation import permutation_p_value
from quebra.analyzers.min_p import min_p, row_p_values


@pytest.mark.unit
def test_two_members_at_b2_match_hand_enumeration() -> None:
    """Oracle: by hand. Rows (observed, perm 1, perm 2): a = 9, 3, 7 and b = 4, 6, 1.

    Row p-values, the share of rows at least as large: a = 1/3, 1, 2/3 and b = 2/3, 1/3, 1.
    Row min-p = 1/3, 1/3, 2/3. Two of the three rows have min-p <= the observed 1/3, so the
    panel's p is 2/3, twice the smallest member p: here the panel pays the full Bonferroni
    cap for its two members.
    """
    stats = {"a": [9.0, 3.0, 7.0], "b": [4.0, 6.0, 1.0]}
    np.testing.assert_array_equal(
        row_p_values(np.column_stack(list(stats.values()))),
        [[1 / 3, 2 / 3], [1.0, 1 / 3], [2 / 3, 1.0]],
    )
    result = min_p(stats)
    assert result.member_p_values == {"a": 1 / 3, "b": 2 / 3}
    assert result.p_value == 2 / 3
    assert result.argmin_labels == ("a",)
    assert result.n_permutations == 2


@pytest.mark.unit
def test_a_duplicated_member_costs_nothing() -> None:
    """Oracle: analytic. A member and an increasing copy of it rank every row identically, so
    every row's min-p is the member's own p, and the panel's p equals the member's p
    exactly. Shifting one member's permutation rows against the other's would break the
    pairing and raise the panel's p toward twice the member's."""
    rng = np.random.default_rng(7)
    null = rng.normal(size=999)
    observed = float(np.quantile(null, 0.97))
    a = np.r_[observed, null]
    result = min_p({"a": a, "a_copy": 2.0 * a + 5.0})
    assert result.member_p_values["a"] == result.member_p_values["a_copy"]
    assert result.p_value == result.member_p_values["a"]
    assert result.argmin_labels == ("a", "a_copy")


@pytest.mark.unit
def test_member_p_values_are_the_checks_p_values_ties_included() -> None:
    """Oracle: `permutation_p_value`, the checks' p, which counts ties against. A tie with
    the observed value is evidence for the null: a = 5, 5, 3 gives 2/3, not 1/3."""
    assert min_p({"a": [5.0, 5.0, 3.0]}).member_p_values["a"] == 2 / 3
    rng = np.random.default_rng(11)
    stats = {f"m{j}": rng.integers(0, 4, size=200).astype(float) for j in range(3)}
    result = min_p(stats)
    for label, values in stats.items():
        assert result.member_p_values[label] == permutation_p_value(
            values[0], values[1:]
        )


@pytest.mark.unit
@pytest.mark.parametrize(
    "stats, match",
    [
        ({}, "at least one member"),
        ({"a": [1.0, 2.0, 3.0], "b": [1.0, 2.0]}, "same number of rows"),
        ({"a": [1.0]}, "at least one permutation"),
        ({"a": [1.0, np.nan, 3.0]}, "non-finite"),
        ({"a": [1.0, np.inf, 3.0]}, "non-finite"),
        ({"a": [1.0, -np.inf, 3.0]}, "non-finite"),
        ({"a": [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]}, "one value per row"),
    ],
    ids=["empty", "unequal", "no-permutation", "nan", "inf", "-inf", "two-dimensional"],
)
def test_a_malformed_panel_raises(stats, match) -> None:
    """Oracle: specification. None of these has a calibrated answer."""
    with pytest.raises(ValueError, match=match):
        min_p(stats)


@pytest.mark.properties
@given(
    data=st.integers(min_value=1, max_value=4).flatmap(
        lambda n_members: st.integers(min_value=2, max_value=30).flatmap(
            lambda n_rows: st.lists(
                st.lists(
                    st.integers(min_value=0, max_value=5),
                    min_size=n_rows,
                    max_size=n_rows,
                ),
                min_size=n_members,
                max_size=n_members,
            )
        )
    )
)
def test_the_panel_p_lies_between_its_member_bounds(data) -> None:
    """Oracle: analytic bounds. The argmin member's own rows already put the panel's p at
    or above the smallest member p, and the union of the members' rows caps it at the
    member count times that p (Bonferroni). It is also a multiple of 1 / (B + 1)."""
    stats = {f"m{j}": [float(x) for x in column] for j, column in enumerate(data)}
    result = min_p(stats)
    smallest = min(result.member_p_values.values())
    n_rows = result.n_permutations + 1
    assert smallest <= result.p_value <= min(1.0, len(stats) * smallest) + 1e-12
    assert abs(result.p_value * n_rows - round(result.p_value * n_rows)) < 1e-9


@pytest.mark.policy
def test_min_p_is_a_leaf() -> None:
    """Oracle: the module docstring's contract. Its source imports nothing from quebra, at
    module level or inside a function, and importing it in a fresh process loads no other
    quebra module."""
    tree = ast.parse(Path(min_p_module.__file__).read_text())
    imported = {
        name
        for node in ast.walk(tree)
        for name in (
            [alias.name for alias in node.names]
            if isinstance(node, ast.Import)
            else ["." * node.level + (node.module or "")]
            if isinstance(node, ast.ImportFrom)
            else []
        )
    }
    assert imported <= {"__future__", "collections.abc", "dataclasses", "numpy"}
    probe = (
        "import sys, quebra.analyzers.min_p; "
        "print(' '.join(sorted(m for m in sys.modules if m.split('.')[0] == 'quebra')))"
    )
    loaded = subprocess.run(
        [sys.executable, "-c", probe], capture_output=True, text=True, check=True
    ).stdout.split()
    assert loaded == ["quebra", "quebra.analyzers", "quebra.analyzers.min_p"]


@pytest.mark.properties
@given(
    data=st.integers(min_value=1, max_value=4).flatmap(
        lambda n_members: st.integers(min_value=2, max_value=25).flatmap(
            lambda n_rows: st.lists(
                st.lists(
                    st.integers(min_value=0, max_value=4),
                    min_size=n_rows,
                    max_size=n_rows,
                ),
                min_size=n_members,
                max_size=n_members,
            )
        )
    )
)
def test_no_more_than_k_rows_reach_k_over_n(data) -> None:
    """Oracle: exchangeability, made deterministic. Under the null every row is equally
    likely to be the observed one, so the size claim holds iff, taking each row in turn as
    the observed one, at most k rows get a panel p <= k / (B + 1), for every k."""
    table = np.array(data, dtype=float).T  # rows x members
    n_rows = table.shape[0]
    panel = []
    for r in range(n_rows):
        order = [r, *(i for i in range(n_rows) if i != r)]
        stats = {f"m{j}": table[order, j] for j in range(table.shape[1])}
        panel.append(min_p(stats).p_value)
    for k in range(1, n_rows + 1):
        assert sum(p <= k / n_rows + 1e-12 for p in panel) <= k
