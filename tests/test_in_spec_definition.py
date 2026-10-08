"""One definition of in spec: `analyzers/windows.margin`, with in spec meaning `margin >= 0`.

Three guards, each with a different reach:

- a static scan of `src/quebra` and `jobs/` for a comparison or subtraction that involves a
  threshold value outside `analyzers/windows.py`, which finds code no test runs. It reads
  operators, augmented subtraction, numpy and pandas comparison calls, string subscript
  keys and the value slot of an unpacked ladder entry, inside the loop over a ladder whose
  name ends in `thresholds` or `ladder`. It does not follow a threshold through `.get`,
  `getattr`, an assignment to a new name or a function argument, nor a ladder value used
  after its loop or a ladder under another name;
- a behavioural flip: negate `windows.margin` and every quantity that should follow from
  it must flip, which finds a private classification however it is spelled;
- an agreement case with reads exactly at the threshold, in both directions, against the
  carve, the in-spec fraction, time to first crossing and cumulative time out of spec.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import numpy as np
import pytest

from quebra.analyzers import event_table, kaplan_meier, windows
from quebra.analyzers.within_calibration_compute import (
    build_within_calibration_panel_data,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFINITION = REPO_ROOT / "src" / "quebra" / "analyzers" / "windows.py"
SCAN_ROOTS = ("src/quebra", "jobs")

# An identifier that names a threshold VALUE: thr, thr_val, threshold_value, inf_thr,
# chi_threshold, gap_threshold_s. Labels and ladders (threshold_label, thresholds) do not.
THRESHOLD_VALUE = re.compile(r"(?i)^(\w+_)?(thr|threshold)(_(value|val|s|us|hz|h))?$")
DEFINITION_CALLS = {"margin", "in_spec_mask"}
# Calls that compare or subtract their operands: numpy ufuncs and pandas methods.
COMPARISON_CALLS = {"greater", "greater_equal", "less", "less_equal", "subtract"}
COMPARISON_CALLS |= {"gt", "ge", "lt", "le", "sub"}
# A ladder entry is (label, threshold value, big_values_good), so a loop over a ladder
# names the value in its middle slot, whatever it calls it.
LADDER = re.compile(r"(?i)(thresholds|ladder)$")
_COMPREHENSIONS = ast.ListComp | ast.SetComp | ast.GeneratorExp | ast.DictComp

# Comparisons with a threshold that classify no read of a metric, keyed by (file,
# scope, `ast.unparse` of the comparison): a moved or reformatted line keeps its entry,
# and a different comparison in the same scope does not inherit it. Scope names classes too.
_DIVERGENCE = "ordered['p_signed_diff_mean'].abs() > threshold"
NOT_A_METRIC_THRESHOLD = {
    # which rungs of the nines ladder the data spans; it picks thresholds
    (
        "src/quebra/analyzers/fidelity.py",
        "panel_thresholds",
        "inf_min < inf_thr < inf_max",
    ),
    # divergence thresholds on p-value differences
    (
        "src/quebra/analyzers/instrument_validation.py",
        "_divergence_levels",
        _DIVERGENCE,
    ),
    (
        "src/quebra/analyzers/instrument_validation.py",
        "_divergence_points",
        _DIVERGENCE,
    ),
    (
        "jobs/bench/xi_ties.py",
        "main",
        "sub['p_signed_diff_mean'].abs() > DIVERGENCE_THRESHOLD",
    ),
    # the chi-squared fit-quality filter
    (
        "src/quebra/transforms/filter.py",
        "run",
        "np.asarray(current['chi_squared'], dtype=float) < chi_threshold",
    ),
}


def _names(node: ast.AST) -> list[str]:
    """Identifiers under `node`, not descending into a call to the definition itself."""
    found: list[str] = []
    stack = [node]
    while stack:
        current = stack.pop()
        if isinstance(current, ast.Call):
            func = current.func
            name = (
                func.attr
                if isinstance(func, ast.Attribute)
                else getattr(func, "id", "")
            )
            if name in DEFINITION_CALLS:
                continue
        if isinstance(current, ast.Name):
            found.append(current.id)
        elif isinstance(current, ast.Attribute):
            found.append(current.attr)
        elif isinstance(current, ast.Subscript):
            found.extend(
                key.value
                for key in ast.walk(current.slice)
                if isinstance(key, ast.Constant) and isinstance(key.value, str)
            )
        stack.extend(ast.iter_child_nodes(current))
    return found


def _operands(node: ast.AST) -> list[ast.AST]:
    """The operands of a comparison or subtraction, or nothing for any other node."""
    if isinstance(node, ast.Compare):
        return [node.left, *node.comparators]
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Sub):
        return [node.left, node.right]
    if isinstance(node, ast.AugAssign) and isinstance(node.op, ast.Sub):
        return [node.target, node.value]
    if isinstance(node, ast.Call):
        func = node.func
        name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", "")
        if name in COMPARISON_CALLS:
            receiver = [func.value] if isinstance(func, ast.Attribute) else []
            return [*receiver, *node.args, *(k.value for k in node.keywords)]
    return []


def _ladder_values(node: ast.AST) -> set[str]:
    """Names a loop or comprehension binds to the value slot of a ladder entry."""
    if isinstance(node, ast.For | ast.AsyncFor):
        loops: list[ast.For | ast.AsyncFor | ast.comprehension] = [node]
    elif isinstance(node, _COMPREHENSIONS):
        loops = list(node.generators)
    else:
        return set()
    found: set[str] = set()
    for loop in loops:
        if any(LADDER.search(n) for n in _names(loop.iter)):
            for entry in ast.walk(loop.target):
                if isinstance(entry, ast.Tuple) and len(entry.elts) == 3:
                    if isinstance(entry.elts[1], ast.Name):
                        found.add(entry.elts[1].id)
    return found


def _threshold_comparisons(source: str) -> list[tuple[int, str, str]]:
    """(line, scope, expression) of every comparison or subtraction on a threshold value.

    The scope is the dotted path of enclosing classes and functions.
    """
    hits: list[tuple[int, str, str]] = []

    def visit(node: ast.AST, scope: str, aliases: frozenset[str]) -> None:
        for child in ast.iter_child_nodes(node):
            inner = scope
            if isinstance(child, ast.FunctionDef | ast.AsyncFunctionDef | ast.ClassDef):
                inner = f"{scope}.{child.name}" if scope else child.name
            known = aliases | _ladder_values(child)
            names = [n for op in _operands(child) for n in _names(op)]
            if any(THRESHOLD_VALUE.search(n) or n in known for n in names):
                hits.append((child.lineno, inner, ast.unparse(child)))
            visit(child, inner, known)

    visit(ast.parse(source), "", frozenset())
    return hits


def _offenders(rel: str, source: str, used: set[tuple[str, str, str]]) -> list[str]:
    """Hits in one file that no allowlist entry covers; adds the entries it uses to `used`."""
    offenders: list[str] = []
    for lineno, scope, expression in _threshold_comparisons(source):
        if (rel, scope, expression) in NOT_A_METRIC_THRESHOLD:
            used.add((rel, scope, expression))
            continue
        line = source.splitlines()[lineno - 1].strip()
        offenders.append(f"  {rel}:{lineno} ({scope}): {line}")
    return offenders


# One private classification per spelling the scanner reads.
SPELLINGS = {
    "operator": "def f(v, thr_val):\n    return v < thr_val\n",
    "subtraction": (
        "def f(s, threshold_value):\n    return np.maximum(s - threshold_value, 0)\n"
    ),
    "augmented": "def f(v, thr):\n    v -= thr\n    return v\n",
    "string_key": "def f(w):\n    return w['value'] >= w['threshold_value']\n",
    "numpy_call": "def f(v, thr):\n    return np.greater_equal(v, thr)\n",
    "pandas_method": "def f(s, thr):\n    return s.ge(thr)\n",
    "unpacked_ladder": (
        "def f(v, ladder):\n"
        "    for i, (label, level, good) in enumerate(ladder):\n"
        "        yield v >= level\n"
    ),
    "ladder_comprehension": (
        "def f(v, thresholds):\n    return [v < lv for _, lv, _ in thresholds]\n"
    ),
}


@pytest.mark.policy
def test_no_threshold_comparison_outside_the_definition() -> None:
    """Oracle: AGENTS.md section 5; only `windows.margin` compares a metric with a threshold."""
    # Positive control on the scanner itself: every shape it exists to find, and the one
    # shape it must let through.
    for shape, source in SPELLINGS.items():
        assert _threshold_comparisons(source), f"the scanner misses {shape}"
    assert not _threshold_comparisons(
        "def h(v, thr):\n    return windows.margin(v, thr, True) >= 0\n"
    )
    # An entry exempts one expression in one scope: not a method of the same name, and not
    # a second comparison beside the allowlisted one.
    filter_py = "src/quebra/transforms/filter.py"
    method = "class A:\n    def run(self, v, thr):\n        return v >= thr\n"
    assert _threshold_comparisons(method)[0][1] == "A.run"
    assert _offenders(filter_py, method, set())
    assert _offenders(filter_py, "def run(v, thr):\n    return v >= thr\n", set())

    offenders: list[str] = []
    used: set[tuple[str, str, str]] = set()
    for root in SCAN_ROOTS:
        for path in sorted((REPO_ROOT / root).rglob("*.py")):
            if path == DEFINITION or "__pycache__" in path.parts:
                continue
            rel = path.relative_to(REPO_ROOT).as_posix()
            offenders += _offenders(rel, path.read_text(encoding="utf-8"), used)
    stale = sorted(NOT_A_METRIC_THRESHOLD - used)
    assert not stale, f"allowlist entries that match nothing in the tree: {stale}"
    assert not offenders, (
        "A metric is compared with a threshold outside analyzers/windows.py. Classify "
        "through windows.in_spec_mask or read the carve's in_spec column instead.\n"
        + "\n".join(offenders)
    )


def _panel(t_h: np.ndarray, series: np.ndarray, thresholds):
    carved = windows.run(
        windows.WindowsInputs(
            t_rel_s=t_h * 3600.0,
            values=series,
            thresholds=thresholds,
            dataset_id="definition",
        )
    )
    data = build_within_calibration_panel_data(
        t_h=t_h,
        primary_series=series,
        primary_label="metric",
        thresholds=thresholds,
        meta={"dataset": "unit-test"},
        windows=carved.windows_in_spec,
        reads=carved.reads,
        gap_spans_s=carved.diagnostics["gap_spans_s"],
        kaplan_meier=kaplan_meier.kaplan_meier_set(
            event_table.event_tables_from_carve(carved)
        ),
    )
    return carved, data


@pytest.mark.unit
def test_negating_the_margin_flips_every_quantity_that_follows_from_it(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Oracle: the definition. With margin negated, in spec and out of spec trade places."""
    t_h = np.linspace(0.0, 12.0, 481)
    series = 3.0 + np.sin(t_h + 0.3)
    thresholds = [("3", 3.0, True)]
    assert not np.any(series == 3.0)  # equality is the one place negation is not a swap

    base_carve, base = _panel(t_h, series, thresholds)
    # Bound as a default, so rebinding `original` later cannot reach the patch.
    original = windows.margin
    monkeypatch.setattr(
        windows, "margin", lambda v, thr, good, defn=original: -defn(v, thr, good)
    )
    flip_carve, flip = _panel(t_h, series, thresholds)

    cols = ["t_birth_s", "t_death_s", "censored", "n_reads"]
    for flipped, expected in (
        (flip_carve.windows_in_spec, base_carve.windows_out_of_spec),
        (flip_carve.windows_out_of_spec, base_carve.windows_in_spec),
    ):
        assert (
            flipped[cols]
            .reset_index(drop=True)
            .equals(expected[cols].reset_index(drop=True))
        )

    base_in = base_carve.reads["in_spec"].to_numpy(dtype=bool)
    flip_in = flip_carve.reads["in_spec"].to_numpy(dtype=bool)
    np.testing.assert_array_equal(flip_in, ~base_in)

    assert flip.reliability.occupancy["3"] == pytest.approx(
        1.0 - base.reliability.occupancy["3"], abs=1e-12
    )
    first_in_h = float(t_h[np.flatnonzero(base_in)[0]] - t_h[0])
    assert flip.reliability.ttf_per_threshold["3"] == pytest.approx(first_in_h)
    total_h = float(t_h[-1] - t_h[0])
    assert flip.reliability.cumulative_time_per_threshold["3"][-1] == pytest.approx(
        total_h - base.reliability.cumulative_time_per_threshold["3"][-1]
    )

    swap = {windows.STATE_IN_SPEC: windows.STATE_OUT_OF_SPEC}
    swap.update({v: k for k, v in swap.items()})
    assert flip.reliability.compliance_state_series["3"] == [
        (start, end, swap[state])
        for start, end, state in base.reliability.compliance_state_series["3"]
    ]
    # Damage integrates the excess beyond the threshold; negated, it integrates the
    # excess on the other side.
    assert flip.reliability.cumulative_damage_per_threshold["3"][-1] == pytest.approx(
        np.trapezoid(np.maximum(series - 3.0, 0.0), t_h), rel=1e-12
    )
    assert base.reliability.cumulative_damage_per_threshold["3"][-1] == pytest.approx(
        np.trapezoid(np.maximum(3.0 - series, 0.0), t_h), rel=1e-12
    )


@pytest.mark.unit
@pytest.mark.parametrize(
    ("values", "big_values_good"),
    [([10.0, 5.0, 5.0, 1.0, 1.0], True), ([1.0, 5.0, 5.0, 10.0, 10.0], False)],
    ids=["big_values_good", "small_values_good"],
)
def test_reads_at_the_threshold_are_in_spec_for_every_consumer(
    values: list[float], big_values_good: bool
) -> None:
    """Oracle: a hand trace. Reads 1 and 2 sit exactly on the threshold and are in spec.

    One-hour spacing. Reads 0-2 are in spec and 3-4 out, so the carve holds one window
    born at scan start and dying at the first out-of-spec read (3 h), the in-spec fraction
    charges three of four intervals (0.75), the first crossing is at 3 h, and one hour is
    out of spec. Calling a read at the threshold out of spec moves every one of these.
    """
    t_h = np.arange(5.0)
    series = np.asarray(values)
    thresholds = [("5", 5.0, big_values_good)]
    carved, data = _panel(t_h, series, thresholds)

    w = carved.windows_in_spec
    assert len(w) == 1
    assert w["birth_type"].iloc[0] == windows.BIRTH_SCAN_START
    assert w["death_type"].iloc[0] == windows.DEATH_DOWN_CROSSING
    assert w["duration_s"].iloc[0] == pytest.approx(3 * 3600.0)

    rel = data.reliability
    assert rel.occupancy["5"] == pytest.approx(0.75, abs=1e-12)
    assert rel.ttf_per_threshold["5"] == pytest.approx(3.0)
    assert rel.cumulative_time_per_threshold["5"][-1] == pytest.approx(1.0)
    assert rel.compliance_state_series["5"] == [
        (0.0, 3.0, windows.STATE_IN_SPEC),
        (3.0, 4.0, windows.STATE_OUT_OF_SPEC),
    ]
