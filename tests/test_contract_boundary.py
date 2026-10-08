"""Where `icontract` may be used, and the measured reason for the boundary.

Oracle: `icontract.ViolationError`'s own MRO, and the exception tuples `check_ledger` catches.
Both are read at runtime rather than quoted, so a change in either breaks this file.

THE FINDING. `icontract.ViolationError` subclasses `AssertionError`, NOT `ValueError`.
`analyzers/check_ledger.py` catches `(ValueError, KeyError)` at four sites to convert a check
failure into a REPORTED ROW rather than aborting the run. A contract inside `analyzers/checks/`
would therefore raise straight past those handlers and kill the job - which is exactly the
control flow SPEC 0008 exists to prevent: the band is always computed and always drawn, and a
check outcome annotates it rather than stopping it.

So contracts are allowed where a violation SHOULD abort - a caller handed the estimator the
wrong dtype - and forbidden where the ledger is meant to catch and report. That is a boundary,
not an inconsistency, and it is enforced here rather than left to a comment.
"""

from __future__ import annotations

import ast
from pathlib import Path

import icontract
import numpy as np
import pytest

from quebra.analyzers import kaplan_meier
from quebra.analyzers.kaplan_meier import KNOWN_DEATH_TYPES, KaplanMeierInputs

pytestmark = pytest.mark.policy


SRC = Path(__file__).resolve().parents[1] / "src" / "quebra"
CHECKS = SRC / "analyzers" / "checks"


def _imports_icontract(path: Path) -> bool:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            if any(a.name.split(".")[0] == "icontract" for a in node.names):
                return True
        elif isinstance(node, ast.ImportFrom):
            if (node.module or "").split(".")[0] == "icontract":
                return True
    return False


def test_a_violation_error_is_not_a_value_error():
    """The measured fact the boundary rests on. If icontract ever makes ViolationError a
    ValueError, this fails and the boundary can be revisited."""
    assert issubclass(icontract.ViolationError, AssertionError)
    assert not issubclass(icontract.ViolationError, ValueError)


def test_no_module_under_analyzers_checks_imports_icontract():
    """The boundary itself. A contract here would abort a run that the ledger is supposed to
    turn into a reported row."""
    offenders = sorted(p.name for p in CHECKS.glob("*.py") if _imports_icontract(p))
    assert offenders == [], (
        f"icontract must not be used inside analyzers/checks/: {offenders}. A "
        f"ViolationError escapes check_ledger's `except (ValueError, KeyError)` and kills "
        f"the job instead of producing a 'not computed' row."
    )


def test_the_one_place_it_is_used_is_the_one_the_spec_named():
    """Positive control. Without it the boundary test above would pass on a package that
    dropped contracts entirely, and R8.5b would have decided nothing."""
    users = sorted(
        str(p.relative_to(SRC)) for p in SRC.rglob("*.py") if _imports_icontract(p)
    )
    assert users == ["analyzers/kaplan_meier.py"], users


def test_the_contract_rejects_a_float_flag_array():
    """The invariant R8.0.4 named: a nonzero float would be read as an observed death."""
    inputs = KaplanMeierInputs(
        duration_min=np.array([1.0, 2.0, 3.0]),
        death_observed=np.array([1.0, 0.0, 1.0]),  # floats, not bools
        label="probe",
        threshold_label="thr",
    )
    with pytest.raises(icontract.ViolationError, match="boolean array"):
        kaplan_meier.run(inputs)


def test_the_same_data_as_real_bools_is_accepted():
    """The discriminator: the contract must reject the dtype, not the data."""
    inputs = KaplanMeierInputs(
        duration_min=np.array([1.0, 2.0, 3.0]),
        death_observed=np.array([True, False, True]),
        label="probe",
        threshold_label="thr",
    )
    curve = kaplan_meier.run(inputs)
    assert len(curve.time_min) > 1


def test_an_unknown_death_type_is_refused_by_a_bare_raise_not_a_contract():
    """The other half of the boundary. This check inspects a frame column and names the
    offending values, which a precondition cannot do usefully - so it stays a bare raise,
    and its exception type is one the ledger's handlers understand."""
    import pandas as pd

    windows = pd.DataFrame(
        {
            "threshold_label": ["thr", "thr"],
            "side": ["in_spec", "in_spec"],
            "t_before_birth_s": [0.0, 100.0],
            "t_birth_s": [10.0, 110.0],
            "t_last_s": [60.0, 230.0],
            "t_death_s": [70.0, 230.0],
            "birth_type": ["up_crossing", "up_crossing"],
            "death_type": ["down_crossing", "teleported"],
            "duration_s": [60.0, 120.0],
            "birth_observed": [True, True],
            "censored": [False, True],
        }
    )
    with pytest.raises(ValueError, match="unknown death_type"):
        kaplan_meier.make_inputs_from_windows(
            windows, threshold_label="thr", label="probe"
        )

    assert "teleported" not in KNOWN_DEATH_TYPES
