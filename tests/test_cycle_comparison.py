"""Two devices' calibration cycles compared through the record comparison, descriptively.

Synthetic status histories only: each device's cycles, checks and figure data come from the
functions its cycles job calls (`cycles_recipe`).
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import quebra.analyzers.cycles as cycles
from quebra.analyzers.check_outcome import PERMUTATION_KEYS
from quebra.analyzers.cycle_comparison import compare_cycle_pair
from quebra.analyzers.event_table import CYCLES
from quebra.analyzers.record_comparison import (
    NO_TEST,
    TAU_RULE_MIN_LARGEST_AGE,
    pair_table,
)
from quebra.core.types import Norm

_T0 = 1_700_000_000.0
_REPO = Path(__file__).resolve().parents[1]
# (start, end) of each calibrating spell, in seconds from the record start, and the record
# length. A: cycles 100, 130, 60 and 270 censored at the end; B: 80, 100, 50 and 30.
_A = ([(0, 10), (110, 120), (250, 260), (320, 330)], 600.0)
_B = ([(0, 10), (90, 100), (200, 210), (260, 270)], 300.0)


def _device(name: str, spec) -> tuple[cycles.CycleTable, cycles.CyclesPanelData]:
    calibrating, length_s = spec
    norm = Norm(
        {
            "spell_start_unix_s": np.array([_T0 + a for a, _ in calibrating], float),
            "spell_end_unix_s": np.array([_T0 + b for _, b in calibrating], float),
            "spell_status": np.array(["calibrating"] * len(calibrating)),
            "meta": {
                "dataset_id": name,
                "record_start_unix_s": _T0,
                "record_end_unix_s": _T0 + length_s,
            },
        }
    )
    table = cycles.cycles(norm)
    checks = cycles.cycle_checks(
        table,
        pd.read_csv(_REPO / "jobs/bench/results/size_table.csv"),
        run_set=PERMUTATION_KEYS,
        alpha=0.05,
        min_events_pass=35,
        tie_cutoff_distinct=5,
        lag_max=5,
        n_permutations=99,
        seed=7,
        min_events_per_segment=2,
    )
    return table, cycles.build_cycles_panel_data(
        table, checks, cycles.status_shares(norm)
    )


def _pair():
    (table_a, panel_a), (table_b, panel_b) = _device("A", _A), _device("B", _B)
    pairs = compare_cycle_pair(
        table_a, table_b, panel_a, panel_b, tau_rule=TAU_RULE_MIN_LARGEST_AGE
    )
    return pairs, panel_a, panel_b


@pytest.mark.unit
def test_the_devices_restricted_means_to_the_shared_tau_match_a_hand_table() -> None:
    """Oracle: a hand table. A's cycles are 100, 130, 60 and a censored 270; B's 80, 100,
    50 and a censored 30. tau is the smaller largest age, 100, set by B. A's curve is 1 to
    60 and 3/4 to 100: 60 + 40 * 3/4 = 90. B's is 1 to 50, 2/3 to 80 and 1/3 to 100:
    50 + 20 + 20/3 = 230/3. The point difference is 90 - 230/3 = 40/3.

    Each MTBC runs to its device's largest cycle. B's is 100, the shared tau, so B's MTBC is
    its restricted mean, 230/3, SE included. A's is 270: its curve is 1/2 to 130 and 1/4 to
    270, so its MTBC is 90 + 30 / 2 + 140 / 4 = 140, above its restricted mean."""
    pairs, panel_a, panel_b = _pair()
    assert pairs.ladder == ("",) and list(pairs.items) == [("", CYCLES)]
    pair = pairs.get("", CYCLES)
    assert (pair.record_a, pair.record_b) == ("A", "B")
    assert pair.tau.rmst_tau_s == 100.0 and pair.tau.tau_set_by == "B"
    assert pair.summary_a.rmst_s == pytest.approx(90.0, abs=1e-12)
    assert pair.summary_b.rmst_s == pytest.approx(230.0 / 3.0, abs=1e-12)
    assert pair.rmst_difference_s == pytest.approx(40.0 / 3.0, abs=1e-12)
    assert (panel_b.mtbc.tau_s, panel_a.mtbc.tau_s) == (100.0, 270.0)
    assert panel_b.mtbc.rmst_s == pytest.approx(230.0 / 3.0, abs=1e-12)
    assert panel_b.mtbc.rmst_se_s == pytest.approx(pair.summary_b.rmst_se_s, abs=1e-12)
    assert panel_a.mtbc.rmst_s == pytest.approx(140.0, abs=1e-12)
    counts = (pair.n_windows_a, pair.n_deaths_a, pair.n_windows_b, pair.n_deaths_b)
    assert counts == (4, 3, 4, 3)


@pytest.mark.unit
def test_the_comparison_tests_nothing_and_says_so_in_its_table() -> None:
    """Oracle: specification. Descriptive only: no log-rank and no SE for the difference,
    each None with the reason in its table row, while each device keeps its own interval
    and the check outcome its own cycles figure draws."""
    pairs, panel_a, panel_b = _pair()
    pair = pairs.get("", CYCLES)
    assert not pair.tested
    assert pair.logrank_z is None and pair.logrank_p is None
    assert pair.rmst_difference_se_s is None
    assert pair.summary_a.rmst_se_s > 0.0 and pair.summary_b.rmst_se_s > 0.0
    assert pair.checks_a == panel_a.checks and pair.checks_b == panel_b.checks
    assert pair.summary_a.checks == panel_a.checks
    assert pair.summary_b.checks == panel_b.checks
    table = pair_table(pairs, comparison="A2", evidence="two devices").set_index(
        ["record", "quantity"]
    )
    for quantity in ("rmst_difference_se_s", "logrank_z", "logrank_chi2", "logrank_p"):
        row = table.loc[("A vs B", quantity)]
        assert pd.isna(row["value"]) and row["note"] == NO_TEST, quantity
    assert table.loc[("A vs B", "rmst_difference_s"), "value"] == pytest.approx(40 / 3)
    outcome = table.loc[("A vs B", "check_outcome"), "note"]
    assert outcome == f"A: {panel_a.checks.summary()}; B: {panel_b.checks.summary()}"


@pytest.mark.unit
def test_a_cycles_figure_of_the_other_device_is_refused() -> None:
    """Oracle: specification. The check outcome attached to a device must be its own, on
    either side of the pair."""
    (table_a, panel_a), (table_b, panel_b) = _device("A", _A), _device("B", _B)
    with pytest.raises(
        ValueError, match="cycles figure of 'B' given for the cycles of 'A'"
    ):
        compare_cycle_pair(
            table_a, table_b, panel_b, panel_b, tau_rule=TAU_RULE_MIN_LARGEST_AGE
        )
    with pytest.raises(
        ValueError, match="cycles figure of 'A' given for the cycles of 'B'"
    ):
        compare_cycle_pair(
            table_a, table_b, panel_a, panel_a, tau_rule=TAU_RULE_MIN_LARGEST_AGE
        )


@pytest.mark.integration
def test_the_a2_job_compares_the_two_devices_own_cycles_and_figures(in_repo) -> None:
    """Oracle: specification. The composite builds no cycle: its one comparison step reads
    each device job's materialized cycles and cycles figure data, Tuna-17 first, and its
    table is materialized."""
    from quebra.cli import _module_from_path

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        job = _module_from_path(_REPO / "jobs/composite/compare_a2_tuna17_tuna9.py").job
    steps = {node.fn.__name__: node for node in job.dag.values()}
    assert {fn: node.fn_name for fn, node in steps.items()} == {
        "compare_cycle_pair": "a2_pair",
        "pair_table": "a2_pair_table",
    }
    refs = [
        (r.included.job.name, r.node_name) for r in steps["compare_cycle_pair"].inputs
    ]
    assert refs == [
        ("cycles_tuna17", "cycles_tuna17_cycles"),
        ("cycles_tuna9", "cycles_tuna9_cycles"),
        ("cycles_tuna17", "cycles_tuna17_panel_data"),
        ("cycles_tuna9", "cycles_tuna9_panel_data"),
    ]
    assert steps["compare_cycle_pair"].kwargs == {
        "tau_rule": TAU_RULE_MIN_LARGEST_AGE,
        "conf_level": 0.95,
    }
    assert {getattr(s, "name", "") for s in job.sinks} == {
        "compare_a2_tuna17_tuna9_table"
    }
