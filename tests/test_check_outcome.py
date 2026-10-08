"""The output-builder reshapes a ledger into grids, and never branches on a result.

Oracle: the ledger's own column contract (`check_ledger.LEDGER_COLUMNS`) and its verdict
vocabulary. These are reshape tests, not statistical ones - the builder computes no statistic,
it selects and pivots rows that `check_ledger.run` already scored.

The load-bearing claim is the one Sera restated: **a check outcome must never prevent
execution.** Control flow that depends on what the data happened to say is unpredictable, and
the whole design is that the band draws in every case and the checks annotate it. So the tests
below drive every unhappy result - `fail`, `underpowered`, `not computed`, a missing row, an
empty ledger - and assert that each one BUILDS and RENDERS. Only wiring mistakes raise.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

from quebra.analyzers import check_outcome as sel
from quebra.analyzers.check_outcome import (
    build_check_outcome,
)
from quebra.analyzers.checks.battery import ROW_KEYS, row_key, run_battery
from quebra.analyzers.checks.result import CALIB_PERMUTATION, Segment
from quebra.analyzers.check_ledger import (
    LEDGER_COLUMNS,
    VERDICT_ABSENT,
    VERDICT_FAIL,
    VERDICT_NOT_COMPUTED,
    VERDICT_PASS,
    VERDICT_UNDERPOWERED,
)
from quebra.plots.check_outcome_plot import CheckOutcomePlot

pytestmark = pytest.mark.unit


THRESHOLD = "3.0 µs"
DATASETS = ("rec_a", "rec_b")
CLOCKS = ("in_spec", "calendar")

C1_PERM = ("c1_lewis_robinson", "permutation", "")
C6_PERM = ("c6_exchangeability", "permutation", "")


def _row(dataset, clock, key, verdict, p_value=0.4, notes=""):
    blank = dict.fromkeys(LEDGER_COLUMNS)
    blank.update(
        dataset_id=dataset,
        threshold_label=THRESHOLD,
        clock=clock,
        check_id=key[0],
        calibration=key[1],
        variant=key[2],
        verdict=verdict,
        p_value=p_value,
        notes=notes,
    )
    return blank


def _ledger(*rows) -> pd.DataFrame:
    return pd.DataFrame(list(rows), columns=LEDGER_COLUMNS)


def _every_verdict_ledger() -> pd.DataFrame:
    """One of each verdict, so no test below rests on a single happy outcome."""
    return _ledger(
        _row("rec_a", "in_spec", C1_PERM, VERDICT_PASS, 0.62),
        _row("rec_a", "calendar", C1_PERM, VERDICT_FAIL, 0.004, "rejected"),
        _row("rec_b", "in_spec", C1_PERM, VERDICT_UNDERPOWERED, 0.31),
        _row("rec_b", "calendar", C1_PERM, VERDICT_NOT_COMPUTED, None, "tau singular"),
        _row("rec_a", "in_spec", C6_PERM, VERDICT_PASS, 0.55),
        _row("rec_a", "calendar", C6_PERM, VERDICT_PASS, 0.51),
        _row("rec_b", "in_spec", C6_PERM, VERDICT_FAIL, 0.01, "rejected"),
        _row("rec_b", "calendar", C6_PERM, VERDICT_PASS, 0.44),
    )


def _build(ledger, run_set=(C1_PERM, C6_PERM), display_set=(C1_PERM, C6_PERM)):
    return build_check_outcome(
        ledger,
        threshold_label=THRESHOLD,
        run_set=run_set,
        display_set=display_set,
        datasets=DATASETS,
        clocks=CLOCKS,
        alpha=0.05,
    )


# ------------------------------------------------- a result never prevents execution


@pytest.mark.parametrize(
    "verdict",
    [VERDICT_PASS, VERDICT_FAIL, VERDICT_UNDERPOWERED, VERDICT_NOT_COMPUTED],
)
def test_every_verdict_builds_and_renders(verdict):
    """The rule, one verdict at a time. A `fail` must not raise, must not empty the grid,
    and must reach a drawn figure exactly like a `pass` does."""
    ledger = _ledger(
        _row(
            "rec_a",
            "in_spec",
            C1_PERM,
            verdict,
            None if verdict == VERDICT_NOT_COMPUTED else 0.2,
        )
    )
    outcome = _build(ledger, run_set=(C1_PERM,), display_set=(C1_PERM,))
    assert len(outcome.grids) == 1
    assert outcome.grids[0].verdicts[0][0] == verdict
    figure = CheckOutcomePlot("probe").build_matplotlib(outcome)
    assert figure is not None


def test_absent_is_a_different_string_from_the_ledgers_own_not_computed():
    """R8.3: a cell nobody has a row for must not read the same as one the check declined to
    answer. Sharing one string would let a total lookup failure render as an honest grid."""
    assert VERDICT_ABSENT != VERDICT_NOT_COMPUTED
    from quebra.plots import theme

    assert theme.verdict_color(VERDICT_ABSENT) != theme.verdict_color(
        VERDICT_NOT_COMPUTED
    )


def test_cells_are_looked_up_by_dataset_id_while_labels_are_only_shown():
    """The regression test for that bug. The ledger keys rows by `dataset_id`; the figure
    shows something friendlier. Looking up by the label must not be possible to do by
    accident, so the two are separate arguments and only one of them is a key."""
    ledger = _ledger(_row("rec_a", "in_spec", C1_PERM, VERDICT_PASS, 0.7))
    outcome = build_check_outcome(
        ledger,
        threshold_label=THRESHOLD,
        run_set=(C1_PERM,),
        display_set=(C1_PERM,),
        datasets=("rec_a",),
        dataset_labels=("Record A",),
        clocks=CLOCKS,
        alpha=0.05,
    )
    grid = outcome.grids[0]
    assert grid.datasets == ("rec_a",), "the key must stay the ledger's own id"
    assert grid.dataset_labels == ("Record A",), "the label must reach the renderer"
    assert grid.verdicts[0][0] == VERDICT_PASS, (
        "the cell resolved to absent, so the lookup used the label instead of the id"
    )


def test_mismatched_label_count_raises_rather_than_mislabelling_a_row():
    """Positional pairing, so a length mismatch would put one record's name on another's row."""
    with pytest.raises(ValueError, match="dataset_labels has"):
        build_check_outcome(
            _ledger(),
            threshold_label=THRESHOLD,
            run_set=(C1_PERM,),
            display_set=(C1_PERM,),
            datasets=("rec_a", "rec_b"),
            dataset_labels=("only one",),
            clocks=CLOCKS,
            alpha=0.05,
        )


def test_an_entirely_empty_ledger_still_builds_and_renders():
    """The extreme case: nothing ran at all. Every cell reads absent and the figure draws,
    because a band with no check outcome is still a band that must be shown."""
    outcome = _build(_ledger(), run_set=(C1_PERM,), display_set=(C1_PERM,))
    assert outcome.grids[0].verdicts == ((VERDICT_ABSENT, VERDICT_ABSENT),) * len(
        DATASETS
    )
    assert CheckOutcomePlot("probe").build_matplotlib(outcome) is not None


# ------------------------------------------------- the reshape itself


def test_the_grid_axes_are_datasets_by_clocks_in_the_order_given():
    outcome = _build(_every_verdict_ledger())
    grid = outcome.grids[0]
    assert grid.datasets == DATASETS
    assert grid.clocks == CLOCKS
    assert grid.verdicts[0][0] == VERDICT_PASS  # rec_a / in_spec
    assert grid.verdicts[0][1] == VERDICT_FAIL  # rec_a / calendar
    assert grid.verdicts[1][0] == VERDICT_UNDERPOWERED  # rec_b / in_spec


def test_the_display_set_selects_and_the_run_set_is_still_recorded():
    """Running everything and showing one: the artifact must still say what ran."""
    outcome = _build(
        _every_verdict_ledger(), run_set=(C1_PERM, C6_PERM), display_set=(C6_PERM,)
    )
    assert [g.key for g in outcome.grids] == [C6_PERM]
    assert set(outcome.run_set) == {C1_PERM, C6_PERM}
    assert outcome.display_set == (C6_PERM,)


def test_a_missing_cell_is_absent_not_a_fabricated_verdict():
    """rec_b has no C1 row at all here. The cell must say so rather than inherit rec_a's."""
    ledger = _ledger(_row("rec_a", "in_spec", C1_PERM, VERDICT_PASS, 0.7))
    grid = _build(ledger, run_set=(C1_PERM,), display_set=(C1_PERM,)).grids[0]
    assert grid.verdicts[0][0] == VERDICT_PASS
    assert grid.verdicts[1][0] == VERDICT_ABSENT
    assert grid.p_values[1][0] is None


def test_a_threshold_absent_from_the_ledger_raises_rather_than_rendering_absent():
    """A wiring mistake, not a result. Silently yielding an all-absent grid renders as an
    honest figure about a rung nobody scored - the same mistake `kaplan_meier` raises on."""
    ledger = _ledger(_row("rec_a", "in_spec", C1_PERM, VERDICT_PASS, 0.7))
    with pytest.raises(KeyError, match="is not in the ledger"):
        build_check_outcome(
            ledger,
            threshold_label="9.0 µs",
            run_set=(C1_PERM,),
            display_set=(C1_PERM,),
            datasets=DATASETS,
            clocks=CLOCKS,
            alpha=0.05,
        )


def test_rows_at_another_threshold_are_not_mixed_in():
    ledger = _ledger(
        _row("rec_a", "in_spec", C1_PERM, VERDICT_PASS, 0.7),
        {
            **_row("rec_a", "calendar", C1_PERM, VERDICT_FAIL, 0.01),
            "threshold_label": "9.0 µs",
        },
    )
    grid = _build(ledger, run_set=(C1_PERM,), display_set=(C1_PERM,)).grids[0]
    assert grid.verdicts[0][1] == VERDICT_ABSENT, "a row from another rung leaked in"


def test_asked_for_but_unanswered_is_reported_separately_from_never_asked():
    """The distinction R8.3 insisted on, now visible on the artifact."""
    ledger = _ledger(_row("rec_a", "in_spec", C1_PERM, VERDICT_PASS, 0.7))
    outcome = _build(ledger, run_set=(C1_PERM, C6_PERM), display_set=(C1_PERM,))
    assert outcome.asked_but_unanswered == (C6_PERM,)

    # The discriminator: a row nobody asked for is not reported as unanswered.
    only_c1 = _build(ledger, run_set=(C1_PERM,), display_set=(C1_PERM,))
    assert only_c1.asked_but_unanswered == ()


# ------------------------------------------------- only wiring mistakes raise


def test_an_empty_display_set_builds_and_renders_a_stated_placeholder():
    """Showing none of them is legitimate configuration, not an error."""
    outcome = _build(_every_verdict_ledger(), display_set=())
    assert outcome.grids == []
    assert CheckOutcomePlot("probe").build_matplotlib(outcome) is not None


def test_check_thresholds_refuses_a_ladder_the_artifact_was_not_built_at():
    outcome = _build(_every_verdict_ledger())
    outcome.check_thresholds([THRESHOLD])
    with pytest.raises(ValueError, match="not on the requested ladder"):
        outcome.check_thresholds(["9.0 µs"])


# -------------------- the run-set and display-set (R8.3)


N_PERM = 99
SEED = 5


def _perm_segment() -> Segment:
    """One well-posed segment: tau clears T_N, gaps positive, enough of them for every row."""
    x = np.array([1.0, 2.0, 1.5, 2.5, 1.2, 1.8, 2.2, 1.1])
    return Segment(x=x, tau=20.0)


def _run_battery_on(**flags):
    return run_battery(
        [_perm_segment()], n_perm=N_PERM, rng=np.random.default_rng(SEED), **flags
    )


# --------------------------------------------------------------- the vocabulary


def test_tau_dependent_keys_match_what_the_battery_actually_drops():
    """The oracle. Derived by running the battery both ways, so the constant cannot drift
    away from the flag it describes."""
    with_tau = {row_key(r) for r in _run_battery_on()}
    without = {row_key(r) for r in _run_battery_on(include_tau_checks=False)}
    assert sel.TAU_DEPENDENT_KEYS == frozenset(with_tau - without)


# --------------------------------------------------------------- normalisation


def test_normalise_imposes_a_canonical_order_so_identity_does_not_depend_on_typing_order():
    """Two jobs naming the same rows must produce the same step kwarg, or they would get
    different run identities for the same analysis."""
    a = sel.normalise([ROW_KEYS[3], ROW_KEYS[0]], field="run_set")
    b = sel.normalise([ROW_KEYS[0], ROW_KEYS[3]], field="run_set")
    assert a == b
    assert list(a) == [k for k in sel.ALL_KEYS if k in {ROW_KEYS[0], ROW_KEYS[3]}]


def test_a_two_part_key_is_refused_rather_than_padded():
    """Addressing by check id alone is the granularity this design rejects; silently
    accepting a short key would reintroduce it."""
    with pytest.raises(ValueError, match="triple"):
        sel.normalise([("c1_lewis_robinson", CALIB_PERMUTATION)], field="run_set")


# --------------------------------------------------------------- the two selections


def test_the_two_are_independent():
    """Changing the display-set must not change the run-set, and vice versa."""
    run_a, display_a = sel.resolve(sel.ALL_KEYS, [ROW_KEYS[1]])
    run_b, display_b = sel.resolve(sel.ALL_KEYS, [ROW_KEYS[1], ROW_KEYS[7]])
    assert run_a == run_b
    assert display_a != display_b

    run_c, display_c = sel.resolve(sel.PERMUTATION_KEYS, [ROW_KEYS[1]])
    assert display_c == display_a
    assert run_c != run_a


def test_a_display_row_outside_the_run_set_raises():
    """It would render as a blank that reads 'no answer' when the truth is 'never asked'."""
    with pytest.raises(ValueError, match="outside the run_set"):
        sel.resolve([ROW_KEYS[1]], [ROW_KEYS[0]])


def test_an_empty_run_set_raises():
    with pytest.raises(ValueError, match="run_set is empty"):
        sel.resolve([], [])


# --------------------------------------------------------------- derived flags


def test_the_derived_flags_really_drive_the_battery():
    """Ties the derived flags to observed behaviour rather than to their own definition."""
    run = tuple(k for k in sel.ALL_KEYS if k not in sel.TAU_DEPENDENT_KEYS)
    flags = sel.battery_flags(run)
    produced = {
        row_key(r)
        for r in _run_battery_on(include_tau_checks=flags["include_tau_checks"])
    }
    assert not (produced & sel.TAU_DEPENDENT_KEYS)


# --------------------------------------------------------------- filtering


def test_a_row_that_answered_nothing_is_distinguishable_from_one_never_asked_for():
    """R8.3's last acceptance, on the real drop. Asking for the tau rows and then running
    without them must report them MISSING, not silently return a shorter list."""
    results = _run_battery_on(include_tau_checks=False)
    asked_for_everything = sel.missing_rows(results, ROW_KEYS, row_key)
    assert set(asked_for_everything) == set(sel.TAU_DEPENDENT_KEYS)

    # The discriminator: rows never asked for are not reported missing.
    asked_for_what_ran = sel.missing_rows(
        results, [k for k in ROW_KEYS if k not in sel.TAU_DEPENDENT_KEYS], row_key
    )
    assert asked_for_what_ran == ()


# -------------------- the band carries its check outcome (AGENTS.md section 5)


def test_a_band_artifact_always_carries_its_check_fields():
    """The positive control the phase's central claim was missing.

    Oracle: `AGENTS.md` section 5, "a band is never reported without its check outcome
    attached". Before this, the claim was asserted only by a docstring - a test named
    "nothing but failures renders" never touched a survival curve at all, so a job that
    suppressed a band on a `fail` would have kept the suite green.

    Attachment is structural: the fields exist on every `KaplanMeierComparison`, and an
    unassessed band says so rather than saying nothing.
    """
    from quebra.analyzers.kaplan_meier import KaplanMeierComparison

    bare = KaplanMeierComparison(curves=[], threshold_label=THRESHOLD)
    assert bare.checks_asked == ()
    assert "NOT ASSESSED" in bare.check_summary(), (
        "a band with no checks must SAY it was not assessed; silence reads as 'fine'"
    )


def test_an_attached_outcome_names_what_was_asked_and_what_did_not_answer():
    """Section 5's second half: naming both is the part that stops an all-grey grid from
    satisfying the rule while saying nothing."""
    from quebra.analyzers.kaplan_meier import KaplanMeierComparison

    attached = KaplanMeierComparison(
        curves=[],
        threshold_label=THRESHOLD,
        assumption_id="a1_renewal_durations",
        checks_asked=("c1 permutation", "c2 permutation"),
        checks_unanswered=("c2 permutation",),
        check_verdicts=(("C1", "rec_a", VERDICT_PASS), ("C1", "rec_b", VERDICT_FAIL)),
    )
    summary = attached.check_summary()
    assert "2 checks asked" in summary
    assert "1 fail" in summary and "1 pass" in summary
    assert "c2 permutation" in summary, "an unanswered check must be named, not dropped"
    assert attached.assumption_id == "a1_renewal_durations"


def test_a_failing_check_does_not_empty_the_band_artifact():
    """The rule the whole phase exists to serve: a `fail` annotates, it never suppresses.
    Asserted on the band artifact, not on the grid figure - that was the gap."""
    from quebra.analyzers.kaplan_meier import KaplanMeierComparison

    all_failed = KaplanMeierComparison(
        curves=[],
        threshold_label=THRESHOLD,
        checks_asked=("c1 permutation",),
        check_verdicts=(("C1", "rec_a", VERDICT_FAIL), ("C1", "rec_b", VERDICT_FAIL)),
    )
    assert all_failed.checks_asked, "the outcome must survive a total rejection"
    assert "2 fail" in all_failed.check_summary()


# ------------ the run-set restricts what the LEDGER computes, not only what it names


def test_the_battery_switches_are_derived_from_the_run_set_and_reach_the_ledger():
    """The run-set must RESTRICT, not merely declare.

    Oracle: `CheckLedgerInputs`' own fields. Before this the builder accepted only
    `include_c3`, so a run-set excluding the asymptotic rows still computed them and the
    `run_set` written onto the artifact was false.
    """
    from quebra.analyzers import check_ledger

    class _W:
        windows_in_spec = None
        meta: dict = {}
        diagnostics: dict = {}

    inputs = check_ledger.make_inputs_from_windows(
        _W(),
        None,
        thresholds=[(THRESHOLD, 3.0e-6, True)],
        alpha=0.05,
        min_events_pass=35,
        tie_cutoff_distinct=5,
        lag_max=5,
        n_permutations=19,
        seed=3,
        **sel.battery_flags(sel.PERMUTATION_KEYS),
    )
    assert inputs.include_c2_asymptotic is False, (
        "PERMUTATION_KEYS excludes C2-asymptotic, so the ledger must not compute it"
    )
    assert inputs.include_c3 is False
    assert inputs.include_tau_checks is True  # C1/C2 permutation rows are tau rows


def test_the_switches_restrict_but_are_coarser_than_the_run_set():
    """MEASURED, so the caveat in `check_ledger.py` is a number rather than a hedge.

    `PERMUTATION_KEYS` declares 6 rows. The battery computes 9 unrestricted and 8 under the
    derived switches: restriction is real, and the 2-row overshoot is `include_tau_checks`
    being all-or-nothing over five rows, so C1-asymptotic and CvM-asymptotic ride along with
    C1's and C2's permutation rows. The run-set names AT MOST what ran.

    If this fails because the numbers moved, the caveat needs re-measuring, not deleting.
    """
    unrestricted = {row_key(r) for r in _run_battery_on()}
    restricted = {
        row_key(r)
        for r in _run_battery_on(
            **{
                k: v
                for k, v in sel.battery_flags(sel.PERMUTATION_KEYS).items()
                if k != "include_c3"
            }
        )
    }
    assert len(unrestricted) == 9
    assert len(restricted) == 8, (
        "the run-set no longer restricts what the battery computes"
    )
    assert unrestricted - restricted == {("c2_anderson_darling", "asymptotic", "")}
    assert len(restricted) - len(sel.PERMUTATION_KEYS) == 2, (
        "the coarseness gap moved; re-measure the caveat in check_ledger.py"
    )
