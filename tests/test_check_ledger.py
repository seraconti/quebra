"""The ledger's verdict rule is the decision this whole increment exists to make.

A p-value alone must never produce a `pass`. These tests pin each of the three conditions
independently, so a regression that drops one of them fails here rather than in a figure.

The file also pins what the ledger records about the R that produced its C3 rows, and that
it still builds every row when R is absent, unusable, or was never asked for.
"""

from __future__ import annotations

import pathlib
import sys
from itertools import product

import numpy as np
import pandas as pd
import pytest

from quebra.analyzers import windows
from quebra.analyzers.check_ledger import (
    VERDICT_FAIL,
    VERDICT_INCONCLUSIVE,
    VERDICT_NOT_COMPUTED,
    VERDICT_PASS,
    VERDICT_TIES,
    _verdict,
    stream_for,
)
from quebra.analyzers.checks import battery
from quebra.analyzers.checks.result import (
    CALIB_PERMUTATION,
    CALIB_R_COPULA,
    CLOCK_CALENDAR,
    CLOCK_SIDE,
    CheckResult,
)
from quebra.analyzers.windows import SIDE_IN_SPEC

pytestmark = pytest.mark.unit


def _result(check: str = "c6_exchangeability", p: float | None = 0.4) -> CheckResult:
    return CheckResult(
        check=check,
        statistic=1.0,
        p_value=p,
        calibration=CALIB_PERMUTATION,
        clock="in_spec",
        n_events=100,
        n_segments=1,
    )


BASE = dict(
    alpha=0.05,
    n_events=100,
    min_events_pass=35,
    n_distinct=100,
    tie_cutoff=5,
    bench_accepted=True,
)


def test_all_three_conditions_met_is_the_only_route_to_pass():
    """Oracle: the verdict rule, at the event floor itself: n_events == min_events_pass
    meets the floor."""
    verdict, reason = _verdict(_result(), **BASE)
    assert verdict == VERDICT_PASS and reason == ""
    verdict, reason = _verdict(_result(), **{**BASE, "n_events": 35})
    assert verdict == VERDICT_PASS and reason == ""


@pytest.mark.parametrize(
    "override, expected_reason",
    [
        ({"n_events": 34}, "n_events 34 < 35"),
        (
            {"bench_accepted": False},
            "bench found this check miscalibrated at this event count",
        ),
        ({"bench_accepted": None}, "no bench cell for this check at this event count"),
    ],
)
def test_each_condition_alone_makes_a_non_rejection_inconclusive(
    override, expected_reason
):
    """Oracle: the verdict rule. A non-rejection with any one of the three conditions unmet
    (too few events, a size the bench found wrong, no bench cell) is `inconclusive`, never
    `pass`, and its reason is that condition alone, one below the event floor included."""
    verdict, reason = _verdict(_result(), **{**BASE, **override})
    assert verdict == VERDICT_INCONCLUSIVE
    assert reason == expected_reason


def test_a_rejection_is_a_fail_and_carries_the_calibration_caveat():
    verdict, reason = _verdict(_result(p=0.01), **BASE)
    assert verdict == VERDICT_FAIL and "rejected" in reason
    verdict, reason = _verdict(_result(p=0.01), **{**BASE, "bench_accepted": False})
    assert verdict == VERDICT_FAIL
    assert "miscalibrated" in reason, (
        "a rejection from a miscalibrated check must say so, or the ledger overstates it"
    )


def test_ties_block_a_rank_check_but_not_a_duration_check():
    """C5/C6/C3 are rank or copula statistics; C1/C2 use the durations directly."""
    thin = {**BASE, "n_distinct": 3}
    assert _verdict(_result("c6_exchangeability"), **thin)[0] == VERDICT_TIES
    assert _verdict(_result("c5_rank_autocorr"), **thin)[0] == VERDICT_TIES
    assert _verdict(_result("c1_lewis_robinson"), **thin)[0] == VERDICT_PASS


def test_a_missing_p_value_is_not_computed_and_keeps_its_reason():
    result = _result(p=None)
    object.__setattr__(result, "notes", "R unavailable")
    verdict, reason = _verdict(result, **BASE)
    assert verdict == VERDICT_NOT_COMPUTED and reason == "R unavailable"


def test_the_precedence_is_not_computed_then_ties_then_rejection():
    """Order matters: a rank check with 3 distinct values has no interpretable p at all,
    so ties must be checked before the p-value is read."""
    verdict, _ = _verdict(_result(p=0.001), **{**BASE, "n_distinct": 3})
    assert verdict == VERDICT_TIES, "ties outrank a rejection for a rank statistic"


def test_the_per_threshold_seed_is_stable_across_processes():
    """Oracle: the shipped `check_ledger.stream_for`, evaluated under three different
    PYTHONHASHSEED values.

    `hash()` is randomised per process, so a seed built from it would make every
    permutation p-value differ between runs while provenance stayed identical.

    This test calls the SHIPPED derivation. An earlier version re-typed the crc32
    expression as a string literal and so asserted only that `zlib.crc32` is a function
    of its argument - it stayed green when the module was mutated to use `hash()`, which
    is the whole defect it existed to catch.
    """
    import os
    import subprocess
    import sys

    code = (
        "from quebra.analyzers.check_ledger import stream_for;"
        "print(stream_for('3 \\u00b5s', 'in_spec'))"
    )
    seen = set()
    for hashseed in ("0", "1", "12345"):
        env = {**os.environ, "PYTHONHASHSEED": hashseed}
        out = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, env=env
        )
        assert out.returncode == 0, out.stderr
        seen.add(out.stdout.strip())
    assert len(seen) == 1, f"per-threshold stream is not stable: {seen}"
    # In-process agreement too, so a mutation cannot pass by being uniformly wrong across
    # subprocesses while disagreeing with the value the ledger actually uses.
    assert seen == {str(stream_for("3 µs", "in_spec"))}


def test_the_ledger_derives_its_stream_through_stream_for(monkeypatch):
    """Oracle: the call site at `check_ledger.run`, not the helper it calls.

    The test above pins `stream_for` itself. That is not enough on its own: extracting the
    expression narrowed the mutable surface rather than covering it, so re-inlining
    `hash(...)` at the call site would leave `stream_for` correct, unused, and the suite
    green. This asserts the routing, so both halves of the seed derivation are pinned.

    Each row's clock is paired with the stream it drew, against the names written out here
    (issue #19 addendum, decision 1): the side clock's stream and bench cells are named
    `in_spec`, as before the rename, so the in-spec bag keeps its p-values and bench cells.
    Reading the names back from `BENCH_CLOCK` would pass with the mapping swapped.
    """
    import quebra.analyzers.check_ledger as ledger_module
    from quebra.analyzers.checks.result import BENCH_CLOCK

    assert BENCH_CLOCK == {"side": "in_spec", "calendar": "calendar"}

    calls: list[tuple[str, str]] = []
    row_clock: list[str] = []
    pairs: set[tuple[str, str]] = set()
    score_segments = ledger_module.score_segments

    def spy(label: str, clock: str) -> int:
        calls.append((label, clock))
        pairs.add((row_clock[-1], clock))
        return stream_for(label, clock)

    def scoring(*args, **kwargs):
        row_clock.append(kwargs["clock"])
        return score_segments(*args, **kwargs)

    monkeypatch.setattr(ledger_module, "stream_for", spy)
    monkeypatch.setattr(ledger_module, "score_segments", scoring)

    # Varied durations: a fixed-period carve is one degenerate segment, which the ledger
    # drops before it draws a stream.
    ledger_module.run(
        ledger_module.CheckLedgerInputs(
            windows=_varied_carve(),
            bench_size_table=_bench_table(),
            thresholds=[("3 µs", 3.0, True)],
            n_permutations=49,
            seed=7,
            include_c3=False,
        )
    )

    assert calls, "check_ledger.run no longer routes its rng stream through stream_for"
    assert pairs == {("side", "in_spec"), ("calendar", "calendar")}
    assert {label for label, _ in calls} == {"3 µs"}


def test_a_ledger_holds_each_bag_clock_and_check_once_at_every_threshold():
    """Oracle: the ledger's contract. Every (threshold, bag, clock) cell holds each check's
    row exactly once, for both bags of the carve; a record that failed reads `not computed`
    on both bags and both clocks at every threshold, under its own dataset id."""
    import quebra.analyzers.check_ledger as ledger_module
    from quebra.analyzers.checks.result import CLOCKS
    from quebra.analyzers.event_table import SIDES

    ladder = [("3 µs", 3.0, True), ("4 µs", 4.0, True)]
    rows = ledger_module.run(
        ledger_module.make_inputs_from_windows(
            _varied_result(ladder),
            _bench_table(),
            thresholds=ladder,
            alpha=0.05,
            min_events_pass=35,
            tie_cutoff_distinct=5,
            lag_max=5,
            n_permutations=19,
            seed=7,
        )
    ).rows
    cells = rows.groupby(["threshold_label", "side", "clock", "check_id"]).size()
    assert (cells == 1).all(), cells[cells != 1]
    expected = set(product(("3 µs", "4 µs"), SIDES, CLOCKS))
    assert set(cells.index.droplevel("check_id")) == expected
    assert "(all)" not in set(rows["check_id"]), "premise: every cell was scored"
    # Both bags read the bench, on both clocks: the side clock's cells are keyed `in_spec`.
    benched = rows.groupby(["side", "clock"])["bench_size_cell"].apply(
        lambda cells: bool((cells != "").all())
    )
    assert benched.all(), benched

    failed = ledger_module.run(
        ledger_module.CheckLedgerInputs(
            windows=pd.DataFrame({"threshold_label": []}),
            bench_size_table=_bench_table(),
            thresholds=ladder,
            dataset_id="broken",
            record_failure="load failed: E",
        )
    ).rows
    assert sorted(zip(failed["threshold_label"], failed["side"], failed["clock"])) == (
        sorted(expected)
    )
    assert set(failed["dataset_id"]) == {"broken"}
    assert set(failed["verdict"]) == {VERDICT_NOT_COMPUTED}


def test_tie_stats_counts_what_it_says():
    from quebra.analyzers.check_ledger import _tie_stats

    distinct, tied = _tie_stats(np.array([1.0, 1.0, 2.0, 3.0]))
    assert distinct == 3
    assert tied == pytest.approx(0.5)
    assert _tie_stats(np.array([]))[0] == 0


def test_ledger_refuses_a_partial_artifact():
    """`check_thresholds` is the construction-time completeness contract."""
    from quebra.analyzers.check_ledger import CheckLedger

    ledger = CheckLedger(rows=pd.DataFrame({"threshold_label": ["3 µs"]}))
    ledger.check_thresholds(["3 µs"])
    with pytest.raises(ValueError, match="incomplete CheckLedger"):
        ledger.check_thresholds(["3 µs", "4 µs"])


def _varied_carve() -> pd.DataFrame:
    """A carve with ten distinct durations, above the tie cutoff of five.

    On the side clock its segments end at an event (`tau == T_N`), so C1, C2 and CvM read
    their first n - 1 events there.

    A metric that alternates on a fixed period gives one repeated duration, and the battery
    declines a constant rank vector, so a fixture built that way would report a blank row
    and could not tell a lost check from a degenerate cell.
    """
    return _varied_result().windows_in_spec


def _varied_result(thresholds=(("3 µs", 3.0, True),)) -> windows.WindowsResult:
    """`_varied_carve`'s whole carve: both bags, on `thresholds`."""
    rng = np.random.default_rng(3)
    values: list[float] = []
    for _ in range(60):
        values += [5.0] * int(rng.integers(2, 12))
        values += [1.0] * int(rng.integers(2, 8))
    series = np.asarray(values, dtype=float)
    return windows.run(
        windows.WindowsInputs(
            t_rel_s=np.arange(len(series), dtype=float) * 60.0,
            values=series,
            thresholds=list(thresholds),
            dataset_id="unit",
        )
    )


def _bench_table() -> pd.DataFrame:
    # `Path(__file__)`, not `repo_root()`: the latter is defined by the working directory.
    return pd.read_csv(
        pathlib.Path(__file__).resolve().parents[1]
        / "jobs"
        / "bench"
        / "results"
        / "size_table.csv"
    )


def test_a_ledger_that_never_asked_for_c3_records_not_asked_and_never_probes(
    monkeypatch,
):
    """Oracle: AGENTS.md section 3 - a silent fallback yields a wrong-but-plausible result.

    A run-set without C3 builds no C3 row, so an R version on that artifact would describe
    work that did not happen. `jobs/active/km_with_checks_6d2s.py` ships exactly this shape:
    its run-set is the permutation keys, which set `include_c3=False`. The probe must not
    run at all there, which is asserted by making it raise if anything calls it.
    """
    import quebra.analyzers.check_ledger as ledger_module
    from quebra.analyzers.checks import c3_serial_copula as c3

    def _forbidden() -> str:
        raise AssertionError("rscript_path() was called for a ledger that excluded C3")

    monkeypatch.setattr(c3, "rscript_path", _forbidden)

    ledger = ledger_module.run(
        ledger_module.CheckLedgerInputs(
            windows=_varied_carve(),
            bench_size_table=_bench_table(),
            thresholds=[("3 \u00b5s", 3.0, True)],
            n_permutations=19,
            seed=7,
            include_c3=False,
        )
    )

    assert ledger.r_version == "not asked"
    assert ledger.r_executable == "not asked"
    assert ledger.r_library_paths == ()
    assert not (ledger.rows["check_id"] == c3.CHECK_NAME).any(), (
        "a ledger with include_c3=False must carry no C3 row"
    )


def test_with_no_rscript_the_r_fields_read_absent_and_every_row_is_still_built(
    monkeypatch, tmp_path
):
    """Oracle: AGENTS.md section 3 - a check outcome never stops execution.

    PATH is emptied rather than `rscript_path` patched, so this exercises the same
    `shutil.which` miss a machine without R would produce. It must NOT carry the `r`
    marker: it asserts the R-ABSENT path and has to run where there is no R.
    """
    import quebra.analyzers.check_ledger as ledger_module
    from quebra.analyzers.checks import c3_serial_copula as c3

    monkeypatch.setenv("PATH", str(tmp_path))
    assert c3.rscript_path() is None, "the fixture did not actually hide Rscript"

    ledger = ledger_module.run(
        ledger_module.CheckLedgerInputs(
            windows=_varied_carve(),
            bench_size_table=_bench_table(),
            thresholds=[("3 µs", 3.0, True)],
            n_permutations=19,
            seed=7,
            include_c3=True,
        )
    )

    assert ledger.r_version == "absent"
    assert ledger.r_executable == "absent"
    assert ledger.r_library_paths == ()

    in_spec = ledger.rows[
        (ledger.rows["side"] == SIDE_IN_SPEC) & (ledger.rows["clock"] == CLOCK_SIDE)
    ]
    built = set(zip(in_spec["check_id"], in_spec["calibration"], in_spec["variant"]))
    expected = set(battery.ROW_KEYS) | {(c3.CHECK_NAME, CALIB_R_COPULA, "")}
    # The side clock emits every battery row and C3's. Asserting the SET rather than a
    # count is what makes a dropped row name itself instead of moving a number.
    assert built == expected, (
        f"R absence changed which rows the in-spec cell carries: {built ^ expected}"
    )
    c3_rows = ledger.rows[ledger.rows["check_id"] == c3.CHECK_NAME]
    assert set(c3_rows["clock"]) == {CLOCK_SIDE, CLOCK_CALENDAR}
    assert set(c3_rows["verdict"]) == {VERDICT_NOT_COMPUTED}
    assert set(c3_rows["notes"]) == {"R unavailable"}


def test_c3_n_null_sim_reaches_the_artifact_as_the_caller_set_it(monkeypatch, tmp_path):
    """The simulation count a C3 p-value was drawn against survives materialization.

    Driven with a value that is neither the module default nor any job's, so a field that
    silently fell back to `c3.N_NULL_SIM` would fail here. R is hidden: this pins the
    bookkeeping, and running C3 for real is not needed to do it.
    """
    import pickle

    import quebra.analyzers.check_ledger as ledger_module
    from quebra.analyzers.checks import c3_serial_copula as c3

    monkeypatch.setenv("PATH", str(tmp_path))
    assert c3.rscript_path() is None, "the fixture did not actually hide Rscript"

    ledger = ledger_module.run(
        ledger_module.CheckLedgerInputs(
            windows=_varied_carve(),
            bench_size_table=_bench_table(),
            thresholds=[("3 µs", 3.0, True)],
            n_permutations=19,
            seed=7,
            c3_n_null_sim=137,
        )
    )
    assert ledger.c3_n_null_sim == 137
    assert c3.N_NULL_SIM != 137, "the fixture must not be the module default"
    # `pickle` is what `job.materialize` writes, and the round trip runs the artifact
    # guard, so this is the field on the materialized object rather than in memory only.
    assert pickle.loads(pickle.dumps(ledger)).c3_n_null_sim == 137


def test_the_recorded_n_is_the_n_handed_to_the_bridge(monkeypatch, tmp_path):
    """Oracle: the kwargs `c3.run` actually receives, captured at the call.

    The artifact field and the bridge argument are two separate reads of
    `inputs.c3_n_null_sim`, so the ledger can record an N it did not use. Pinning only the
    field leaves that undetectable: a call site reverted to the module default stays green.
    """
    import quebra.analyzers.check_ledger as ledger_module
    from quebra.analyzers.checks import c3_serial_copula as c3

    # PATH is emptied for the same reason as its siblings: this is a `unit` test in the
    # fast selection, and one that shelled out to R would assert something different on a machine
    # with R than on one without. The kwargs are captured before dispatch, so the bridge
    # never needs to run for the assertion to hold.
    monkeypatch.setenv("PATH", str(tmp_path))
    assert c3.rscript_path() is None, "the fixture did not actually hide Rscript"

    seen: dict[str, object] = {}
    real_run = c3.run

    def _capture(segments, **kwargs):
        seen.update(kwargs)
        return real_run(segments, **kwargs)

    monkeypatch.setattr(c3, "run", _capture)

    ledger = ledger_module.run(
        ledger_module.CheckLedgerInputs(
            windows=_varied_carve(),
            bench_size_table=_bench_table(),
            thresholds=[("3 \u00b5s", 3.0, True)],
            n_permutations=19,
            seed=7,
            include_c3=True,
            c3_n_null_sim=137,
        )
    )

    assert seen.get("n_null_sim") == 137, (
        f"the bridge was handed {seen.get('n_null_sim')}, not the 137 the inputs carried"
    )
    assert ledger.c3_n_null_sim == 137, "and the artifact must record the same N"


def test_a_version_probe_that_raises_is_recorded_and_does_not_stop_the_ledger(
    monkeypatch, tmp_path
):
    """Oracle: AGENTS.md section 3 - a check outcome never stops execution.

    `r_version()` calls a subprocess, so it can raise where the interpreter hangs
    (`TimeoutExpired`), vanishes mid-call (`OSError`), or emits bytes that are not valid
    UTF-8 (`ValueError`). None is a `RuntimeError`. The probe runs above the threshold loop
    and outside every check's own error handling, so an escape kills the ledger before a
    single row exists. This pins the guard that stops it; without the test the guard can be
    deleted with the suite green.
    """
    import quebra.analyzers.check_ledger as ledger_module
    from quebra.analyzers.checks import c3_serial_copula as c3

    stub = tmp_path / "Rscript"
    stub.write_text("#!/bin/sh\nexit 0\n")
    stub.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    assert c3.rscript_path() is not None, "the stub interpreter must be visible"

    def _raises(*_args, **_kwargs):
        raise OSError("interpreter vanished mid-call")

    monkeypatch.setattr(c3, "r_version", _raises)

    ledger = ledger_module.run(
        ledger_module.CheckLedgerInputs(
            windows=_varied_carve(),
            bench_size_table=_bench_table(),
            thresholds=[("3 \u00b5s", 3.0, True)],
            n_permutations=19,
            seed=7,
            include_c3=True,
        )
    )

    assert ledger.r_version == "probe failed (OSError)", ledger.r_version
    assert ledger.r_executable == str(stub)
    assert len(ledger.rows) > 0, "the ledger must still build every row"


def test_an_unusable_r_is_recorded_on_the_ledger_and_does_not_stop_it(
    monkeypatch, tmp_path
):
    """Oracle: AGENTS.md section 3, again - for R present but not usable.

    The common case is `copula` not installed, and it is not the same case as no R at all:
    `rscript_path` finds an executable, so both probes run and both fail. Stood in for by a
    script that exits non-zero on every invocation, which is what such an R does; no `r`
    marker, because no R is involved.
    """
    import quebra.analyzers.check_ledger as ledger_module
    from quebra.analyzers.checks import c3_serial_copula as c3

    broken = tmp_path / "Rscript"
    broken.write_text(
        f"#!{sys.executable}\n"
        "import sys\n"
        'sys.stderr.write(\'Error in packageVersion("copula") : '
        "there is no package called \\'copula\\'\\nExecution halted\\n')\n"
        "sys.exit(1)\n"
    )
    broken.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))

    ledger = ledger_module.run(
        ledger_module.CheckLedgerInputs(
            windows=_varied_carve(),
            bench_size_table=_bench_table(),
            thresholds=[("3 µs", 3.0, True)],
            n_permutations=19,
            seed=7,
            include_c3=True,
        )
    )

    assert ledger.r_executable == str(broken)
    assert ledger.r_version.startswith("probe failed"), ledger.r_version
    assert "copula" in ledger.r_version, "the record must say what failed"
    assert ledger.r_library_paths == ("probe failed (RuntimeError)",), (
        "a failed library probe must leave a visible sentinel: an empty tuple is also what "
        "`not asked` and `absent` carry, so silence here would hide the failure"
    )
    in_spec = ledger.rows[
        (ledger.rows["side"] == SIDE_IN_SPEC) & (ledger.rows["clock"] == CLOCK_SIDE)
    ]
    built = set(zip(in_spec["check_id"], in_spec["calibration"], in_spec["variant"]))
    expected = set(battery.ROW_KEYS) | {(c3.CHECK_NAME, CALIB_R_COPULA, "")}
    assert built == expected, (
        f"an unusable R changed which rows the in-spec cell carries: {built ^ expected}"
    )
    assert set(ledger.rows[ledger.rows["check_id"] == c3.CHECK_NAME]["verdict"]) == {
        VERDICT_NOT_COMPUTED
    }


# The checks a ledger can run that the bench has never calibrated. DERIVED below from the
# bench tables rather than trusted: if C3 is ever benched, or a new uncalibrated check is
# registered, the test fails and this constant has to move with it.
UNCALIBRATED_CHECKS = {"c3_serial_copula"}


def test_every_ledger_check_is_benched_or_declared_uncalibrated():
    """Oracle: jobs/bench/results/size_table.csv, read rather than assumed.

    A check with no bench cell has no measured size, so neither its rejections nor its
    non-rejections carry a level. The ledger already demotes the non-rejection to
    `inconclusive`; the rejection keeps `fail` and must say so in its note. Nothing pinned
    WHICH checks are in that state, so a later pass could bench C3, or register a seventh
    check with no calibration, and the distinction would quietly stop matching the tables.
    """
    from quebra.analyzers.checks import c3_serial_copula as c3

    size = pd.read_csv(
        pathlib.Path(__file__).resolve().parents[1]
        / "jobs"
        / "bench"
        / "results"
        / "size_table.csv"
    )
    benched = set(size["check"].unique())
    # What a ledger can run: the battery's rows, and C3 when a caller asks for it.
    runnable = {key[0] for key in battery.ROW_KEYS} | {c3.CHECK_NAME}

    assert UNCALIBRATED_CHECKS == runnable - benched, (
        "the ledger's checks with no bench cell have changed. If a check was benched, "
        "remove it from UNCALIBRATED_CHECKS here and check that "
        "analyzers/instrument_validation.py's tier rows moved with it. If a new "
        "uncalibrated check was registered, it needs a row in the bench before its "
        "verdicts can be read as levelled."
    )
    # No second assertion here. `UNCALIBRATED_CHECKS <= runnable` is true by construction
    # once the equality above holds, and "not in ROW_KEYS" is already pinned by
    # tests/test_independence_survey.py, whose subject that is.


def test_an_uncalibrated_rejection_says_no_level_stands_behind_it():
    """Oracle: AGENTS.md section 3 - a check outcome never stops execution, and section 4.

    The twin of the `inconclusive` demotion on the non-rejection path. Without it a
    rejection from a check the bench has never run reads exactly like a calibrated one, and
    the survey draws 60 such cells today.
    """
    from quebra.analyzers.check_ledger import VERDICT_FAIL, _verdict
    from quebra.analyzers.checks.result import CALIB_R_COPULA, CheckResult

    rejected = CheckResult(
        check="c3_serial_copula",
        statistic=9.9,
        p_value=0.001,
        calibration=CALIB_R_COPULA,
        clock="in_spec",
        n_events=80,
        n_segments=1,
        n_censored_dropped=0,
        notes="",
    )
    verdict, note = _verdict(
        rejected,
        alpha=0.05,
        n_events=80,
        min_events_pass=0,
        n_distinct=80,
        tie_cutoff=5,
        bench_accepted=None,
    )
    assert verdict == VERDICT_FAIL, (
        "the statistic did land in the tail; do not hide that"
    )
    assert "UNCALIBRATED" in note and "no level stands behind" in note, note


def test_score_segments_with_nothing_to_score_writes_not_computed():
    """Oracle: the ledger's verdict vocabulary; an empty segment list is a row, not a raise.

    The shared scoring half is reached by callers whose segments were all dropped
    (`checks.battery.drop_unusable_segments`). It must still write the cell, with its
    reason, so the band beside it is drawn; the reason names both counts once a segment
    was degenerate.
    """
    from quebra.analyzers.check_ledger import score_segments

    rows = score_segments(
        [],
        n_dropped=3,
        n_degenerate=0,
        n_used=4,
        label="cycles",
        clock="calendar",
        dataset_id="unit",
        acceptance=pd.DataFrame(),
        grid=[20, 35],
        alpha=0.05,
        min_events_pass=35,
        tie_cutoff_distinct=5,
        lag_max=5,
        n_permutations=99,
        seed=0,
        min_events_per_segment=2,
        include_tau_checks=True,
        include_c2_asymptotic=False,
        include_c3=False,
        c3_n_null_sim=10,
    )
    assert len(rows) == 1
    (row,) = rows
    assert row["verdict"] == VERDICT_NOT_COMPUTED
    assert row["notes"] == "no segment reached 2 events (3 dropped)"
    assert (row["threshold_label"], row["clock"], row["dataset_id"]) == (
        "cycles",
        "calendar",
        "unit",
    )
    assert row["n_windows_used"] == 4
    (row,) = score_segments(
        [],
        n_dropped=3,
        n_degenerate=2,
        n_used=4,
        label="cycles",
        clock="calendar",
        dataset_id="unit",
        acceptance=pd.DataFrame(),
        grid=[20, 35],
        alpha=0.05,
        min_events_pass=35,
        tie_cutoff_distinct=5,
        lag_max=5,
        n_permutations=99,
        seed=0,
        min_events_per_segment=2,
        include_tau_checks=True,
        include_c2_asymptotic=False,
        include_c3=False,
        c3_n_null_sim=10,
    )
    assert row["verdict"] == VERDICT_NOT_COMPUTED
    assert row["notes"] == (
        "no segment left: 3 below 2 events, 2 with every duration identical"
    )


def test_a_segment_of_identical_durations_is_dropped_and_the_others_answer():
    """Oracle: a carve built so that the middle of three stretches, split by read gaps,
    alternates on a fixed period: every in-spec window in it lasts three reads and every
    gap between its births five, so its durations are identical on both clocks. The ledger
    drops that one segment, names it in every row's notes, and every check answers on the
    other two, with the rows the same carve gives without the middle stretch: same
    segments, same permutations.
    """
    import quebra.analyzers.check_ledger as ledger_module

    rng = np.random.default_rng(5)

    def varied() -> list[float]:
        values: list[float] = []
        for _ in range(20):
            values += [5.0] * int(rng.integers(2, 12))
            values += [1.0] * int(rng.integers(2, 8))
        return values

    stretches = [varied(), [1.0] + ([5.0] * 3 + [1.0] * 2) * 10, varied()]

    def ledger_rows(kept: list[int]) -> pd.DataFrame:
        t_rel_s: list[float] = []
        values: list[float] = []
        for index in kept:
            t_rel_s += [index * 1.0e6 + 60.0 * i for i in range(len(stretches[index]))]
            values += stretches[index]
        carved = windows.run(
            windows.WindowsInputs(
                t_rel_s=np.asarray(t_rel_s),
                values=np.asarray(values),
                thresholds=[("3 µs", 3.0, True)],
                dataset_id="unit",
            )
        )
        gaps = carved.diagnostics["gap_spans_s"]
        assert len(gaps) == len(kept) - 1, "premise: a read gap between stretches"
        return ledger_module.run(
            ledger_module.CheckLedgerInputs(
                windows=carved.windows_in_spec,
                bench_size_table=_bench_table(),
                thresholds=[("3 µs", 3.0, True)],
                dataset_id="unit",
                gap_spans_s=gaps,
                n_permutations=99,
                seed=7,
            )
        ).rows.pipe(lambda rows: rows[rows["side"] == SIDE_IN_SPEC])

    rows = ledger_rows([0, 1, 2])
    assert "(all)" not in set(rows["check_id"])
    columns = ["clock", "check_id", "calibration", "statistic", "p_value", "verdict"]
    pd.testing.assert_frame_equal(rows[columns], ledger_rows([0, 2])[columns])
    for clock in (CLOCK_SIDE, CLOCK_CALENDAR):
        cell = rows[rows["clock"] == clock]
        assert set(cell["n_segments"]) == {2}, clock
        assert (
            cell["notes"]
            .str.contains(
                "1 segment(s) with every duration identical dropped", regex=False
            )
            .all()
        ), clock
        answering = [key[0] for key in battery.ROW_KEYS]
        answered = cell[cell["check_id"].isin(answering)]
        assert len(answered) == len(answering), clock
        assert answered["p_value"].notna().all(), answered[["check_id", "notes"]]


def test_cycle_shaped_segments_ending_at_an_event_answer_and_count_what_they_read():
    """Oracle: the failure-truncation rule (Rigdon and Basu 2000, Sect. 4.3) on the scoring
    half the cycles path calls (`cycles.cycle_checks` -> `score_segments`). The segments are shaped as
    `cycles.battery_segments` builds them: tau is the cumulative sum of the complete cycles,
    plus a censored tail only where one exists. Two of three end at an event, so C1, C2
    and CvM answer, their rows count 43 - 2 = 41 events, look the bench up and judge the
    floor by that, and say how many segments ended at an event; C5 and C6 count all 43."""
    from quebra.analyzers.calibration_summary import (
        bench_acceptance_at_n,
        nearest_bracketing_n,
    )
    from quebra.analyzers.check_ledger import EVENT_TIME_CHECKS, score_segments
    from quebra.analyzers.checks.result import Segment

    rng = np.random.default_rng(11)
    segments = []
    for size, tail in ((14, 0.0), (17, 0.0), (12, 2.5)):
        x = rng.weibull(1.5, size=size) * 10.0
        tau = float(np.cumsum(x)[-1]) + tail
        segments.append(Segment(x=x, tau=tau, n_censored_dropped=int(tail > 0)))
    assert [s.ends_at_event for s in segments] == [True, True, False]
    acceptance = bench_acceptance_at_n(_bench_table(), alpha=0.05)
    grid = sorted({int(n) for n in acceptance["n_target"]})
    rows = pd.DataFrame(
        score_segments(
            segments,
            n_dropped=0,
            n_degenerate=0,
            n_used=46,
            label="cycles",
            clock=CLOCK_CALENDAR,
            dataset_id="unit",
            acceptance=acceptance,
            grid=grid,
            alpha=0.05,
            min_events_pass=35,
            tie_cutoff_distinct=5,
            lag_max=5,
            n_permutations=99,
            seed=0,
            min_events_per_segment=2,
            include_tau_checks=True,
            include_c2_asymptotic=False,
            include_c3=False,
            c3_n_null_sim=10,
        )
    ).set_index("check_id")
    for check in EVENT_TIME_CHECKS:
        row = rows.loc[check]
        assert row["verdict"] != VERDICT_NOT_COMPUTED, (check, row["notes"])
        assert row["n_events"] == 41, check
        assert row["bench_n_target"] == nearest_bracketing_n(grid, 41), check
        assert (
            "2 of 3 segment(s) end at an event: first n - 1 events read"
            in (row["notes"])
        ), check
    for check in ("c5_rank_autocorr", "c6_exchangeability"):
        assert rows.loc[check, "n_events"] == 43, check
        assert "end at an event" not in rows.loc[check, "notes"], check


def test_an_event_ended_row_meets_the_floor_and_reads_the_bench_by_the_events_it_read():
    """Oracle: the verdict rule at the 35-event floor, on the events read. One segment of 35
    gaps that ends at an event reads 34, so C1, C2 and CvM do not reach the floor and a
    non-rejection reads `inconclusive` for it. With CvM's bench cells taken out, C1's and
    C2's notes say their cell was measured without the n - 1 form and CvM's says nothing of
    a cell it did not read."""
    from quebra.analyzers.calibration_summary import bench_acceptance_at_n
    from quebra.analyzers.check_ledger import EVENT_TIME_CHECKS, score_segments
    from quebra.analyzers.checks import cvm_cramer_von_mises as cvm
    from quebra.analyzers.checks.result import Segment

    x = np.random.default_rng(5).exponential(10.0, size=35)
    segment = Segment(x=x, tau=float(np.cumsum(x)[-1]))
    assert segment.ends_at_event and segment.n_events_read == 34
    acceptance = bench_acceptance_at_n(_bench_table(), alpha=0.05)
    grid = sorted({int(n) for n in acceptance["n_target"]})
    rows = pd.DataFrame(
        score_segments(
            [segment],
            n_dropped=0,
            n_degenerate=0,
            n_used=35,
            label="3 µs",
            clock=CLOCK_SIDE,
            dataset_id="unit",
            acceptance=acceptance[acceptance["check"] != cvm.CHECK_NAME],
            grid=grid,
            alpha=0.05,
            min_events_pass=35,
            tie_cutoff_distinct=5,
            lag_max=5,
            n_permutations=99,
            seed=0,
            min_events_per_segment=2,
            include_tau_checks=True,
            include_c2_asymptotic=False,
            include_c3=False,
            c3_n_null_sim=10,
        )
    ).set_index("check_id")
    clause = "the bench cell was measured without the n - 1 form"
    for check in EVENT_TIME_CHECKS:
        row = rows.loc[check]
        assert row["p_value"] > 0.05, f"premise: {check} does not reject"
        assert row["n_events"] == 34, check
        assert row["verdict"] == VERDICT_INCONCLUSIVE, (check, row["notes"])
        assert "n_events 34 < 35" in row["notes"], (check, row["notes"])
        assert (clause in row["notes"]) == (check != cvm.CHECK_NAME), (
            check,
            row["notes"],
        )
    assert "no bench cell" in rows.loc[cvm.CHECK_NAME, "notes"]


def test_c1_keeps_its_answer_when_c2_cannot_give_one(monkeypatch):
    """Oracle: the issue's rule ("the retry no longer drops C1"). With C2 forced to raise,
    its row reads `not computed` with its own reason in every cell of both bags and clocks,
    and C1 and every other check still answer."""
    import quebra.analyzers.check_ledger as ledger_module
    import quebra.analyzers.checks.c2_anderson_darling as c2

    def refuse(*_args, **_kwargs):
        raise ValueError("forced C2 refusal")

    monkeypatch.setattr(c2, "statistic", refuse)
    rows = ledger_module.run(
        ledger_module.make_inputs_from_windows(
            _varied_result(),
            _bench_table(),
            thresholds=[("3 µs", 3.0, True)],
            alpha=0.05,
            min_events_pass=35,
            tie_cutoff_distinct=5,
            lag_max=5,
            n_permutations=19,
            seed=7,
        )
    ).rows
    assert len(rows.groupby(["side", "clock"])) == 4, "premise: both bags, both clocks"
    for (side, clock), cell in rows.groupby(["side", "clock"]):
        by_check = cell.set_index("check_id")
        c2_row = by_check.loc["c2_anderson_darling"]
        assert c2_row["verdict"] == VERDICT_NOT_COMPUTED, (side, clock)
        assert "forced C2 refusal" in c2_row["notes"], (side, clock)
        others = by_check.drop(index="c2_anderson_darling")
        assert others["p_value"].notna().all(), (side, clock, others["notes"].to_dict())


def test_the_p_value_ladder_draws_cvm_from_its_permutation_row():
    """Oracle: the ledger's own CvM row. A ledger carries one CvM row per cell, its
    permutation row, so the in-spec bag's figure draws a CvM line through that row's
    p-value on the side clock."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    import quebra.analyzers.check_ledger as ledger_module
    from quebra.ledger_recipe import ledger_bag
    from quebra.panels.check_ledger import CheckLedgerPanel

    ledger = ledger_module.run(
        ledger_module.CheckLedgerInputs(
            windows=_varied_carve(),
            bench_size_table=_bench_table(),
            thresholds=[("3 µs", 3.0, True)],
            n_permutations=19,
            seed=7,
        )
    )
    rows = ledger.rows
    cvm = rows[
        (rows["check_id"] == "cvm_cramer_von_mises")
        & (rows["side"] == SIDE_IN_SPEC)
        & (rows["clock"] == CLOCK_SIDE)
    ]
    assert list(cvm["calibration"]) == [CALIB_PERMUTATION]
    assert cvm["p_value"].notna().all(), "premise: CvM answers on the side clock"
    fig = CheckLedgerPanel(name="ladder_probe").build_matplotlib(
        ledger_bag(ledger, side=SIDE_IN_SPEC)
    )
    try:
        (ladder,) = [ax for ax in fig.axes if ax.get_ylabel() == "p-value"]
        lines = [line for line in ladder.get_lines() if line.get_label() == "CvM"]
        assert len(lines) == 1, [line.get_label() for line in ladder.get_lines()]
        assert list(lines[0].get_ydata()) == list(cvm["p_value"].astype(float))
    finally:
        plt.close(fig)


def test_the_ledger_figure_refuses_rows_of_both_bags():
    """Oracle: the panel's contract, one bag per figure. A whole ledger holds both bags,
    and drawing it would put two bags' rows in one table cell; `ledger_bag` selects one."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    import quebra.analyzers.check_ledger as ledger_module
    from quebra.ledger_recipe import ledger_bag
    from quebra.panels.check_ledger import CheckLedgerPanel

    ledger = ledger_module.run(
        ledger_module.make_inputs_from_windows(
            _varied_result(),
            _bench_table(),
            thresholds=[("3 µs", 3.0, True)],
            alpha=0.05,
            min_events_pass=35,
            tie_cutoff_distinct=5,
            lag_max=5,
            n_permutations=19,
            seed=7,
        )
    )
    with pytest.raises(ValueError, match="draws one bag"):
        CheckLedgerPanel(name="both_probe").build_matplotlib(ledger)
    plt.close(
        CheckLedgerPanel(name="one_probe").build_matplotlib(
            ledger_bag(ledger, side=SIDE_IN_SPEC)
        )
    )


def test_the_ledger_figure_is_laid_out_by_tight_layout():
    """Oracle: matplotlib's tight_layout, which warns when it skips an axes, as it skips
    every axes of a grid whose spacing is set locally. The figure holds both tables and the
    ladder with its twin event-count axis."""
    import warnings

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    import quebra.analyzers.check_ledger as ledger_module
    from quebra.ledger_recipe import ledger_bag
    from quebra.panels.check_ledger import CheckLedgerPanel

    ledger = ledger_module.run(
        ledger_module.CheckLedgerInputs(
            windows=_varied_carve(),
            bench_size_table=_bench_table(),
            thresholds=[("3 µs", 3.0, True)],
            n_permutations=19,
            seed=7,
        )
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fig = CheckLedgerPanel(name="layout_probe").build_matplotlib(
            ledger_bag(ledger, side=SIDE_IN_SPEC)
        )
        fig.canvas.draw()
    plt.close(fig)
    labels = [ax.get_ylabel() for ax in fig.axes]
    assert labels.count("p-value") == labels.count("Events (windows)") == 1, labels
    assert not [w for w in caught if "layout" in str(w.message)], [
        str(w.message) for w in caught
    ]
