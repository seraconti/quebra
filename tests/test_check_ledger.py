"""The ledger's verdict rule is the decision this whole increment exists to make.

A p-value alone must never produce a `pass`. These tests pin each of the three conditions
independently, so a regression that drops one of them fails here rather than in a figure.

The file also pins what the ledger records about the R that produced its C3 rows, and that
it still builds every row when R is absent, unusable, or was never asked for.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

from quebra.analyzers import windows
from quebra.analyzers.check_ledger import (
    VERDICT_FAIL,
    VERDICT_NOT_COMPUTED,
    VERDICT_PASS,
    VERDICT_TIES,
    VERDICT_UNDERPOWERED,
    _verdict,
    stream_for,
)
from quebra.analyzers.checks import battery
from quebra.analyzers.checks.result import (
    CALIB_PERMUTATION,
    CALIB_R_COPULA,
    CheckResult,
)

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
    verdict, reason = _verdict(_result(), **BASE)
    assert verdict == VERDICT_PASS and reason == ""


@pytest.mark.parametrize(
    "override, expected_fragment",
    [
        ({"n_events": 20}, "n_events 20 < 35"),
        ({"bench_accepted": False}, "miscalibrated"),
        ({"bench_accepted": None}, "no bench cell"),
    ],
)
def test_each_condition_alone_blocks_a_pass(override, expected_fragment):
    """A non-rejection with any one condition unmet is `underpowered`, never `pass`."""
    verdict, reason = _verdict(_result(), **{**BASE, **override})
    assert verdict == VERDICT_UNDERPOWERED
    assert expected_fragment in reason


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
    """
    import quebra.analyzers.check_ledger as ledger_module

    calls: list[tuple[str, str]] = []

    def spy(label: str, clock: str) -> int:
        calls.append((label, clock))
        return stream_for(label, clock)

    monkeypatch.setattr(ledger_module, "stream_for", spy)

    t = np.arange(80, dtype=float) * 60.0
    values = np.where((np.arange(80) // 7) % 2 == 0, 5.0, 1.0)
    carved = windows.run(
        windows.WindowsInputs(
            t_rel_s=t,
            values=values,
            thresholds=[("3 µs", 3.0, True)],
            dataset_id="unit",
        )
    )
    # `Path(__file__)`, not `repo_root()`: the latter is defined by the working directory,
    # so it breaks when the suite runs from elsewhere. Same idiom as test_checks_cvm.py.
    bench_csv = (
        pathlib.Path(__file__).resolve().parents[1]
        / "jobs"
        / "bench"
        / "results"
        / "size_table.csv"
    )
    bench = pd.read_csv(bench_csv)
    ledger_module.run(
        ledger_module.CheckLedgerInputs(
            windows=carved.windows,
            bench_size_table=bench,
            thresholds=[("3 µs", 3.0, True)],
            n_permutations=49,
            seed=7,
            include_c3=False,
        )
    )

    assert calls, "check_ledger.run no longer routes its rng stream through stream_for"
    assert {clock for _, clock in calls} == {"in_spec", "calendar"}
    assert {label for label, _ in calls} == {"3 µs"}


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

    It does NOT produce every battery row. On the in-spec clock `tau == T_N`, so C1 and C2
    are singular there and read `not computed`, and CvM asymptotic is not re-emitted. Tests
    asserting completeness must name those exclusions rather than pin a count.

    A metric that alternates on a fixed period gives one repeated duration, and the battery
    declines a constant rank vector, so a fixture built that way would report a blank row
    and could not tell a lost check from a degenerate cell.
    """
    rng = np.random.default_rng(3)
    values: list[float] = []
    for _ in range(60):
        values += [5.0] * int(rng.integers(2, 12))
        values += [1.0] * int(rng.integers(2, 8))
    series = np.asarray(values, dtype=float)
    carved = windows.run(
        windows.WindowsInputs(
            t_rel_s=np.arange(len(series), dtype=float) * 60.0,
            values=series,
            thresholds=[("3 µs", 3.0, True)],
            dataset_id="unit",
        )
    )
    return carved.windows


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

    in_spec = ledger.rows[ledger.rows["clock"] == "in_spec"]
    built = set(zip(in_spec["check_id"], in_spec["calibration"], in_spec["variant"]))
    expected = set(battery.ROW_KEYS) | {(c3.CHECK_NAME, CALIB_R_COPULA, "")}
    # The in-spec clock re-emits every battery row except CvM asymptotic, which it drops
    # because `tau == T_N` makes it singular there. Asserting the SET rather than a count
    # is what makes a dropped row name itself instead of moving a number.
    assert built == expected - {("cvm_cramer_von_mises", "asymptotic", "")}, (
        f"R absence changed which rows the in-spec cell carries: {built ^ expected}"
    )
    c3_rows = ledger.rows[ledger.rows["check_id"] == c3.CHECK_NAME]
    assert set(c3_rows["clock"]) == {"in_spec", "calendar"}
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
    # fast gate, and one that shelled out to R would assert something different on a machine
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
    in_spec = ledger.rows[ledger.rows["clock"] == "in_spec"]
    built = set(zip(in_spec["check_id"], in_spec["calibration"], in_spec["variant"]))
    expected = set(battery.ROW_KEYS) | {(c3.CHECK_NAME, CALIB_R_COPULA, "")}
    assert built == expected - {("cvm_cramer_von_mises", "asymptotic", "")}, (
        f"an unusable R changed which rows the in-spec cell carries: {built ^ expected}"
    )
    assert set(ledger.rows[ledger.rows["check_id"] == c3.CHECK_NAME]["verdict"]) == {
        VERDICT_NOT_COMPUTED
    }
