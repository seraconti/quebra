"""C3 - copula-based serial independence test, over a file bridge to R.

Wraps `copula::serialIndepTest` (Genest & Remillard's empirical-copula test based on the
Mobius decomposition of the serial independence hypothesis). The copula 1.1.7 help page
cites Genest and Remillard, "Tests of independence and randomness based on the empirical
copula process", Test 13:335-369 (2004). UNVERIFIED: the paper was not opened, so the
Mobius-decomposition description carries no page locator. It is here because it tests a
strictly stronger null than C5/C6 - full serial independence at all lags jointly, not just
zero rank autocorrelation - and no maintained Python implementation was located.

**Bridge, not bindings.** CSV out, `Rscript`, CSV back. No `rpy2`: rpy2 pins an ABI against
a specific R build and turns "R is missing" into an import-time failure of the whole
package, whereas a subprocess turns it into a per-call `None`. The pipeline must import
cleanly on a machine without R.

**C3 runs and is UNCALIBRATED.** Smoke test on iid exponential input at `lag.max=5`,
`n.sim=1000`, `seed=1`: n=50 gives statistic 0.0057908 / p 0.9575, n=150 gives 0.0071326 /
p 0.9036, n=355 gives 0.0076331 / p 0.8656.

THE DATA RECIPE IS PART OF THE NUMBER. Those three series come from ONE
`numpy.random.default_rng(0)` drawing `exponential(1.0)` sequentially in the order n=50,
then n=150, then n=355 - not re-seeded per n, and not prefixes of one long draw. A statistic
without its recipe cannot be reproduced, which is the defect this paragraph exists to avoid.
SO IS THE SEED, for a different reason: the R side simulates its own null, so the statistic
is seed-invariant and reproduces exactly while the p-value does not - at the module default
`seed=0` the first two series give 0.9476 and 0.9116.

Failing to reject data that satisfies the null is the expected outcome and is a SMOKE TEST,
NOT calibration: C3 still has no size or power evidence, because the bench never ran it and
`battery.ROW_KEYS` has no C3 row.

When R is missing, every entry point below returns `p_value=None` with
`notes="R unavailable"` rather than raising - a missing optional interpreter is a fact about
the environment, not a defect in the data, and the other five checks must still run. The
promotion report scores five checks, not six. Do not read a `not computed` row as a pass.

**No multi-process form.** `serialIndepTest` takes one series. Concatenating segments would
manufacture lag pairs across read gaps - the exact relation the carve refuses to assert -
and combining per-segment p-values needs a dependence-free combination rule that is not
established for this statistic. So `m > 1` returns `None` with that reason recorded rather
than a number whose meaning nobody can defend.
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from quebra.analyzers.checks._rank_serial import MAX_LAG_CAP
from quebra.analyzers.checks.result import (
    CALIB_R_COPULA,
    CLOCK_IN_SPEC,
    CheckResult,
    Segment,
    concatenated_gaps,
    segment_sizes,
    validate_segment,
)

CHECK_NAME = "c3_serial_copula"

R_SCRIPT = Path(__file__).with_name("c3_serial_copula.R")

# The R side simulates its own null distribution; this is how many draws it uses. Kept
# here rather than in the R file so the Python result can record it.
#
# It is the DOMINANT cost: the simulation is what makes C3 grow steeply with n, so a survey
# over hundreds of cells can trade it down. `run` takes it as a parameter for that reason -
# at N = 200 the p-value resolution is 1/201 = 0.005, still far below any alpha this project
# reads, while the run is several times cheaper. The value used is recorded in `notes` on
# every result, because a C3 p-value means something different at a different N.
N_NULL_SIM = 1000

_UNAVAILABLE = "R unavailable"


def rscript_path() -> str | None:
    """Absolute path to `Rscript`, or None. Every other probe here goes through it."""
    return shutil.which("Rscript")


def r_library_paths() -> tuple[str, ...]:
    """The library paths a NON-vanilla `Rscript` would search, asked of R itself.

    Hard-coding `~/R/library` would bind the tool to one machine. Asking R keeps the answer
    correct wherever it runs, and `--vanilla` below then searches exactly what an ordinary
    `Rscript` would - the isolation `--vanilla` buys is from the user PROFILE, not from the
    installed packages.
    """
    executable = rscript_path()
    if executable is None:
        return ()
    completed = subprocess.run(
        [executable, "-e", "cat(paste(.libPaths(), collapse='\\n'))"],
        capture_output=True,
        text=True,
        errors="replace",
        timeout=60,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"could not ask Rscript for .libPaths(): exited {completed.returncode}\n"
            f"stderr: {completed.stderr.strip()}"
        )
    # Returns `()` for an empty answer rather than a sentinel. The other consumer of this
    # function joins the tuple into `R_LIBS` for the real bridge call, so a diagnostic string
    # here would become a filesystem search path. Telling an empty answer apart from a probe
    # that never ran is the LEDGER's problem, and `check_ledger._r_provenance` does it.
    return tuple(line.strip() for line in completed.stdout.splitlines() if line.strip())


def r_version() -> str | None:
    """The R build and the installed `copula` version, or None when there is no `Rscript`.

    Both halves are one answer. The statistic comes from `copula::serialIndepTest` and the
    p-value from `copula`'s own simulated null, so an R build alone does not identify what
    produced a C3 row. Asked without `--vanilla`, for the same reason `r_library_paths`
    is: the package reported must be the one the bridge will load.

    A non-zero exit is REPORTED in the returned string, not raised. That is the one place
    this differs from `r_library_paths`, and the reason is the caller: this is a provenance
    field built once per ledger, outside any check's own error handling, and its common
    failure is R present with `copula` not installed. That must annotate the record rather
    than stop a ledger whose other five checks are pure Python. The string names the
    failure, so no reader can mistake it for a version.
    """
    executable = rscript_path()
    if executable is None:
        return None
    completed = subprocess.run(
        [
            executable,
            "-e",
            'cat(R.version.string, "copula", as.character(packageVersion("copula")))',
        ],
        capture_output=True,
        text=True,
        errors="replace",
        timeout=60,
        check=False,
    )
    if completed.returncode != 0:
        detail = " ".join(
            (completed.stderr.strip() or completed.stdout.strip()).split()
        )
        return f"probe failed (exit {completed.returncode}): {detail}"[:200]
    reported = " ".join(completed.stdout.split())
    # An empty string is the CheckLedger field default, so a probe that ran and said
    # nothing would be indistinguishable from one that never ran.
    if not reported:
        return "probe returned nothing"
    return reported


def _unavailable_result(
    segments: list[Segment], clock: str, note: str, statistic: float = float("nan")
) -> CheckResult:
    return CheckResult(
        check=CHECK_NAME,
        statistic=statistic,
        p_value=None,
        calibration=CALIB_R_COPULA,
        clock=clock,
        n_events=int(sum(len(s.x) for s in segments)),
        n_segments=len(segments),
        n_censored_dropped=int(sum(s.n_censored_dropped for s in segments)),
        notes=note,
    )


def _invoke_rscript(
    x: np.ndarray, max_lag: int, seed: int, timeout_s: float, n_null_sim: int
) -> tuple[float, float]:
    """Run the R script on `x` and return `(statistic, p_value)`.

    It is written to fail loudly rather than plausibly - a non-zero exit, a missing output
    file or an unparseable row all raise, so a run either produces a real number or a
    traceback, never a silently wrong p-value.

    **`--vanilla` needs `R_LIBS` handed to it explicitly.** `--vanilla` implies
    `--no-environ`, which drops `R_LIBS_USER` and therefore the user library from
    `.libPaths()`: measured here, plain `Rscript` sees `/home/sera/R/library` and
    `Rscript --vanilla` does not, so `copula` was invisible and every call raised. Keeping
    `--vanilla` is right for a reproducibility tool - it refuses the user profile and any
    side effect hiding in it - so the fix is to pass the search path as an explicit
    environment variable rather than to drop the flag. The set searched is recorded as well
    as made explicit: `check_ledger.CheckLedger.r_library_paths` carries it on the
    materialized artifact, so a run against a different library set is distinguishable
    after the fact.
    """
    executable = rscript_path()
    if executable is None:
        raise RuntimeError("_invoke_rscript called without Rscript on PATH")
    env = dict(os.environ)
    env["R_LIBS"] = os.pathsep.join(r_library_paths())
    with tempfile.TemporaryDirectory(prefix="qre_c3_") as workdir:
        work = Path(workdir)
        in_path = work / "durations.csv"
        out_path = work / "result.csv"
        with in_path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["duration_s"])
            writer.writerows([[repr(float(v))] for v in x])
        completed = subprocess.run(
            [
                executable,
                "--vanilla",
                str(R_SCRIPT),
                str(in_path),
                str(out_path),
                str(int(max_lag)),
                str(int(n_null_sim)),
                str(int(seed)),
            ],
            capture_output=True,
            text=True,
            timeout=timeout_s,
            check=False,
            env=env,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"Rscript exited {completed.returncode} for {CHECK_NAME}.\n"
                f"stdout: {completed.stdout.strip()}\nstderr: {completed.stderr.strip()}"
            )
        if not out_path.exists():
            raise RuntimeError(
                f"Rscript exited 0 but wrote no result file for {CHECK_NAME}. "
                f"stdout: {completed.stdout.strip()}"
            )
        with out_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise RuntimeError(f"expected exactly one result row, got {len(rows)}")
    row = rows[0]
    try:
        return float(row["statistic"]), float(row["p_value"])
    except (KeyError, TypeError, ValueError) as exc:
        raise RuntimeError(f"unparseable result row from R: {row!r}") from exc


def run(
    segments: list[Segment],
    *,
    clock: str = CLOCK_IN_SPEC,
    max_lag: int = MAX_LAG_CAP,
    seed: int = 0,
    timeout_s: float = 900.0,
    n_null_sim: int = N_NULL_SIM,
) -> CheckResult:
    """Serial independence via `copula::serialIndepTest`, over the R file bridge.

    `timeout_s` was 120.0 and that was too small to be honest. Wall clock at
    `N_NULL_SIM = 1000` on iid exponential input, `seed=1`, under R 4.5.3 with `copula`
    1.1.7 on x86_64: SECONDS at n = 50, UNDER TEN SECONDS at n = 150, and OF ORDER A MINUTE
    at n = 355. No interval is quoted, and that is deliberate: the cost depends on machine
    load, which no run here controlled, and every interval this docstring has quoted was
    falsified by the next run.

    Whether that crosses a 120 s ceiling is a property of the MACHINE and its load, not of
    the statistic. The measurement that set this default recorded 130.2 s at n = 355 and did
    cross it, and such a row silently became `declined:` once the ledger stopped raising -
    a timeout wearing the clothes of a check outcome.

    NO GROWTH EXPONENT IS CLAIMED. An earlier draft of this docstring said "roughly n^2.8";
    that number is not in the data, and three points do not determine a power law here. The
    two segments disagree: n = 50 to 150 roughly triples the cost over a factor 3 in n,
    while n = 150 to 355 multiplies it by roughly eight over a factor 2.4. The cost also
    depends on the machine, on its load, and on `N_NULL_SIM`. Measure it for the n you
    actually have.

    900 s covers every n the T2* ladder has produced so far, and the evidence above n = 355
    is contradictory rather than thin. One unverified observation records n = 682 still
    running at 580 s. The shipped ledger in
    `output/check_ledger_q1_070423_4fc898_20260813_170015/` records a COMPLETION at the same
    n: C3 at n_events = 682 on the in-spec clock, statistic 1.511, p 0.0005, with no row
    timing out. The two are compatible rather than contradictory, and neither pins the
    margin at larger n. It still fails fast rather than hanging a job forever. A row that
    exceeds it becomes `not computed` with the timeout recorded, which is the correct
    outcome; raise `timeout_s` at the call site if that row is wanted.
    """
    if not segments:
        raise ValueError("C3 needs at least one segment")
    for segment in segments:
        validate_segment(segment, require_strict_tau=False)

    if rscript_path() is None:
        return _unavailable_result(segments, clock, _UNAVAILABLE)
    if len(segments) > 1:
        return _unavailable_result(
            segments,
            clock,
            f"skipped: serialIndepTest has no multi-process form (m={len(segments)})",
        )

    x = concatenated_gaps(segments)
    n = int(sum(segment_sizes(segments)))
    if n <= max_lag + 1:
        return _unavailable_result(
            segments, clock, f"skipped: n={n} too short for lag.max={max_lag}"
        )

    statistic, p_value = _invoke_rscript(x, max_lag, seed, timeout_s, n_null_sim)
    return CheckResult(
        check=CHECK_NAME,
        statistic=statistic,
        p_value=p_value,
        calibration=CALIB_R_COPULA,
        clock=clock,
        n_events=n,
        n_segments=1,
        n_censored_dropped=int(sum(s.n_censored_dropped for s in segments)),
        notes=f"lag.max={max_lag} N={n_null_sim} seed={seed}",
    )
