"""Is the independence assumption safe on this device, or only on the record we looked at?

Runs the check ledger over EVERY 6D2S dataset and reshapes the result into one grid per
instrument - datasets down, the T2* threshold ladder across, both clocks stacked. One
figure per entry of `independence_survey.SURVEY_KEYS`, which is `battery.ROW_KEYS` plus C3.

**Why this exists.** `jobs/composite/check_ledger_q1.py` answers the question for one qubit
on two dates. A single record rejecting iid could be that record; thirty-four records
rejecting it is a property of the instrument or the device. Neither figure can be read off
the other.

**In jobs/composite/, not jobs/active/.** `main.py` globs `jobs/active` for `run --all`, and
this job carves thirty-four datasets. Same reason the other composites live here, even
though this one loads its datasets directly rather than through `job.include`: there is no
per-dataset job to include for thirty-two of them.

**C3 is ON, and is the one instrument here with no bench cell.** Every other grid can be
read against a measured size; C3's cannot, so a green C3 cell means "did not reject" and
nothing stronger. It is included anyway because a survey that omits an instrument answers a
different question from the one it appears to answer - but the figure says what it lacks.

**CvM is included.** It is in `battery.ROW_KEYS`, so it arrives
through `run_battery` like the rest. It earns its place on the IN-SPEC clock: `tau == T_N`
silences C1 and C2 there, and CvM's integrand carries no `1/(s(1-s))` weight, so its
permutation row is defined where theirs are not.

**The carve is repeated here, and that is a real cost.** The active T2* jobs wire the same
chain; `check_ledger_q1.py` avoids duplicating it by `include`-ing them. That is not
available at thirty-four datasets, so the parameters below must stay equal to those jobs'
or the survey and the panels describe different windows. `tests/test_independence_survey.py`
asserts that equality rather than trusting this comment.

Run:  quebra run jobs/composite/independence_survey.py
"""

from __future__ import annotations

import re

import pandas as pd

import quebra.analyzers.check_ledger as check_ledger
import quebra.analyzers.t2star as t2star
import quebra.analyzers.windows as windows
from quebra.analyzers.independence_survey import build_independence_survey
from quebra.analyzers.t2star import T2StarResult
from quebra.analyzers.windows import WindowsResult
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.recipes import RAMSEY_CONFIG, _filter_step, _final_stage
from quebra.plots.independence_survey_plot import (
    SURVEY_PLOTS,
    IndependenceSurveyOverviewPlot,
)
from quebra.schemas.track912 import track912Schema

# The logical name and category. `include` resolves JOB_ID, so this
# file can move without breaking any composite; recategorising costs one string edit.
JOB_ID = "independence_survey"
JOB_FAMILY = "independence"
# Not swept by a bare `run --all`: it re-runs sub-jobs and/or is long. Selectable
# with `--family independence` or by path. This is the declaration that replaced the old
# "jobs/composite/ is not swept" directory rule.
JOB_SWEEP = False

PREFIX = "independence_survey"

# The same ladder every T2* job uses. Repeated rather than imported so the ladder this
# survey scored is visible in this file and lands on the provenance label.
# `k / 1e6`, NOT `k * 1e-6`. The two differ: 1e-6 is not exactly representable, so
# `5 * 1e-6` is 4.9999999999999996e-06 while the literal `5e-6` the T2* jobs write is
# 5e-06 - two of the ten rungs (5 and 10) land on a different double. Dividing by the
# exactly-representable 1e6 reproduces the literal on all ten.
# The gap is about 4e-22 and no measured value will fall in it, so this has never changed
# a window. It is fixed because the survey CLAIMS to score the same ladder as the panel,
# and a claim that is true only to 15 significant figures is not the claim being made.
THRESHOLDS: list[tuple[str, float, bool]] = [
    (f"{k} µs", k / 1e6, True) for k in range(1, 11)
]

# Carve parameters. MUST equal jobs/active/t2star_*.py - see the module docstring.
GAP_MULT = 10.0
K = 1.0
USE_UNCERTAINTY = True

# Verdict parameters, declared here because choosing them is deciding what counts as
# evidence. Equal to jobs/composite/check_ledger_q1.py so the two ledgers are comparable -
# EXCEPT `SEED`, which is deliberately different and must be: an identical permutation seed
# across two jobs would make their Monte Carlo error identical rather than independent, so a
# rung that happened to sit near alpha would land the same way in both and look like
# corroboration. The seed is not a comparability parameter; the five below are.
ALPHA = 0.05
MIN_EVENTS_PASS = 35
TIE_CUTOFF_DISTINCT = 5
LAG_MAX = 5
N_PERMUTATIONS = 999
SEED = 20260814  # NOT equal to check_ledger_q1's, on purpose - see above.

# C3 IS included. It is the only out-of-process check and by far the most expensive, and
# it is here because a survey that silently omits an instrument answers a different
# question from the one it appears to answer.  The cost is paid down through the R
# side's own simulation count rather than through a timeout. Per call at N = 1000 the
# cost runs from seconds at n = 50 to about a minute at n = 355;
# `analyzers/checks/c3_serial_copula.py` carries the measured spread. Over this survey's
# populated cells (415 at the time this was written, not re-counted since) that SUGGESTS
# a quarter of an hour if every cell were as small as n = 50 and about six hours if
# every cell were as large as n = 355. It is not a bracket: those are the smallest and
# largest n MEASURED, not the smallest and largest cell, and no cell has been checked
# against them. The actual total has not been measured. A per-cell timeout would buy it
# back by discarding the LARGEST cells - exactly the ones with the most events and the
# most to say - so it is the wrong lever. Dropping N to 200 keeps every cell and costs
# resolution instead: the p-value floor becomes 1/201 = 0.005, still an order below the
# 0.05 this ledger reads. The value is recorded on the artifact as
# `CheckLedger.c3_n_null_sim`, not on the row, whose note the verdict reason replaces -
# and it has to be recorded somewhere, because a C3 p-value at N = 200 is not the same
# object as one at N = 1000.
INCLUDE_C3 = True
C3_N_NULL_SIM = 200

# Every 6D2S record, in filename order. Listed explicitly rather than globbed: a glob makes
# the run identity depend on the contents of a directory outside the repo, so adding a file
# would silently change what this figure describes. Probed across all thirty-four
# carve cleanly, total 32 s.
DATASET_FILES: tuple[str, ...] = (
    "040423_6D2S_qubit1.pickle",
    "050423_6D2S_qubit1.pickle",
    "070423_6D2S_qubit1.pickle",
    "070623_6D2S_qubit5.pickle",
    "070623_6D2S_qubit6.pickle",
    "070723_6D2S_qubit4.pickle",
    "070723_6D2S_qubit5.pickle",
    "070723_6D2S_qubit6.pickle",
    "080623_6D2S_qubit5.pickle",
    "080623_6D2S_qubit6.pickle",
    "090623_6D2S_qubit5.pickle",
    "090623_6D2S_qubit6.pickle",
    "100423_6D2S_qubit1.pickle",
    "100723_6D2S_qubit4.pickle",
    "110723_6D2S_qubit3.pickle",
    "120423_6D2S_qubit1.pickle",
    "120423_6D2S_qubit1_after.pickle",
    "120423_6D2S_qubit1_before.pickle",
    "120423_6D2S_qubit5.pickle",
    "120423_6D2S_qubit6.pickle",
    "210423_6D2S_qubit1.pickle",
    "210423_6D2S_qubit2.pickle",
    "210423_6D2S_qubit3.pickle",
    "220423_6D2S_qubit1.pickle",
    "270623_6D2S_qubit2.pickle",
    "270623_6D2S_qubit2_after.pickle",
    "270623_6D2S_qubit2_before.pickle",
    "280623_6D2S_qubit1.pickle",
    "280623_6D2S_qubit2.pickle",
    "280623_6D2S_qubit2_after.pickle",
    "280623_6D2S_qubit2_before.pickle",
    "290623_6D2S_qubit4.pickle",
    "290623_6D2S_qubit5.pickle",
    "290623_6D2S_qubit6.pickle",
)

BENCH_SIZE_TABLE = Dataset(path="jobs/bench/results/size_table.csv", schema=None)


def _label(filename: str) -> str:
    """`070423_6D2S_qubit1.pickle` -> `q1 070423`, with any suffix kept.

    Qubit first so the rows group by device element when sorted, which is the comparison
    a reader makes: the same qubit on different days sits together.
    """
    stem = filename[: -len(".pickle")]
    match = re.match(r"(\d{6})_6D2S_qubit(\d+)(?:_(\w+))?$", stem)
    if not match:
        raise ValueError(f"unrecognised dataset filename: {filename!r}")
    date, qubit, suffix = match.groups()
    return f"q{qubit} {date}" + (f" {suffix}" if suffix else "")


def _qubit(filename: str) -> int:
    return int(re.search(r"qubit(\d+)", filename).group(1))


def _t2star_run(norm: object) -> T2StarResult:
    return t2star.run(t2star.make_inputs_from_norm(norm))  # type: ignore[arg-type]


def _windows_run(
    result: T2StarResult, gap_mult: float, k: float, use_uncertainty: bool
) -> WindowsResult:
    return windows.run(
        windows.make_inputs_from_frame(
            result.frame,
            time_col="t_rel_s",
            value_col="t2star_s",
            sigma_col="t2star_error_s" if use_uncertainty else None,
            thresholds=THRESHOLDS,
            dataset_id=str(result.meta.get("dataset_id", "")),
            gap_mult=gap_mult,
            k=k,
            use_uncertainty=use_uncertainty,
        )
    )


def _ledger_run(
    window_result: WindowsResult,
    bench_size_table: pd.DataFrame,
    *,
    alpha: float,
    min_events_pass: int,
    tie_cutoff_distinct: int,
    lag_max: int,
    n_permutations: int,
    seed: int,
    include_c3: bool,
    c3_n_null_sim: int,
) -> check_ledger.CheckLedger:
    return check_ledger.run(
        check_ledger.make_inputs_from_windows(
            window_result,
            bench_size_table,
            thresholds=THRESHOLDS,
            alpha=alpha,
            min_events_pass=min_events_pass,
            tie_cutoff_distinct=tie_cutoff_distinct,
            lag_max=lag_max,
            n_permutations=n_permutations,
            seed=seed,
            include_c3=include_c3,
            c3_n_null_sim=c3_n_null_sim,
        )
    )


job = Job(PREFIX)
_bench = job.load_df(BENCH_SIZE_TABLE)

_ledgers = []
for _filename in DATASET_FILES:
    _name = _label(_filename).replace(" ", "_")
    _ds = Dataset(
        path=f"data/real_private/6D2S/{_filename}",
        schema=track912Schema,
        qubit=_qubit(_filename),
        device="6D2S",
        extra={"run_name": _filename[: -len(".pickle")]},
    )
    _node = job.load(_ds)
    _filtered = job.step(_filter_step(RAMSEY_CONFIG), _node, name=f"filter_{_name}")
    _final = job.step(_final_stage, _filtered, name=f"final_{_name}")
    _result = job.step(_t2star_run, _final, name=f"t2star_{_name}")
    _wins = job.step(
        _windows_run,
        _result,
        name=f"windows_{_name}",
        gap_mult=GAP_MULT,
        k=K,
        use_uncertainty=USE_UNCERTAINTY,
    )
    _ledgers.append(
        job.step(
            _ledger_run,
            _wins,
            _bench,
            name=f"ledger_{_name}",
            alpha=ALPHA,
            min_events_pass=MIN_EVENTS_PASS,
            tie_cutoff_distinct=TIE_CUTOFF_DISTINCT,
            lag_max=LAG_MAX,
            n_permutations=N_PERMUTATIONS,
            seed=SEED,
            include_c3=INCLUDE_C3,
            c3_n_null_sim=C3_N_NULL_SIM,
        )
    )

_survey = job.step(
    build_independence_survey,
    *_ledgers,
    name="independence_survey",
    dataset_labels=tuple(_label(f) for f in DATASET_FILES),
    alpha=ALPHA,
)

job.materialize(_survey, name=f"{PREFIX}_grids")

# The single-image overview first: it is the one a reader opens to orient, and the
# per-instrument figures are where a p-value is actually read.
job.figure(
    IndependenceSurveyOverviewPlot,
    _survey,
    targets=["static", "academic"],
    title=f"{PREFIX} overview all checks",
)

for _PlotClass in SURVEY_PLOTS:
    job.figure(
        _PlotClass,
        _survey,
        targets=["static", "academic"],
        title=f"{PREFIX} {_PlotClass.__name__}",
    )
