# SPEC 0007 - The R boundary

**Phase:** 6
**Budget:** 12-19 h. Above `quebraplan.md`'s 6-14 h, and the reason is that the phase's premise
inverted. Phase 6 is written for a machine with no R. R 4.5.3 is installed here with `copula`
1.1.7, so 6.1 and 6.3 are removed, 6.4 shrinks to a local test, and three items the plan does
not contain are added: four tracked documents that assert R is absent, a ledger field that
records a filesystem path where it claims a version, and a composite job that now shells out to
R without saying so.
**Depends on:** SPEC 0006 complete and committed (`278a2fd`).
**Feeds:** `quebraplan.md` 9.1, whose ADR on the subprocess R bridge rather than `rpy2` consumes
R7.6 directly, and 9.4, whose disclosure section consumes the dated search in R7.0.
**Blocks:** nothing. `quebraplan.md` §1 records Phase 6 as blocking nothing and that remains
true. Registering C3 in `battery.ROW_KEYS` is blocked BY this spec's Not done section, which is
the opposite relation and is stated there.
**Reference:** `spec/quebraplan.md` §2, Phase 6, items 6.1 through 6.5, and §7 item 2.

**Numbering.** 0007 is reserved for the R boundary by four sites: `specinstallabity02.md:68`,
`spectests06.md:19`, `spectests06.md` in its Not done section, and `specvalidity08.md:13`. This
spec takes that number so all four stay correct.

Rationale lives in `spec/quebraplan.md`. This file contains requirements and acceptance criteria
only, plus the measurements in R7.0 that change what the phase should do. Do not re-derive or
re-argue the decisions below.

The governing constraint: **R is present on this machine and four tracked documents say it is
absent. The work of this phase is not porting and not packaging. It is making the artifacts a
reviewer reads true, and recording which R produced a number.**

---

## Preconditions

Baseline on `dev` at `278a2fd`, measured rather than assumed:

- `pytest --collect-only -q` gives **592 tests**. Every requirement states its collect delta.
- `pytest -m r --collect-only -q` collects nothing and exits 5: 592 deselected. The `r` marker is
  declared at `pyproject.toml:110` and has **zero members**.
- `scripts/fast-selector.txt` is `not slow and not heavy and not r`, so the fast gate already
  excludes the tier before it has anything in it.
- `Rscript` is on PATH. `R version 4.5.3 (2026-03-11)`, `copula` 1.1.7, `randtests` 1.0.2, all
  four packages resolving only from `/home/sera/R/library`, which `--vanilla` drops and
  `c3_serial_copula.py` restores through an explicit `R_LIBS`.
- `.github/workflows/ci.yml` is one job over Python 3.11, 3.12 and 3.13. There is no R job and
  its header says an R job is out of scope there.
- `pyproject.toml:59-60` carries `r = []` under `[project.optional-dependencies]`, commented
  `# placeholder for now`.

---

## R7.0 - Premises in the plan that do not survive contact with the code

### R7.0.1 The phase was written for a machine with no R, and that inverted

Phase 6 reads throughout as preparation for an absent interpreter. R 4.5.3 is installed here with
all four packages the project uses. **Nine tracked sites still assert or imply the opposite**, one
of them generated into a tracked artifact. Enumerated, because a count without a list is what let
an earlier draft of this spec ship an acceptance grep that passed over five of them:

1. `jobs/bench/report.py:536-540`, emitted into `jobs/bench/results/promotion_report.md:5`:
   "`Rscript` is absent on this machine" and "no code path past `_invoke_rscript` has ever
   executed". Both false.
2. `src/quebra/analyzers/checks/battery.py:19`: "the bench calls it separately when R exists (it
   does not here)".
3. `tests/test_checks_c3_bridge.py:3-6`, module docstring: "R is absent on this machine" and "The
   code past `_invoke_rscript` has never run".
4. `tests/test_checks_c3_bridge.py:66-67` says "the promotion report's four-check scope" while
   `report.py:536` says five checks. They cannot both be right.
5. `src/quebra/analyzers/checks/c3_serial_copula.R:10-11`: "Historic: R was absent on the machine
   this was written on, so this script had never run."
6. `src/quebra/analyzers/independence_survey.py:16-20`: "Seven grids, because `battery.ROW_KEYS`
   has seven entries" (it has nine) and "C3 is NOT here ... the survey turns it off to stay
   affordable", while the composite sets `INCLUDE_C3 = True` and `:42-51` adds `C3_KEY` to
   `SURVEY_KEYS`.
7. `src/quebra/plots/independence_survey_plot.py:5`: "seven classes in `ROW_KEYS` order". The code
   at `:180-181` was corrected; the docstring was not.
8. `docs/iid_checks/C3_serial_copula.md:63-68`: `not computed` is "every row in every ledger
   produced so far".
9. `docs/iid_checks/LIMITATIONS.md:46` carries the `n^2.8` growth exponent that
   `c3_serial_copula.py:210-215` retracted as "not in the data".

Two sites already tell the truth and are NOT in this list, which is why the count is nine rather
than a rounder number: `docs/iid_checks/C3_serial_copula.md:3-8` opens "R is installed and C3 has
executed" (it contradicts itself at `:63-68`, which is site 8), and `docs/iid_checks/LIMITATIONS.md:42`
already says "Where R is present C3 executes".

The first work item is therefore document correction, not a port and not packaging.

### R7.0.2 §7 item 2 is resolved and the resolution is not recorded anywhere

`quebraplan.md` §7 item 2 reads "**Does a Python Genest-Rémillard serial independence test
already exist?** If so, Phase 6 collapses entirely. Nobody checked." Searched 2026-09-14: the
implementations located are R only, `copula::serialIndepTest` and `copula::multSerialIndepTest`.
Neither `scipy.stats` nor `statsmodels` provides one. The phase does not collapse.

The file is unedited at `278a2fd` and still reads "Nobody checked". A spec built on the opposite
answer that leaves the question standing is the stale-claim class `AGENTS.md` §4 exists to
prevent, so R7.1 amends it in place, on the precedent of `spectests06.md` and
`specvalidity08.md` filing their findings into §7.

### R7.0.3 6.1 would produce a statistic with no consumer

6.1 asks for `randtests::bartels.rank.test` ported to NumPy/SciPy, "This one genuinely is short",
venue check "Reduces install friction". There is no friction to reduce. Bartels exists in exactly
three places: `jobs/rscripts/reference_values.R:145-147` and `:156-158`, four rows of the
committed fixture `jobs/reference/r_reference_values.csv`, and a paragraph at
`tests/test_r_cross_implementation.py:219-223` that declines to use it:

> Cross-checking C5's own statistic against `randtests::bartels.rank.test` would be an
> apples-to-oranges comparison and is NOT done here [...] The Bartels values are in the fixture
> as context, not as a target.

No Python module computes a Bartels statistic and no pipeline step consumes one. Removing
`randtests` would also not remove the R requirement, because `copula` is the package that gates
the bridge. **Declined**, with the re-open trigger recorded in Not done.

### R7.0.4 `pip install quebra[r]` cannot deliver what 6.2 asks of it

6.2 keeps `copula::serialIndepTest` "behind `pip install quebra[r]` plus an `Rscript` subprocess
bridge", citing core/extras splits as standard practice and `rpy2` as shipping the pattern. The
extra is empty and must stay empty: `copula`, `XICOR`, `energy` and `randtests` are CRAN
packages, `Rscript` is an interpreter, and none is pip-installable. The bridge's own Python
imports are `csv`, `os`, `shutil`, `subprocess`, `tempfile`, `pathlib` and `numpy`, and numpy is
already a hard dependency. The `rpy2` analogy does not transfer, because `rpy2` IS a Python
package, which is why an extra works for it, and `c3_serial_copula.py:8-11` refuses `rpy2` by
name. So `pip install quebra[r]` is a no-op and will remain one. **Reshaped** by R7.4.

### R7.0.5 6.2's justification names the wrong check and the wrong role

6.2 calls the copula test "your Check 2 validity gate". Wrong twice. C2 in this codebase is
`c2_anderson_darling`, an unrelated check with its own `ROW_KEYS` rows. And no check is a gate:
`AGENTS.md` §5 states that the band is always computed and always drawn and that promoting a
check to a gate is an open question. C3 is also not load-bearing today. `battery.ROW_KEYS` has
nine entries and none is C3, `jobs/bench/results/size_table.csv` and `power_table.csv` carry zero
C3 rows, and `jobs/bench/report.py:536-540` says outright that C3 carries no evidence and that
its silence is not a pass.

The decision 6.2 reaches - keep the R call, do not reimplement - survives and is already shipped.
Its stated reason does not.

### R7.0.6 6.3's actionable raise would violate a domain invariant

6.3 asks the bridge to "check for R at **call** time and raise an actionable error naming the
extra and the `Rscript` PATH requirement". It deliberately does not raise:
`c3_serial_copula.py:227-228` returns `CheckResult(p_value=None, notes="R unavailable")`.

**The obvious argument for that decline is wrong and is not used here.** A raise would NOT kill the
run. `check_ledger.py:481` is `except (ValueError, KeyError, RuntimeError, subprocess.TimeoutExpired)`,
and the comment above it anticipates exactly this case: "With R installed the call executes for
real, and an uncaught RuntimeError kills the whole ledger job instead of writing the `not computed`
row this except clause exists to write." A `RuntimeError` is caught, written as a `declined:` row
and scored `VERDICT_NOT_COMPUTED` at `:266-267`. The `icontract` analogy fails on its own terms:
that ban exists because `ViolationError` subclasses `AssertionError` and escapes a narrower clause,
whereas this clause names `RuntimeError` explicitly.

**The decline stands on a different argument.** Because that except clause is broad, a raise would
make "R is absent" and "the bridge blew up" arrive at the ledger as the same `declined:` row,
indistinguishable after the fact. That is a loss of information the `not computed` path currently
preserves, and it is already a recorded review finding. Naming the extra would also misdirect, per
R7.0.4.

The import-cleanly half of 6.3 is already satisfied and R7.5 adds the guard it never had.
**Declined** as to the raise.

### R7.0.7 There are no R tests for a CI job to run

6.4 wants R tests skipping when R is absent and **required** in a dedicated
`r-lib/actions/setup-r` job, budget 3-5 minutes. The `r` marker has zero members, so a `setup-r`
job would install R and run nothing. The one R-aware test skips in the wrong direction:
`tests/test_checks_c3_bridge.py:58-69` skips when `rscript_path()` is **not** None, so on this
machine it is inert and prints an unactioned instruction every run.

The 3-5 minute budget has no measurement behind it and the two branches differ by roughly an
order of magnitude: CRAN `copula` 1.1-7 is `NeedsCompilation: yes` and imports `gsl`, whose
`SystemRequirements` is the GNU Scientific Library >= 2.5. No Posit binary for that version was
located. **Reshaped**: the marker and the tests land here, the CI job is deferred in Not done
with the measurement that would settle it. State the consequence plainly, because 6.4 warned about
it: after this phase the bridge's only pinning tests run in no automated gate. They run on a
developer machine that has R, and nowhere else.

### R7.0.8 The phase is about provenance, not install friction

The live costs are not packaging. `check_ledger.py:323` writes an `Rscript` PATH into a field
named `r_version` that `check_ledger.py:19-21` documents as holding the version, so two runs
against `copula` 1.1.7 and 1.2.0 are indistinguishable in provenance. `c3_serial_copula.py:138-141`
records that the `R_LIBS` search path "is NOT written to any provenance artifact", and all four
packages here live only in a user-specific directory. For a tool whose first principle is
provenance, an optional out-of-process dependency whose version and library path never enter the
run record is the gap that matters. R7.2 closes it.

### R7.0.9 Turning R on silently changed what a shipped job produces

`jobs/composite/independence_survey.py:111` ships `INCLUDE_C3 = True`. With R absent every C3 cell
produced `None` and drew as `not computed`. With R present the next survey run shells out for
every eligible cell and produces real, **uncalibrated** p-values in the same cells. Nothing in the
artifact, the provenance record or the check outcome distinguishes the two states after the fact.
The survey artifact also hardcodes `"c3_excluded": True` into `meta` regardless of what ran, and
nothing reads the key back. R7.5 measures the transition and deletes the false key.

---

## R7.1 - Make the nine tracked sites true, and amend `quebraplan.md` §7 item 2

**What.** Correct all nine sites enumerated in R7.0.1. `jobs/bench/report.py:535-540` calls
`c3.rscript_path()` and emits what it finds rather than a hardcoded string; regenerate
`jobs/bench/results/promotion_report.md` with `python jobs/bench/report.py`, which reads the two
committed CSVs and does not re-run the bench. Amend `spec/quebraplan.md` §7 item 2 in place to
record the 2026-09-14 search, its result, and that Phase 6 proceeds.

**Why.** Every one of the nine is a tracked claim that is now false, and one is generated into a
tracked artifact. `AGENTS.md` §4 requires a measured number in prose to come from the artifact it
cites in the state it ships.

**Acceptance.** An ENUMERATED site list, not a grep sweep. A grep-only acceptance is what let an
earlier draft pass over five of the nine, which is the positive-control failure `AGENTS.md` §4
names. Each of the nine is shown corrected at the checkpoint. The three greps
(`"is absent on this machine"`, `four-check`, `n\^2.8`) are kept as a REGRESSION guard on top, not
as the criterion. `python jobs/bench/report.py` exits 0 and its diff against the committed report
is shown.

**Budget.** Collect +0 to +2. Files 10.

**Settled in this pass, 2026-09-15.** The three per-call C3 figures were in dispute: one spot
check had measured 1.58 s at n=50 against a quoted 3.9 s. Re-measured, and the outcome splits.

The **statistics and p-values are correct** and are kept: n=50 gives 0.0057908 / 0.9575, n=150
gives 0.0071326 / 0.9036, n=355 gives 0.0076331 / 0.8656, matching the shipped text to every digit
it quotes. What was missing is the **data recipe**, which is the actual defect: the shipped text
says "on iid exponential input at `seed=1`", but `seed=1` is the C3 SIMULATION seed and the
data-generation seed was never recorded, so no reader could reproduce the statistic. The recipe,
recovered by testing candidates against the shipped values, is one `numpy.random.default_rng(0)`
drawing `exponential(1.0)` sequentially in the order n=50, n=150, n=355. It is recorded at every
site that quotes the numbers.

The **timings are re-quoted**, measured 2026-09-15 under R 4.5.3 and `copula` 1.1.7 at
`n_null_sim=1000`: 2.1 s at n=50, 7.1 s at n=150, 54.6 s at n=355, against the shipped 3.9 s,
14.7 s and 130.2 s. Each site now names the machine context, because a wall-clock number without
its machine is not checkable. No growth exponent is quoted anywhere:
`c3_serial_copula.py:210-215` retracted "n^2.8" as "not in the data" and
`docs/iid_checks/LIMITATIONS.md:46` still carried it.

With this settled, later requirements may quote a C3 runtime.

---

## R7.2 - Make the `r` tier safe to populate

**What.** Add a `pytest_itemcollected` hook in `conftest.py` recording every collected item
**pre-deselection**, and have `tests/test_marker_discipline.py`'s totality guard read that record
rather than `request.session.items`, which is post-deselection. Add a session-scoped
`requires_rscript` fixture so each new `r` test does not re-implement the probe. Correct the marker
description at `pyproject.toml:110`, `Makefile:62` and `quebraplan.md:235`, which carry three copies
of it.

**AMENDED during implementation, on a measurement.** This requirement first named
`pytest_collection_modifyitems` with `tryfirst=True`. That form LOSES a hook-ordering race: a
competing `tryfirst` deselector registered later runs first and removes the item before the recorder
sees it. Measured twice independently, by the implementer and then by a verifier reproducing it from
scratch: with such a competitor present, `itemcollected` sees 596 items and `tryfirst` sees 595, and
the planted item is visible only to the former. `pytest_itemcollected` fires per item inside
collection itself, structurally before any `modifyitems` implementation runs, so it cannot depend on
ordering between plugins. The SPEC is corrected here rather than the code, so a later session does
not restore the losing form.

**Why.** Measured: a planted `r`-marked test with no tier marker fails the totality guard on a full
run and **passes** under the repo's own fast selector, which `.github/workflows/ci.yml:44-45` makes
the only selector CI runs. Today that is harmless because the tier is empty. The moment it gains a
member the tier axis stops being total under CI. The hook form is chosen over the one-line
alternative (deleting ` and not r` from `scripts/fast-selector.txt`, as `Makefile:6-7` does for
`real`) because the record does not depend on which selector ran.

**What it does NOT buy, measured.** Totality under `-m "unit"`. The guard is itself a `policy` test,
so that selector deselects the guard and no hook can help. Nor does it survive a path restriction,
which narrows what is collected at all. The full run is what makes the axis total, and the shipped
docstring says so.

**Acceptance.** MUTATION, run and reported: plant an `r`-marked, tier-less test; the guard goes red
under `pytest -m "not slow and not heavy and not r"` after the change and green before it. Show both.
A second test pins the record as PRE-deselection by running a child pytest under a narrowing
selector and asserting the record exceeds what survived: without it, a replacement returning the
post-deselection list would leave the suite green.

**Budget.** Collect +1. Files 4.

**Not part of this.** `Makefile:69-78`'s `test-real` target line already exists, so the coupling it
documents is already closed. Verifying it is a read, not an edit, and it is not counted in the file
budget.

---

## R7.3 - The ledger must record which R produced the p-value, and at what N

**What.** Replace `CheckLedger.r_version` (`check_ledger.py:130`, `:323`) with three fields:
`r_version` holding what `Rscript -e 'cat(R.version.string, as.character(packageVersion("copula")))'`
returns, `r_executable` holding the path the field holds today, and `r_library_paths` holding the
tuple `c3.r_library_paths()` already computes. Add `c3_n_null_sim` as a fourth field, promoted from
`CheckLedgerInputs:185`. Add the version probe beside `r_library_paths` in `c3_serial_copula.py`.
Correct `check_ledger.py:19-21`.

**Why.** R7.0.8. `c3_n_null_sim` is the more load-bearing of the two gaps and was nearly missed:
`lag_max` and `seed` are already on the artifact, `c3_n_null_sim` is not, and `check_ledger.py:518-531`
overwrites the C3 `notes` string (which carried `N`) with the verdict reason, so `N` is dropped from
the row as well. `jobs/composite/independence_survey.py:112` sets `C3_N_NULL_SIM = 200` against a
module default of 1000, and `c3_serial_copula.py` states that a p-value at one `N` is not the same
object as one at another. A shipped survey therefore records p-values whose `N` appears nowhere.
It is included here rather than deferred because R7.3 already accepts the artifact break below, so
the fourth field is free.

**Acceptance.** `r_version` matches `Rscript -e 'cat(R.version.string)'` and names the `copula`
version. `c3_n_null_sim` on a materialized ledger equals the value the job passed. With `Rscript`
stripped from PATH the three R fields read absent and nothing raises, and a test asserts the ledger
still builds all nine rows in that state; the oracle is `AGENTS.md` §5, that a check outcome never
stops execution.

**Budget.** Collect +2 to +4. Files 4.

**Known risk, and it is accepted rather than avoided.** `CheckLedger` subclasses
`StaleArtifactGuard` (`check_ledger.py:112-113`), and `core/_artifact_guard.py:43-51` raises
`ValueError` on unpickle when a field is missing, derived from `dataclasses.fields()`. **Adding
these fields makes every already-materialized `CheckLedger` pickle in `output/` unreadable**, so
every composite run with `--reuse-deps` over an existing ledger will raise until the sub-job is
re-run. That is the decision taken here: the fields are worth the regeneration. It is a
provenance-continuity cost, not a caching one, because `core/closure.py:13-16` records that reuse
already keys on a repo-wide `git_commit` and any commit invalidates every cached artifact anyway.

**What this risk is NOT.** An earlier draft warned about the Mermaid label. That warning was
unfounded and its acceptance could not fail: the label is built from step kwargs
(`provenance.py:128-133`), and a dataclass field on the RESULT cannot reach it. Separately, identity
`code` is the job file text plus a sha256 of every `quebra` module in the static import closure
(`core/job.py:239-250`, `core/closure.py:191-198`), so editing `check_ledger.py` moves the identity
of every ledger job by construction. So do R7.1's `battery.py` edit and R7.6's
`independence_survey.py` edit. Identity movement is expected in this phase, not a stop condition.

---

## R7.4 - Three `r`-marked tests

**What.**
(a) **Write** a test putting the bridge against `jobs/reference/r_reference_values.csv` on **both**
statistic and p-value at seed 707, N 1000, `lag.max` 5, in `tests/test_r_cross_implementation.py`,
whose index is that oracle. `analyzers/instrument_validation.py:558-575` is a `TierRow` carrying
`TIER_ABSENT` that NAMES this comparison and says it "needs R at test time, which the suite refuses
to require". It is a declaration that the test is absent, not a test: no test in the suite invokes
the bridge, and `test_r_cross_implementation.py:231-245` reads the fixture without R. Promote that
`TierRow` when the test lands.
(b) **Add** the fixture's ten meta rows checked against the local R, BESIDE the eight `>= 0`
assertions at `tests/test_r_cross_implementation.py:69-72` rather than replacing them. The `>= 0`
assertions are vacuous, but `r_value()` (`:50-57`) raises when a row is missing or duplicated, so
the existing test does assert one real property unconditionally: that the fixture carries version
metadata at all. Moving that behind `r` would delete it from every CI leg, which `AGENTS.md` §7
governs and this spec does not attempt.
(c) The bridge's timeout path, driven directly as `c3.run(..., timeout_s=0.5)`.

**Why.** Naming an oracle is not detecting anything. (a) must assert both quantities because the
statistic is seed- and N-invariant, so a statistic-only assertion cannot detect a swapped
`seed`/`n_sim` argument order.

**Acceptance.** `pytest -m r -q` selects 3 and exits 0 with `Rscript` present, and reports 3 skipped
with PATH stripped. MUTATION 1: swap the `n_sim`/`seed` positional arguments at
`c3_serial_copula.py:163-165`; test (a) goes red while a statistic-only version stays green. Show
both. MUTATION 2: delete `env["R_LIBS"]`; test (a) goes red.

**Budget.** Collect +3 to +5. Files 5.

**Why (c) does not go through the ledger.** `check_ledger.py:463-473` calls `c3.run()` with no
`timeout_s` and `CheckLedgerInputs:155-185` has no timeout field, so the ledger's `declined:` branch
cannot be reached by a real timeout without new plumbing. The three available routes were
monkeypatching, which `AGENTS.md` §7 forbids for this tier, plumbing `timeout_s` through
`CheckLedgerInputs`, which is unbudgeted scope, and a 15-minute test. A direct bridge test is
honest about covering the bridge's timeout and not the ledger's handling of it. The ledger's
`declined:` branch therefore remains untested and is recorded in Not done.

**Stop condition, three-way.** If the fixture-reproduction test does not reproduce the fixture, the
cause is one of: a defect in the shipped bridge, a stale fixture, or R or `copula` having moved
under it. Read the fixture's `meta` rows against the local R FIRST, which is what test (b) exists
to make cheap, and only then decide. Do not widen a tolerance.

---

## R7.5 - Settle the `r` extra and publish the real install route

**What.** Either delete the `r` extra from `pyproject.toml` or keep it empty with a comment of
three lines or fewer stating that R and CRAN packages are not pip-installable and that the name
exists to match the `r` pytest marker. Either way `# placeholder for now` goes. Add the two actual
install commands to `README.md`. Do not change `_UNAVAILABLE = "R unavailable"`.

**Why.** R7.0.4. A placeholder that looks like an install route is what makes 6.2's phrasing read
as implementable. `specinstallabity02.md:68` created this artifact for this phase.

**Acceptance.** `git grep -n "placeholder for now"` returns nothing. `README.md` names the two
commands. A test asserts `import quebra` succeeds with `Rscript` stripped from PATH, which is the
guard 6.3's surviving half never had. The install half is NOT re-asserted here: `scripts/acceptance.sh`
already builds a wheel, installs it outside the repo and runs the suite from a non-checkout
directory, and a pytest test shelling out to `pip install` would be slow, network-touching and
environment-mutating, which `AGENTS.md` §1 rules out.

**Budget.** Collect +0 to +1. Files 3.

**Decision owed to Sera.** Delete versus keep-empty is a packaging call.

---

## R7.6 - Record that R appearing changed what the survey produces

**What.** Delete the false `"c3_excluded": True` at `analyzers/independence_survey.py:304`; nothing
reads it and it has been a recorded review finding. Run the composite survey with R present and
report, as numbers: wall clock, how many C3 cells produced a p-value, how many read `not computed`,
and the reason distribution across multi-process, too-short and declined.

**Why.** R7.0.9. A shipped job changed what it produces because of an environment change, and
nobody has measured the difference.

**Acceptance.** `git grep -n c3_excluded src/ jobs/ tests/` returns nothing. The transition is
reported as numbers, not as a direction.

**Budget.** Collect +1 to +2. Files 3, plus a bounded run.

**Wall-clock ceiling, because this is otherwise unbounded.** The run uses the job's own
`C3_N_NULL_SIM = 200`, NOT the module default of 1000; the "6.5 h over 415 cells" figure in the
job's comment is for `N = 1000` and is not this configuration's cost. Cap the requirement at **2 h
of wall clock**. If it overruns, fall back to five datasets rather than all 34 and say so in the
report. A full-survey run at an unmeasured cost is not inside a 12-19 h phase, and `AGENTS.md` §9
makes budget a STOP rather than a report.

**Stop condition.** If the regenerated C3 grid shows rejections, whether an uncalibrated rejection
may appear in a thesis figure is a scientific call and Sera's, not the implementer's.

---

## R7.7 - The delegation claim

**What.** Write the delegation paragraph at the top of `docs/iid_checks/C3_serial_copula.md`, four
sentences: the statistic is computed by `copula::serialIndepTest` and not by this package; the
bridge is a subprocess and a CSV round trip rather than `rpy2`, and the package imports without R;
the bridge is pinned against the reference implementation's own output by an `r`-marked test naming
the fixture; C3 is UNCALIBRATED, with no bench cell, no size and no power evidence. File the
bench-runner defects recorded in Not done into `quebraplan.md` §7 as numbered entries.

**Why.** `quebraplan.md` 6.5. It cannot be written before R7.1 through R7.6 land, because the
repository's own tracked artifacts currently contradict it.

**Acceptance.** Every one of the paragraph's four claims is traceable to a command in R7.1 through
R7.6's acceptance. No claim in it asserts calibration, power or level.

**Budget.** Collect +0. Files 3.

---

## Mutation ledger for checkpoints 7.2 and 7.3

Every test these two checkpoints add or change, with the mutation it must catch and the observed
result. Run 2026-09-27 against the settled tree. A test with no row here is not evidence.

| # | Mutation applied to the source | Test that must catch it | Result |
|---|---|---|---|
| M1 | `_r_provenance(inputs.include_c3)` becomes `_r_provenance(True)` | `test_a_ledger_that_never_asked_for_c3_records_not_asked_and_never_probes` | RED |
| M2 | `not asked` collapses to `absent` | same | RED |
| M3 | the `try/except` around `c3.r_version()` is deleted | `test_a_version_probe_that_raises_is_recorded_and_does_not_stop_the_ledger` | RED |
| M4 | the library sentinel becomes a silent `()` | `test_an_unusable_r_is_recorded_on_the_ledger_and_does_not_stop_it` | RED |
| M5 | the bridge is handed `c3.N_NULL_SIM` instead of the inputs' value | `test_the_recorded_n_is_the_n_handed_to_the_bridge` | RED |
| M6 | the artifact records `c3.N_NULL_SIM` instead of the inputs' value | `test_c3_n_null_sim_reaches_the_artifact_as_the_caller_set_it` | RED |
| M7 | the fixture returns the post-deselection list | `test_the_record_is_pre_deselection_not_post` | RED |
| M8 | the rejected `modifyitems(tryfirst=True)` hook form is restored | `test_the_record_survives_a_competing_tryfirst_deselector` | RED |
| M9 | `pytest_itemcollected` becomes a no-op | `test_every_collected_test_carries_exactly_one_tier_marker` | RED |
| M10 | `r_version`'s empty-stdout sentinel is deleted | `test_a_probe_that_answers_nothing_is_distinguishable_from_one_that_never_ran` | RED |
| M11 | the ledger's empty-library sentinel is deleted | same | RED |
| M12 | `absent` collapses to `not asked` | `test_with_no_rscript_the_r_fields_read_absent_and_every_row_is_still_built` | RED |

M12 covers the test R7.3's own Acceptance clause names ("a test asserts the ledger still builds
all nine rows in that state"). It was missing from the first version of this table, which claimed
to list every test these checkpoints touch: a totality claim that was not total.

**M10 and M11 were GREEN when first measured**, and the whole suite passed with either sentinel
deleted. Both sentinels had been added in response to a review finding with no test, which is the
defect this ledger exists to make visible. The test named against them was written afterwards.

**The harness that found them was itself defective on its first run.** It passed a malformed node
id (`tests/test_check_ledger.py::`, trailing colons) for those two rows, which pytest treats as a
usage error and exits non-zero, so both scored as caught without running anything. A mutation
harness reporting a non-zero exit as a caught mutant cannot fail either. Re-run against the whole
suite, both survived.

---

## Checkpoints

The marker-tier work comes BEFORE the ledger work, because R7.3's acceptance is an R-dependent
assertion and no `r`-marked test may land until R7.2 has made the tier axis total.

```
CHECKPOINT 7.0 - this spec lands. No behaviour change.
CHECKPOINT 7.1 - the nine tracked sites say what is true; quebraplan §7 item 2 records the
                 resolved search; the promotion report is regenerated. (R7.1)
CHECKPOINT 7.2 - the r tier can gain members; the planted-test mutation is reported. (R7.2)
CHECKPOINT 7.3 - the ledger records R, copula, the library path and c3_n_null_sim. Existing
                 CheckLedger pickles are knowingly invalidated. (R7.3)
CHECKPOINT 7.4 - three r-marked tests and two mutations run. (R7.4)
CHECKPOINT 7.5 - the r extra is settled and the install route is published. (R7.5)
CHECKPOINT 7.6 - the survey transition measured, inside the 2 h cap; the false meta key gone. (R7.6)
CHECKPOINT 7.7 - the delegation paragraph and the quebraplan §7 entries. (R7.7)
```

The amendment to `quebraplan.md` §7 item 2 belongs to CHECKPOINT 7.1 only. An earlier draft assigned
it to both 7.0 and 7.1.

---

## Not done

- **Porting `randtests::bartels.rank.test`** (6.1). Declined per R7.0.3. Re-open trigger: a
  consumer appears that needs a Bartels statistic in Python.
- **An `r-lib/actions/setup-r` CI job** (6.4). Deferred, not declined. The cost is unmeasured and
  brackets by roughly 20x between a Posit binary and a source build reaching `gsl` and the GNU
  Scientific Library at the system level. What would settle it: one timed run of `setup-r` plus
  `install.packages("copula")` on a clean ubuntu-latest runner. Consequence stated in R7.0.7: until
  it exists, the bridge's only pinning tests run in no automated gate.
- **The ledger's `declined:` branch.** `check_ledger.py:481-493` catches `RuntimeError` and
  `TimeoutExpired`, and with R present it is the only thing between a 900 s timeout and a dead
  survey job. R7.4(c) tests the bridge's timeout, not the ledger's handling of it, because the
  ledger does not pass `timeout_s`. Closing this properly means plumbing `timeout_s` through
  `CheckLedgerInputs`, which is a separate change.
- **Registering C3 in `battery.ROW_KEYS` and benching it.** Blocked, and the blockers are defects
  rather than cost. `jobs/bench/runner.py:207-213` increments `seen` unconditionally but `hits` only
  when `p_value is not None and < ALPHA`, and C3 is the only check in the tree that can return
  `None`. So a registered C3 would report a rejection rate of 0 on every multi-segment cell, and
  `analyzers/calibration_summary.py:441-449` scores on `|z|` against a Bonferroni threshold, giving
  `|z|` near 10 at a few thousand replicates: C3 would be published as **miscalibrated**, not as
  conservative, on cells where it never ran. `runner.py:169`, `:181` and `:200` catch only
  `(ValueError, KeyError)` while the bridge raises `RuntimeError` and `TimeoutExpired`, so one
  bridge failure aborts a multi-hour run. `report.py:658-673` asserts in generated prose that every
  row of a replicate reads the same permutation set, which a C3 row makes false. All three must land
  before C3 enters `ROW_KEYS`, not after. R7.7 files them into `quebraplan.md` §7.
- **Caching the simulated null.** `c3_serial_copula.R:44` documents a Python-side cache that does
  not exist. R7.1 deletes the claim rather than building the thing, because a cache introduces
  dependence between replicates within a bench cell and `runner.py:243`'s binomial standard error
  would no longer be the error of the reported size.
- **Changing `_UNAVAILABLE`.** Pinned at `tests/test_check_ledger.py:90-94` and it reaches drawn
  figure cells. The actionable text goes in README prose instead.
- **Whether the `m > 1` decline and the n-too-short decline are the right scientific calls.**
  `c3_serial_copula.py:29-34` states the reasoning; this phase verified only that they return
  `None` with a reason.
- **A test for the C3 tie verdict.** `check_ledger.py:85` lists `c3_serial_copula` in
  `TIE_SENSITIVE_CHECKS` but `tests/test_check_ledger.py:83-87`, whose docstring names C3, asserts
  only C6, C5 and C1.
- **`spec/spec01hygiene.md:6,8` and `spec/specinstallabity02.md:7,9` cite `spec/PLAN.md`, which
  does not exist.** Not this phase's errand. Recorded so the next reader knows it is uncaught
  rather than clean: the stale-reference guard excludes `spec/` deliberately.

---

## Open questions this phase does not settle

1. What a full `independence_survey` run costs with C3 live at `N_NULL_SIM = 200`. R7.6 caps the
   measurement at 2 h rather than answering it.
2. Whether the regenerated C3 grid shows rejections.
3. Whether a GitHub runner would resolve the same R package versions. The simulated p-value is the
   part most at risk from a `copula` version change, and R7.4's exact p-value assertion is what
   would detect it.
4. Whether an uncalibrated check belongs in a published figure at all. R7.6's acceptance makes it
   concrete for the first time.
