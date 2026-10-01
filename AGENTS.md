# AGENTS.md

QUEBRA is a reliability-statistics toolkit for qubit quality time series. It carves
threshold-excursion durations out of a metric record and applies survival and repairable-systems
statistics to them. Correctness, reproducibility and provenance come first. Python 3.11+.

The pipeline is a lazy DAG: nothing runs until a sink (figure or materialize) resolves.

**Implemented today:** Kaplan-Meier, MTBF and MTBC, window carving with an explicit gap policy,
an independence check battery, and the metric analyzers.
**Not implemented:** Nelson-Aalen, log-rank, RMST, MCF. `grep -rli` finds no module for any of
them. Do not describe them as shipping.

This is the only agent-facing file, and `CLAUDE.md` is a symlink to it so it loads every
session. Read `spec/quebraplan.md` only when told which phase to work on.

---

## 1. Hard rules

**Never run `git add`, `git commit`, `git push`, `git tag`, `git reset`, `git checkout`.**
Sera commits. You stop at checkpoints and say so. These are denied at user scope, so an attempt
will fail. Do not work around it by invoking git through Python or a shell script.

**Never write into `data/`, `output/`, `output_backup/` or `output_backup2/`.** You may read
all of them. You may write code that reads and writes them at runtime. You may not edit a
dataset or a materialised artifact directly. Machinery, never evidence. `data/` splits three
ways by redistribution status: `real_private/` holds the embargoed records, inside the
checkout; `real_public/` is for records that may be shared and `simulated/` for payloads that
are regenerated.

**Never modify `.claude/settings.json`, `.claude/hooks/`, or your own permissions.**

**Do not self-review deterministic properties.** If a property can be checked by running a
command, run it and report the exit code. Do not assert that lint passes, that types check, or
that tests pass without having run them.

**Ask before installing anything.** No `pip install`, no new dependency, without approval.

---

## 2. Commands

```
pytest                                        # full suite
ruff check --fix . && ruff format .           # lint and format
quebra run jobs/active/<job>.py               # run one job
quebra run --all                              # every job declaring JOB_SWEEP = True
quebra run --all --family t2star              # every job in one family, sweep opt-out included
quebra inspect [job_file]                     # print the graph without running it
```

`make lint`, `make test` and `make clean` wrap the first two. `make types` runs mypy over
`core/`; `make arch` runs the import-linter contract. `make check` is all four, and `make
deps` runs deptry.

`make check` uses the tools installed here. CI installs the newest version every `>=`
admits into a fresh non-editable environment, so `check` passing is not evidence about CI.
`make check-ci` runs the same steps against a clean resolve in a throwaway environment.
That is the one to run before a push.

---

## 3. Principles

Each carries its reason, so it holds as the code grows. Judge against these, not against current
state.

**Steps are pure compute.** The DAG caches and re-runs steps by content hash, so a step must be a
deterministic function of its inputs. Input and output therefore live at the edges: disk access
in loaders, file writes in render targets, drawing in plot and panel classes. A step never reads
disk and never imports matplotlib or plotly.

Panels and plots **are** the render layer and they draw only. Every data-derived quantity lives
in the typed artifact, computed by an output-builder such as
`analyzers/within_calibration_compute.build_within_calibration_panel_data`. The renderer keeps only
functions of axes and theme: adaptive limits, decade-guide ticks, colours. Cumulative time and
its relatives are data, not render-local views.

**Physical quantities carry unit suffixes** (`_hz`, `_rel_s`, `_unix_s`, `_local_dt`, `_utc_dt`).
Unit and clock mismatch is the costly bug class here; see `docs/TIME_SEMANTICS.md`. Canonical
frequency keys: `rabi_hz`, `qubit_frequency_hz`, `raw_frequency_hz`, `delta_hz`.

**Errors are raised, not swallowed.** This is a reproducibility tool, so a silent fallback yields
a wrong-but-plausible result. `Job.load()` raises rather than defaulting `run_start` to
`t_raw[0]`, which would be roughly 1970.

**`icontract` preconditions only where a violation SHOULD abort.** It is a runtime dependency
and it is used in exactly one module, `analyzers/kaplan_meier.py`. It is FORBIDDEN under
`analyzers/checks/`: `icontract.ViolationError` subclasses `AssertionError`, not `ValueError`,
so it escapes the `except (ValueError, KeyError)` that `analyzers/check_ledger.py` uses to turn
a check failure into a reported row, and would kill a run that is supposed to draw a band with
a `not computed` cell beside it. A check outcome never stops execution. Enforced by
`tests/test_contract_boundary.py`, not by this paragraph.

**A file earns its existence from a consumer or a boundary, never from a line count.** Split on
a layer the two halves must not cross (drawing or disk on one side, a pure step on the other), on
a fan-in of three or more sibling importers, or on a piece whose docstring declares it
replaceable. Otherwise one file per consumer: a module serving one figure keeps its vocabulary,
typed artifacts, helpers and output-builder together, as `analyzers/independence_survey.py` does.
A contract split across six files is a contract no reader can see, so fragmentation costs more
than length. Judge size against the largest single-consumer module already in the same package,
not against a quota.

**Results are typed dataclasses, not raw dicts**, each defined in the file with the step that
returns it. Analyzers use `make_inputs_from_norm(...)` then `run(inputs)`.

**Formats and render targets are added by registering** (`loaders/registry.py`,
`plots/targets.py`), not by ad-hoc code elsewhere. Schema validation goes through pandera. One
decorated function adds a format or a target with zero other changes.

**Provenance and outputs are append-only.** Every figure sink and every materialize emits
`.prov.json` and `.prov.md`, and `output/` is never deleted. Step kwargs that affect output
appear on the Mermaid label.

**No class hierarchies for job families.** Categorisation is by directory and a discriminator,
never by subtype polymorphism. Do not introduce `ValidationJob` / `SurveyJob` base classes.

---

## 4. Claims discipline

Each of these has recurred across sessions. The instance is kept so the rule sticks.

**Never call a mean "the cell".** If a filter leaves a factor varying, the number is a pooled
mean: name the factor and carry the spread. Recurred three times: `report.size_by_n`,
`censoring_effect`, then `calibration_summary.size_vs_n` AFTER the first two were fixed. At n=20
the two Weibull shapes give 0.0650 and 0.0865; a single line labelled "the primary null cell" hid
a spread wider than the effect being drawn.

**A measured number in a docstring must come from the artifact it cites, in the state it ships.**
Quoting a scratchpad run is the same defect as inventing it. Instances: "0.0607" (no such cell;
the table says 0.0757), "44-90 null cells" (44-70), "z_crit 3.3-3.5" (3.25-3.38), a test pinning
0.061 borrowed from a different generator (measured 0.0634). Compute it into the prose, or cite
the table and stop.

**Fix every site of a class, not the one that failed.** A guard at one comparison and not its
twin is not a fix. Instances: `TAU_MARGIN` applied to `tau <= T_N` but not `tau < T_N` (refused
valid records 23% of the time); gap truncation applied to interior segments but not the final
one; the negative-variance split applied to `gamma_hat` but not `gamma_hat_batch`.

**A positive control must fail when the thing it guards is broken.** Assert against the REAL
output, not a copy of the reference. Instances: the carve control compared `reference` to
`reference` and passed with a deliberately broken carve; the AST walk missed `importlib` and
relative imports that its own docstring named; `check_agreement` aligned on a key missing `clock`
and INVENTED 36 pairings.

**State the aggregation and the multiplicity before reading a verdict off it.** Max-over-grid
flattered power by 7x; a flat +/-0.01 size band rejected all seven rows because the max of 44-90
deviations is roughly 3 SE by chance. Say which aggregate gates, and correct for how many cells
it saw.

**When a claim cannot be cheaply verified, write it as the open question it is.** "I have not
measured this" is cheaper than the review that finds it false.

---

## 5. Domain invariants

Violating these is a scientific error, not a style problem. None of them fails a test.

**In-spec means the metric is at or above the threshold.** For T2\*, in-spec is
`T2* >= threshold`. Never invert this.

**Never resample a metric time series onto a uniform grid.** Resampling destroys and invents
threshold crossings; measured loss is roughly 25 to 45 percent of real crossings at working
thresholds. If a function needs regular spacing, it is the wrong function.

**Statistical licensing.** The Kaplan-Meier confidence band rests on the renewal assumption
(`analyzers/assumptions.a1_renewal_durations`), and so would a k-sample log-rank comparison when
one exists - it does not today, and `kaplan_meier.compare()` returns a DISTANCE, not a test.

**A band is never reported without its check outcome attached**, and the attached outcome must
name which checks were asked for and which produced no answer. Attachment is a FIELD on the band
artifact (`KaplanMeierComparison.checks_asked` / `checks_unanswered` / `check_verdicts`), not two
sibling files in one directory: a reader who opens the band alone must still see what was
checked. Empty means NOT ASSESSED, and `check_summary()` says so - which is a different claim
from "assessed and nothing rejected". A grid of `not computed` cells satisfies the letter of
"attached" and says nothing, which is why the naming half is part of the rule.

The band is ALWAYS computed and ALWAYS drawn. A failed check annotates it and never suppresses
it: control flow that depends on what the data happened to say is unpredictable, and a reader is
better served by a band they are told not to trust than by a missing one. Promoting a check to an
actual gate is an open question (`spec/quebraplan.md` section 7), not current behaviour.

**Locked vocabulary.** These are not synonyms and must never be substituted.

- `within-calibration`: metric series, KM and NA estimators, threshold excursion windows.
- `across-calibration`: calibration event records, MCF, repair effectiveness.
- Do not use `repairable` / `non_repairable` as OUR vocabulary. The literature's own term is
  a separate matter: `repairable system` is
  standard usage from Ascher and Feingold and from Rigdon and Basu, and it stays in prose that
  cites that field, because rewriting it there would make the sentence false.
  `panels/across_calibration.py` carries the canonical note on why our tiers are named after the
  calibration boundary instead.
- `in-spec fraction` for the within-calibration quantity. `availability` is reserved for the
  across-calibration systems tier.
- Load-bearing terms, never reworded: **window, read, bag, check, band, scan clock, window age,
  birth type, run-set, display-set, check outcome**.
- `check outcome`, never "trust annotation" or similar. The code's word is `verdict`, and both
  the ledger and the survey exist to stop a non-rejection reading as reassurance; "trust" invites
  exactly that reading.

---

## 6. The codebase

### Layout (real)

The packages live under `src/quebra/`. Paths below are relative to that, except
`jobs/`, `tests/`, `scripts/` and `docs/`, which stay at the repository root: `jobs/` is the
researcher's analysis configuration rather than library code, and keeping it out of the
wheel is what stops `output/` being written into site-packages.

```
core/        dataset.py, job.py, runner.py, types.py (Norm, Measurement, CalibrationEvent)
provenance.py
loaders/registry.py            built-ins: .csv  .yaml/.yml  .h5/.hdf5  .pkl/.pickle
schemas/     base.py, track912.py (912-day Ramsey), calibration_log.py
transforms/  filter.py, interpolate.py, lookup_prior.py (public check_unix_s lives here)
analyzers/   allan.py, fidelity.py, t2star.py, tlf.py, mtbf.py, psd.py (stub)
             windows.py  in-spec window carving: gap policy, censoring, read table
             signal_band.py, distinguish_band.py, reliability_band.py  the three nested
             band contracts of the within-calibration panel; shape_stats.py
             calibration_summary.py  reshapes the bench tables into the four typed
             artifacts the calibration figures draw; also re-simulates the P-P curve
             checks/  C1 Lewis-Robinson, C2 Anderson-Darling, CvM Cramer-von Mises,
             C5 rank autocorrelation, C6 exchangeability, and C3 copula-via-R (needs
             Rscript; not in ROW_KEYS, so the battery never runs it).
             battery.py runs the five permutation checks off ONE shared permutation set.
             These are steps: pure compute, and the pipeline may import them.
             Permutation calibration needs an explicit rng. block_permutations raises on
             None, because defaulting to OS entropy made p-values irreproducible while the
             run identity stayed unchanged.
panels/      within_calibration.py, across_calibration.py   generic render components;
             adapters feed them. _within_calibration_render.py is the
             functions-of-axes half. The artifact itself is built in
             analyzers/within_calibration_compute.py and typed in
             analyzers/within_calibration_data.py.
plots/       base.py, targets.py, theme.py, *_plot.py   targets: static, academic, poster
jobs/active/ ramsey_*.py, ramsey_2x2_*.py, t2star_*.py, mtbf_*.py, check_calibration.py,
             check_ledger_6d2s_q*.py (checks only: per-record ledgers + per-qubit survey)
jobs/composite/ compare_*.py, independence_survey.py   job.include + .ref across datasets;
             declare JOB_SWEEP = False
src/quebra/recipes.py          RAMSEY_CONFIG + configure_ramsey_job orchestrator, and
                               wire_t2star_carve, the one T2* carve.
                               Reusable library code: keeping it out of jobs/ is what
                               lets the CLI avoid putting the caller's directory on
                               sys.path, which is forbidden.
src/quebra/ledger_recipe.py    configure_check_ledger_job and the survey wiring. Apart
                               from recipes.py so the checks stay out of the T2* jobs'
                               identity closure.
src/quebra/cli.py              the CLI, installed as the `quebra` console script.
                               Anchors output/ and the --all glob on the WORKING
                               DIRECTORY, never on __file__, so an installed copy
                               cannot write into site-packages.
scripts/acceptance.sh          clean-venv acceptance: builds the wheel, installs it
                               outside the repo, and runs the suite from a directory
                               that is not the repository.
tests/                         tracked
```

`jobs/reference/` holds TRACKED external validation data consumed by both the suite and a
figure: the published load-haul-dump record (Kvaloy and Lindqvist, Technometrics 2020,
Section 8.1) and the R
reference values written by `jobs/rscripts/reference_values.R`. It sits beside the jobs that
declare it as a Dataset. `reference_values.R` is never run by the test suite; the fixtures are
committed so pytest works without R.

`jobs/bench/` holds a TRACKED calibration study for `analyzers/checks/`: five generating arms,
the real carve via `analyzers/windows.py` primitives, `results/` (`size_table.csv`,
`power_table.csv`) and `promotion_report.md`. It is a study, not a pipeline layer. NOTHING
outside `jobs/bench/` may import it; `tests/test_bench_isolation.py` asserts this across every
pipeline package. Figures needing its numbers declare `jobs/bench/results/*.csv` as a Dataset and
`job.load_df` it, so the dependency runs through provenance instead of around it.

Gitignored: `output/`, `output_backup/`, `output_backup2/`. Generated, append-only.
`scripts/` holds `acceptance.sh`, `promote_run.py`, `make_data_manifest.py`,
`make_fixtures.py` and `make_job_manifest.py`.

Note: `provenance.get_git_commit` and `is_tree_clean` both anchor on the working directory, so
they describe the project you are running in, not the one this package is installed into.

### Job shape

See any `jobs/active/ramsey_*.py`.

- Load the main dataset with `job.load(ds) -> norm`; load the companion with
  `job.load_df(ds) -> DataFrame`.
- Enrich with `lookup_prior(main, comp, fields=[...], aliases={...})`. Source columns
  (`frequency`, `Rabi_frequency`) stay as-is and map to `qubit_frequency_hz` / `rabi_hz` via
  aliases.
- Either call `configure_ramsey_job(job, enriched, profile=..., include_fidelity=..., ...)` or
  wire steps by hand with `job.step(fn, *inputs, name=...)`.
- End in `job.figure(Panel, panel_node, targets=["static","academic"], title=...)`.

### Data

- Main 912-day Ramsey: `data/real_private/6D2S/{DDMMYY}_6D2S_qubit{N}.pickle`
- Companion calibration pickles: `data/real_private/companion/qubit{N}.pickle`
- Dataset pickles are read-only inputs.

---

## 7. Tests

`tests/` is flat and stays flat. `spec/quebraplan.md` 5.1's six-directory layout was measured
against the suite and **declined** in `spec/spectests06.md`: directories cross-cut the
oracle-and-subject index below, and several files legitimately hold more than one tier. The tier
is carried by a **marker**, not a directory.

**Two orthogonal marker axes.** `unit`, `properties`, `statistical`, `integration`, `validation`,
`regression` and `policy` classify the QUESTION a test answers; every collected test carries
exactly one, enforced by a guard. `slow`, `heavy`, `real` and `r` classify its COST or
REQUIREMENT and are independent: a test may be `statistical` and `slow`.

**Oracle rule, effective now.** Any test asserting a statistical result must name its oracle in
the test name or the first line of the docstring: an analytic value, a reference implementation,
or a simulation truth. A test that cannot name an oracle is a change detector, not evidence.
Do not write tests that assert what the code currently returns.

**Naming an oracle is not detecting anything.** A test can cite a source and still be unable to
fail: because it re-types the value it claims to check, because it divides by the quantity it
claims to pin, or because both sides of an identity descend from the helper being mutated. All
three shapes were found in this suite. When a test is the evidence for a claim, break the code it
guards and confirm it goes red.

**Record the mutation and its result in the checkpoint banner.** An unrecorded mutation is
indistinguishable from one that was not run, and the failure it guards against is silent: a test
that cannot fail passes forever and reads as evidence. Two instances, both in one phase. A test
written to catch a fixture being swapped probed the stash behind the fixture instead, so the swap
left the suite green. A test written to prove a collection hook beats a competing deselector built
the competitor as a `-p` plugin, which registers before the root conftest and therefore cannot win
the race, so it passed whichever hook shipped. Both were mutation-verified only after review asked.
A sentinel added in response to a review finding went in twice with no test at all, and deleting
either left the whole suite green.

**Test files are indexed by oracle and subject, not by source module.** Two test files sharing
both an oracle and a subject are one file, at any length: splitting them duplicates the fixture
and leaves neither able to show which one is the evidence.

**Deleting a test.** Volume of deletion is the wrong metric: a deleted test that was catching
something is an undetectable regression, and it is the one operation whose damage is invisible to
every gate here. Before removing one, all three must hold.

1. **Name the property it asserts, and name the test that still asserts it.** If no other test
   does, it is not redundant - it is the only evidence.
2. **Show the removal is safe by mutation, not by reading.** Break the code the test guards and
   demonstrate the REMAINING suite goes red. If the suite stays green the deletion is refused and
   the mutation result is reported, because a test whose guarded code can break with the suite
   still green has found a second gap rather than proved itself redundant.
3. **Record it in the checkpoint banner**: the property on `Not done`, the mutation result on
   `Known risk`.

Tests requiring R **skip** when `Rscript` is absent. They never pass with mocked values.

---

## 8. Docs hygiene

This repo has drifted here before. Hold the line.

- `docs/` holds reference docs: `TIME_SEMANTICS`, `PANEL_CONTRACT`, `FIGURE_STANDARD`, `JOBS`,
  `WRITING_A_JOB`, `WRITING_A_SCHEMA`, `GOLD_STANDARD` (every published number the checks
  reproduce, with its locator), and `iid_checks/` (one page per check plus limitations).
  Architecture rationale lives there, not in this file. Refresh docs; do not narrate evolving
  architecture here.
- Every `.md` file is tracked. `.gitignore` carries no blanket `*.md` or `.*` rule, so `.md`
  needs no `!` exception. The three `!` rules it does carry are data-manifest carve-outs, not
  doc ones.
- `FIGURE_STANDARD` binds figures you add or edit. The existing panels are not yet conformant and
  that doc says so.
- Do not create new long `.md` docs unprompted.
- Every claim in a doc must match the code. If unsure, verify against the code; do not assert.
- Do not document deferred or speculative design as if normative. Label it deferred, and keep it
  out of this file.
- No status or progress tables in agent-facing docs; they go stale. Track status in commits and
  issues.

---

## 9. Checkpoint protocol

Work is cut into numbered checkpoints defined in the phase spec. At a checkpoint, stop and print
exactly this, then wait:

```
CHECKPOINT <n.n> - <one line: what this checkpoint achieved>

  Changed:      <paths>  (<count> files)
  Gates run:    ruff <exit> | pytest <exit> (<n> passed) | <other> <exit>
  Budget:       collect +<actual> of +<low>..+<high> stated | files +<actual> of +<stated>
  Not done:     <what a reader might assume was done but was not>
  Known risk:   <what could break, especially cached identities>
  Suggested:    <conventional-commit message, one line>

  Review `git diff` and commit if you see fit. I will not proceed until you say so.
```

**`Budget` is a STOP, not a report.** Every requirement in a phase spec states an expected
collect delta and file count. Exceeding either by more than 2x halts the phase: print the banner,
say which requirement overran and why, and wait. Do not carry the overrun into the next
checkpoint.

The line exists because the numbers already did and were absorbed anyway. SPEC 0008's
per-requirement envelopes summed to +19 to +34; the phase landed at +99, every envelope exceeded
2 to 4x, at every checkpoint, while the spec's own "any difference explained rather than absorbed"
went unenforced. A budget nothing halts on is a budget that is not being kept.

**Never review a tree that is being edited.** Freeze it, run the review once, apply the findings,
then re-review if you must. A reviewer measuring a moving tree reports gate results that describe
no state that ever existed, and it will keep finding defects indefinitely because it is reading the
previous round's fixes. Measured: one review watched the tree change four times underneath it, and
two of its findings were about tests the previous round had added an hour earlier.

`Not done` and `Known risk` are mandatory and must not read "none" unless that is literally true.
They are what makes the diff review fast. Do not proceed past a checkpoint on your own
initiative, even when the next step seems obvious.

---

## 10. Writing rules

These apply to every docstring, comment, spec, ADR and doc page you write.

- **No em dashes anywhere.** Spaced hyphens.
- **No measured wall-clock number as a claim about what the code costs.** State the order of
  magnitude and say it is load-dependent. A timing is a property of a machine under a load nobody
  recorded, so a quoted interval is falsified by the next run: three successive bands in one
  docstring each were.
  Two things this does NOT cover. A dated historical measurement that EXPLAINS A CONSTANT is a
  design record, not a claim about current cost, and it stays: `c3_serial_copula.run()` keeps the
  130.2 s that set its 900 s timeout, because deleting it would leave the constant unexplained.
  And analytic values, counts and versions are reproducible; they belong in prose with their recipe.
- No "surfacing", "brings into view", "data-driven", "delve", "leverage" as a verb, "firstly" as
  an orphaned ordinal, or "excellent" as hyperbole.
- No overclaiming. "To our knowledge" is used deliberately and sparingly, not as a hedge.
- Short declarative paragraphs. Colon expansions over dense subordinate clauses.
- Never invent a section number, equation number, figure number, or citation. If you do not know
  the locator, write "no source located".
- **A module docstring states what the module is for and what contract it holds, not how the
  author got there.** "Arm C turned out to be a null arm" and "looking a row up by its display
  label silently missed every cell" are debugging history: they date on the first refactor and
  teach nothing about the code in front of the reader. If a past bug still constrains the design,
  state the constraint. This binds code you write or edit; instances already committed, such as
  `jobs/bench/arms.py`, are not a cleanup errand.
- **No spec or phase identifier in a source docstring.** `SPEC 0008 R8.x` in a header is a phase
  artifact left in a permanent file, and it outlives the phase by years. Cite a spec mid-file
  where one decision needs its source, as the existing `quebraplan` references do.
- **Section 8's ban on normative deferred design covers docstrings.** Write what the module does
  now, not what it is positioned to become.
- Statistical docstrings carry a `Validity assumptions` section with four fields per assumption:
  assumption, diagnostic, consequence of violation, reference.
- Specs and docs in Markdown, not LaTeX. What makes a spec work is numbered requirements,
  explicit acceptance criteria and a done-when clause, not the markup.

---

## 11. Scope and workflow

Do not build deferred or speculative work unprompted: composite job machinery beyond what exists,
a real `psd.py` (currently a stub), or large new subsystems. If you ARE asked to, implementing it
is correct and welcome. These are scope calls, not standing prohibitions.

- New functionality: Plan mode (Shift+Tab). The plan is approved by Sera before any edit.
- Implement on Sonnet by default, Opus when Sonnet is not enough.
- Before a checkpoint: `/review` (Opus reviewer, fresh context). ruff and tests green.
- Cold second opinion on a plan or architecture: `/plan-critique`.
- Check a doc against the code: `/doccheck`.