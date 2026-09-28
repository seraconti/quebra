# QUEBRA hardening plan

Cross-adjudicated from five independent survey reports, verified against the live JOSS
review checklist (retrieved 20 August 2026) and the SlopCodeBench preprint.

Evidence tags follow the survey convention: `[DOCUMENTED]` retrieved from a primary source,
`[MEASURED]` a numerical result from a study, `[PRACTICE]` widespread convention without
controlled evidence, `[INFERRED]` my own reasoning.

---

## 0. The two facts that set the order

**Fact 1 - the clock.** The JOSS review checklist now contains, verbatim:

> Open development: Was the software developed openly from early stages? For projects with
> recent public repositories, is there at least six months of public development history with
> evidence of releases, public issues/pull requests, and ideally external engagement?

`[DOCUMENTED]` https://joss.readthedocs.io/en/latest/review_checklist.html

**The repository is already public, so the clock is already running.** Check the actual date
(`gh repo view seraconti/qretool --json createdAt`, or the first push to the public remote) and
write it down, because it sets your earliest JOSS submission date and therefore whether JOSS is
a thesis deliverable or a post-thesis one. What is *not* yet running is the second half of the
same checkbox: "evidence of releases, public issues/pull requests". Tags and issues are the part
you still control, and they are addressed in Phase 3C rather than now, for the reason you gave -
issues opened against an architecture that is about to change describe a repository that will
not exist in a month.

**Plan the thesis and the TQE methods paper to stand without JOSS regardless.** JOSS acceptance
is not on your critical path; reproducibility of the methods paper is.

**Fact 2 - installability is the only thing a reviewer can independently verify.** The
Functionality checkbox is "Does installation proceed as outlined in the documentation?"
`[DOCUMENTED]` A reviewer runs your documented steps on their machine. `core/paths.py`
resolving data as "one directory above the git root" fails that literally. Every other
recommendation in every report is unverifiable until this is fixed.

The ordering below is **dependency order plus clock order**, not importance order.

---

## 1. Phase budget and sequencing

| Phase | What | Hours | Blocks what | Can it wait? |
|---|---|---|---|---|
| **0** | Hygiene week: gitignore, pyc, session logs, agent governance, rename, licence | **6-10** | everything | **No.** Nothing later is tracked or enforceable without it. |
| **1** | Installability: pyproject, src layout, importlib.resources | **8-16** | 1B, 2, 8 | No. Sphinx autodoc needs an importable package. |
| **1B** | Data and outputs layout | **6-10** | 5, 7.5 | No. The test tiers are defined in terms of these directories. |
| **2** | Boundaries: import-linter, deptry | **3-6** | 3 | No, but it is cheap. Adopt with the ratchet so it does not block the refactor. |
| **3** | Hash narrowing, logical job IDs, job parameterisation, `v0.2.0` + issues | **12-24** | 5 (the `jobs/` tier), 9.2 | Partially. 3.1 and 3.2 before 3.3, always. |
| **4** | Validity contract: assumption records, contracts, test-naming traceability | **10-18** | 5.3, 8.2 | Yes, but it is the distinctive part and the thesis wants it. |
| **5** | Test architecture | **30-55** | 9 | No. Largest block, highest scientific value. |
| **6** | R boundary | **6-14** | nothing | **Yes** - and check the open question first; it may vanish. |
| **7** | Presentation: plot contract, styles, deterministic PDF | **10-20** | thesis figures, poster | Partly. Do it before you regenerate thesis figures, not after. |
| **8** | Documentation and diagrams | **20-36** | 9 | No, if JOSS matters. |
| **9** | Paper artifacts: ADRs, comparison table, AI disclosure | **10-18** | JOSS submission | No. Also feeds the thesis directly. |

**Total 121-227 hours**, or roughly 15 to 28 focused working days. That is a semester-shaped
programme, not a sprint, and it runs in parallel with thesis writing rather than instead of it.

**The one hard ordering constraint:** 0 → 1 → 1B → 2 → 3.1 → 3.2 → 3.3. Everything after that
reorders freely against your other deadlines. Phase 6 is the only phase that might disappear
entirely.

---

## 2. Ordered intervention sequence

### Phase 0 - Hygiene week (6-10 h, do it as one block)

This is not deferred hygiene and it is not "the boring part before the real work". Every other
phase in this document is untracked or unenforceable until it lands. Budget one focused block.

| # | Change | Hours | Why | Evidence | Standard? | Venue check |
|---|---|---|---|---|---|---|
| 0.1 | Fix `.gitignore`: stop excluding `.*`, `*.md`, `tests/`. Add the new exclusions from Phase 1B in the same commit. | 0.5 | From outside, the repo currently appears to have no tests, no docs, no CI config. Every subsequent item is untracked until this is fixed. Do it once, correctly, with the data and outputs rules already in place, so you are not rewriting `.gitignore` three times. | `[DOCUMENTED]` JOSS checks "Automated tests" and "Community guidelines" | Universal | JOSS: two checkboxes fail outright |
| 0.2 | Untrack any committed `.pyc` | 0.5 | Tracked `.pyc` files permanently break the clean-working-tree gate in your reuse policy. Already identified as a prerequisite by the plan critic. | `[INFERRED]` from your identity scheme | n/a | n/a - internal correctness |
| 0.3 | Preserve agent session logs to a tracked directory now | 0.5 | Claude Code `.jsonl` transcripts are subject to a documented deletion sweep. Once gone, the AI usage disclosure section and the thesis methods chapter become unevidenced assertion. Nothing later can recreate them. | `[DOCUMENTED]` JOSS AI usage disclosure is a required section | Emerging | JOSS: required section. IEEE TQE: AI disclosure mandated. |
| 0.4 | Agent governance: permissions, deny rules, commit protocol (see §3) | 2-3 | Without this, every subsequent phase costs you a keystroke per tool call, or you turn permissions off and lose the safety. Do it before the work, not after. | `[DOCUMENTED]` Claude Code permission model | `[PRACTICE]` | n/a |
| 0.5 | Rename to QUEBRA, fix the `repairable`/`non_repairable` to `within-calibration`/`across-calibration` drift in the same pass, tag `v0.1.0` | 2-4 | One-way door, and it must precede the camera-ready URL. GitHub redirects old paths. A rename is the only cheap moment to fix load-bearing vocabulary, and your own rule says these are not synonyms. The tag also starts satisfying "evidence of releases". | `[INFERRED]`; `[DOCUMENTED]` JOSS Open development checkbox mentions releases | n/a | EQTC: name and URL must agree |
| 0.6 | Decide the licence deliberately | 0.5 | GPL-3.0 constrains vendor use (OrangeQS, IBM). No report raised it. One-way door once anyone else depends on it. | `[INFERRED]` | n/a | JOSS: any OSI licence passes; the constraint is your adoption goal |

Nothing here touches architecture, so nothing here is invalidated by Phases 1 to 3.

---

### Phase 1 - Installability (8-16 h)

| # | Change | Why | Evidence | Standard? | Venue check |
|---|---|---|---|---|---|
| 1.1 | `pyproject.toml` with PEP 621 metadata, hatchling backend | Standard project metadata mechanism; every downstream tool (`deptry`, `import-linter`, `ruff`, `pytest`) configures from it. | `[DOCUMENTED]` PEP 621 | Universal in scientific Python | JOSS: "clearly-stated list of dependencies, ideally handled with an automated package management solution" |
| 1.2 | `src/quebra/` layout, real namespace, `__init__.py` everywhere | Implicit namespace packages work until setuptools auto-discovery has to choose. `core`, `loaders`, `plots`, `schemas`, `transforms`, `jobs/composite` currently lack them. | `[PRACTICE]` scientific Python ecosystem default | Universal | JOSS: install must succeed from a clean environment |
| 1.3 | Replace `core/paths.py` resolution with `importlib.resources` | Path resolution assuming a git checkout is the single most likely concrete review failure. | `[DOCUMENTED]` stdlib since 3.9 | Universal | JOSS: Functionality checkbox |
| 1.4 | Acceptance command: `pip install .` in a fresh venv, then import smoke test, then fast tests | Catches exactly the current failure mode. Make this the definition of done for Phase 1. | `[INFERRED]` | `[PRACTICE]` | JOSS: this is literally what the reviewer does |

Do this as a behaviour-preserving commit series. Do not fold the scientific refactor into it.

---

### Phase 1B - Data and outputs layout (6-10 h)

No report answered this. It matters because it decides what a reviewer can run, what a
collaborator can hand you under embargo, and what your identity scheme hashes.

**The two path mechanisms are different and must not be conflated.** `importlib.resources` is for
*packaged* resources: the tiny synthetic fixtures that ship inside the wheel. It is the wrong tool
for datasets. Datasets live outside the package and are found through a resolved data root.
`[INFERRED]`

```
data/                        <- outside src/, never packaged
  real_public/               <- committed if small; pooch registry if not
  simulated/                 <- generated; only seeds + manifests committed
  real_private/              <- gitignored entirely
    MANIFEST.toml            <- COMMITTED: filename, sha256, provenance, embargo status
output/                     <- gitignored entirely
  <job-id>/<identity>/
    run.json                 <- COMMITTED for published figures only
src/quebra/
  _fixtures/                 <- tiny, packaged, reached via importlib.resources
```

| # | Change | Hours | Why | Evidence | Standard? | Venue check |
|---|---|---|---|---|---|---|
| 1B.1 | Three-way `data/` split: `real_public/`, `simulated/`, `real_private/` | 2 | The split is not cosmetic, it is a *licence and redistribution* boundary. It lets you accept non-disclosable data from a collaborator, compute on it, and publish results, without ever risking a commit. It also makes the boundary legible in the C4 Context diagram (8.5). | `[DOCUMENTED]` ACM artifact guidance permits proprietary artifacts to remain unavailable provided the dependency is documented and proxies are included | `[PRACTICE]` | ACM badging: proxies required. JOSS: original-data results presented in a submission must be reproducible by reviewers, so anything shown in the paper must have a public or synthetic path. |
| 1B.2 | `data/real_private/MANIFEST.toml`, committed, listing filename, sha256, provenance, embargo status | 1 | This is the single highest-value small artifact in the whole plan. A reviewer sees exactly what was used without receiving it; you can verify a collaborator sent the same file twice; and the sha256 is already what your identity scheme folds in, so the manifest is a human-readable view of a hash you compute anyway. | `[INFERRED]`; aligned with `[DOCUMENTED]` ACM guidance on documenting proprietary dependencies | `[PRACTICE]` | ACM/ICSE: this is what "document the dependency" means concretely |
| 1B.3 | Data root resolution order: `QUEBRA_DATA_ROOT` env var, then `quebra.toml` at repo root, then a `platformdirs` user data dir | 2-3 | Replaces "one directory above the git root". Env var first so CI and a reviewer's machine can both work without editing files. Not `importlib.resources`, which cannot reach outside the wheel. | `[INFERRED]` | `[PRACTICE]` | JOSS: Functionality checkbox - a reviewer must be able to point the tool at their own data |
| 1B.4 | Presence detection, not failure: a loader for a missing private dataset raises a named `DataUnavailable` carrying what was expected and where it was looked for | 1 | Anyone without your pickles should get a legible message and a working simulated path, not a stack trace. Same pattern as the R degradation in 6.3. | `[INFERRED]` | `[PRACTICE]` | JOSS: reviewer will hit this on day one |
| 1B.5 | `output/` gitignored entirely; commit `run.json` only for figures that appear in a paper | 1-2 | Correct, outputs are too big. But an ungitignored `output/` with nothing committed means published figures have no auditable trail. The manifest is kilobytes and carries node IDs, input hashes, artifact IDs, software version, timings and cache status. Commit those and the figure is reproducible without the artifacts. | `[INFERRED]` | `[PRACTICE]` | JOSS: reproducible execution of paper figures |
| 1B.6 | `src/quebra/_fixtures/` for the tiny packaged synthetic fixtures, reached via `importlib.resources` | 1 | These ship in the wheel so that `pip install quebra && pytest --pyargs quebra` works for a reviewer with no repository checkout. | `[DOCUMENTED]` stdlib | `[PRACTICE]` | JOSS: reviewer verification path |

---

### Phase 2 - Boundaries (3-6 h)

| # | Change | Why | Evidence | Standard? | Venue check |
|---|---|---|---|---|---|
| 2.1 | `import-linter` `layers` contract in `pyproject.toml`, running in CI | `core/job.py` imports `loaders`, `transforms`, `plots.base` and root `provenance`, so the compute core sits above every layer it is meant to sit beneath. A diagram does not fix this; a failing build does. | `[DOCUMENTED]` import-linter v2.13, released 2026-07-03 | `[PRACTICE]`, growing | JOSS: no explicit item, but it is the evidence behind the required "Software design" section |
| 2.2 | Adopt against the violating codebase using `ignore_imports` as a ratchet | `unmatched_ignore_imports_alerting` defaults to `error`, which lets you enumerate current violations, freeze them, and remove entries as you fix them. Without this you cannot turn the contract on until the refactor is finished. | `[DOCUMENTED]` import-linter option surface | `[PRACTICE]` | n/a |
| 2.3 | Add `acyclic_siblings` once `jobs/` splits into `validation/`, `survey/`, `comparison/` | Catches cycles among the new sibling packages, which is the failure mode the reorganisation introduces. | `[DOCUMENTED]` | New | n/a |
| 2.4 | `deptry` in CI | Detects missing, unused and misclassified dependencies against `pyproject.toml`. Cheap. | `[DOCUMENTED]` | `[PRACTICE]` | JOSS: dependency list must be accurate |

`import-linter` over `tach` and `pytest-archon`: native TOML config, precise dependency-path
reporting, and the rule you need is dependency direction rather than a public-interface system.
All three reports that discussed it agreed.

---

### Phase 3 - Identity before parameterisation (12-24 h)

**This is the phase the majority of the reports got wrong.** Four of five said "collapse the
63 job files" as a standalone Adopt at 4-12 hours. It is not standalone.

| # | Change | Why | Evidence | Standard? | Venue check |
|---|---|---|---|---|---|
| 3.1 | Narrow the hash: per-step sources rather than the module file, resolved parameter row rather than the table | Today, editing `analyzers/km.py` already invalidates all 63 jobs. File count is not what protects your cache. Parameterise before narrowing and the blast radius gets *worse*, not better. | `[INFERRED]`, and the only report that raised it inspected the repo | n/a | n/a - this is internal correctness |
| 3.2 | Replace repo-relative path strings with stable logical job IDs plus a registry | `job.include("jobs/active/ramsey_q1_100423.py")` means every folder reorganisation is a breaking change and every cached identity moves. Logical ID determines the graph node; content hashes determine identity; filesystem location is never identity. | `[DOCUMENTED]` R `targets` demonstrates exactly this separation (target names are logical, storage paths are implementation detail) | `[PRACTICE]` | JOSS: supports the "Software design" section |
| 3.3 | Collapse 63 near-identical jobs to one family definition plus a reviewed manifest | Copy-pasted `_T2STAR_THRESHOLDS` ladders are a table pretending to be source code. | `[DOCUMENTED]` `targets` static branching is designed for precisely "target instances are variations of a shared definition" | `[PRACTICE]`: Hydra, Snakemake wildcards, papermill, `targets` | JOSS: reviewer will notice the duplication |
| 3.4 | Emit a generated, human-readable job manifest | 63 explicit files are trivially greppable and diffable; one parameterised definition hides per-job intent. Generate the readable view so reviewability survives the abstraction. | `[INFERRED]` | `[PRACTICE]` | JOSS: reproducibility of reported results |
| 3.5 | Declarative categorisation: module-level `JOB_FAMILY = "survey"` constant plus `pkgutil` auto-discovery | You explicitly do not want class hierarchies. A constant plus auto-discovery is zero framework, and recategorising costs one string edit. Entry points and `pluggy` solve an extensibility problem you do not have. | `[PRACTICE]` | Widespread | n/a |

**Phase 3C - open the issues (1-2 h, at the end of Phase 3, not before).** Your reason for
holding them is correct: an issue written against `jobs/active/ramsey_q1_100423.py` describes a
file that will not exist. Once logical IDs and the job manifest are in place, the architecture
those issues reference is stable, and the four MSc subprojects can be written as issues that a
student could actually pick up. Open them then, plus tag `v0.2.0`. This is what starts satisfying
the second half of the Open development checkbox - "evidence of releases, public issues/pull
requests, and ideally external engagement" - and the Collaborative effort checkbox, which
lists as not acceptable a single author with no evidence of community engagement.
`[DOCUMENTED]` Waiting three weeks costs nothing on a six-month clock that is already running;
opening them now costs you accuracy.

---

### Phase 4 - The validity contract (10-18 h)

This is the artifact that survives from your university SE project, and it is the one that is
genuinely distinctive.

| # | Change | Why | Evidence | Standard? | Venue check |
|---|---|---|---|---|---|
| 4.1 | An assumption record per estimator: assumption, diagnostic, consequence of violation, whether it raises / warns / reports | Your requirements *are* statistical validity conditions. The KM band and the log-rank merge verdict are two consequences of one assumption (ABGK multiplicative intensity model) and are revoked together. That is a domain-level object, not a runtime invariant. | `[INFERRED]`, aligned with the one report that drew the distinction | `[PRACTICE]`: `lifelines` documents assumptions in prose and warns at runtime; `pingouin` runs assumption checks inline; `statsmodels` exposes diagnostics as separate callables | JOSS: not a checklist item, but it is what makes a reviewer trust a statistics package |
| 4.2 | `icontract` **only** for hard computational invariants (positivity of durations, binary event codes, non-empty risk set) | Two of the five reports would have had you encode graded validity conditions as preconditions that raise. Check 3 is a claim selector with a worst-verdict rule. A precondition cannot express that. | `[DOCUMENTED]` icontract; `[INFERRED]` for the boundary | `[PRACTICE]`, minority | Not checked anywhere |
| 4.3 | Estimator returns or exposes the diagnostic outcome rather than silently deciding the assumption holds | Failing branch means no licensed band and no licensed merge verdict. The caller must be able to see that. | `[INFERRED]` | `[PRACTICE]` | n/a |
| 4.4 | Traceability by test naming, not a matrix: `test_A3_noninformative_censoring_*`, assumption IDs cited in docstrings | Grep becomes the trace query. The measured traceability benefit (roughly 24% faster, 50% more correct) came from 71 subjects consuming pre-built links for third-party code on paper. No study has measured traceability at n=1. | `[MEASURED]` Mäder & Egyed; `[MEASURED]` Rempel & Mäder TSE 2017 (24 medium-to-large OSS projects) | `[PRACTICE]` | JOSS: not checked. Regulated industry: mandatory, and you are correctly far from it. |

---

### Phase 5 - Test architecture (30-55 h, the largest single block)

**Start by naming the problem honestly.** Your current suite was agent-written, and
agent-written tests have a characteristic failure mode: they assert *what the code currently
does* rather than what it should do. Those are change-detector tests. They go green forever, they
break on every legitimate refactor, and they provide zero evidence of statistical correctness.
The measured result behind this is that **test-oracle strength correlates far more strongly with
suite effectiveness than coverage does** `[MEASURED]`.

So the entry rule for the rewrite is: **every test in the `statistical/` and `jobs/` tiers must
name its oracle** - an analytic value, a reference implementation output, or a simulation truth -
in the test name or the first docstring line. Any test that cannot name an oracle either belongs
in `unit/` or should be deleted. That single rule is what converts an agent-generated suite into
a defensible one, and it is cheap to audit with grep.

#### 5.1 Six tiers, with a distinct question each

| Tier | Directory | The question it answers | Oracle | Data | Runs where |
|---|---|---|---|---|---|
| **Unit** | `tests/unit/` | Does this function do the small thing it claims? | the specification of the function | none | CI, every push |
| **Properties** | `tests/properties/` | Do the mathematical invariants hold under arbitrary inputs? | algebra | Hypothesis-generated | CI, every push |
| **Statistical** | `tests/statistical/` | Is the estimator *correct*? | analytic case, R/`lifelines` reference, or simulation truth | generated at test time | CI (coarse grid) + local (full grid) |
| **Integration** | `tests/integration/` | Does the machinery wire up? Step to job to artifact, identity, cache hit and miss, cycle detection, clean-tree gate | the runner contract | tiny packaged fixtures | CI, every push |
| **Validation jobs** | `tests/jobs/` | Does a real end-to-end job produce the right *verdict*? | known truth from the generator | simulated, known parameters | CI (small) + local (real) |
| **Regression** | `tests/regression/` | Did a refactor change a published number? | golden values from a pinned tagged run | `real_private/` | **local only** |

The distinction you asked about - validation jobs versus integration testing - is the fourth and
fifth rows, and they are genuinely different. Integration asks *does the DAG resolve, cache and
materialise correctly*, and it should pass even if every estimator is wrong. Validation jobs ask
*does the whole pipeline reach the correct scientific conclusion on data whose answer we know*,
and it should fail if the estimator is wrong even when the wiring is perfect. Conflating them is
why agent-written suites feel comprehensive and prove nothing. `[INFERRED]`

#### 5.2 The CI / local partition

Markers, declared in `pyproject.toml` so `--strict-markers` catches typos:

```toml
[tool.pytest.ini_options]
addopts = "--strict-markers"
markers = [
  "slow: >30 s",
  "heavy: full simulation grids or large memory",
  "real: requires gitignored data under data/real_private/",
  "r: needs Rscript and the copula package, which is all c3_serial_copula.R calls",
]
```

| Job | Selector | Budget | Why there |
|---|---|---|---|
| CI - fast | `-m "not slow and not heavy and not real and not r"` | under 3 min | Runs on every push. This is what a JOSS reviewer runs. It must pass with **zero** private data and **zero** R. |
| CI - R | `-m "r"`, on `r-lib/actions/setup-r` | 5-8 min | Gated as a separate job so the core matrix stays fast. Must be *required*, not skipped, or the bridge rots silently. |
| CI - nightly | `-m "slow or heavy"` | 20-40 min | Full simulation grids. Scheduled, not per-push, because coverage of the parameter space matters more than latency. |
| Local | `-m "real or regression"` | as long as it takes | The only tier that touches `data/real_private/`. Never on Actions, both because the data cannot leave and because a public runner log is an exfiltration surface. |

Two rules that keep the partition honest: **`real` and `r` tests skip when the prerequisite is
absent, never pass with mocked values** `[DOCUMENTED]`; and **CI must never see
`data/real_private/`**, which holds because the tree is gitignored and the fast-CI selector
excludes the `real` marker. A public runner log is an exfiltration surface, so the `real` tier
never gets a CI job at all rather than getting one that skips.

#### 5.3 The work

| # | Change | Hours | Why | Evidence | Standard? | Venue check |
|---|---|---|---|---|---|---|
| 5.1 | Audit and triage the existing agent-written suite against the oracle rule: keep, relabel, or delete | 4-8 | Do this *first*. Rewriting on top of change-detector tests means you carry them forever. Expect to delete more than you keep, and expect that to feel bad and be correct. | `[MEASURED]` oracle strength over coverage | `[PRACTICE]` | JOSS: "Automated tests" is judged on meaningfulness, not count |
| 5.2 | Synthetic generator: 1/f, RTN, Ornstein-Uhlenbeck, with known ground truth and **independently exposed** latent, observation, censoring/repair processes and seed | 12-20 | Two functions at once: the only test oracle that is not "agreement with the current code", and the only redistributable CI fixture. Independent exposure is what lets an estimator failure be attributed to a specific source of misspecification instead of "something is off". | `[MEASURED]` canonical refs below | `[PRACTICE]` | JOSS: this is what makes the tool usable by someone without your pickle files, which is the substantial-scholarly-effort argument |
| 5.3 | Five-layer estimator validation, in this order: analytic special cases, reference implementation cross-check, simulation with known truth over a varied grid, invariants, then figures | 8-14 | A green unit suite can coexist with a scientifically wrong estimator. This hierarchy is more defensible than any coverage number. | `[PRACTICE]`: `lifelines`, `statsmodels`, `scikit-survival` and `arviz` all lean on reference comparison plus analytic special cases | Dominant real pattern | JOSS: Automated tests. TQE: reviewers will read this as the methods contribution |
| 5.4 | `hypothesis` for algebraic and invariance properties only: survival monotone non-increasing, cumulative hazard non-decreasing, permutation invariance, shape and type preservation, degenerate-case equivalence | 3-6 | Do not ask Hypothesis to discover statistical validity. Use it to stress the deterministic mathematics *around* the estimator. | `[DOCUMENTED]` hypothesis | `[PRACTICE]`, growing | Not checked |
| 5.5 | Runner-contract integration tests: acyclicity, deterministic identity, no reuse under code or content mismatch, no reuse under a dirty tree, valid topological order, cycle detection reports the full chain | 4-6 | These are exactly the invariants one report suggested modelling in Alloy. As executable tests they are cheaper, they run in CI, and they cannot drift from the code. This is also the evidence behind the ADR in 9.1. | `[INFERRED]` | `[PRACTICE]` | JOSS: supports the required Software design section |
| 5.6 | Coverage reported, not enforced. Roughly 90% on estimators and identity hashing, roughly 70% overall, 0% expected on `jobs/active/` | 2-3 | The four reports proposed 90%, 85-90%, 70-80% and "do not enforce". All four numbers are invented. Inozemtseva & Holmes found the coverage-effectiveness correlation collapses once suite size is controlled; the paper was named ICSE Most Influential N-10 in 2024. JOSS specifies no percentage anywhere. Report it in CI so regressions are visible; do not gate on it. | `[MEASURED]` Inozemtseva & Holmes, ICSE 2014, DOI 10.1145/2568225.2568271 | Contested | JOSS: automated tests checked, no percentage |

Canonical noise-model references, verified as real and correctly attributed in the JOSS-focused
report, usable directly in the thesis: Kasdin (1995) and Timmer & König (1995) for 1/f;
Uhlenbeck & Ornstein (1930) and Gillespie (1996) *Exact numerical simulation of the
Ornstein-Uhlenbeck process* for OU. Prefer the exact Gaussian transition over a generic
Euler-Maruyama integrator when the transition law is known. For RTN, parameterise switching
rates and amplitudes directly and validate empirical state occupancy and switching-time
distributions before the signal is used to validate anything else.

---

### Phase 6 - The R boundary (6-14 h)

| # | Change | Why | Evidence | Standard? | Venue check |
|---|---|---|---|---|---|
| 6.1 | Port `randtests::bartels.rank.test` to NumPy/SciPy, validated against R output as oracle | It is a rank-difference statistic with a known asymptotic normal approximation. This one genuinely is short. | `[INFERRED]` | `[PRACTICE]` | Reduces install friction |
| 6.2 | **Keep** `copula::serialIndepTest` behind an `Rscript` subprocess bridge. RESHAPED: the `pip install quebra[r]` half of this row is impossible and the extra was deleted, see item 11 of section 7 and `specrboundary07.md` R7.0.4 | It is a Cramér-von Mises functional of the empirical copula process with a multiplier-bootstrap null. It is also your Check 2 validity gate. Reimplementing your own load-bearing validity test and validating it against the implementation you just deleted is a scientific risk, not an engineering win. | `[INFERRED]`; the most epistemically careful report reached the same conclusion | `[PRACTICE]`: core/extras splits via `[project.optional-dependencies]` are standard, and `rpy2` itself ships this pattern | JOSS: extras splits are accepted; `rpy2` draws reviewer friction |
| 6.3 | Import the package successfully without R; check for R at **call** time and raise an actionable error naming the extra and the `Rscript` PATH requirement | Graceful degradation at call time, never at import time. | `[DOCUMENTED]` | Universal | JOSS: core install must succeed cleanly |
| 6.4 | R-dependent tests `pytest.importorskip`/skip when R is absent, and **required** in a dedicated `r-lib/actions/setup-r` CI job | They must not silently pass with mocked values. Budget 3-5 minutes per run; gate the job so the core matrix stays fast. | `[DOCUMENTED]` r-lib/actions | `[PRACTICE]` | JOSS: reviewer runs the core matrix |
| 6.5 | Say in the paper that the copula test is delegated to the reference implementation | This reads as more careful than a home rolled null distribution, and it strengthens the DAG-runner argument: the R step really is a process boundary, so a file-based engine costs nothing *there*, which makes the rest of the in-memory argument credible. | `[INFERRED]` | n/a | TQE/JOSS: honest scoping |

**Open question that could delete this whole phase:** nobody checked whether a Python
implementation of the Genest-Rémillard serial independence test already exists. Check PyPI and
JOSS before spending anything here.

---

### Phase 7 - Presentation layer (10-20 h)

| # | Change | Why | Evidence | Standard? | Venue check |
|---|---|---|---|---|---|
| 7.1 | Plot contract: accept an optional `ax`, draw into it, return the `Axes`, never call `show()` or `savefig()` inside library code | Makes every plot independently testable and composable into arbitrary layouts, which is exactly your panel-as-step goal. | `[PRACTICE]`: `seaborn` axes-level functions, `lifelines` `plot(ax=...)`, `arviz`, `statsmodels` | Universal in scientific Python | JOSS: not checked; thesis and poster both depend on it |
| 7.2 | Compose with `plt.subplot_mosaic`; a thin composer assembles multi-panel figures | Named-axes mosaics, no duplicated plotting logic. | `[DOCUMENTED]` matplotlib | Standard | n/a |
| 7.3 | Two `.mplstyle` targets, `academic` and `poster`, applied via `plt.style.context(...)`; plotting code stays theme-agnostic | Replaces `plots/theme.py`, which is a qubit colour lookup rather than a theme. Your `plots/targets.py` already has the `["static", "academic"]` machinery for "poster" to slot into. | `[PRACTICE]`, matplotlib rcParams | Standard | Poster and thesis both need it; newtx matching is the constraint |
| 7.4 | `pdf.fonttype: 42`, `ps.fonttype: 42`, deterministic metadata, `svg.hashsalt` | Forces font embedding, avoids Overleaf compilation errors, makes figure files diffable and hashable. | `[PRACTICE]` | Standard | Thesis: hard requirement for newtx consistency |
| 7.5 | Write the plotted numerical vectors beside every published figure as CSV or NPZ | A PDF alone does not preserve numerical provenance. Your figure sinks already write both pickle and PDF; add the tabular form. | `[INFERRED]` | `[PRACTICE]` | JOSS: reproducibility of reported results |
| 7.6 | **Skip** pixel-baseline figure tests entirely, or restrict to 2-3 flagship figures | `pytest-mpl`'s own docs tell you to encode matplotlib and FreeType versions into the hash-library filename. Your theming rewrite invalidates every baseline anyway. Assert on `axes.shape` and the underlying arrays instead, as `arviz` does. | `[DOCUMENTED]` pytest-mpl | Contested | Not checked anywhere |

---

### Phase 8 - Documentation and diagrams (20-36 h)

Estimates are inflated over the reports' figures because you review every comma. AI drafting
roughly halves the writing time and roughly doubles the reviewing time; net saving is real but
smaller than it looks. Budget accordingly and do not let a deadline force you to ship unreviewed
generated prose, because generated documentation that is subtly wrong is worse than none - it
tells a reviewer the software does something it does not.

The three-stream partition from §3.1 applies here: `docs/` is the only stream that renders.
`spec/` and `docs/adr/` go in `exclude_patterns` and in the hatch `exclude` list.

| # | Change | Why | Evidence | Standard? | Venue check |
|---|---|---|---|---|---|
| 8.1 | Sphinx + `numpydoc` + `napoleon` + Read the Docs, structured by Diátaxis | `autodoc` is unmatched for docstring extraction, and it is the scientific Python default (NumPy, SciPy, astropy, lifelines). `numpydoc` over Google style because it is the ecosystem standard and has direct tooling consequences. Four of five reports agree; the Quarto dissent is wrong for API reference. | `[DOCUMENTED]` Sphinx autodoc, numpydoc | Ecosystem default | JOSS: "Functionality documentation" and "Example usage" checkboxes |
| 8.2 | Diátaxis mapping: tutorial = first qubit analysis walkthrough; **how-to = authoring a job** (first-class, as pytest documents writing plugins); reference = estimator signatures and assumptions; explanation = the statistical rationale | "How to write a job" is the load-bearing document for your stated goal of independent step use. | `[DOCUMENTED]` Diátaxis; adopted by NumPy and the Scientific Python development guide | Standard | JOSS: Example usage |
| 8.3 | Documentation as tests: `pytest --doctest-modules`, `sybil` for examples in Markdown/rST, `sphinx -W --nitpicky`, `interrogate` docstring-coverage gate | Highest-return testing available to you, because examples are part of the public scientific contract. Autodoc also depends on the package being importable, which is why Phase 1 gates this. | `[DOCUMENTED]` sybil, Sphinx | `[PRACTICE]` | JOSS: examples must actually work |
| 8.4 | Graphviz/DOT as canonical diagram source; Mermaid for GitHub rendering only | DOT is the only route to genuine font unification with newtx (via `dot2tex` to PGF/TikZ), `sphinx.ext.graphviz` is bundled and emits embedded PDF for LaTeX builds, and every source-analysis and workflow tool already emits DOT. Mermaid needs a headless Chromium via `mmdc` and has no font-matching story. | `[DOCUMENTED]` Graphviz font handling, Structurizr export formats | `[PRACTICE]` | Thesis: newtx matching |
| 8.5 | Draw C4 Level 1 (Context) and Level 2 (Container); **generate** Level 3 (Component) from the real import graph via `grimp` or `import-linter drawgraph`; skip Level 4 | C4's own documentation says "you don't need to use all 4 levels; only those that add value" and that context plus container suffice for most teams. The Container diagram is the highest-value single drawing: it makes the unpublishable-dataset boundary, the optional R subprocess and the path-resolution constraint legible as architecture rather than embarrassment. A class diagram of a package whose primary abstraction is the pure function would be a fabrication. | `[DOCUMENTED]` c4model.com | `[PRACTICE]` | JOSS: "Software design" required section |
| 8.6 | Emit `run.dot` and `run.json` beside each run's artifacts | The graph is data your runner already holds. Snakemake (`--dag`), Nextflow (`-with-dag`) and DVC (`dvc dag`) all do this. Lower risk than deriving a static architecture diagram from source, and it is the figure your thesis wants. | `[DOCUMENTED]` | `[PRACTICE]` | JOSS: provenance evidence |
| 8.7 | CI diagram-freshness check: regenerate, then `git diff --exit-code` | Cheap generic idiom for generated artifacts (Kubernetes `hack/verify-codegen.sh`). Note honestly: the reports searched and found no repository doing this *for diagrams specifically*. Treat as a cheap idiom, not an established standard. | `[PRACTICE]` for generated artifacts; `[INFERRED]` for diagrams | Not established | Not checked |
| 8.8 | Community guidelines: `CONTRIBUTING.md`, issue templates, support pathway, `CITATION.cff`, code of conduct | Currently absent; README is 225 characters. | `[DOCUMENTED]` | Standard | JOSS: explicit checkbox, three sub-items |

---

### Phase 9 - Paper artifacts (10-18 h)

| # | Change | Why | Evidence | Standard? | Venue check |
|---|---|---|---|---|---|
| 9.1 | 5-8 ADRs in `docs/adr/`: bespoke DAG runner, content-hashed identity folding git cleanliness, directories rather than class hierarchies for job families, subprocess R bridge rather than `rpy2`, package layout, flat namespace | These are exactly the questions an examiner asks. Note the honest caveat: the ADR literature measures template ergonomics and self-reported satisfaction, never defect rates. Adopt them because they are a byproduct-generating device for prose you must write anyway, not because they are proven to work. | `[PRACTICE]` Nygard / MADR | Widespread | JOSS: feeds the required "Software design" section directly |
| 9.2 | One-page workflow-engine comparison table: input interface, caching identity, transitive identity, provenance, in-memory composition, and why each missing property matters | JOSS explicitly accepts reimplementation "provided that they meet the criteria and cite prior similar work", and the State of the field rubric marks as **not acceptable**: "Ignores existing similar tools or fails to justify why contributing to existing projects wasn't appropriate". Your differentiator is not "a DAG"; it is typed in-memory step values plus content-derived identity plus transitive composite identity plus git-cleanliness gating. R `targets` is the strongest conceptual comparator even though it is R; Hamilton is the nearest Python analogue and lacks the code-hash-plus-clean-tree identity model. | `[DOCUMENTED]` JOSS State of the field rubric; `[DOCUMENTED]` `targets` design spec | `[PRACTICE]` | JOSS: required section, with an explicit "bad" rating for omitting it |
| 9.3 | Do **not** claim a performance benefit for the runner | 63 qubit-nights is not a big-data problem. The honest argument is type safety and object fidelity. An unmeasured speed claim trips "Performance: if there are any performance claims, have they been confirmed?" | `[DOCUMENTED]` JOSS checklist | n/a | JOSS: a checkbox you would otherwise fail for free |
| 9.4 | AI usage disclosure section, drafted from the preserved development ledger: agent/tool/version, date, task ID, spec or issue referenced, files changed, test results, human review action, commit hash | Required section. Non-disclosure is treated as an ethical breach with sanctions up to desk rejection and post-publication withdrawal. Retaining every prompt verbatim is unnecessary unless the thesis is a study of prompt engineering. | `[DOCUMENTED]` JOSS checklist and policies | Emerging, becoming universal | JOSS: required section. IEEE TQE: AI disclosure mandated. |
| 9.5 | Methods wording: "transitioned from conversation-driven development to a specification-backed agent workflow", not "specification-driven" | Claiming the spec pipeline improved code quality is an unsupported assertion a sharp examiner will catch. See §5 below. | `[MEASURED]` | n/a | Thesis defence |

---

## 3. Agent governance: permissions, commits, and where instructions live

### 3.1 The four places text can live, and the rule for each

Your instinct that `.claude/` clashes with the repo docs is right, but the diagnosis is slightly
off. The clash is not disclosure, it is **category**. `.claude/` currently holds prose that is not
configuration. Separate the four streams and the ugliness goes away.

| Stream | Lives in | Audience | Format | Ships in the wheel? | On the docs site? | Tracked in git? |
|---|---|---|---|---|---|---|
| **Package documentation** | `docs/` | someone using QUEBRA | MyST/rST, built by Sphinx | no | **yes** | yes |
| **Specification and decisions** | `spec/`, `docs/adr/` | you, the agent, the examiner | plain Markdown | no | **no** (`exclude_patterns`) | yes |
| **Agent operating instructions** | `AGENTS.md` at root, `.claude/` for config only | the agent | Markdown + JSON | no | no | yes |
| **Thesis and papers** | separate repo or Overleaf | examiners, reviewers | LaTeX | no | no | **no** - cites a git tag |

Four rulings that follow:

**Stop writing LaTeX for the agent.** The `.tex` worked, but not because it was LaTeX. It worked
because it had numbered requirements, explicit acceptance criteria, and a done-when clause per
item. That structure reproduces exactly in Markdown, the agent parses Markdown more reliably,
and you stop paying tokens for markup generation. Keep the structure, drop the backslashes.
`[INFERRED]`

**`.claude/` being public is fine, and is in fact load-bearing.** It is the primary evidence for
the JOSS AI usage disclosure section and for your thesis methods chapter. Kiro explicitly
recommends committing specification files alongside code so they become part of development
history `[DOCUMENTED]`. The fix for "it is ugly" is to make `.claude/` contain only
configuration - `settings.json`, `hooks/`, `agents/`, `commands/` - and move every prose-shaped
file to `spec/`.

**The versioning clash is solved by exclusion, not by hiding.** `spec/` and `docs/adr/` are
excluded from the Sphinx build via `exclude_patterns` and from the wheel via the hatch `exclude`
list. They are in git, they never appear on the published docs site, they never reach a user.
There is no version to clash because they are never rendered alongside the API reference.
`[INFERRED]`

**The thesis never imports from the repo.** It cites a git tag and, later, a Zenodo DOI. Figures
cross that boundary as PDF plus the companion data file from 7.5, never as source. This is the
only rule that keeps two documents with different release cadences from corrupting each other.
`[INFERRED]`

### 3.2 Permissions: how to stop pressing yes without handing over sudo

Verified against the current Claude Code permission model `[DOCUMENTED]`:

- There are three rule lists - `allow`, `ask`, `deny` - plus `defaultMode`.
- Rules are evaluated **deny, then ask, then allow**, and the first match wins. Rule specificity
  does not change the order, so a broad deny cannot carry allowlist exceptions.
- Settings load from four scopes: enterprise managed, project `.claude/settings.json`, local
  `.claude/settings.local.json` (gitignored), and user `~/.claude/settings.json`.
- **Permission rules merge across scopes rather than overriding.** A deny defined at any scope
  stays in force and a less authoritative scope cannot loosen it.
- A deny blocks the tool **even in `bypassPermissions` mode**.

**The answer to "the agent must not rewrite its own permissions" follows directly from that
list, and it is not a rule in the project file.** Put the hard denies in your **user-scope**
`~/.claude/settings.json`. The agent has write access to the repository; it does not have write
access to your home directory unless you add it to `additionalDirectories`. User-scope denies
merge in and no project file can loosen them. Project-scope denies on `.claude/**` are a useful
second layer, but they are the layer the agent could in principle edit, so they are not the
mechanism you rely on. `[DOCUMENTED]` + `[INFERRED]`

Second point worth internalising: **spoken constraints do not survive context compaction.**
"Don't commit" said in conversation disappears when the context is compacted. Hard limits belong
in `deny`. `[DOCUMENTED]`

**User scope (`~/.claude/settings.json`) - the layer the agent cannot reach:**

```json
{
  "permissions": {
    "deny": [
      "Bash(sudo:*)",
      "Bash(rm -rf:*)",
      "Bash(git commit:*)",
      "Bash(git add:*)",
      "Bash(git push:*)",
      "Bash(git tag:*)",
      "Bash(git reset:*)",
      "Bash(git checkout:*)",
      "Bash(pip install:*)",
      "Bash(curl:*)",
      "Bash(wget:*)",
      "Read(./.env)",
      "Read(./.env.*)",
      "Edit(**/.claude/settings.json)",
      "Edit(**/.claude/hooks/**)"
    ]
  }
}
```

**Project scope (`.claude/settings.json`, committed) - the working comfort layer:**

```json
{
  "permissions": {
    "defaultMode": "acceptEdits",
    "allow": [
      "Read(./**)",
      "Edit(src/quebra/**)",
      "Edit(tests/**)",
      "Edit(docs/**)",
      "Edit(spec/**)",
      "Bash(git status)",
      "Bash(git diff:*)",
      "Bash(git log:*)",
      "Bash(git mv:*)",
      "Bash(mv:*)",
      "Bash(mkdir:*)",
      "Bash(ruff:*)",
      "Bash(mypy:*)",
      "Bash(lint-imports:*)",
      "Bash(deptry:*)",
      "Bash(pytest:*)",
      "Bash(sphinx-build:*)"
    ],
    "ask": [
      "Bash(Rscript:*)",
      "Edit(pyproject.toml)",
      "Edit(.gitignore)",
      "Edit(.github/workflows/**)"
    ],
    "deny": [
      "Edit(data/**)",
      "Edit(output/**)"
    ]
  }
}
```

`defaultMode: "acceptEdits"` is the Cowork-like feel you want: file edits proceed silently, shell
stays gated. **Never use `--dangerously-skip-permissions`** on a machine with your credentials on
it; it is designed to short-circuit the settings file. `[DOCUMENTED]`

Note the two data rules. The agent **can read** everything under `data/`, including
`real_private/`, because it cannot debug a loader against data it cannot see. What it cannot do is
`Edit(data/**)` or `Edit(output/**)`: it may change the machinery, never the evidence. A dataset
or a materialised artifact modified by an agent mid-session is a silent scientific error with no
diff to catch it. The exfiltration risk is handled one layer up, by the `git push` and `curl`
denies in user scope, not by blinding the agent. `[INFERRED]`

### 3.3 Commit protocol: checkpoints, not commits

The agent never commits. That is enforced by the four `git` denies above, not by an instruction
in `AGENTS.md`, because an instruction is advisory and a deny is not.

Every phase below is cut into numbered checkpoints. At a checkpoint the agent stops and prints a
fixed banner:

```
CHECKPOINT 3.2 - logical job IDs replace path strings

  Changed:      src/quebra/core/registry.py (new), core/job.py, jobs/composite/*.py  (14 files)
  Gates run:    ruff 0 | mypy 0 | lint-imports 0 | pytest -m "not heavy" 0 (218 passed)
  Not done:     jobs/active/* still carry inline threshold ladders (checkpoint 3.3)
  Known risk:   cached identities for all composites will differ after this commit
  Suggested:    refactor(core): resolve sub-jobs by logical ID via registry

  Review `git diff` and commit if you see fit. I will not proceed until you say so.
```

Two details make this work rather than become theatre:

1. **A `Stop` hook runs the gate suite**, so the "Gates run" line is a result and not a claim.
   Hooks are the only mechanically-enforcing element in the Claude Code extension model; the
   `CLAUDE.md` text is static context. `[DOCUMENTED]`
2. **The banner names what was not done.** Agent-written summaries default to describing success.
   Requiring an explicit "not done" and "known risk" line is what makes the diff review fast,
   because you know where to look.

Put the banner format in `AGENTS.md` as a literal template. Symlink `CLAUDE.md -> AGENTS.md`;
Claude Code reads `CLAUDE.md`, and the symlink is the documented workaround for keeping one
canonical file that other agent runtimes can also read. `[DOCUMENTED]`

---

## 4. Toolkit and framework list

Adopt. Everything here is either an ecosystem default or has a named venue checkbox behind it.

| Function | Tool | Why this one |
|---|---|---|
| Packaging | `pyproject.toml` (PEP 621) + hatchling | Standard metadata; every other tool configures from it |
| Resource paths | `importlib.resources` | stdlib; the only correct answer for installed packages |
| Architecture enforcement | `import-linter` v2.13 | Native TOML, precise dependency-path reporting, `ignore_imports` ratchet for adopting against a violating codebase. Preferred over `tach` (rich public-interface system you do not need) and `pytest-archon` (rules as pytest tests, a style preference) |
| Lint and format | `ruff` (already on PostToolUse) | Cheapest deterministic gate. Honest note: no empirical evidence links linters to defect reduction in research software |
| Types | `mypy` or `pyright`, **on `core/` only** | Typed in-memory results are your core design claim; typing the module that carries them is coherent. Whole-repo typing at this scale is a churn generator |
| Dependency hygiene | `deptry` | Missing / unused / misclassified deps against `pyproject.toml` |
| Contracts | `icontract` | Hard computational invariants only. `deal` is the alternative; PEP 316 was rejected, so there is no stdlib answer |
| Tests | `pytest`, `pytest-cov`, `hypothesis` | Property-based tests for invariants, not for statistical validity |
| Test data | Committed synthetic fixtures + your own generator; `pooch` only if you later publish a benchmark | `pooch` is a registry for public downloads; you have none yet |
| Docs | Sphinx + `numpydoc` + `napoleon` + `pydata-sphinx-theme` + Read the Docs | Ecosystem default; unmatched autodoc |
| Doc honesty | `pytest --doctest-modules`, `sybil`, `sphinx -W --nitpicky`, `interrogate` | Turns examples into acceptance criteria |
| Doc structure | Diátaxis | Adopted by NumPy and the Scientific Python development guide |
| Diagrams (canonical) | Graphviz / DOT, `dot2tex` for newtx matching | Only route to genuine LaTeX font unification |
| Diagrams (GitHub) | Mermaid | Renders natively in issues and READMEs; nothing else |
| Diagram generation | `grimp`, `pydeps`, `import-linter drawgraph` | Filter aggressively; 105 modules unfiltered is unreadable |
| Architecture description | C4 levels 1-2 drawn, level 3 generated + ADRs | See 8.5 and 9.1 |
| Plots | matplotlib `subplot_mosaic`, `.mplstyle` targets, optional `SciencePlots` | See Phase 7 |
| CI | GitHub Actions; `r-lib/actions/setup-r` for the R job only | Gate the R job so the core matrix stays fast |
| Agent config | `AGENTS.md` at root, `CLAUDE.md` symlinked to it, `.claude/` hooks | Hooks are the only mechanically-enforcing element; `CLAUDE.md` and skills are advisory context |
| Decisions | `docs/adr/` (Nygard format) | See 9.1 |
| Provenance | Your existing runner + emitted `run.dot` / `run.json` | Keep; it is the defensible novelty |

**Explicitly do not adopt**, with the reason:

| Not adopting | Why |
|---|---|
| GitHub Spec Kit, Amazon Kiro, any formal spec pipeline | Every artifact is prompt text with zero enforcement. Evidence is null-to-negative (§3). ThoughtWorks places spec-driven development at Assess with the caution that the workflows "remain elaborate and opinionated" |
| Requirements traceability matrix | Measured benefit is for unfamiliar developers maintaining someone else's code with links pre-built for them. Demote to test naming (4.4) |
| Use cases, use case diagrams, activity diagrams, sequence diagrams (except optionally one), domain class diagram, deployment view, runtime view | They answer questions a library does not have: who are the actors, where does it deploy, what interactive state exists |
| Alloy | Your invariants (acyclicity, deterministic identity, no reuse under mismatch, valid topological order) are cheaper and more convincing as executable property tests |
| `rpy2` | Compiles against the R C-API; the least installable option. Every report except one says avoid |
| Containerisation as the default install path | Not a JOSS requirement. Build one optional reproducibility image *if* you route the methods paper through an artifact-badging venue |
| `pytest-mpl` baselines across all figures | FreeType and matplotlib version fragility; your theming rewrite invalidates them anyway |
| SBOM / AI-provenance metadata | Neither SPDX 3.0.1 nor CycloneDX 1.7 has a field meaning "an AI wrote this source file" |
| ISO/IEC 42001, ISO/IEC 5338 | They govern AI *systems*. Citing them for a Kaplan-Meier package is a category error a reviewer will notice |
| Rewriting onto Snakemake / Prefect / Dagster | Months of work to lose your defining feature. File-based engines would force serialisation of typed in-memory results, which is the exact constraint the runner exists to avoid |

---

## 5. The spec-first question, settled

All five reports independently concluded: skip the formal spec pipeline. That consensus is
evidence-backed rather than models agreeing with each other.

- **SlopCodeBench** (Orlanski et al., arXiv:2603.24755, March 2026, revised May 2026; verified
  live): the base "just-solve" prompt had the **best strict performance**, with plan-first
  averaging a **3.6 percentage point drop**; quality-aware prompting reduced *initial* erosion
  and verbosity but did **not** slow the degradation rate, improve pass rates, or reduce cost,
  and increased cost per checkpoint by **12.1% on average**. `[MEASURED]` Not yet peer reviewed;
  treat the point estimates as indicative.
- **Fu et al., "The First Prompt Counts the Most!"** (Proc. ACM SE, FSE 2025, DOI
  10.1145/3728947): over 95% of successfully implemented functionalities are achieved in the
  first round, suggesting LLMs struggle to use iteratively supplemented requirements.
  `[MEASURED]`
- **ClarifyGPT** (FSE 2024, DOI 10.1145/3660810): requirement *clarification* improved GPT-4
  Pass@1 in the evaluated setting. `[MEASURED]` This tests clarification, not the workflow.

**The defensible synthesis:** explicit requirements, decomposition and clarification improve
coding-agent *inputs*. The end-to-end claim that spec-first workflows reduce defects is
unresolved. Adopt `AGENTS.md`, `docs/adr/` and a light `spec/` for consequential changes because
they solve your **durability and context-loss** problem and because they are a thesis and
methods-paper asset. Do not claim they improved the code.

Spend the freed hours on deterministic gates instead. Each of `ruff`, `mypy`, `import-linter`,
`pytest`, `deptry`, `sphinx -W`, `sybil` runs in CI without a model and converts "please review"
into "the build is red." That is the actual answer to your credit-consumption problem: move
objective checks off the model, and ask the agent for semantic review only after they are green.

---

## 6. What the venues actually check, side by side

| Practice | Regulated ceiling | Academic RSE | Your target venues | Verdict |
|---|---|---|---|---|
| Formal requirements spec | ISO/IEC/IEEE 29148; DO-178C, IEC 62304 mandate it `[DOCUMENTED]` | Advocated, rarely formal | JOSS checks Statement of need, not an SRS `[DOCUMENTED]` | Adapt |
| Traceability matrix | Mandatory and audited | Associated with defect rate at team scale `[MEASURED]` | Not checked `[DOCUMENTED]` | Skip |
| Architecture description | ISO/IEC/IEEE 42010 | One package-level view is usually enough | JOSS **Software design** now a required section `[DOCUMENTED]` | Adopt-lite |
| ADRs | Change Control Boards | Advocated | Not explicit; resolves the "why your own runner" query | Adopt-lite |
| Coverage floor | MC/DC mandatory at DAL A | Weakly evidenced `[MEASURED]` | "Automated tests" checked, **no percentage** `[DOCUMENTED]` | Report, do not enforce |
| CI | Mandatory | Standard | Reviewer expectation | Adopt |
| Documentation | Formal manuals | Strongly advocated | JOSS explicit: install, examples, API, community `[DOCUMENTED]` | Adopt |
| Containerisation | Controlled images | Advocated | JOSS not required; ACM/ICSE strongly encourages, and warns installation over ~30 min is unlikely to be accepted `[DOCUMENTED]` | Adopt-lite, later |
| AI-use disclosure | Emerging, SBOM-adjacent | Emerging | JOSS **required section**; IEEE mandates disclosure `[DOCUMENTED]` | Adopt |
| Provenance capture | Cryptographic audit logs | FAIR4RS calls for it | Artifact evaluation considers reproducibility | Adopt (you already have it) |
| Public history | n/a | n/a | JOSS: **six months** `[DOCUMENTED]` | Hard gate |
| Community engagement | n/a | n/a | JOSS: "single author with no evidence of community engagement" is listed as not acceptable `[DOCUMENTED]` | Hard-ish gate |
| Artifact badging | n/a | ACM badges; uptake measured at 9.8% in an SE corpus `[MEASURED]` | IEEE TQE and ACM TQC: **encouraged, opt-in, not required** `[DOCUMENTED]` | Skip unless the methods paper routes through a badging venue |

The disagreement across the three columns is sharp and worth stating in the thesis: regulated
industry asks whether every safety requirement traces through design and verification; RSE asks
whether the software stays reusable and trustworthy; JOSS asks whether a real, installable,
tested, documented research software project exists and whether its contribution is clear. The
regulated column exists in the thesis only so an examiner can see how far from it you
legitimately sit.

---

## 7. Open questions, ordered by value

1. **When the six-month clock starts** - public repo date, or first commit? A pre-submission
   editorial query to JOSS resolves it and determines whether JOSS is a thesis deliverable.
2. **RESOLVED 2026-09-14. Does a Python Genest-Rémillard serial independence test already
   exist?** No maintained one was located. Searched: the implementations found are R only,
   `copula::serialIndepTest` and `copula::multSerialIndepTest`; neither `scipy.stats` nor
   `statsmodels` provides the test; `statsmodels` is not installed here. **Phase 6 does not
   collapse**, and SPEC 0007 proceeds on that answer. Two caveats bound this. A search is a
   negative result with a date on it, not a proof, so it is worth one minute of re-checking
   before anyone reimplements the statistic. And the answer does not change Phase 6's cost much
   either way, because the phase turned out to be about provenance and the truth of tracked
   documents rather than about the statistic: see `spec/specrboundary07.md` R7.0.
3. **What passes as a "research impact statement"** for a tool with one dataset. The section is
   new; read three to five JOSS papers accepted after January 2026.
4. **Licence.** GPL-3.0 versus permissive given vendor-adoption goals. One-way door, untouched by
   every report.
5. **Whether Metriq integration or the Unitary microgrant generates the external-engagement
   evidence** the new checklist requires.
6. **Whether any repository does CI diagram-freshness checking** (8.7). Searched, not found.
   Cheap idiom regardless.
7. **Does read-level LONG-MEMORY noise reach the durations, where AR(1) does not?**
   `jobs/bench/arms.py:397` measures that read-level AR(1) induces no positive
   duration-level dependence at any rho, reaching only -0.086 at rho = 0.99, because the
   level crossings of a stationary Gaussian process regenerate. That is a SHORT-memory
   argument and it does not obviously carry to 1/f. If it does not, a 1/f arm would produce
   duration dependence the battery could be scored against. It is a bench arm and a
   finding, not a test oracle: its induced duration dependence would have unknown
   magnitude, whereas Arm E already supplies duration-level dependence with analytic ground
   truth. SPEC 0006 R6.4 declined to build it for that reason.

8. **The shared median-spacing defect in `analyzers/tlf.py:134` and
   `analyzers/allan.py:20-25`.** Both take a median of read spacings and treat it as the
   sample period: `tlf` multiplies run lengths by it to get dwells in seconds, `allan`
   passes `rate = 1/dt` to `allantools`. On a gapped record neither quantity is elapsed
   time. `allan` at least filters to positive diffs; `tlf` does not. Fixing either changes
   a shipped figure, so it is a deliberate decision rather than a passing repair.

9. **`allan.py:141-142` divides `fractional_adev` by `carrier_hz`, the mean of
   `qubit_frequency_hz`, while mode `"fractional"` at `:63` divides by `f0_hz`, the mean of
   `delta_hz`.** Two different denominators under one name, which is the unit-suffix bug
   class `AGENTS.md` section 3 calls the costly one here. Recorded rather than fixed
   because the choice changes a published figure.

10. **Which check verdicts, if any, should GATE a Kaplan-Meier band rather than annotate it.**
   SPEC 0008 built the annotation and deliberately built no gate: collapsing a ledger to one
   licence is not decidable without the device survey analysed, and every candidate rule
   measured so far revokes on every real record for reasons unrelated to serial independence
   (C3 is `not computed` without R; in-spec C1/C2 are singular at roughly 73% of replicates).
   Deciding it needs the survey read, a stated rule for which clock and which rows license, and
   a stated mapping for `underpowered` and `not computed`. It is the ADR Phase 9.1 will have to
   write either way, so it is named here rather than left implicit in a passing test.

---

11. **A container image is the only honest way to make the R dependency installable, and
   nothing installs it today.** Raised while settling SPEC 0007 R7.5. C3 needs an R
   interpreter plus the CRAN `copula` package, and neither is a Python distribution, so
   `pip install quebra[r]` cannot ever work: the empty extra that implied otherwise was
   deleted. `rpy2` does not help - it is pip-installable but links against an R you must
   already have, and `analyzers/checks/c3_serial_copula.py` refuses it by name because that
   turns a missing R into an import-time failure of the whole package. So today the only
   install route is prose in the README telling a reader to run `apt`/`dnf`/`brew` and then
   `install.packages("copula")`, and the only thing verifying that route is a human reading it.
   An image would make the R side reproducible, pin the `copula` version the p-values depend
   on (which `CheckLedger.r_version` now records but nobody can currently reproduce), and give
   the `r` tier somewhere to run other than one developer laptop.
   **What this already collides with.** Section 6 rates containerisation "Adopt-lite, later"
   and section 4 scopes it as one optional image *if* the methods paper routes through an
   artifact-badging venue. It also notes ACM/ICSE warn that an install over roughly 30 minutes
   is unlikely to be accepted, and building `copula` from source is the part at risk there.
   **What would settle it:** one timed build of an image carrying R plus `copula`, which is the
   same measurement the deferred `r-lib/actions/setup-r` CI job in `specrboundary07.md` needs.
   Do them together or not at all. Not scheduled, and not part of Phase 6.

   **The scope grew during Phase 7, and the case is now stronger than the R install route.**
   Two findings there share one root cause: nothing pins the versions a figure was produced
   under, and nothing records them.
   `pyproject.toml` declares LOWER BOUNDS ONLY (`matplotlib>=3.9`, and the same shape for every
   runtime dependency), so a clean tree at a given commit does NOT fix what is installed: two
   people on one commit can hold matplotlib 3.9 and 3.10.9, and CI resolves newest on every
   leg. Identity folds the job file, the module hashes and a clean-tree check; it does not fold
   the environment, and `core/identity.py:17` says so deliberately.
   Phase 7 then removed matplotlib's version string from every written PDF, because leaving it
   in made figure bytes move on each release even when the drawing did not. That was right for
   the bytes, and it deleted the only record of which renderer produced a figure, because
   `provenance.build_prov_record` carries no package version at all.

   **The proposal, as ONE piece of work.** An image fixes the versions; provenance records
   them. Neither half is worth much alone: recording an unpinned version documents drift
   without preventing it, and pinning without recording leaves a reader unable to tell which
   pin applied.

   **The obstacle to name up front.** `check_ledger.py:20` states that the provenance record's
   "schema is closed", which is why the R version went onto the materialized `CheckLedger`
   artifact rather than into the `.prov.json`. Putting package versions into provenance
   therefore REOPENS a decision Phase 6 made, and it should be reopened explicitly rather than
   by accident. A figure has no typed artifact of its own to carry them, which is the asymmetry
   that makes the R workaround unavailable here.

## 7B. Remaining statistical tools and checks

What is implemented, what is calibrated, and what is neither. Calibration here means a bench
cell: a measured size, and where it matters a measured power, against a generating arm with
known truth. A check with no bench cell produces numbers whose LEVEL is unmeasured, so neither
its rejections nor its non-rejections are evidence at alpha.

Which instruments exist is `AGENTS.md`'s "Implemented today / Not implemented" list and
`battery.ROW_KEYS`; neither is restated here, because a second copy is what goes stale. What
this section adds is the CALIBRATION status, which is decidable from the bench tables: the
five checks in `battery.ROW_KEYS` have size and where relevant power rows, and C3 has none.
`tests/test_check_ledger.py::test_every_surveyed_check_is_benched_or_declared_uncalibrated`
derives that split rather than asserting it, so it cannot drift.

### C3 is the one shipped check with no calibration, and it now draws numbers

Measured 2026-09-27, the full survey with `INCLUDE_C3 = True` and `C3_N_NULL_SIM = 200`, 34
datasets, 30 minutes of wall clock: **248 of 680 C3 cells produced a p-value, and 60 of those
read `fail`** (46 in-spec, 14 calendar). Before R was installed on this machine every one of
those cells read `not computed`. `battery.ROW_KEYS` has nine entries and none is C3; both bench
tables carry zero C3 rows.

A guard now exists so this cannot be lost: `check_ledger._verdict` names the uncalibrated state
in the note of every such rejection, and `tests/test_check_ledger.py` derives the uncalibrated
set from `size_table.csv` and fails if it stops matching. Benching C3 is what removes them.

**What is NOT settled, and needs thought before anyone builds it.**

1. **Which calibration C3 should even have.** The other five are scored on a permutation or an
   asymptotic null that the bench can resample. C3's null is simulated ON THE R SIDE from
   `(n, lag.max, N, seed)` by `copula::serialIndepTestSim`, so it shares nothing with the
   `PermutationSet` that pairs every other row's Monte Carlo noise. A bench cell for C3 is
   therefore not the same object as a bench cell for C1, and `report.py`'s generated prose
   asserting that every row of a replicate reads the same permutation set becomes false the
   moment a C3 row exists.
2. **Whether the R bridge changes what a calibration means.** The p-value depends on the
   installed `copula` version and on `N`. `CheckLedger` now records both, but a bench cell
   measured under one `copula` and read under another is an open question nobody has posed.
3. **Cost, and the multi-process hole.** C3 returns `None` for `m > 1`, which is most of the
   bench's censoring grid at n >= 50, so a C3 bench can only be calibrated at censoring 0 for
   the larger n. That removes cells from inside the promotion report's own envelope and changes
   the Bonferroni threshold for the C3 row relative to every other row.
4. **A blocking defect, before any of the above.** `jobs/bench/runner.py:207-213` counts a
   `None` p-value into `seen` but never into `hits`, and C3 is the only check that can return
   `None`. A registered C3 would report a rejection rate of 0 on every multi-segment cell and
   `calibration_summary.py` would publish it as MISCALIBRATED on cells where it never ran. Two
   more in the same file: `runner.py:169`, `:181` and `:200` catch only `(ValueError, KeyError)`
   while the bridge raises `RuntimeError` and `TimeoutExpired`, so one bridge failure aborts a
   multi-hour run; and `report.py:665-673` asserts the shared-permutation-set pairing that item
   1 makes false.

**Scope.** This is one pass, not four: the unimplemented estimators and C3's calibration share a
question (what is the oracle, and what does the bench cell mean) and the same blocking defect
list. It is not Phase 6, which stops at making the R boundary honest about what it does and does
not know.

### The check path reads a record as gap-free when its gap list is omitted

Found at the Phase 7 review and not fixed there. The panel builders now require the carve's
`gap_spans_s` (`specpresentation09.md` R9.4c). The independence-check path does not:
`CheckLedgerInputs.gap_spans_s` defaults to `None` (`analyzers/check_ledger.py:176`),
`make_inputs_from_windows` reads it through `getattr(window_result, "diagnostics", {}) or {}` and
`.get("gap_spans_s")` (`:225-231`), and `_segment_starts` and `segments_from_windows` take it as
optional (`analyzers/checks/_multiprocess.py:160`, `:202`). Omitted, segmentation falls back to
birth types alone, which the module's own docstring says cannot see a gap flanked by out-of-spec
reads: two renewal processes separated by unobserved hours are scored as one segment. What
happens next depends on the clock. On the in-spec clock, a censored window left mid-segment
raises, and the ledger reports that cell as not computed. On the calendar clock the merge is
silent. On both, a merge with no censored window between the two processes is silent.

Production is not affected today: the ledger is built from a `WindowsResult`, and `windows.run`
always records the key. The defect is the API, the same one R9.4c closed for the panel.

**Why it is here rather than fixed.** `jobs/bench/arms.py:314` calls `segments_from_windows`
without a gap list, on a regular synthetic grid that has no gaps, so requiring the argument
changes the bench's calling code, and the bench is the calibration evidence.

**Open question.** Require it everywhere and have the bench pass `[]` explicitly, which should
change no bench number because its grid has no gaps, but that has to be shown by re-running a
bench cell rather than asserted; or keep the fallback and name the callers allowed to use it. No
source located for which is right; it is a design decision.

---

## 8. Three things worth being told bluntly

**One.** You are over-engineering the development *process* and under-engineering the package
*boundary*. `core/job.py` importing `loaders`, `transforms` and `plots` is more dangerous than
the absence of a spec pipeline. A package survives without a beautiful architecture document; it
does not survive arbitrary dependency inversion once anyone else imports it.

**Two.** You have internalised "a reviewer may ask why I wrote my own runner" as a major threat
and "no `pyproject.toml`, tests and docs gitignored" as minor. It is backwards. The runner has a
real, articulable justification that an ADR plus one table fully neutralises. The `.gitignore`
means that from the outside the repo appears to have no tests and no documentation at all, and
the public-history clock has not started.

**Three.** The clock item cannot be bought back with effort, and it is the only item on the list
with that property. Everything else in this plan is 90 to 150 hours of work you control. Phase 0
is one day and a supervisor conversation, and it determines whether JOSS is reachable at all.
Do it this week.

---

## 9. Possible next ideas that came up while implementing

Not requirements, not scheduled, and not costed. Each is here because building something
else made it visible, and each is written so it can be reasoned about separately from the
session that raised it. Nothing below has been decided, and none of it is a commitment.

### 9.1 An assumption record for non-informative censoring

**Where it came from.** SPEC 0008 R8.2 built assumption records and shipped exactly one,
`a1_renewal_durations`, diagnosed by the six checks. Writing it made the neighbouring gap
obvious: the Kaplan-Meier band leans at least as hard on the censoring being non-informative,
and no shipped check tests that.

**What the assumption says.** A window right-censored by a read gap or by the scan ending
must carry no information about how much longer it would have lasted. If long-lived windows
are preferentially the ones still running when observation stops, censoring is informative
and the survival curve overstates survival.

**Why it is not simply a missing check.** The question is currently handled OUTSIDE this
framework, as a property of the measurement scheme rather than of the record, and that is a
deliberate position rather than an oversight. Writing it into an `Assumption` record would
mean stating the justification in the package, and the justification is a scheme-level
argument that has not been written down yet.

**What would have to be settled first.**

- On what basis censoring is taken to be non-informative here: which feature of the scan
  scheme supports it, and under what conditions it would stop holding.
- Whether anything data-side can be said at all, or only bounded. The standard difficulty is
  that right-censoring shows only which of the lifetime and the censoring time came first, so
  a mechanism that is informative and one that is not can produce the same observed
  distribution. Confirm this against a source before relying on it; it is stated here as the
  reason to think carefully, not as a settled result.
- Whether the useful artifact is a diagnostic at all, or a sensitivity analysis. Item 5.2's
  generator exposes the censoring process independently, so a known informative mechanism
  could be imposed and the resulting bias in the curve measured. That reports the COST of a
  violation rather than detecting one, and it may be the more honest deliverable.

**What is already in place if it is taken up.** `analyzers/assumptions.py` is shaped so a new
record is one literal plus one line in `_RECORDS`: `diagnostic_checks` lives on the record, so
there is no second registry to update, the id shape is enforced, and the `undetected`
disposition already exists for an assumption the package cannot see fail. A record with an
empty `diagnostic_checks` is constructible today and is refused only if it also claims to
report or raise something.

**Related, and separable.** `spec/specvalidity08.md` R8.1 measured a second thing worth its
own decision: a zero-duration window ending a block makes `tau == T_N` exactly on the calendar
clock, which drops C1 and C2 for a whole record. That is a carve and truncation question, not
a censoring one, and it is recorded there rather than here.

### 9.2 `independence_survey.py` has quietly become the vocabulary module

**Where it came from.** A naming review during SPEC 0008 CHECKPOINT 8.4b. The review's verdict
on the new code was that only a word was wrong, but it flagged this as a real structural note
about EXISTING code, and it is recorded rather than acted on because it touches shipped
modules.

**What it is.** `analyzers/independence_survey.py` is named as, and documented as, the
output-builder for one figure family. It also defines the canonical row-key vocabulary that
several unrelated things now import:

| Constant | line | what it is |
|---|---|---|
| `C3_KEY` | 48 | the one row key `battery.ROW_KEYS` omits |
| `SURVEY_KEYS` | 51 | `ROW_KEYS` + `C3_KEY`, the spanning vocabulary |
| `CHECK_LABELS` | 55 | row key to display label |
| `CHECK_NULL` | 79 | row key to the null it tests, in prose |
| `VERDICT_ORDER` | 120 | verdict draw order |

**Honest attribution: SPEC 0008 caused this, it did not find it.** Before that phase the only
importer of the vocabulary was `plots/independence_survey_plot.py`, that module's own figure,
which is a normal builder-to-plot pairing and not a smell. Adding
`analyzers/check_selection.py` (imports `C3_KEY`, `SURVEY_KEYS`) and
`analyzers/check_outcome.py` (imports `CHECK_LABELS`, `CHECK_NULL`) is what made a figure
module the shared home for a vocabulary two other subsystems depend on.

**Size, measured.** The five constants span roughly 75 of the file's 333 lines. Importers
today: two analyzers, one plot, one job, two test modules. A move is mechanical - no logic
changes - but it touches every one of those import sites.

**Is it already planned? No.** Nothing in sections 1 to 8 covers it. The nearest items are
2.3 and 3.2, and both are about reorganising `jobs/`, not `analyzers/`. So if it is worth
doing it needs its own home; it will not arrive as a side effect of a later phase.

**Options, none chosen.**

1. **Leave it.** The coupling is real but static, the arch contract is intra-layer so nothing
   is violated, and the names are correct wherever they are read. Cost: a reader looking for
   the row-key vocabulary has to know to open a figure builder.
2. **Move the five constants into `analyzers/check_selection.py`**, which already owns `RowKey`
   and `ALL_KEYS` and is named for the job. `independence_survey` then imports them back. Small
   and mechanical, but it moves shipped constants and every cached identity that reaches them.
3. **A dedicated `analyzers/checks/vocabulary.py`.** Cleanest by name, and it puts the row-key
   vocabulary beside the checks it names rather than beside either consumer. Largest diff.

**When.** Option 2 or 3 is cheapest immediately after a phase that already moves identities,
and most expensive just before a figure is promoted for publication. It is not urgent: nothing
is wrong today, and the note exists so the decision is made deliberately rather than by a
future agent noticing the same thing and refactoring unasked.
