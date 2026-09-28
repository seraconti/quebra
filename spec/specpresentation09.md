# SPEC 0009 - The presentation layer

**Phase:** 7
**Budget:** 10-16 h, below `quebraplan.md`'s 10-20 h. R9.0 declines 7.1 and 7.2 outright and finds
7.4 and 7.5 built on measurements that do not hold, which frees the hours the plan allocated to a
render-contract rewrite. They go instead to three defects the plan does not contain: figure
artifacts that no longer load, a curve drawn across read gaps, and a font the poster target needs
that nothing declares. The figure lands below the plan's ceiling rather than at it because a cold
review moved two items out: the durability fix is a module-alias remap rather than a new
serialisation format across fifteen artifact types, and the vocabulary sweep is deferred to its own
spec. Both are recorded in Not done with the measurement that moved them.
**Depends on:** SPEC 0006 complete and committed (`278a2fd`). Independent of SPEC 0007; the two
touch no common file.
**Feeds:** `quebraplan.md` 8.2, whose reference material describes the figures, and 9.4.
**Blocks:** regenerating thesis and poster figures. `quebraplan.md` §1 says "Do it before you
regenerate thesis figures, not after", and R9.1, R9.2, R9.4a and R9.4b are the requirements that
statement is about: each changes what a regenerated figure contains or how it renders elsewhere.
**Reference:** `spec/quebraplan.md` §2, Phase 7, items 7.1 through 7.6.

**Numbering.** 0009 rather than 0007, because 0007 is reserved for the R boundary by four sites
(`specinstallabity02.md:68`, `spectests06.md` twice, `specvalidity08.md:13`) and is taken by
`specrboundary07.md`. 0008 is the validity contract. 0009 is the next free number, stated here so
a reader seeing a Phase 7 spec numbered 0009 does not assume a missing document.

Rationale lives in `spec/quebraplan.md`. This file contains requirements and acceptance criteria
only, plus the measurements in R9.0 that change what the phase should do.

The governing constraint: **a figure is evidence only while the numbers behind it can still be
read and the drawing does not assert something the data does not. Both are currently false, and
neither is what Phase 7 was scoped to fix.**

---

## Preconditions

Baseline on `dev` at `278a2fd`, measured rather than assumed:

- `pytest --collect-only -q` gives **592 tests**.
- matplotlib **3.10.9** in the project venv; 3.10.8 is a user-site copy. An earlier draft said
  3.10.8, having measured under the wrong interpreter. Measurements below are venv measurements.
  `pyproject.toml:33` declares `matplotlib>=3.9` with no upper pin, and CI resolves newest across
  three Python versions.
- **20 concrete renderers**: 16 plot classes under `plots/` and 4 panel classes under `panels/`,
  each implementing `build_matplotlib(result, style="default") -> plt.Figure`. Nothing accepts an
  `ax`. **Five draw into a single Axes** (`km_survival_plot.py:114`, `mtbc_hist_plot.py:64`,
  `panels/comparison.py:51`, `instrument_validation_plot.py:76`, `calibration_plot.py:260`); the
  other 15 create multiple axes, 7 of them with a data-dependent count.
- **3 of those 20** are invoked by any test. 17 are never constructed by the suite.
- `plots/theme.py` is 469 lines: `RCPARAMS` at `:98`, `apply_rcparams` at `:367`, `style_context`
  at `:383`, `_TARGET_FOR_STYLE` at `:195`, `scaled_text` / `_TEXT_SCALE` at `:451-454`, a dead
  commented poster block at `:138-159` and the live one at `:160-190`.
- `plots/targets.py` registers four targets: `static` (PDF 300 dpi), `academic` (PDF 600 dpi),
  `poster` (PNG 600 dpi), `interactive` (plotly HTML). **No renderer implements a plotly backend**:
  11 `build_plotly` overrides all raise, plus the base at `plots/base.py:18-19`, and five renderers
  do not override at all.
- Zero repository hits for `pytest-mpl`, `pytest_mpl`, `mpl_image_compare`, `subplot_mosaic`,
  `pdf.fonttype`, `ps.fonttype`, `svg.hashsalt`, and savefig `metadata=`.
- `tests/test_style_baseline.py` pins hardcoded style literals at `BASELINE = 17` and asserts it
  **both ways**: `total <= BASELINE` at `:55` and `total == BASELINE` at `:66`. Measured: `panels/`
  10 and `plots/` minus `theme.py` 7, summing to exactly 17.
- **Figure artifacts that no longer load.** Measured three ways, because the loosest denominator
  flatters the case: over every run directory, 356 of 388 top-level pickles fail (91.8%); over the
  latest run of each of the 51 job names, 40 of 46 (87.0%); over the latest run of each of the
  **13 jobs that still exist**, **21 of 27 (77.8%)**. The defect is therefore not an artifact of
  development churn. An earlier draft quoted 664 of 696 (95.4%) from a recursive count.
- `output/` holds 244 run directories across 51 job names and roughly 320 MB.

---

## R9.0 - Premises in the plan that do not survive contact with the code

### R9.0.1 7.1's axes-level contract does not fit these renderers

7.1 asks every plot to "accept an optional `ax`, draw into it, return the `Axes`, never call
`show()` or `savefig()` inside library code".

The `savefig` half is **already satisfied**: `savefig` appears only in the three registered render
targets in `plots/targets.py`, which are the sink layer, and `plt.show` appears nowhere.

The `ax` half does not fit. **15 of the 20 implementations create multiple axes in one call and
cannot return an `Axes` at all**, 7 of those with a data-dependent axes count. The other five draw
into a single Axes but all five also do figure-level work (`tight_layout`, and in two cases
`fig.patch.set_facecolor`), so converting them means splitting that out, and no consumer exists for
the conversion. The cited exemplars do not have this shape: seaborn's axes-level functions and
`lifelines`' `plot(ax=...)` each draw one artist family into one axes. These are composites, and the
correct analogue is seaborn's **figure-level** functions, which do not take `ax`.

An earlier draft of this section said "exactly one class converts cleanly". That was wrong by five
and is corrected above. The decline stands on the 15 composites.

The premise that the drawing code is not composable is also already false at the panel level: 19
signatures take `ax: plt.Axes`, 18 of them in `panels/`, and the panels are internally sets of
`_draw_*(ax, data)` functions, which is what `AGENTS.md` §3 requires. **Declined.**

### R9.0.2 7.2's mosaic cannot express the layout

7.2 wants `plt.subplot_mosaic` and "a thin composer". Measured against the venv's matplotlib
3.10.9: `subplot_mosaic` applies `width_ratios`, `height_ratios` and `gridspec_kw` to the **outer**
layout only, its own docstring saying "In the case of nested layouts, this argument applies only to
the outer layout"; it builds nested levels with a bare `subgridspec(rows, cols)` forwarding no
kwargs; and `sharex`/`sharey` are figure-wide booleans. `panels/within_calibration.py` needs all of
it: `:98` outer `height_ratios`, `:105-107` a second level with `width_ratios` and `wspace`,
`:122-124` a third with `height_ratios` and `hspace`, `:109` `sharey`, `:130` `sharex`, and a
data-dependent outer row count at `:89-91`. **Declined.**

### R9.0.3 7.3's premise about `theme.py` is false

7.3 says `.mplstyle` files replace `plots/theme.py`, "which is a qubit colour lookup rather than a
theme", and that `plots/targets.py` "already has the `["static", "academic"]` machinery for
'poster' to slot into". Both halves are false. `theme.py` is a theme by every structural measure in
the Preconditions, and the qubit colour lookup is one of roughly twenty things in it. `poster` is
already a registered fourth target, used by `jobs/active/mtbc_q6.py:71` and
`jobs/active/km_poster_6d2s.py:143`.

What is actually missing is smaller: the live poster block at `theme.py:161` hardcodes
`"font.family": "Roboto"` with no fallback stack and no check that the face exists, and `:138-159`
is a dead commented duplicate. **Reshaped** into R9.2.

### R9.0.4 7.4 bundles four settings under one rationale, and it holds for none of them as stated

7.4 asks for `pdf.fonttype: 42`, `ps.fonttype: 42`, deterministic metadata and `svg.hashsalt`,
because they "force font embedding, avoid Overleaf compilation errors, make figure files diffable
and hashable". Measured:

- Output is **already deterministic** under both `fonttype` 3 and 42 across processes under
  randomized `PYTHONHASHSEED`, with a stable font subset tag. For PDF the sole variable field is
  `/CreationDate`.
- **`svg.hashsalt` is inert.** No SVG target exists.
- **Fonts are already embedded and subset** under the default `fonttype` 3. What 42 changes is the
  Unicode mapping, Type 3 with a custom encoding becoming CID TrueType Identity-H. That is a text
  extraction property, not a determinism one.
- **The "diffable and hashable" claim survives only within one matplotlib version.** Every PDF
  carries `Creator: Matplotlib v{version}` and `Producer: ... v{version}`, and the poster PNG
  carries a `Software` chunk naming the version. Against an unpinned `matplotlib>=3.9` with CI
  resolving newest, figure bytes move on every matplotlib release regardless of `CreationDate`.

**Reshaped** into R9.1, which treats determinism and font portability as two things rather than one.

### R9.0.5 7.5's premise that the numbers are already written is false in three ways

7.5 says "your figure sinks already write both pickle and PDF; add the tabular form".

First, **the pickles do not preserve the numbers**: 21 of 27 fail on the tightest denominator. A
pickle stores the fully qualified module path of its class, so a rename or a move makes it
unreadable. `StaleArtifactGuard` in `core/_artifact_guard.py` additionally raises by design
when a pickle of one of its subclasses lacks a field the class has since gained.

Second, **"beside every published figure" has no addressee**. A published figure is one promoted by
`scripts/promote_run.py`, and promotion copies only `provenance/*.prov.json` and `*.prov.md`. It
never copies the PDF and never copies the pickle. Writing a file into the run directory changes what
`published/` contains by zero bytes, and `published/` has never been created.

Third, **the sink's unit is the artifact node, not the figure**.
`jobs/composite/independence_survey.py:282-295` renders 11 figures off a single node.

**Reshaped** into R9.3, which repairs durability rather than adding a second format beside a broken
one. The tabular form is separated and deferred.

### R9.0.6 7.6 recommends a state that already holds, and misses the coupling that matters

There are no pixel baselines to invalidate: `pytest-mpl` is not installed, declared or imported. The
recommendation costs nothing and is adopted by default.

What 7.6 does not mention is `tests/test_style_baseline.py`'s two-sided pin, which couples any work
that moves a literal to a test edit, and its blind spot on the dict form. **Reshaped** into R9.7.

---

## R9.1 - Make figure output byte-reproducible within a matplotlib version

**What.** Pass `metadata={"CreationDate": None}` at the two **PDF** `savefig` calls in
`plots/targets.py:52` and `:63`. Pass `metadata={"Software": None}` at the **PNG** call at `:80`.
Pin `Creator` and `Producer` to constants on the PDF calls.

**Why.** R9.0.4. The keys differ by format and an earlier draft prescribed `CreationDate` for all
three: PNG has no `CreationDate`, so on the poster that would have been a **silent no-op**.
matplotlib seeds PNG metadata with a `Software` key naming its own version, which is the poster's
real variable field. Pinning `Creator`/`Producer` is what extends byte-equality across matplotlib
versions rather than only within one.

**Acceptance.** A test renders each of the three registered matplotlib targets twice in separate
processes under different `PYTHONHASHSEED` and asserts identical sha256 per target. MUTATION: revert
one of the three metadata arguments; the test goes red naming that target. Run it and report both.
The poster must be included, because it is the target the wrong key would have silently skipped.

**Budget.** Collect +1 to +2. Files 2.

**Known risk.** None to figure content. `bbox_inches="tight"` and dpi are untouched.

---

## R9.2 - The poster target must not depend silently on an undeclared font

**What.** `theme.py:161` sets `"font.family": "Roboto"`, a user font under `~/.local/share/fonts`
that is not in `pyproject.toml`, cannot be pip-installed and is not provisioned by CI, so the poster
renders with different typography and different bytes on any other machine, silently. Give the
poster block an explicit fallback stack and make a missing face **loud**. Delete the dead commented
block at `theme.py:138-159`.

**Why.** `AGENTS.md` §3: errors are raised, not swallowed, because a silent fallback yields a
wrong-but-plausible result. A poster whose typography depends on which laptop rendered it is that.

**Acceptance.** A test asserts the configured poster family resolves to itself. Two traps must be
avoided and the acceptance names both. It must capture the `matplotlib.font_manager` **logger**, not
warnings, because the miss is `_log.warning` and `pytest.warns` would never fire. And `findfont` is
memoised behind an `lru_cache`, so a log-capture test is order-dependent and can be inert if an
earlier test warmed the cache: clear the cache in the fixture, or call
`findfont(..., fallback_to_default=False)` and assert on the raised error instead. MUTATION: set the
family to a name no system has; the test goes red. Run it after a full-suite run, not alone, to show
the cache clause works.

**Budget.** Collect +1 to +2. Files 2.

**Decision owed to Sera.** Vendor the font, declare it as an environment requirement, or drop to a
guaranteed face. Making the fallback loud is the requirement; which of the three follows is a call
about how the poster is produced.

---

## R9.3 - A figure's numbers must survive a rename

**What.** Give the unpickler a module-alias map so artifacts written under an old module path
resolve to the current class: a custom `Unpickler.find_class`, or pickling under a stable registry
name rather than the defining module path. One module, plus its registration at the read sites.

**Why.** R9.0.5, and it is chosen over writing a second file beside the pickle for three measured
reasons. It **repairs the existing artifacts**, where a new sidecar format would leave them dead
and help only future runs.

**AMENDED during implementation, on a measurement.** This paragraph first said the remap repairs
"all 21 of the live-job failures". Measured on the latest run of every job file under
`jobs/active/` and `jobs/composite/`, every `*.pkl` directly in it (the compare job's
`subjobs_output/` holds 4 more, not counted): 29 pickles, 8 load without the alias,
24 load through it, and **20 load complete, so the alias repairs 12**. Of the 16 it makes
loadable, 4 are stale. `KaplanMeierComparison` in both `km_poster_6d2s` pickles lacks the four
check-outcome fields and would read as "NOT ASSESSED"; `TLFResult` in both `ramsey_q1_100423` TLF
pickles lacks `fit_failed` and would crash `plots/tlf_plot.py:88`. Neither class carries
`StaleArtifactGuard`, so `load_artifact` checks every dataclass it can reach and refuses all four.

The 5 it cannot load are not a gap in the alias table. Every one references a class that was
RENAMED AND RESTRUCTURED, not moved: `NonRepairablePanelData`, `RepairablePanelData` and
`CompareNonRepairableData`, from the vocabulary AGENTS.md section 5 retired, the first of them
before the within-calibration panel was split into three bands. An alias maps a moved package to
the same class; mapping these onto their successors would load old field values into a new shape.
They stay unreadable with an error that says so, and the tests pin that the alias refuses to
rename a class. Re-running those jobs is the only honest repair.

The denominator moved from the Preconditions' 27 (6 loading) to 29 (8 loading) because of one
run. At `278a2fd` the latest `independence_survey` run was `a21fd1_20260830_110425`, which holds
no pickle; the run written since, `d6fe38_20260927_095456`, adds 2 that load without the alias.
That accounts for both moves, and the 21 failures are the same 21: the same rule with runs cut at
`278a2fd` gives 27, 6 and 21.

**It does not repair the composite transport**, which an earlier draft gave as the second
reason. `core/runner.py` reads a sub-job's pickle back in `_locate_artifact`, now through
`load_artifact`, but its reuse gate admits only a run at the current commit on a clean tree, and
such a run is written under current paths, so no pre-move artifact reaches that read in
production. The runner uses the loader so that every read gets the completeness check; the test
that drives it fakes a provenance record and proves the wiring only. The choice therefore rests on
the first reason, artifacts opened by hand, and the third: it is one module against roughly
fifteen artifact types, one of which is a plain `dict` (`plots/tlf_plot.py:12-18`) with nowhere
to hang a method.

**Acceptance.** A test writes an artifact, renames the defining module, and loads it back through
the remap. MUTATION: remove the alias entry; the load fails. Run both. Separately, re-run the
Preconditions measurement and report how many of the 27 live-job pickles load after the change; the
number is the requirement's whole justification and must be stated, not asserted. (29 at
re-measurement; see the amendment.)

**Budget.** Collect +1 to +2. Files 2 to 3. Measured: collect +4, files 3, at the 2x line.

**Known risk.** The alias table must be extended on every future package move. That is one line
per move, and omitting it fails loudly at load rather than silently, which is the same failure mode
as today and no worse. A module renamed inside a package cannot be expressed in the table: the
lookup key is only the first component of the pickled path.

---

## R9.4a - Nothing is drawn across a read gap

**What.** `panels/within_calibration.py:903` and `:955` draw the cumulative-time and
cumulative-damage curves over the full `t_h` array, while `:245` and `:259` correctly break the
primary trace with `render.observed_slices`. Break both cumulative curves the same way.

**Why.** One site of the class is fixed and two are not, which is the "fix every site of a class"
defect `AGENTS.md` §4 records. The sharpest evidence is internal:
`analyzers/within_calibration_compute.py:199-200`, inside `_observed_dt_h`, states "The figure
already refuses to draw across a gap; this makes the numbers agree with it." That sentence is false
for these two curves, so the code carries a docstring asserting a property it does not have.

**Scope of the claim, stated precisely.** `_observed_dt_h` zeroes gap intervals, so the cumulative
arrays are flat across a gap by construction: the drawn segment is horizontal and both endpoints are
true artifact values. This is ink over unobserved time, not an invented rate. It is worth fixing and
it is not the same severity as the primary trace would be.

**Acceptance.** A test builds panel data with a known gap and asserts the drawn path breaks at the
gap for both curves, against the same `observed_slices` the primary trace uses. MUTATION: restore
the full-array draw at one site; the test goes red naming that curve.

**Budget.** Collect +1 to +2. Files 3. Measured: collect +4, files 2, at the 2x line.

**Known risk.** This changes a shipped figure, which is why `quebraplan.md` §1 puts this phase
before regenerating thesis figures.

**FOUND during implementation.** The primary trace, called the already-correct site above, is only
correct when its caller supplies the carve's gap list. That trap is its own requirement, R9.4c.

## R9.4b - The across-calibration median line must not bridge empty bins

**What.** `panels/across_calibration.py:207-209` draws the 14-day median line across every empty
bin, because `panels/_across_calibration_compute.py:33-37` omits empty bins from `centers` with
`if len(vals) == 0: continue` rather than emitting `NaN`. Emit `NaN`.

**Why.** Separated from R9.4a deliberately. This is a **compute** change to
`AcrossCalibrationPanelData.binned_interval_stats`, not a render change, so it has a different
acceptance and a different blast radius. It also must NOT be justified by the read-gap rule:
`docs/FIGURE_STANDARD.md:128` defines `gap` as a within-calibration inter-read interval, and
`AGENTS.md` §5 forbids substituting the two tiers' vocabularies. The rule here is the plainer one,
that a line joining two points asserts something about the interval between them.

**Acceptance.** On the artifact, not the axes: an empty bin yields `NaN` in
`binned_interval_stats`. `tests/test_across_calibration_builder.py:57-59` and `:76-79` currently
pass under a NaN-emitting variant, so a new assertion is required rather than an existing one
tightened.

**Scope addition, recorded at the checkpoint.** The same function dropped the latest event
whenever the span is an exact multiple of `bin_days`: `np.digitize` puts a value on the last
edge one past the final bin. The final bin is now closed on the right, as in `np.histogram`,
with its own test. On `mtbf_q1` (2894 intervals between 2895 events, 903.007 d) the drop does not
fire, one bin is empty (centre day 77), and the 64 populated bins are identical to the shipped
artifact.

**Fixed at review on Sera's decision.** A line needs two consecutive finite points, so a
populated bin with no populated neighbour, interior between two empty bins or first or last beside
one, drew nothing on the median or p90 line (its IQR fill is a hairline). The renderer marks
exactly those bins on both lines. A record of two or more bins with no empty bin draws as before;
a one-bin record gains the markers, since its single point drew nothing before either. `mtbf_q1`
has no lone bin.

**Budget.** Collect +1 to +2. Files 3. Measured: collect +3, files 3.

## R9.4c - The gap list cannot be omitted on the way to the panel builder

**Added at the review of 9.4-9.6, on Sera's decision.** Found while implementing R9.4a, whose
"already-correct" primary trace is only correct when its caller supplies the gaps.

**What.** `build_within_calibration_panel_data`, `t2star.make_panel_data` and
`fidelity.make_panel_data` took `gap_spans_s` as an OPTIONAL parameter defaulting to empty, and
`recipes.py:196` and `:380` read it with `diagnostics.get("gap_spans_s")`. Below them, the four
private helpers that hand the list to `_observed_dt_h` took `gap_spans_h=None`, and
`SignalBand.gap_spans_h` defaulted to empty. Make it required at every one of these, have the
builder refuse `None`, and index the key in `recipes.py`.

**Why.** Omit the list and the panel has no gaps at all, silently, so:

- every curve on the scan clock, the primary trace included, draws straight across the hole; and
- cumulative time, cumulative damage and occupancy count the unobserved hours as observed, because
  `_observed_dt_h` is handed the same empty list. That is a wrong NUMBER, not just ink.

The window and read tables cannot stand in for the list: the read table has no gap column, and
the window table marks a gap only where an in-spec read borders it, so a gap during an
out-of-spec stretch leaves no trace there. Measured on the shipped two-gap fixture
(`_with_read_gaps`) with the list omitted: the primary trace is one line from 0 to 20 h while the
carve logged `gaps=2`; cumulative time out of spec at 4 µs reads 20.0 h against 14.0 h; occupancy
at 3 µs reads 0.2915 against 0.4165.
Production was not affected: both adapters threaded the value through from the carve, and
`windows.run` always populates it. The defect is the API, in a codebase whose rule is that errors
are raised, not swallowed.

**Acceptance.** A test fails if any of these regains a default, if the builder accepts `None`, or if
either adapter softens `None` before the builder sees it. MUTATION: re-add the default at each
site in turn, remove the `None` check, and make each adapter forward `gap_spans_s or []`; each
goes red. A second test pins the NUMBER on the two-gap fixture against a hand sum: cumulative
time out of spec at 4 µs ends at the observed 14 h, and every state timeline marks both holes
unobserved. MUTATION: hand the reliability band or the distinguish band an empty list, or stop
`_observed_dt_h` zeroing gap intervals; each goes red. Reverting `recipes.py` to `.get` is an equivalent mutant
once the builder refuses `None`, and is recorded as such.

**Budget.** Recorded as measured, not forecast, because the requirement was added after the work:
collect +2. Files 9 (`within_calibration_compute.py`, `t2star.py`, `fidelity.py`, `recipes.py`,
`signal_band.py`, `docs/PANEL_CONTRACT.md`, and the three test files whose gap-free fixtures now
pass the list explicitly).

**Not done.** The same pattern in the independence-check path, `analyzers/check_ledger.py:176`
and `:231` and `analyzers/checks/_multiprocess.py:160` and `:202`. `jobs/bench/arms.py:314` calls
`segments_from_windows` without a gap list on a gap-free synthetic grid, so changing it is a
check-machinery decision. Filed in `quebraplan.md` section 7B.

---

## R9.5 - `FIGURE_STANDARD` describes the real backlog, and a panel states what it dropped

**What.** Three parts. The vocabulary sweep is NOT one of them; see Not done.

(a) `docs/FIGURE_STANDARD.md:6-9` is stale in both directions. It names three vocabulary violations
and there are more; and it says "no panel yet reports how much data it dropped" when **seven now
do**, several citing the document by name (`within_calibration.py:394-396`, `:458-459`;
`check_ledger.py:14`, `:193`; `check_outcome_plot.py:92-104`; `independence_survey_plot.py:276`;
`instrument_validation_plot.py:120-129`). Correct the document to the measured state and point its
vocabulary backlog at the spec that owns it.

(b) `panels/within_calibration.py:821-878` draws a survival estimator that discards censored
windows and prints no count. `analyzers/reliability_band.py:214-219` has the counts, and `:240-245`
`print`s the total to **stdout**, which `docs/FIGURE_STANDARD.md:89-91` explicitly rules out
("Not in the caption, not in the log"). Put the count in the panel and delete the `print`, which is
a side effect inside a step.

(c) `plots/km_survival_plot.py:22-24` and `:173-178` both state that the legend carries `n` as the
last surviving trace of how much data backs each curve; `:180-182` emits only the label and the
median. Either emit `n` or strike the claim.

**Why.** A binding document that misdescribes its own backlog stops being checkable; the panel that
has the censored counts in hand is the one panel that does not print them; and a declared mitigation
that does not exist is worse than an undeclared deviation.

**Acceptance.** `FIGURE_STANDARD.md` names the seven reporting panels and no longer claims none
report. The censored count is visible in the rendered panel. `git grep -n "print(" src/quebra/analyzers/reliability_band.py`
returns nothing. The km legend either carries `n` or no longer claims to.

**Budget.** Collect +1 to +2. Files 4.

---

## R9.6 - A rendered caption may not carry a hardcoded measured number

**What.** `plots/instrument_validation_plot.py:226` prints `36 complete gaps, time censored at
2000 h, 1 censored gap dropped` as a string literal, while `InstrumentValidationData.dropped`
(`analyzers/instrument_validation.py:101`) ships as `{}` and is never read. Compute the caption from
the artifact, and populate `dropped`.

**Why.** `AGENTS.md` §4: a measured number in prose must come from the artifact it cites, in the
state it ships. A figure caption is the one place a stale number is guaranteed to outlive its
source, because nothing recomputes it and a reader cannot tell.

**Acceptance.** `git grep -n "36 complete gaps"` returns nothing. MUTATION: change the generator so
the gap count differs; the caption follows and a test asserts it does.

**Budget.** Collect +1. Files 3.

---

## R9.7 - Close the style ratchet's blind spot, and pin the number last

**What.** Two parts, deliberately at opposite ends of the phase.

(a) **First checkpoint:** widen the pattern in `tests/test_style_baseline.py` to catch the dict form
`**{"fontsize": ...}`, and add that form to the `:7-9` list of what the scope does not catch.

**AMENDED during implementation, on a measurement.** This requirement said to match the bare key
and accept a floor of 19. Measured instead: matching only a NUMERIC LITERAL in either spelling
scores **15**. The bare-key pattern counted two `fontsize=theme.X` references as hardcoding, which
is a false positive against the very behaviour the ratchet exists to encourage, and the two dict
sites it would newly have caught are theme references too. So the Not done alternative below is
what shipped: it is both the stronger gate and the lower floor, and it removes the two units of
slack this requirement's Known risk warned about rather than accepting them. Mutation-verified in
three directions: a planted `fontsize=9` is caught, a planted `**{"fontsize": 9}` is caught, and a
planted `fontsize=theme.CAPTION["fontsize"]` is not counted.

(b) **Last checkpoint:** re-pin `BASELINE` to the measured total, once every other requirement that
moves a literal has landed.

**Why.** The docstring currently misleads a reader into taking 17 for the true count. Measured, the
widened pattern gives **19**: two further sites, `plots/check_outcome_plot.py:138` and
`plots/independence_survey_plot.py:143`.

**The ordering was backwards in an earlier draft and the correction matters.** `:66` asserts
equality, so every later change to a literal forces another re-pin. What must go first is the
**pattern**; what must go last is the **number**. Pinning the number first guarantees re-pinning it
again at the end.

**The Known risk this carried is now spent.** It warned that counting the two
`theme.ON_FILL_TEXT["fontsize"] - 1` sites would raise the floor to 19 and permanently grant two
units of slack, because the test sums one total with no per-site attribution. Taking the numeric
literal alternative removes that: the floor is 15 and every counted site is a real literal. The
final count still goes in the checkpoint banner.

**Acceptance.** The widened pattern is run before and after and both counts are reported. At the
final checkpoint `BASELINE` equals the measured total and both assertions pass.

**Budget.** Collect +0. Files 2.

---

## R9.8 - Remove render code that nothing reaches

**What.** Three deletions.

(a) `plots/interpolation_stage_plot.py` in full. `InterpolationStagePlot` (`build_matplotlib` at
`:25`, 165 lines) is imported by **no job and no test**, and `plots/__init__.py` is empty. This is
the largest dead module in the render layer and an earlier draft of this spec proposed relabelling
two of its axis labels, which would have spent budget on a file that should not exist.

(b) `plots/fidelity_helpers.py:8-41`, `make_fidelity_figure`, 34 of that module's 47 lines, zero
consumers. The surviving 4-line `apply_common_style` is imported by `within_calibration.py:40`,
`across_calibration.py:43` and `comparison.py:18`, and moves into `theme.py`, after which the module
goes. The gain is one fewer module: the import direction is `panels` to `plots.fidelity_helpers`
and becomes `panels` to `plots.theme`, which all three already have, so no layer crossing is removed.

(c) Unregister the `interactive` target. No renderer implements a plotly backend and no job requests
it. Note that a job asking for it already fails with "`{cls}` has no plotly backend", so only the
test is missing, not the behaviour.

**Why.** `AGENTS.md` §3: a file earns its existence from a consumer or a boundary. A registered
target that cannot render is machinery a reader will reasonably assume works.

**Acceptance.** `git grep -n InterpolationStagePlot` and `git grep -n make_fidelity_figure` return
nothing. `sorted(RENDER_TARGETS)` matches what `targets.py` documents. A test asserts a job
declaring an unregistered target fails with a message naming the reason.

**Budget.** Collect +0 to +1. Files 6.

**Follow-through this requirement should claim.** `plotly>=6.0` (`pyproject.toml:34`) is a hard
runtime dependency existing only for the dead target and its `-> go.Figure` annotations.
Unregistering `interactive` is what makes dropping it an accurate-dependency change rather than a
capability removal. Whether to drop it is a decision owed to Sera.

---

## Checkpoints

R9.7(a) is first because the ratchet pattern must be correct before other requirements move
literals. R9.7(b) is last because the pinned number must be pinned once, at the end.

```
CHECKPOINT 9.0 - this spec lands. No behaviour change.
CHECKPOINT 9.1 - the style ratchet catches the dict form; its docstring is honest. (R9.7a)
CHECKPOINT 9.2 - PDF and PNG output is byte-reproducible; the poster is included. (R9.1)
CHECKPOINT 9.3 - the poster font no longer fails silently; the dead theme block is gone. (R9.2)
CHECKPOINT 9.4 - figure artifacts survive a rename; the live-job pickles are re-measured. (R9.3)
CHECKPOINT 9.5 - no cumulative curve is drawn across a read gap. One shipped figure changes. (R9.4a)
                 The gap list cannot be omitted on the way to the builder. (R9.4c, added at review)
CHECKPOINT 9.6 - the across-calibration median line breaks at empty bins. (R9.4b)
CHECKPOINT 9.7 - FIGURE_STANDARD is accurate; the censored count is in the panel; the km legend
                 claim is true. (R9.5)
CHECKPOINT 9.8 - the hardcoded caption is computed from its artifact. (R9.6)
CHECKPOINT 9.9 - dead render code is gone. (R9.8)
CHECKPOINT 9.10 - BASELINE re-pinned to the final measured count. (R9.7b)
```

---

## Not done

- **7.1's axes-level plot contract.** Declined per R9.0.1. Re-open trigger: a composer that needs to
  place two existing figures on one canvas, which is the only consumer the conversion would have.
- **7.2's `subplot_mosaic` composer.** Declined per R9.0.2.
- **Moving `RCPARAMS` into `.mplstyle` files.** Declined. A form change with no consumer, and
  `theme.py:378-388` documents an ordering constraint any rewrite must preserve and could silently
  break. `scaled_text` and the palette cannot live in an `.mplstyle` at all.
- **The vocabulary sweep, and a policy test over the locked vocabulary.** DEFERRED TO ITS OWN SPEC,
  not declined, because it is larger than it looks and does not fit this phase's budget. Three
  reasons, all of which that spec must address. The enumeration in an earlier draft summed to 18,
  not the 21 it claimed. The tick-label class cannot be fixed in the render layer at all: threshold
  labels come from `analyzers/t2star.py:117-128` and `recipes.py:47-48`, which are **step kwargs**
  reaching the provenance label and **dict keys** in every band artifact
  (`analyzers/reliability_band.py:203-228`), so editing them invalidates every cached t2star and KM
  identity. And a grep-for-forbidden-terms test would force wrong rewordings, because
  `docs/FIGURE_STANDARD.md:43-60` reads `window | not interval` as "do not say interval when you
  mean a window", not as a global ban: in the across-calibration tier an inter-event interval
  genuinely is an interval and is not a window, and `AGENTS.md` §5 forbids substituting the tiers.
  `plots/mtbc_hist_plot.py:116`'s "Number of intervals" has no faithful replacement. Note also that
  a test asserting both "returns zero" and "ratchets rather than asserting a fixed list" is only
  satisfiable at `BASELINE = 0`.
- ~~Requiring a numeric literal in the style ratchet~~ ADOPTED at CHECKPOINT 9.1 instead of the
  bare-key widening, on the measurement in R9.7. It is the stronger gate and it is one more regex
  branch, and it lands the floor at 15 rather than 19.
- **`pdf.fonttype: 42` and `ps.fonttype: 42`.** Deferred, not declined, and separated from R9.1
  because R9.0.4 shows they buy text extractability rather than determinism. What would settle it:
  one Overleaf compile of a thesis figure at the default and at 42. `svg.hashsalt` is declined
  outright: no SVG target exists.
- **Pixel-baseline figure tests.** Declined, per 7.6 and R9.0.6. Note that R9.1 makes the poster
  hashable, so the objection that a poster hash gate would break on every matplotlib release is
  answered by pinning rather than by avoiding hashes.
- **A tabular form beside each figure (`to_frames()` or CSV).** Deferred, not declined. R9.3 makes
  the existing artifacts readable, which is the measured defect; a tabular form is a separate
  question about what a reader receives, and it should be scoped against the group's existing
  Zenodo deposit layout (a sibling of this checkout, holding per-figure directories of source data
  plus a plotting notebook, whose pickles are plain `pandas` objects and all load) rather than
  invented here. That layout's path must be named in any spec that cites it, so the comparison can
  be re-measured.
- **Testing the 17 renderers no test constructs.** Out of scope and larger than this phase.
  Recorded because R9.0.6 measured it: the real figure-test gap is not the absence of pixel
  baselines, it is that most figures have no test of any kind.

---

## Open questions this phase does not settle

1. What the plotted vectors ARE, as a contract. R9.3 makes the artifact readable, which is what the
   renderer is handed, not what matplotlib drew. For panels those coincide closely; for any plot
   that computes something render-local they do not. Filing this into `quebraplan.md` §7 is owed.
2. Whether `matplotlib` needs an upper pin. Three findings touch it: the PNG `Software` chunk, the
   PDF `Creator`/`Producer` fields, and `subplot_mosaic` behaviour measured at one version.
3. Whether `plotly` can be dropped from the runtime dependencies once `interactive` is unregistered.
