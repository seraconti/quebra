# Panel Contract

Defines the interface between metric modules and `WithinCalibrationPanel`.

---

## Architecture

`WithinCalibrationPanel` is a **generic** rendering component. It knows nothing
about fidelity, Allan deviation, or any other specific metric. All
domain-specific knowledge lives in adapter functions beside their analyzer
(`analyzers/t2star.py::make_panel_data`, `analyzers/fidelity.py::make_panel_data`)
that convert a typed result into `WithinCalibrationPanelData`.

Those adapters are called by a job STEP, never at draw time. `FidelityPlot` and
`T2StarPlot` must not build panel data inside `build_matplotlib`, where no DAG node
could supply the window tables; both are gone and jobs use `WithinCalibrationPanel`
directly off a panel-data step.

The panel does **not** import from `analyzers.*` for domain logic. It does import
`analyzers.windows` constants for the birth/death vocabulary, which is shared, not
domain-specific.

---

## Three bands

`WithinCalibrationPanelData` is composed of one contract per band, each produced by its own
step, each complete on its own, each replaceable without touching the others:

| band | step | question it answers |
|---|---|---|
| `signal` | `analyzers/signal_band.py` | what was measured, before any threshold |
| `distinguish` | `analyzers/distinguish_band.py` | can a reader tell in from out at all |
| `reliability` | `analyzers/reliability_band.py` | what follows from the 2-state carve |

The outer class owns only what all three share - the ladder, `primary_label`, `traces`,
the render flags - plus `meta`.

**Two guard requirements, both silent if forgotten.** `StaleArtifactGuard` derives its
key set from `dataclasses.fields(cls)`, so the outer class alone would validate nothing
but the four band names: **every band inherits the guard**, and no band uses
`slots=True` (the guard's non-dict branch fires on the slots tuple even for a valid
load). And because the per-threshold maps moved off the class that owns `thresholds`,
each band exposes `check_thresholds(labels)` which the outer `__post_init__` calls -
without it the construction-time completeness contract silently drops to nothing.

**Both compliance timelines are drawn, never merged.** Band 2 renders the 4-state view
(uncertainty included), band 3 the 2-state view the carve actually used, on a shared
x-axis and vertically aligned. If they look the same, resolvability costs nothing. If
they differ, that difference is the finding.

**`estimator` is a field.** The reliability band declares which estimator produced its
survival curve and the axis label is derived from it, so the label cannot go stale.
`cumulative_hazard`, `band_lower` and `band_upper` are reserved and are `None` here. The
estimator exists - `analyzers/kaplan_meier.py` - but this panel does not consume it; wiring
it in changes only that band.

**Occupancy is fraction of OBSERVED time**, gap intervals excluded. There is one
definition: `reliability.occupancy` and the renderer's >=5% timeline cull read the same
number, and `_cumulative_time_out_of_spec` uses the same denominator. Counting gap time
credited unobserved hours to whichever state held at the left edge.

**Not drawn**: the detail/zoom view and the 30-minute median/IQR/p90 view, with
`binned_stats_per_trace` and `adaptive_ylim`. Extra `traces` are now overlaid on the
signal axis rather than in a subplot of their own.

**Shape statistics are reported only where they have support.** Complete windows only
(`up_crossing` birth, uncensored), then `shape_min_reads`; a threshold with fewer surviving
windows than `SHAPE_SUPPORT_FLOOR` is named as unsupported rather than drawn. On the
shipped T2* ladder that is one threshold. Every median ships with its defined-count, and
`rho2` is never rendered without Chatterjee's xi and the falling-limb rho beside it -
a symmetric excursion drives Spearman to zero by construction.

---

## WithinCalibrationPanelData - what it holds

Three band artifacts plus render flags. The per-read arrays a reader might expect here
(`t_h`, `values`, `sigma`, `gap_spans_h`) live on `SignalBand`; polarity is per threshold, the
third element of each `thresholds` tuple.

```python
signal: SignalBand              # analyzers/signal_band.py
distinguish: DistinguishBand    # analyzers/distinguish_band.py
reliability: ReliabilityBand    # analyzers/reliability_band.py
meta: dict[str, object]
thresholds: list[tuple[str, float, bool]]   # (label, value, big_values_good)
primary_label: str
traces: list[tuple[str, np.ndarray]] | None
use_log_scale: bool
color: object
include_cumulative_time: bool
include_cumulative_damage: bool
include_ttf: bool
```

Each band is replaceable without touching the others, which is why they are three modules and
not one. `analyzers/within_calibration_data.py` is the typed contract; read it rather than a
field list here, which is what went stale.

---

## Panel-internal computations

The following are computed by `build_within_calibration_panel_data` from
`(t_h, values, thresholds, damage_fn)`, polarity per threshold, plus the window and read
tables. They are NOT separate metric modules. No analyzer module should reimplement
them.

**Window carving is not one of them.** It lives in `analyzers/windows.py` and reaches
the builder as two required DataFrames, so the gap policy, censoring and per-read state
are the same facts in the artifact, the figure and any downstream analysis. The builder
raises rather than carving a second time.

| Computation | Method | Integration rule | Notes |
|---|---|---|---|
| Threshold compliance timeline | `distinguish_band._timeline_segments` | - | Gantt bars over per-read `state`; 2 or 4 states depending on `use_uncertainty` |
| Window survival | `_window_survival` | - | Empirical survival, **censored windows dropped** |
| CV, initial value, range | `_draw_summary` | - | Summary text |
| Per-threshold window stats | `_analyze_threshold_windows` + `_carve_counts` | - | Above/below counts, mean, p90, plus `n_windows`, `n_censored`, `n_endurance_bags`, `n_gaps` |
| **Cumulative time out of spec** | `_cumulative_time_out_of_spec` | Left-Riemann | Step-function indicator; result in hours |
| **Cumulative damage** | `_cumulative_damage` | Trapezoidal | Continuous damage_rate; result in primary_unit · h |
| **TTF (first crossing time)** | `_ttf` | - | Scalar per threshold; shown in summary text |

### Why left-Riemann for time out of spec, trapezoidal for damage

The out-of-spec indicator is a step function: `{0, 1}`. Trapezoidal
integration would interpolate between 0 and 1 at transitions, which is
physically wrong - a moment is either in-spec or out-of-spec. Left-Riemann
correctly assigns the state at the left edge of each interval.

The damage rate (excess after applying `damage_fn`) is continuous assuming
`damage_fn` is continuous (the linear default is). Trapezoidal integration
is second-order accurate and appropriate for continuous integrands.

### Polarity convention (big_values_good per threshold)

Each threshold tuple carries a `big_values_good: bool` as its third element.
Polarity is per-threshold - different thresholds in the same panel can have
different polarities.

For `big_values_good=False` (e.g. infidelity - lower is better):
  - `out_of_spec[i] = values[i] > threshold_value`
  - `excess[i] = max(values[i] - threshold_value, 0)`
  - TTF: first `i` where `values[i] > threshold_value`

For `big_values_good=True` (e.g. T2* - higher is better):
  - `out_of_spec[i] = values[i] < threshold_value`
  - `excess[i] = max(threshold_value - values[i], 0)`
  - TTF: first `i` where `values[i] < threshold_value`

### Default damage_fn

`None` → identity on excess (`lambda x: x`). Produces linear damage in the
excess above/below the threshold. The y-axis label is
`"Cumulative damage ({primary_label} · h)"`. Callers using nonlinear
`damage_fn` should note that the label is not automatically updated - they
may supply a descriptive `primary_label` that includes units if needed.

### Empty thresholds

If `thresholds == []`, the new subplots render with "No thresholds defined"
and axis turned off. The panel never crashes on empty thresholds.

---

## Window statistics: two definitions on one panel

The summary block's above/below stats (`_analyze_threshold_windows`) and the survival
curve (the window table) answer different questions and do not have to agree:

- **above/below** is direction-agnostic run-length bookkeeping on the raw series. It has
  no gap policy and no censoring; a run spanning a read gap is one run.
- **the window table** applies the gap policy, labels births and deaths, and marks
  censored windows. The survival curve uses only uncensored windows.

Both are correct for what they measure. Do not read the summary `count` as the number of
windows in the window table.

`ReliabilityBand.survival_curve_min` is an empirical survival function over **uncensored**
window lengths only: `S(x) = #{w >= x} / n`. `analyzers/kaplan_meier.py` uses the censored
windows rather than discarding them; this field does not, and the two are not interchangeable.

---

## Color scheme (D4)

Threshold colors come from `plots/theme.py::threshold_color`, which samples
`THRESHOLD_CMAP` across the ladder. They are used in four places:
- On the primary axis: one dashed horizontal line per threshold
- On the survival subplot: one curve per threshold
- On cumulative time and cumulative damage subplots: one curve per threshold

Threshold index `i` maps to `threshold_color(i, len(pd_.thresholds))` in all
four locations, so colors correspond visually.

Sampling across `len(thresholds)` replaced a fixed 8-entry list indexed
`i % 8`, which gave two thresholds the same color on any ladder longer than
8 - the shipped T2* ladder has 10.

Timeline bars are colored by per-read spec state through `plots/theme.py::state_color`,
not by the threshold color: `in_spec`, `out_of_spec`, and - when the carve ran with
`use_uncertainty=True` - `in_spec_uncertain` / `out_of_spec_uncertain`, which are the
crisp colors washed toward white. A read is uncertain when
`abs(value - threshold) < k * sigma`. Uncertainty is an annotation only: it never moves
a window boundary.

---

## Layout

Built in `panels/within_calibration.py`'s `build_matplotlib` from a row list, so the row set
depends on the include flags and the height is dynamic. Fixed width 16 inches.

Rows, in order: `signal`, `distinguish_timeline`, `distinguish_detail`,
`reliability_timeline`, `survival`, `cumulative` (flagged), `summary`. Band 1's row subdivides
into primary plus value / sigma / relative-error histograms; band 2's detail row is a shape
heatmap with xi and rho sub-panels.

No diagram is reproduced here. The row list and its height ratios are one literal in
`build_matplotlib`; a copy in this file goes stale silently, and the previous copy did.

---

## Adding a new metric

To add a new metric that uses `WithinCalibrationPanel`:

1. Write a compute function or analyzer that returns a typed result.
2. Write an adapter function `make_<metric>_panel_data(result) -> WithinCalibrationPanelData`.
   Set `t_h`, `values`, `gap_spans_h` and `sigma` on the `SignalBand`, and `thresholds` (whose
   third element carries polarity) on the panel data. `gap_spans_h` has no default: pass the
   carve's gaps, or `[]` for a record with none.
3. Write a plot class that calls the adapter and delegates to `WithinCalibrationPanel`.

**No metric needs to modify the panel.** The panel's internal computations
(including cumulative time, cumulative damage, and MTTR) cover the standard
within-calibration degradation analysis surface. Only add panel-internal
computation if a new view is fundamentally about `(t, values,
thresholds)` and cannot be expressed as a step result or adapter field.

### What is panel-internal forever (R1 decision)

Cumulative time out of spec, cumulative damage, and MTTR are **not**
candidates for separate analyzer modules. They are panel-internal views
derived from `(t, values, thresholds, damage_fn)`. Any
metric that supplies these inputs automatically gets all three views.

---

## Fidelity adapter

`analyzers/fidelity.py::make_panel_data(result: FidelityResult, ...) -> WithinCalibrationPanelData`

Decisions made by the adapter:
- `values` = infidelity (clipped to ≥ 1e-16 for log scale)
- `primary_label` = "Infidelity"
- `thresholds` = "nines" thresholds in infidelity units (e.g. 0.01 = 99%
  fidelity, 0.001 = 99.9%) - only thresholds that the data actually
  crosses (`inf_min < thr < inf_max`) are included
- `use_log_scale` = True
- `big_values_good` = False per threshold (above infidelity threshold = out of spec)
- For longrun profile: `traces` includes f−fₘₑₐₙ and f−f₀ infidelity traces

---

## File size note

`panels/within_calibration.py` is the largest module in the tree by a wide margin, and it is
past what one reader holds in a single pass. This repo sets no line quota - a file is judged
against the largest single-consumer module in its own package - so the case against this one is
the reading cost, not a number. A line count written here goes stale silently, so it is not
repeated; `wc -l` is the source.

The split is partly done: `analyzers/within_calibration_compute.py`
builds the artifact, `analyzers/within_calibration_data.py` is the typed contract, and
`panels/_within_calibration_render.py` holds the functions-of-axes half. What remains in the
main file is the
drawing sequence, which is still the largest module in the tree. A further split does not
change any public API and can be done at any time.


---

## How a check verdict reaches a panel

`docs/iid_checks/iid_checks_basics.md` said the durable half of that page folds here once the
licence wiring was decided. SPEC 0008 decided it, and the decision was that there is no licence
and no gate: a verdict ANNOTATES a figure and never suppresses one.

The path, and where each decision is made:

| stage | module | decides |
|---|---|---|
| run | `analyzers/check_ledger.run` | scores each check into a verdict from the p-value, the event count and the bench cell |
| select | `analyzers/check_selection` | the **run-set** (which rows are computed) and the **display-set** (which reach a figure), independently, over the full `(check_id, calibration, variant)` key |
| build | `analyzers/check_outcome.build_check_outcome` | reshapes the selected rows into one `OutcomeGrid` per displayed key |
| draw | `plots/check_outcome_plot` | colour from `theme.verdict_color`, and nothing else |

Four rules a panel author has to keep.

**A verdict never gates a draw.** The estimate is always computed and always drawn. A `fail`
changes a cell's colour and the caption; it does not remove a band, a curve or a panel. Control
flow that branches on what the data said is unpredictable, and a reader is better served by an
estimate they are told not to trust than by a missing one.

**Absent is not the same as not computed.** `check_ledger.VERDICT_NOT_COMPUTED` means the check
ran and declined to answer; `VERDICT_ABSENT` means there is no row for that cell at all. They
carry different tones because a reader who cannot tell them apart cannot tell "we asked and got
nothing" from "we never asked". A builder that collapses them is wrong even when the figure
looks fine.

**The selection is a step kwarg, never a sink-loop filter.** Tuples of strings are what
`core/closure.py` will render into a run identity, so a run-set and a display-set reach
`.prov.json` and the Mermaid label. A filter applied while declaring `job.figure` sinks reaches
neither, and the figure then cannot say what it was asked to show.

**The renderer states, it does not conclude.** Which checks appear is resolved in the
output-builder. The panel says what the verdict is and what the check's null was; it does not
tell the reader what to infer. `panels/_within_calibration_render.shape_annotation` is the older
instance of the same rule.
