# Window Semantics (Project Convention)

What "in spec", a side, a window, its length interval and the event table mean, as the
code defines them. Every definition here names the function that holds it; where this page
and the code disagree, the code is right and this page is stale.

## In spec

- **Margin**: `analyzers/windows.margin(values, threshold_value, big_values_good)` is
  `value - threshold_value` when big values are good (T2\*) and `threshold_value - value`
  when small values are good (infidelity). It raises on a non-finite threshold or value.
- **In spec** means `margin >= 0` (`windows.in_spec_mask`). A read exactly at the threshold
  is in spec, in both directions.
- This is the one definition. The carve classifies each read once, in its own units, and
  writes the result to the read table's `in_spec` column. Every other consumer reads that
  column. `tests/test_in_spec_definition.py` holds this two ways. A static scan of
  `src/quebra` and `jobs/` fails on a threshold comparison outside `windows.py` and a short
  exemption list, within the reach its docstring states. A flip test negates the margin and
  fails when one of the quantities it lists does not flip: the window tables, the read
  table's `in_spec`, the in-spec fraction, time to first crossing, cumulative time,
  compliance states and damage.
- **Side**: in spec (`SIDE_IN_SPEC = "in_spec"`) or out of spec (`"out_of_spec"`). The
  strings are the read states, so a side and the state of its reads are spelled one way.
- The uncertain read states (`|margin| < k * sigma`) annotate reads and never move a window
  boundary.

## Windows

A **window** is a maximal run of consecutive same-side reads inside one gap-free stretch. A
gap is a read step longer than `gap_mult` times the median positive spacing
(`windows.spacing`).

`windows.run` returns two tables with one schema:

- **`windows_in_spec`**: the in-spec windows.
- **`windows_out_of_spec`**: the out-of-spec windows.

The **read table** keeps one row per read and threshold; its `window_index` refers to the
in-spec table. A **bag** is the collection of windows.

For a window with reads at indices `s .. e-1`:

| Column | Meaning |
|---|---|
| `t_birth_s` | first read in the window, `t[s]` |
| `t_before_birth_s` | last read before it, `t[s-1]`; NaN when the birth was not observed |
| `t_last_s` | last read in the window, `t[e-1]` |
| `t_death_s` | first read past the death crossing, `t[e]`; `t_last_s` for a censored window |
| `duration_s` | `t_death_s - t_birth_s`, the Kaplan-Meier age |
| `birth_observed` | the window was born by a crossing, not at scan start or gap resume |
| `censored` | the death was not a crossing (the window ended at a gap or at scan end) |
| `extreme_margin` | smallest margin over its reads: closest approach in spec, minus the depth out of spec |

**Crossings are named relative to spec.** An `up_crossing` enters spec and a
`down_crossing` leaves it. An in-spec window is born by an up_crossing and dies by a
down_crossing; an out-of-spec window the other way round. Consumers read `birth_observed`
and `censored`, never the names.

**Windows whose birth was not observed** have no known age: they are excluded from survival
estimates and counted (`n_unobserved_birth_dropped`). The legacy field name for that count
is `n_endurance_bags`.

**Length interval.** The true crossings lie between reads, and the intervals are closed: an
observed birth lies in `[t[s-1], t[s]]` and an observed death in `[t[e-1], t[e]]`. So the
true length lies in the closed interval `[age_lo, age_hi]`:

- **Death:** `age_lo = t[e-1] - t[s]` and `age_hi = t[e] - t[s-1]`.
- **Censored window:** `age_lo = t[e-1] - t[s]` and `age_hi = inf`.

Duplicate timestamps are common (`windows.spacing`). A duplicate at one crossing pins that
crossing exactly, and the length still ranges over the whole interval. A duplicate at both
crossings makes the interval a single point, `age_lo == age_hi`: the length is known
exactly. `EventTable` holds the same closed interval and checks that the
Kaplan-Meier age lies in it.

One reader departs from the closed reading. `survival.turnbull` reads a non-degenerate
interval half-open, `(age_lo, age_hi]`, as R's `Surv(type = "interval2")` does, and a point
interval as an exact length. Its docstring states what that costs, as a validity
assumption.

On a regular grid the Kaplan-Meier age sits at the middle of a death's interval.

**Censoring convention.** A death's age runs to the first read past the crossing; a censored
window's runs to its last read. Censored ages therefore sit about half a spacing earlier than
the deaths' convention. This is deliberate; do not "fix" one without the other.

## The tiling identity

`windows.run` raises unless two things hold, at every threshold:

- **Tiling.** In-spec windows, out-of-spec windows and the gap spans tile the record
  `[t_0, t_N]`. Consecutive pieces share an endpoint by float equality, in read-index
  order.
- **Accounting.** On each side, the summed durations equal the observed time charged to that
  side's reads. Each interval is charged to the read that opens it; a gap interval is
  charged to nobody. Equivalently, windows times mean length equals the side's fraction
  times observed time.

The two sides of the accounting identity share no code. A fraction taken from the windows'
own durations could never fail.

## The event table

`analyzers/event_table.EventTable` is what every survival estimator reads. There is one per
record, threshold and side. It holds:

- per window: `age_s`, a boolean `event`, `age_lo_s` and `age_hi_s`;
- the risk table: distinct ages, deaths, censorings, and the number at risk;
- `n_unobserved_birth_dropped`.

At a tied age deaths come first, so a window censored at age `a` is still at risk at `a`.
`event_table.from_carve` builds one from a carve, validating the label against the carve's
ladder. An empty side is an empty table, not an error. The constructor also takes any other
right-censored durations directly and validates them; no caller outside `event_table.py`
builds one that way today.

## What an estimate rests on

Each result in `analyzers/survival.py` names its standing in a field, from
`analyzers/check_attachment`: `standing`, or `point_standing` and `interval_standing` on
the curve summaries. Nothing else carries the field, the Kaplan-Meier curve and the in-spec
quantities included; the table gives their standing by the same rule.

| Standing | Meaning | Results |
|---|---|---|
| survives a rejection | true of this record whatever the independence checks say | in-spec fraction, mean window length, placement bracket |
| descriptive, `n_censored` shown | exact when nothing is censored; with censoring it leans on one hazard shared by every window at risk | Kaplan-Meier points, quantile points, restricted mean |
| rests on exchangeability | carries the check outcome | every band, interval and SE; the Nelson-Aalen shape; Turnbull |

**The check outcome** attached to a result names the `a1_renewal_durations` assumption, the
checks asked, and those that gave no answer. Empty means NOT ASSESSED.

**Whether a verdict transfers is UNRESOLVED.** A verdict computed on the battery's sample
(complete in-spec durations on the in-spec clock, or births on the calendar clock) is not
automatically a verdict about the Kaplan-Meier sample. `analyzers/assumptions.py` records
that open question.

**Single-interval Turnbull is biased here.** On a read grid, the chance that a length shows
`n` reads is triangular across its interval, so the interval is informative about the length
and Turnbull's likelihood is the wrong one. The bias does not shrink with more windows.
`survival.turnbull` ships with that stated. The grid-placement question is answered by
`survival.placement_bracket`, which needs no such assumption.
