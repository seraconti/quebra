# Figure Standard

Binds every figure added or edited from here on. `plots/theme.py` owns the colour and
typography half of this; the rules below are the half a module cannot enforce.

**The existing panels do not all conform yet.** Two backlogs, stated as they stand:

- **Dropped data.** These report in the figure what they leave out:
  `panels/within_calibration.py` (fit-error values above the view, relative errors
  excluded, and censored windows each survival curve excluded), `panels/check_ledger.py`
  (thresholds with no usable windows), `plots/check_outcome_plot.py` (checks shown of those
  run, and checks asked for that produced no answer), `plots/independence_survey_plot.py`
  (cells not computed, and that C3 is drawn with no bench cell) and
  `plots/instrument_validation_plot.py` (the tier that is not drawn, and the censored gap the
  load-haul-dump comparison drops). Two record a deviation in their own docstrings and carry
  the count only on the artifact: `plots/km_survival_plot.py` (windows whose birth was not
  observed) and `plots/mtbc_hist_plot.py` (non-positive intervals). `plots/tlf_plot.py` drops
  non-finite values, and `plots/allan_plot.py` non-finite and non-positive fractional-deviation
  points, with no note.
- **Vocabulary.** The backlog is larger than a few axis labels, and part of it cannot be fixed
  in the render layer: threshold labels are step kwargs and artifact keys. It is deferred to
  its own spec; see `spec/specpresentation09.md`, Not done.

Those are targets, not a description of every panel. Conform a panel when you next touch it;
do not sweep them.

This document uses spaced hyphens throughout, per its own rule.

---

## Titles

**A title names the object. It does not argue.** A figure title is a label, not a
finding. The reader draws the conclusion; the title tells them what they are looking at.

- Good: `In-spec window survival, T2* >= 3 us`
- Bad: `Windows tell us more than the raw signal does`

The bad one states a claim the axes cannot be checked against, and it stops being true
the moment the data changes.

**A title is one line.** No sentence-length titles, no rhetorical titles, no question
marks.

- Good: `Interval between calibration checks, qubit 2`
- Bad: `Interval between calibration checks for qubit 2, showing the long tail that the mean does not capture`

If the qualification will not fit, it belongs in the caption.

---

## Axis labels

**Use the locked vocabulary. Never synonymise it.** These words have one meaning each
across the whole tool. A synonym in an axis label makes two figures look like they
measure different things when they do not.

| term | not |
|---|---|
| window | interval, period, span |
| read | sample, point, measurement |
| bag | group, batch |
| check | test, probe |
| band | region, envelope, ribbon |
| scan clock | wall time, elapsed time |
| window age | time since start, duration so far |
| birth type | origin, start kind |
| within-calibration | intra-cal, same-cal |
| across-calibration | inter-cal, cross-cal |
| gap | hole, missing stretch |
| record | row, entry |
| dataset | file, run, series |
| run-set | checks enabled, active checks |
| display-set | shown checks, visible checks |
| check outcome | trust annotation, confidence badge, health indicator |

**Units go in the axis label, in parentheses. Never in the tick labels.**

- Good: `Window age (h)` with ticks `0  12  24  36`
- Bad: `Window age` with ticks `0 h  12 h  24 h  36 h`

Repeating the unit on every tick costs space and breaks alignment between panels that
share an axis.

---

## Captions

**A caption states what is drawn and what was excluded. It does not state a conclusion
the figure does not show.**

- Good: `Survival of in-spec windows for qubit 2, 912-day Ramsey record. Windows that
  begin inside a calibration gap are excluded (n = 14).`
- Bad: `Survival of in-spec windows for qubit 2, showing that drift dominates after the
  first week.`

The bad one asserts a cause the figure cannot separate from any other, and it names a
time scale that is not marked on either axis.

---

## Dropped data

**Every panel that drops data says how much it dropped, in the panel.** Not in the
caption, not in the log, not only in provenance - in the panel, where anyone reading
the figure will see it.

- Good: an annotation reading `excluded: 14 of 302 windows (gap-born)`
- Bad: a filtered series with no count anywhere on the axes

A reader who cannot see the exclusion cannot judge it. This is the same reason
`Job.load()` raises instead of defaulting: a silent drop produces a wrong-but-plausible
figure.

---

## Punctuation

**Spaced hyphens throughout. Never em dashes.** Applies to every string that reaches a
figure - titles, axis labels, legend entries, annotations, captions.

The rest of the repo uses em dashes in prose and docstrings. Do not retro-edit it; this
rule governs rendered figures and this document.

---

## Convention note: interval-censored durations

The crossing is observed at the first out-of-spec read, so a complete duration is the
upper end of an interval of width equal to the read spacing.

---

## The vocabulary now has code behind it

Four of the locked terms stopped being conventions and became columns in
`analyzers/windows.py`. Use them to mean exactly what the carve means:

| term | definition in code |
|---|---|
| window | one row of the window table: an in-spec run with a birth, a death and a duration |
| read | one row of the read table: a single observation, with its state and margin |
| gap | an inter-read interval greater than `gap_mult` times the median positive spacing |
| band | reserved; `BAND_STYLE` in `plots/theme.py` carries its fill and edge alphas |
| run-set | which check rows are COMPUTED: a step kwarg over `analyzers/check_selection.ALL_KEYS` |
| display-set | which check rows reach a figure: an independent subset of the run-set, also a step kwarg |
| check outcome | the verdicts a figure shows beside an estimate. `analyzers/check_outcome.CheckOutcome`. Never called a trust annotation: the code's word is `verdict`, and "trust" invites reading a non-rejection as reassurance |

Two consequences for figures:

- A window that did not die of an observed crossing is **censored**, and a figure that
  counts it as a completed lifetime is wrong. Say which windows were excluded.
- Nothing is drawn across a gap. The trace breaks; the timeline does not bridge it.

A read whose error bar overlaps the threshold is **uncertain**, drawn in the washed tone
of its crisp state. Uncertain is a statement about resolution, not about compliance -
the window boundaries are identical with and without it.
