# Tellability Semantics

Tellability asks whether a record, read on its schedule with its error bars, can tell when a
threshold crossing happened. This page defines the quantities that answer it, in the code's
terms (`analyzers/tellability.py`, `analyzers/spectrum.py`, `analyzers/min_p.py`). It defines;
it does not argue. The rationale for each choice belongs to the methodology chapter.

Every number here is computed on the reads as observed. Nothing is resampled onto a grid.

## Rules

- **No verdict is derived from the data.** No row is turned into a pass or a fail and no cut is
  taken from the record. Every number ships, with its status, and the reader reads the panel.
  The names that do ship describe and do not judge: the panel names the member that attains
  its p, and crossings are classed clean, merged or unbracketed.
- **No tellability number changes what any step computes.** Nothing downstream reads one.
- **One record, one gap rule.** Tellability segments the reads with the carve's rule. It
  refuses windows carved under another `gap_mult`, at another confidence level `k`, or without
  one of its thresholds, compared by label, value and side, and windows carved on other reads,
  of another record or of part of this one, compared by read count and median spacing.
- **A row's numbers agree with its statuses when it is built.** A row whose computed field is
  missing or not finite, or whose not-computed field carries a value, cannot be constructed.

## Status words

Each field group of the result carries one status:

- `computed`: every field of the group is set.
- `not computed: <reason>`: every field is `None`. The reason is one of a fixed set:

  | Reason | Groups that raise it |
  |---|---|
  | `sigma propagation not derived` | every sigma group, as infidelity's sigma status |
  | `fewer than <n> reads` | window, spectrum, G10 halves and self-consistency (2), dependence (10), clipped and run length (1) |
  | `no segment with at least <n> reads` | spacing, bounds and rho_k (2); second differences and the sigma bound (3) |
  | `span shorter than <n> median spacings` | spectrum (2, and 20 for the top tenth of the band), G10 and its halves (20) |
  | `no down-crossing at this threshold` | count exponent, crossing counts |
  | `no in-to-out transition within a segment` | run length |
  | `fewer than two rungs with a down-crossing` | count exponent |
  | `zero variance` | values of zero range (spectrum, G10 halves, second differences, self-consistency, sigma bound, dependence), a sigma of zero range (dependence), an in-spec indicator that never changes (rho_k), no positive spacing (spacing), an edge level that is not positive (noise floor) |
  | `no floor for this metric` | clipped |
  | `no read with a finite sigma` | every sigma group, when sigma is declared and no read carries a finite one |

  When the spectrum is not computed, the noise floor and G10 carry its reason; the noise floor
  carries the sigma status first. Zero range is tested on the values themselves, never as a
  standard deviation above 0: equal values can give a deviation of rounding.
- `not declared`: the record carries no sigma, so the rows that read sigma have nothing to read.

A metric's sigma status is `declared` (given, and read), `not declared` (absent), `not
computed: sigma propagation not derived`, which is infidelity's (its sigma is not derived yet),
or `not computed: no read with a finite sigma`, which the step sets once when a declared sigma
is finite on no read. The rows that read sigma (the noise floor, the dependence panel, sigma's
self-consistency and its bound, the crossing counts and the in-spec fraction bounds) carry that
status whenever it is not `declared`, before any other reason they might have. Tellability reads
the metric's sigma whatever the carve's `use_uncertainty` switch says.

Every result also carries, with no status: the number of reads, of segments, of reads dropped
for a non-finite time or value, and of reads kept with no finite sigma (`None` when no sigma is
declared); and the rows of the file the load dropped before any step ran,
`n_rows_dropped_by_schema` and `n_rows_dropped_invalid` (`None` when the load did not record
them; see `docs/WRITING_A_SCHEMA.md`). The median spacing and the gap threshold belong to the
spacing group.

## Reads, segments, repeated times

A **segment** is a maximal gap-free run of reads under the carve's gap rule
(`windows.observed_reads`): a gap is a step longer than `gap_mult` times the median positive
spacing, and non-finite reads are dropped first. `observed_reads` raises when time steps
backwards or no step between reads is positive. No pair, difference or crossing spans a gap,
and a spacing across a gap is never observed time. The rows that use no pair at all (the
spectral window, the spectrum and G10, `sd(q)`, the variance `v` and the dependence panel)
take every read, across gaps.

Reads at one time are legal. The spectrum, G10 and the second differences collapse them to one
read at their mean value, with the sigma of that mean, `sqrt(sum sigma_i^2) / m`, which is
unknown when any read of the group has no sigma; the noise floor (Q1.2) and the sigma bound
(Q3.3) take that collapsed sigma, so both sides of each ratio see the same reads. Segment sizes
in the second differences' reason count reads after this collapse. The spectral window and the
other rows count every read. The Allan deviation and the telegraph dwells raise on a repeated
time inside a segment, which they cannot space.

## A read with no sigma

A read with no sigma is never sure. The carve, when it draws uncertainty, marks it uncertain on
the side its value falls (`in_spec_uncertain` or `out_of_spec_uncertain`); the crossing
placement (Q3.4) never takes it as confident; the in-spec fraction bounds count the interval it
holds as unsettled. One rule, `windows.confident_mask`. It changes no read in the current
records: the T2* and Ramsey loads (`track912Schema`) drop every row with a missing fit error
before any step runs, and the count of dropped rows is on the artifact
(`n_rows_dropped_by_schema`).

## Spectral quantities (`analyzers/spectrum.py`)

Notation: the reads are at times `t_1 < ... < t_N`; `Delta_bar` is the median positive spacing;
`T` is the record's span, gaps included.

- **First alias frequency** `f_s = 1 / Delta_bar`. **Band edge** `f_e = 1 / (2 Delta_bar)`.
- **Frequency grid**: the independent frequencies `k / T`, `k = 1, 2, ...`, up to `f_e`.
- **Spectral window** `W(f) = |sum_n exp(-2 pi i f t_n)|^2 / N^2`. It depends on the read times
  only. It is 1 at every alias of a perfect lattice and about `1 / N` for scattered times, so
  `W(f_s)` near 0 means nothing above the band edge can fold onto it.
- **Spectrum** `S(f) = 2 Delta_bar P(f)`, with `P` the floating-mean Lomb-Scargle power of the
  values (`scipy.signal.lombscargle`), unweighted. For white noise of variance `sigma^2`, `P`
  averages `sigma^2`, so `S` is the one-sided density `2 sigma^2 Delta_bar`, which integrates to
  `sigma^2` over `[0, f_e]`. For a coloured process on irregular reads, `S` is the density seen
  through the spectral window, which leaks power between frequencies; it is not the density
  itself. Values of zero range have no spectrum. Units: the value's unit squared per hertz.
- **Derivative variance** `L(f) = integral from 0 to f of (2 pi f')^2 S(f') df'`, by the
  trapezoid rule on the grid, with `S` held at its first value below `1 / T`.
- **G10** `= L(f_top) / L(f_top / 10)`, with `f_top` the last grid frequency, within `1 / T`
  below `f_e`, and `L(f_top / 10)` interpolated on the grid, so the decade is exact. G10 needs
  ten grid frequencies, a span of at least 20 median spacings. For the continuous integral a
  flat spectrum gives 1000 and `S ∝ f^-2` gives 10.

  G10 is read as a one-way indicator, and no exponent is derived from it. On irregular reads,
  leakage moves G10 up for a spectrum that falls across the band, so a small G10 means a red
  spectrum. A large G10 can be a red spectrum read through leakage; on sparse, gappy schedules
  a red spectrum can read close to white. A near-white spectrum with power near the band edge
  can read below its truth, but stays large. `tests/test_spectrum.py` pins both directions on
  a sparse, gapped schedule.

## Per metric

- **Spacing.** `Delta_eff = sum(delta^2) / sum(delta)` over the within-segment spacings: the
  mean length of the interval a crossing lands in (Sera's L0 specification, section 8.3).
  Beside it, the mean spacing and the CV, with the `1/n` variance, so that `Delta_eff = mean x
  (1 + CV^2)`.
- **Q1.1** `W(f_s)`.
- **Q1.2** `floor / S(f_e)`, with `floor = 2 sigma_bar^2 Delta_bar`, `sigma_bar^2 =
  mean(sigma_i^2)` over the reads with a finite sigma (repeated times collapsed), and `S(f_e)`
  the mean of `S` over the top 10% of the band, `[0.9 f_e, f_e]`. On irregular reads the edge
  level holds power leaked from lower frequencies, so a ratio below 1 can be leakage, not
  signal.
- **Q2.1** `sigma_hat_D2 / sd(q)`. `sigma_hat_D2` is the Gasser-Sroka-Jennen-Steinmetz noise
  scale within segments: for a read with a neighbour on both sides, the pseudo-residual
  `e = a q_(i-1) + b q_(i+1) - q_i`, with `a = (t_(i+1) - t_i) / (t_(i+1) - t_(i-1))` and `b =
  1 - a`; `sigma_hat_D2^2` is the mean of `e^2 / (a^2 + b^2 + 1)` and `sigma_hat_D2` its
  square root. On equal spacing `sigma_hat_D2^2 = (1/6) E[(second difference)^2]`.
- **Q2.3** G10 on the whole record, and G10 of each half of the reads, each half on its own
  span but up to the whole record's band edge. A half whose values have zero range is not
  computed.
- **Q3.0, the dependence panel.** Is the reported sigma independent of the value? Five members,
  each oriented so that larger is stronger dependence: `|Pearson r|`, `|Spearman r|`,
  Chatterjee's xi with the value as the explanatory variable, the distance correlation, and
  the decile spread `(largest - smallest mean sigma over the ten deciles of q) / mean sigma`.
  Each null row shuffles sigma against q; every member is recomputed on the same rows; `min_p`
  gives every member's own permutation p (ties counted against, as in the check battery) and
  the panel's single-step min-P p. The seed and the number of permutations are step settings.
  The panel's p is conservative, never below `d / (B + 1)` with `d` at most the member count.
  The panel uses the reads with a finite sigma, refuses a negative one, needs at least 10 of
  them and a value and a sigma that both vary, and splits ties in the deciles by position. Its null shuffles sigma
  against q over the whole record, drift included: a slow drift shared by both makes the
  panel's p too small.
- **Q3.2** `(v - mean(sigma^2)) / v`, with `v` the variance (`ddof = 1`) of every read with a
  finite value and `mean(sigma^2)` over the reads with a finite sigma.
- **Q3.3** `sigma_bar / sigma_hat_D2`, both on the reads with repeated times collapsed.
- Q3.2 and Q3.3 are conditional on sigma being independent of the value it accompanies (L0,
  sections 5.3 and 5.4). They are always computed; when the dependence panel finds dependence,
  they are computed but not interpretable.
- **Clipped fraction**: the share of reads at the metric's floor. Infidelity is clipped at
  `fidelity.INFIDELITY_FLOOR`; T2* and frequency have no floor.

## Per threshold

In spec means `windows.margin >= 0`, the one definition. Every threshold reads the windows the
carve made for it. The frequency of a Ramsey record has no thresholds, so it ships the
per-metric rows only, read from `raw_frequency_hz` with `frequency_error_hz` as its sigma by
default; the artifact records the value key, and the step's settings record both.

- **Q2.2, the count exponent.** Thinning by `m` keeps every `m`-th read of each segment,
  starting again in every segment, from each of the `m` offsets in turn; `m` runs over 1, 2,
  3, 5, 8, 13, 20. `n(m)` is the mean number of down-crossings between kept reads and `P(m)`
  the mean number of kept pairs, so `P(m) = sum over segments of (L_s - min(m, L_s)) / m`. No
  thinning step runs: over the `m` offsets the kept pairs are exactly the within-segment pairs
  at lag `m`, and the code counts those. `gamma = 1 - slope of log(n / P) on log m`, fitted on
  the rungs with a down-crossing and a kept pair; the others are listed in `rungs_excluded`,
  and their `count_ratio` is `None` when they keep no pair. Without a down-crossing at `m = 1`
  the row is `no down-crossing at this threshold`. `count_ratio = [n(m) / P(m)] / [n(1) / P(1)]` is
  normalised per kept pair, with limits 1 (independent reads) and `m` (a smooth signal); this
  is NOT L0's `R_m = n(m) / n(1)`, whose limits are `1 / m` and 1.
- **Q3.4, crossing placement**, at `k = 1` and 2. A read is confident when its sigma is finite
  and `|margin| >= k sigma` (`windows.confident_mask`); a read whose sigma is unknown is never
  confident. Each raw down-crossing inside a segment is bracketed by the nearest confident read
  at or before it and at or after it, in the same segment: CLEAN when that bracket runs from a
  confident in-spec read to a confident out-of-spec one with exactly one crossing inside;
  MERGED when such a bracket holds more; UNBRACKETED otherwise: a segment end or a gap comes
  first, or the bracket stays on one side or runs from out of spec to in spec. The counts
  ship, and the clean fraction at `k = 1`.
- **Q4.1, the run-length ratio.** Two in-spec quantities, kept apart: the **in-spec fraction**
  is time-weighted, each within-segment interval held by its left read
  (`within_calibration_compute._threshold_in_spec_frac`); `read_fraction_in_spec` (`A`) is per
  read. The mean in-spec run, in reads, is read from the carve's window table:
  `sum(n_reads - censored) / sum(1 - censored)`, the reads with an observed successor in their
  segment over the observed in-to-out transitions. `run_length_ratio = mean run x (1 - A)`,
  which is 1 for independent reads.
- **In-spec fraction bounds.** An interval is unsettled when it holds a crossing or its left
  read is not sure at the carve's `k`: its sigma is unknown, or `|margin| < k sigma`. The lower
  bound counts every unsettled interval out of spec, the upper bound in spec; observed time is
  the sum of the within-segment spacings. The carve's in-spec fraction lies between them. They
  bound the fraction of time held by reads in spec; an excursion between two reads on the
  same side lies outside both.
- **rho_k** of the in-spec indicator, `k = 1` to 10, within segments: centred on `A`, scaled by
  `A (1 - A)`, averaged over the within-segment pairs at lag `k`. The array stops at the last
  lag that has a pair.

## The Allan deviation's tau

`analyzers/allan.py` treats every segment as read at one spacing, `tau0`, the whole record's
median positive spacing, so `tau = m tau0` is nominal when the reads inside a segment are
unevenly spaced. Its figure prints the within-segment spacing CV beside it. The CV is shown and
not acted on; whether to act on it is an open question for issue #4.

## Where these depart from Sera's L0 specification and from the article

| # | L0 | Here |
|---|---|---|
| L1 | gamma from `n ∝ m^-gamma`; `R_m = n(m) / n(1)` | gamma from `n / P`; `count_ratio` per kept pair |
| L2 | keep every `m`-th read | per segment, averaged over the `m` offsets |
| L3 | `S` at `f_e` | the mean over the top 10% of the band |
| L4 | the estimator of `S` unspecified | `2 Delta_bar` times the floating-mean Lomb-Scargle power |
| L5 | second differences at equal spacing | Gasser-Sroka-Jennen-Steinmetz weights |
| L6 | mean run, edge runs unstated | edge runs kept; reads with an observed successor |
| L7 | every raw crossing | down-crossings inside a segment; a gap-cut bracket is unbracketed |
| L8 | confident at `|q - u| > k sigma` | sigma finite and `>= k sigma`, the carve's rule |
| L9 | Q3.1 a precondition | the dependence panel ships whole; nothing decides |
| L10 | `sigma_bar` | `sigma_bar^2 = mean(sigma_i^2)` |
| L11 | ADEV on the fractional series, with intervals | per segment and pooled; fractional only with a declared `nu0_hz`; `adev_err_hz` is adev over the square root of the pooled terms, which treats them as independent and so understates the uncertainty |
| L12 | sigma rows for the metrics it tabulates (frequency, T2*, offset, amplitude) | infidelity's sigma rows: `sigma propagation not derived` |
| L16 | Q2.3: `alpha = 3 - log10 G10`, quoted beside G10 | G10 only, read one way; no exponent, since leakage on irregular reads biases it |
| L14 | the article's PSD, by interpolation and a Fourier transform (Supplementary notes, "PSD from frequency calibration data", p. 18) | Lomb-Scargle, no interpolation |
| L15 | the article's filtering | the chi-square guard (the released records are already filtered); no sigma clip |

L13 is the definition of the check battery's panel rows, which do not exist yet; it is written
with them.
