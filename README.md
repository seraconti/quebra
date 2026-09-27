# QUEBRA

QUEBRA is a toolkit for QUantum Engineering Reliability Analysis. The name is Portuguese for break, and that is the point of the tool.



## Why

Characterizing a quantum device means measuring it repeatedly and recording how good it was: coherence time, gate fidelity, frequency offset, readout error. Do that long enough, and you have a record spanning weeks or months. The usual summary is the metric average and its variation.

That summary answers a question about typical quality, but not the questions an operator actually has. How long does the device hold specification before dropping out? When it drops out, how long until it is back? Is this device more dependable than that one, or does it only look that way because it was measured on a quieter week? Those are reliability questions, and reliability engineering has answered them for a century, for pumps, bearings, turbines and semiconductors. The statistics exist. They are written for a different kind of input.

They are written for durations: how long something lasts before it stops working. A characterization record is not durations. It is a sequence of quality readings, and it never says the device failed, because nothing failed. That is the gap.

We close it by defining what counts as failure. Fix a specification limit on the metric,
call the device in-spec while it sits on the good side of that limit, and the record resolves
into alternating intervals: stretches in spec and stretches out of it. The in-spec stretches
have lengths, and those lengths are durations.

The important part is that this is a choice, not a measurement. Nothing failed; we defined
a threshold and named the crossing a failure. Every result is conditional on that threshold,
which is why a threshold in QUEBRA is a swept parameter rather than a constant, and why a
result is reported across a ladder of thresholds rather than at one of them.

Once the record is expressed as durations, the canonical reliability toolbox applies to it
unchanged: the estimators, the handling of incomplete observation, the comparison tests, the
recurrence models. QUEBRA is our implementation of that toolbox for this kind of data, plus the
machinery to get from a raw instrument record to durations without corrupting them on the way.

## What it does

A run moves through four stages, each marking what is implemented today.

1. Read the record as measured. Characterization data is irregularly sampled: runs stop, instruments get retuned, nights end. Resampling onto a uniform grid is the obvious convenience, and it is a trap, because interpolation invents crossings the instrument never reported and hides crossings that fell between grid points. At working thresholds this shifts the crossing count by tens of percent and every downstream number with it. QUEBRA carves from observed reads only and never bridges gaps in observation.

Loaders for .csv, .yaml, .h5 and .pickle with schema validation; schemas for the
912-day Ramsey record and for calibration event logs; time carried in explicitly suffixed
units, so clock and unit mismatches fail loudly rather than quietly.

2. Carve the windows. A window is a maximal run of consecutive reads meeting spec. The policy must say what happens at a gap in observation, at a window still open when the record ends, and at a window already running when observation began. Each is a different kind of incomplete information and is recorded as such rather than guessed at.

analyzers/windows.py: gap policy, censoring, and a per-read state table saying why each read
did what it did.

3. Check what the estimators assume. Reliability estimators rely on assumptions about the durations they are handed, chiefly that those durations do not depend on each other in time. A confidence band and a two-sample comparison are consequences of the same assumption, so a failed check revokes both together rather than one of them.

Lewis-Robinson, an Anderson-Darling type trend test, rank autocorrelation, exchangeability,
and a copula serial-independence test through an R bridge. The permutation checks share one
permutation set and require an explicit seed because defaulting to system entropy made
p-values irreproducible while the run identity stayed unchanged. Each check is scored against
a measured calibration bench, so “the check passed” stays distinguishable from “the check had
no power to answer”.

4. Estimate, and report the license with the number. An estimate whose assumptions were never checked differs from one whose were checked and held, and both differ from one whose checks failed. QUEBRA carries the check outcome, event count it was computed at, and power verdict alongside the estimate.

Kaplan-Meier survival of in-spec windows. MTBF and MTBC over calibration event logs. Metric
analyzers for T2*, fidelity, telegraph noise and Allan deviation. Nelson-Aalen, log-rank,
RMST and the mean cumulative function belong to the methodology and are not implemented yet.

## How it works

### Four ideas carry the design.

**Everything is a step**. A step is a pure function of typed inputs returning a typed result. It
does not read disk, and it does not draw. Input and output live at the edges: loaders read,
render targets write, panels draw. A step that touches the world is not a function of its
arguments, and could not be cached.

**Jobs are declarative**. A job is a Python file that names its datasets and wires steps
together. It describes a graph; it does not execute one. Nothing computes until a sink
resolves, which is either a figure or an explicit request to materialize an artifact.

**Identity is content, not filename**. A run is identified by hashes of its inputs and of the
code that transformed them, folded transitively through the graph. Change a threshold, and it is
a different run. Edit an analyzer and everything downstream of it is a different run. A cached
result is reused when, and only when, the content that produced it is identical.

**Provenance is append-only.** Every figure and every materialized artifact is written beside a
provenance record, in JSON for machines and Markdown for people, naming the graph that produced
it and the identities that fed it. Output directories are never overwritten and never deleted.

Together these mean a figure in a paper traces back to the exact data and the exact code that
made it, and that changing one parameter recomputes only what changed.

### Status

QUEBRA is in active development alongside an MSc thesis, and the documentation is under
construction. Interfaces may change. Parts of the methodology described above are not
implemented yet and are marked as such where they appear.

Feedback is welcome, particularly from anyone who runs long characterization campaigns and
disagrees with how this frames the problem. Open an issue, including for questions and for
merely confusing things. See CONTRIBUTING.md.

### Install

QUEBRA installs as a package. Python 3.11 or newer.

```bash
git clone https://github.com/seraconti/quebra.git
cd quebra

# an isolated environment, so QUEBRA's dependencies stay out of your system Python
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate

pip install .                      # or: pip install -e ".[dev]" to work on it
pip install pytest                 # the test runner is not a runtime dependency
pytest tests/
```

`scripts/acceptance.sh` checks the same claim from a harsher angle: it builds a wheel,
installs THAT into a throwaway virtualenv outside the repository, and runs the suite from a
directory that is not the checkout. If a packaging mistake makes the steps above work only
from a git clone, that script is what catches it.

There is no `requirements.txt`. Dependencies are declared in `pyproject.toml` and derived
from the imports that actually appear under `src/quebra/`.

### The optional R check

One check, C3 (`copula::serialIndepTest`), runs out of process through `Rscript`. Everything
else is pure Python, and the package imports and runs without R: C3 reports `not computed`
and the other five checks still answer.

**There is no `pip install quebra[r]`, and there cannot be.** An extra installs Python
distributions, and what C3 needs is an R interpreter plus a CRAN package, neither of which
pip can provide. `rpy2` is pip-installable but does not ship R either: it links against an R
you must already have, which is also why this project uses a subprocess rather than `rpy2`
(see `analyzers/checks/c3_serial_copula.py`). An empty `[r]` extra used to sit in
`pyproject.toml` and was removed, because a name you can type at pip that installs nothing
reads as an install route.

To enable C3, install R and the `copula` package yourself:

```bash
# Debian/Ubuntu: sudo apt install r-base    Fedora: sudo dnf install R    macOS: brew install r
Rscript -e 'install.packages("copula", repos="https://cloud.r-project.org")'
Rscript -e 'packageVersion("copula")'       # confirm it resolves
```

`randtests`, `XICOR` and `energy` are not used by the pipeline, but `make test-r` DOES need
them: one `r`-marked test checks the committed fixture's recorded versions against the local
R, and it asks for all four. Install them too if you intend to run that gate:

```bash
Rscript -e 'install.packages(c("randtests","XICOR","energy"), repos="https://cloud.r-project.org")'
```

They are otherwise needed only to regenerate the fixtures with
`jobs/rscripts/reference_values.R`, which the test suite never runs.

Tests that need R carry the `r` marker and skip without it. `make test-r` runs them;
`make test` excludes them, so the default gate stays green on a machine with no R. A
container image would remove this step entirely and is recorded as a later option in
`spec/quebraplan.md`, not as something this project ships today.

Running a job

```bash
quebra run jobs/active/t2star_q1_070423.py   # run one job
quebra run --all                             # every job except those opting out
quebra inspect jobs/active/km_poster_6d2s.py # print the graph without running it
```

`quebra` is installed by pip as a console script. It anchors on the working directory: run
directories are written to `./output`, and `--all` discovers every job under `./jobs` and
runs all of them except those declaring `JOB_SWEEP = False`, which is how the composites opt
out. Pass `--output-root` to send them elsewhere.

inspect is the fastest way to understand a job: it prints the step graph, including the
keyword arguments that affect each result, without computing anything.

Each run writes into a directory named by its content identity, holding the figures, the
materialized artifacts, and the provenance record.

### Documentation

Reference documentation is in docs/. Start with WRITING_A_JOB.md to run something, then
WRITING_A_SCHEMA.md to point the tool at your own data; the rest covers time and clock
semantics, the panel contract, and the figure standard. It is not published as a site yet.

Citing

CITATION.cff carries the machine-readable form, which GitHub renders as a “Cite this
repository” button.

Licence

GNU General Public License v3.0 or later. See LICENSE.
