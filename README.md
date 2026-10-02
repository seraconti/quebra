# QUEBRA

[![ci](https://github.com/seraconti/quebra/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/seraconti/quebra/actions/workflows/ci.yml)
[![licence: GPL-3.0-or-later](https://img.shields.io/badge/licence-GPL--3.0--or--later-blue.svg)](LICENSE)

**QUantum Engineering Reliability Analysis.** QUEBRA reads the quality metrics a quantum device
already reports (coherence times, fidelities, frequencies) as time series, declares a threshold,
and turns each record into alternating in-spec and out-of-spec windows. It then applies
reliability statistics to those windows and to calibration logs, checks on each record what
those statistics assume, and traces every figure to the content hashes of the data and code that
produced it.

The name is Portuguese for *break*, which is the point of the tool.

> **Status: v0.1.0.** Research software, developed alongside an MSc thesis. Interfaces may
> change. Install from source: QUEBRA is not on PyPI. What is and is not implemented is listed
> under [What it does](#what-it-does); open work is in the
> [issues](https://github.com/seraconti/quebra/issues) and
> [milestones](https://github.com/seraconti/quebra/milestones).

## Statement of need

Characterizing a quantum device means measuring it repeatedly and recording how good it was.
Run long enough, that gives a record spanning nights or months, usually summarized by an average
and a spread. That summary answers a question about typical quality. It does not answer the
questions an operator has: how long does the device stay in specification, how long does it take
to come back, and is one device more dependable than another or only measured on a quieter week?

Those are reliability questions, and reliability engineering has standard statistics for them.
The statistics are written for durations, and a characterization record contains none: it is a
sequence of quality readings, and nothing in it ever says the device failed.

QUEBRA closes that gap by constructing the failures. A threshold on the metric is declared, the
device is in spec while the metric sits on the good side of it, and each crossing out of spec is
a failure. The record then resolves into in-spec and out-of-spec windows, and their lengths are
durations. This is a choice, not a measurement, so every result is conditional on its threshold:
a threshold in QUEBRA is a swept parameter, always reported with the result.

## What it does

A run moves through four stages. Each lists what ships today and what does not yet.

| Stage | Implemented today | Not implemented yet |
|---|---|---|
| **1. Read the record as measured** | Loaders for `.csv`, `.yaml`, `.h5`, `.pickle` with schema validation; schemas for the 912-day Ramsey record and calibration logs; time in explicitly suffixed units | Typed calibration fields beyond the timestamp (end, target, kind, outcome, trigger) |
| **2. Carve the windows** | In-spec windows from observed reads only, with an explicit gap policy and censoring; per-read states, including unobserved stretches and read states by margin | Out-of-spec windows (time to recovery); read-schedule diagnostics |
| **3. Check what the estimators assume** | Trend: Lewis-Robinson (C1), Anderson-Darling type (C2), Cramer-von Mises type. Dependence: rank autocorrelation (C5), exchangeability (C6), copula serial independence through R (C3). One shared, seeded permutation set; two clocks; a check ledger; a calibration bench with size and power tables | A bench cell for C3, whose verdicts are therefore reported as uncalibrated |
| **4. Estimate, and report the licence with the number** | Kaplan-Meier with a log-log Greenwood band; share of time in spec and cumulative excess per threshold; shape statistics (Chatterjee's xi, Spearman, distance correlation); MTBF and MTBC as mean and spread of intervals; metric analyzers for T2\*, fidelity, telegraph noise and Allan deviation | Kaplan-Meier in the panel's reliability band (it still ships a crude empirical estimate); median and restricted mean; Nelson-Aalen and hazard; log-rank; between-record variance; product-limit on calibration cycles |

Interpolating an irregular record onto a uniform grid invents crossings the instrument never
reported and hides others, so QUEBRA carves from observed reads only and never bridges a gap in
observation. A failed check revokes the confidence band and the comparison that rest on the same
assumption together, and the check outcome, the event count and the power verdict travel with the
estimate.

## How it works

Four ideas carry the design.

- **Everything is a step.** A step is a pure function of typed inputs returning a typed result.
  It does not read disk and does not draw: loaders read, render targets write, panels draw.
- **Jobs are declarative.** A job is a Python file that names its datasets and wires steps
  together. It describes a graph; nothing computes until a sink (a figure, or an explicit request
  for an artifact) resolves.
- **Identity is content, not filename.** A run is identified by hashes of its inputs and of the
  code that transformed them, folded transitively through the graph. Change a threshold and it is
  a different run; a cached result is reused only when the content that produced it is identical.
- **Provenance is append-only.** Every figure and artifact is written beside a provenance record,
  in JSON and Markdown, naming the graph and the identities that fed it.

A figure in a paper therefore traces back to the exact data and code that made it, and changing
one parameter recomputes only what changed.

## Install

Python 3.11 or newer. From source; there is no PyPI release.

```bash
git clone https://github.com/seraconti/quebra.git
cd quebra
python -m venv .venv
source .venv/bin/activate              # Windows: .venv\Scripts\activate
pip install .                          # to work on it: pip install -e ".[dev]"
pip install pytest hypothesis          # the test suite's conftest needs hypothesis
pytest -m "not slow and not heavy and not r" tests/
```

`scripts/acceptance.sh` checks the same claim more harshly: it builds a wheel, installs it into a
throwaway environment outside the repository, and runs the suite from a directory that is not the
checkout. Dependencies are declared in `pyproject.toml`; there is no `requirements.txt`.

### Optional: R, for check C3

One check, C3 (`copula::serialIndepTest`), runs out of process through `Rscript`. Everything
else is pure Python; without R, C3 reports `not computed` and the other checks still answer.
No pip extra can install R, so install it yourself:

```bash
# Debian/Ubuntu: sudo apt install r-base   Fedora: sudo dnf install R   macOS: brew install r
Rscript -e 'install.packages("copula", repos="https://cloud.r-project.org")'
```

`make test-r` runs the tests marked `r`; they also need `randtests`, `XICOR` and `energy`, because
one test compares the committed fixture's recorded package versions with the local R.

## Usage

### Without any private data

A synthetic Ramsey record ships inside the package. This runs it through the real T2\* analyzer,
and is the same snippet `scripts/acceptance.sh` runs from outside the repository:

```python
import quebra.analyzers.t2star as t2star
from quebra._fixtures import fixture_path
from quebra.core.dataset import Dataset
from quebra.core.job import _load_dataset

norm = _load_dataset(Dataset(path=fixture_path("ramsey_synthetic.csv"), qubit=1,
                             extra={"run_start_unix_s": 1.7e9}))
result = t2star.run(t2star.make_inputs_from_norm(norm))
print(len(result.frame), "T2* points")
```

### Running jobs

```bash
quebra inspect jobs/active/km_poster_6d2s.py   # print the step graph without computing
quebra run jobs/active/t2star_q1_070423.py     # run one job
quebra run --all                               # every job, except those declaring JOB_SWEEP = False
```

Runs are written to `./output/<job>_<identity>_<timestamp>/`, holding the figures, the materialized artifacts
and the provenance record; `--output-root` sends them elsewhere. `make promote` copies one run's
provenance (not its artifacts) into the committed tree, so a published figure stays auditable.

The jobs under `jobs/active/` read embargoed records. Without them, a run stops with
`DataUnavailable`, which names the file, every location tried, and whether the file is embargoed
(listed in `data/real_private/MANIFEST.toml`) or the path is simply wrong.

### Data layout

| Directory | Contents | In git |
|---|---|---|
| `data/real_private/` | embargoed records | only `MANIFEST.toml` (filename, sha256, provenance) |
| `data/real_public/` | records that may be redistributed | yes, if small |
| `data/simulated/` | regenerated payloads | only seeds and manifests |

The data root resolves from `--data-root`, then `QUEBRA_DATA_ROOT`, then `data_root` in a
`quebra.toml` at or above the working directory, then a per-user data directory. A root named
explicitly that does not exist is an error, never a silent fallback.

## Documentation

Reference documentation lives in [`docs/`](docs/); it is not published as a site yet.

- [Writing a job](docs/WRITING_A_JOB.md): start here to run something
- [Writing a schema](docs/WRITING_A_SCHEMA.md): point the tool at your own data
- [Time semantics](docs/TIME_SEMANTICS.md), [panel contract](docs/PANEL_CONTRACT.md),
  [figure standard](docs/FIGURE_STANDARD.md), [jobs](docs/JOBS.md),
  [gold standard](docs/GOLD_STANDARD.md)
- [The checks](docs/iid_checks/iid_checks_basics.md), one page per check, the
  [bench](docs/iid_checks/BENCH.md) and their [limitations](docs/iid_checks/LIMITATIONS.md)

## Contributing and support

Feedback is welcome, particularly from anyone who runs long characterization campaigns and
disagrees with how QUEBRA frames the problem. Open an
[issue](https://github.com/seraconti/quebra/issues) for bugs, questions, and anything that is
merely confusing. See [CONTRIBUTING.md](CONTRIBUTING.md) and the
[code of conduct](CODE_OF_CONDUCT.md).

QUEBRA is developed with AI coding agents working from written specifications, with every change
reviewed and committed by the author; the session transcripts are kept in `spec/ledger/`.

## Citing

[`CITATION.cff`](CITATION.cff) is the machine-readable form, behind GitHub's "Cite this
repository" button.

```bibtex
@software{conti_quebra_2026,
  author  = {Conti, Sera},
  title   = {{QUEBRA}: Quantum Engineering Reliability Analysis},
  version = {0.1.0},
  year    = {2026},
  url     = {https://github.com/seraconti/quebra}
}
```

## Licence

GNU General Public License v3.0 or later. See [LICENSE](LICENSE).
