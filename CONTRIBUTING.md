# Contributing to QUEBRA

QUEBRA is a research toolkit built alongside an MSc thesis. Contributions are welcome, and
the three things below are what a contributor most often needs.

## How to report a bug

Open an issue at <https://github.com/seraconti/quebra/issues>.

A useful report for this project includes:

- the job file you ran, or the smallest snippet that reproduces the problem
- the full traceback, not a summary of it
- the run directory name under `output/` if a job produced one, since it encodes the
  content identity of the inputs
- whether R is installed, if the problem involves the copula check

A wrong number matters more than a crash here. If a figure or a statistic looks wrong, say
what you expected and why, and include the provenance file from the run directory.

## How to get support

Open an issue with the `question` label, or write to the maintainer address in
`CITATION.cff`. There is no chat channel and no mailing list.

Response is best effort, be kind <3

## How to contribute code

1. Open an issue first for anything larger than a typo. It is cheaper to disagree about an
   approach in an issue than in a pull request.
2. Fork, branch from `main`, and keep the branch focused on one change.
3. Run these before you push:

   ```bash
   make check     # lint, types, import contract, tests, using your installed tools
   make deps      # dependency declarations
   make check-ci  # the same steps against a clean resolve, which is what CI installs
   ```

   `check` and `check-ci` run the same steps against different inputs. CI installs the
   newest version every `>=` admits into a fresh non-editable environment, so `check`
   passing is not by itself evidence that CI will.

   `make cov` reports coverage. Nothing has a threshold and nothing fails on the number.
   Read it for one thing: a module executing code that no test checks. A high percentage
   beside no oracle is the signal, not a low one.

4. Open a pull request describing what changed and what you ran.

### What the reviewer will look for

**A test that fails without your change.** A test that passes either way documents nothing.
A test that asserts a statistical result must name its oracle: an analytic value, a
reference implementation, or a simulation truth. `tests/` is flat and stays flat: the
six-tier directory layout proposed in `spec/quebraplan.md` was declined in
`spec/spectests06.md`, and the tier is carried by a marker instead. Every test carries
exactly one of `unit`, `properties`, `statistical`, `integration`, `validation`,
`regression`, `policy`, enforced by `tests/test_marker_discipline.py`. The cost markers
`slow`, `heavy`, `real` and `r` are a separate, orthogonal axis.

**Claims that match the code.** A number in a docstring must come from the artifact it
cites. If a value cannot be checked cheaply, write it as the open question it is rather than
asserting it.

**Errors raised, not swallowed.** This is a reproducibility tool, so a silent fallback
produces a wrong-but-plausible result, which is worse than a crash.

**Units on physical quantities.** Suffixes such as `_hz`, `_rel_s`, `_unix_s` are load
bearing. Unit and clock mismatch is the costly bug class in this codebase.

**The locked vocabulary.** `window`, `read`, `bag`, `check`, `band`, `scan clock`,
`window age` and `birth type` have one meaning each and are never substituted. The two
analysis tiers are `within-calibration` and `across-calibration`.

### Style

Spaced hyphens, never em dashes, in every docstring, comment and document. Short
declarative paragraphs. Do not claim a result you have not measured.

## Code of conduct

Participation is governed by `CODE_OF_CONDUCT.md`.
