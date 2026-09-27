# mypy is scoped to the layering root AGENTS.md names, not the whole package. Widening it
# is deferred: the modules outside this tree do not type-check yet.
PKG := src/quebra/core
# One definition, read by this file and by scripts/{check_ci,acceptance}.sh. Spelled in
# three places it drifts, and check_ci.sh:70 claims all three run the same arguments.
# `real` is deliberately IN the fast selector: CI has no private tree, so those tests skip
# there, and excluding them locally would move manifest integrity onto a discipline promise.
#
# Resolved against THIS file, not the cwd, and empty is a hard error: pytest reads `-m ""`
# as "no filter", so a missing file would silently run the excluded tiers and still exit 0.
HERE := $(dir $(lastword $(MAKEFILE_LIST)))
FAST_SELECTOR = $(strip $(shell cat $(HERE)scripts/fast-selector.txt 2>/dev/null))
FAST = -m "$(FAST_SELECTOR)"

## Guard, in the recipe rather than at parse time: `$(error ...)` at the top level fires for
## EVERY target, so a missing selector would block `make clean` and `make lint` too - the
## targets you need to recover with.
require-selector:
	@test -n '$(FAST_SELECTOR)' || { \
	  echo "scripts/fast-selector.txt is missing or empty." >&2; \
	  echo "pytest reads -m \"\" as NO FILTER, so this would silently run the excluded tiers." >&2; \
	  exit 1; }
TIERS := unit properties statistical integration validation regression policy

.PHONY: check check-ci lint types arch deps test require-selector test-all test-r test-real tiers cov promote docs clean $(addprefix tier-,$(TIERS))

## Run before every checkpoint. This is what the checkpoint banner reports.
check: lint types arch test

## Run the same gate against a clean dependency resolve in a throwaway environment, which
## is what the workflow does. `check` uses the installed tools, so it cannot see a failure
## caused by a newer release or by the install being non-editable. Run before pushing.
check-ci:
	bash scripts/check_ci.sh

lint:
	ruff check .
	ruff format --check .

types:
	mypy $(PKG)

## Dependency direction. The contract lives in [tool.importlinter]; 
##
## `--no-cache` needed: `.import_linter_cache` can serve a graph built before a
## module MOVED, and the cache file is then newer than the sources it misdescribes, so
## nothing invalidates it
arch:
	lint-imports --no-cache

deps:
	deptry .

test: require-selector
	pytest $(FAST)

## Simulates a reviewer with no private tree. Not a superset of `test`: `test` now
## includes `real`, so neither selector contains the other.
test-all:
	pytest -m "not real"

## Local only today; the CI job is deferred. Needs Rscript plus copula, and also XICOR,
## energy and randtests: one test checks the committed fixture's recorded versions against
## the local R and asks for all four.
test-r:
	pytest -m "r"

## Local only: `real` needs data/real_private/, and `regression` compares against values
## pinned from a run over it, so neither may execute on a public runner.
##
## The target line is load-bearing. Without it this recipe attaches to `test-r`, and the
## private-data selector rides along with whatever calls the R target. That exposure is
## latent, not live: `pytest -m "r"` exits 5 first and make stops on it. It goes live the
## moment any test carries the `r` marker.
##
## `real` now has 9 members, so this selects rather than exiting 5. `regression` is a
## declared marker with no members by decision: its oracle needs a tagged run, and whether
## values derived from an embargoed record may enter a public history is not settled.
test-real:
	pytest -m "real or regression"

## Per-tier selection. The tier axis says what question a test answers; nothing in `check`
## gates on it, so these are for reading the suite, not for CI.
##
## Three tiers are empty and `tier-<name>` exits 5 for each, deliberately. `properties`
## gains members with the hypothesis work. `validation` is empty because the one place the
## wiring and the estimator can jointly be wrong is the carve-to-estimator handoff, tested
## directly rather than through a job. `regression` is empty because its oracle needs a
## tagged run, and whether values derived from an embargoed record may enter a public
## history is unsettled.
$(addprefix tier-,$(TIERS)): tier-%:
	pytest -m "$*"

## How many tests answer each kind of question.
tiers:
	@for t in $(TIERS); do \
	  printf "%-14s %s\n" "$$t" "$$(pytest --collect-only -q -m "$$t" 2>/dev/null | tail -1)"; \
	done

## Reported, never gated: not in `check`, no `fail_under`. Reads whether a module executes
## code nothing checks. Via `pytest --cov` because `pytest-cov` is declared and `coverage`
## is not. The number moves with the selector and the private tree, so both are printed.
cov: require-selector
	@echo "selector:          $(FAST_SELECTOR)"
	@test -d data/real_private && echo "data/real_private: present" || echo "data/real_private: absent"
	pytest $(FAST) --cov=quebra --cov-report=term

## Commit the provenance of a figure that appears in a publication. output/ is gitignored,
## so this copies the kilobytes that make a figure auditable and none of the megabytes.
##   make promote RUN=output/<run-dir> NOTE=""
promote:
	@test -n "$(RUN)" || { echo "usage: make promote RUN=output/<run-dir> NOTE=\"where it appears\""; exit 2; }
	@test -n "$(NOTE)" || { echo "NOTE is required: say where the figure appears"; exit 2; }
	python3 scripts/promote_run.py "$(RUN)" --note "$(NOTE)" $(PROMOTE_FLAGS)

## future integration of sphinx, currently unsupported
docs:
	sphinx-build -W --nitpicky -b html docs docs/_build/html

clean:
	rm -rf .import_linter_cache
	find . -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete
