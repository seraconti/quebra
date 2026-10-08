#!/usr/bin/env bash
# Clean-environment acceptance: install and run from outside the repository.
#
# The claim this phase makes is "a reviewer runs the documented install steps on a machine
# that has never seen this repository". This script is the only thing that can falsify it:
# every other test run uses the editable install, where `src/` is on the path and a
# packaging mistake is invisible.
#
# The virtualenv is created OUTSIDE the repository on purpose. Inside it, `pip install`
# might resolve the source tree instead of the wheel, and the test would pass while proving
# nothing.
#
# Exit codes are printed per step rather than only propagated, so a failure says which step
# failed without reading the whole log.
set -u

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SANDBOX="$(mktemp -d)"
# Remove the sandbox on exit, including on failure. Without this every run leaks a ~620 MB
# venv into /tmp: twenty accumulated during one working session and filled a 12 GB tmpfs,
# at which point `pip install` inside the NEXT run died with ENOSPC and the script reported
# a packaging failure that was really a disk failure. Keep it with SANDBOX_KEEP=1 when a
# failure needs inspecting.
cleanup() {
  if [ -n "${SANDBOX_KEEP:-}" ]; then
    echo "sandbox kept at ${SANDBOX}"
  else
    rm -rf "${SANDBOX}"
  fi
}
trap cleanup EXIT
VENV="${SANDBOX}/venv"
# Somewhere that is emphatically not the repository, and has no project marker of its own.
OUTSIDE="${SANDBOX}/elsewhere"
mkdir -p "${OUTSIDE}"
STATUS=0

step() {
  local name="$1"; shift
  echo ""
  echo "=== ${name} ==="
  "$@"
  local code=$?
  echo "--- ${name}: exit ${code}"
  [ "${code}" -ne 0 ] && STATUS=1
  return 0
}

echo "repo:  ${REPO}"
echo "venv:  ${VENV}"
echo "cwd:   ${OUTSIDE}"

step "build the wheel" env -C "${REPO}" python3 -m build --outdir "${REPO}/dist"

WHEEL="$(ls -t "${REPO}"/dist/*.whl 2>/dev/null | head -1)"
if [ -z "${WHEEL}" ]; then
  echo "no wheel produced; stopping"
  exit 1
fi
echo "wheel: ${WHEEL}"

step "create the venv"  python3 -m venv "${VENV}"
step "install the wheel" "${VENV}/bin/pip" install --quiet "${WHEEL}"
step "install pytest"    "${VENV}/bin/pip" install --quiet pytest hypothesis
# hypothesis, not only pytest: the repository conftest registers a settings profile at
# import time, so without it collection fails here and nowhere else.

# Imported from a directory that is NOT the repository, so a stray relative path or a
# leftover `sys.path` entry cannot mask a packaging error.
step "import quebra from elsewhere" env -C / "${VENV}/bin/python" -c \
  "import quebra, quebra.core.paths, quebra.cli; print('quebra', quebra.__file__)"

step "the console script exists" env -C / "${VENV}/bin/quebra" --help

# A reviewer with none of our data must still reach a real result. This runs
# the packaged fixture through the actual analyzer from `/`, so it fails if the CSV did not
# ship in the wheel or if the fixture is reached by `__file__` arithmetic rather than
# importlib.resources.
step "a packaged fixture runs the pipeline from elsewhere" env -C / "${VENV}/bin/python" -c "
import quebra.analyzers.t2star as t2star
from quebra._fixtures import fixture_path
from quebra.core.dataset import Dataset
from quebra.core.job import _load_dataset
norm = _load_dataset(Dataset(path=fixture_path('ramsey_synthetic.csv'), qubit=1,
                             extra={'run_start_unix_s': 1.7e9}))
result = t2star.run(t2star.make_inputs_from_norm(norm))
assert len(result.frame) == 200, result.frame.shape
print('fixture ->', len(result.frame), 'T2* points')
"

# The suite lives in the checkout - `tests/` is deliberately not in the wheel - but it runs
# from a DIRECTORY THAT IS NOT THE REPOSITORY, because that is the claim the spec makes:
# "import quebra works from an unrelated working directory, the fast test suite passes
# there". An earlier version of this script ran the step with `env -C "${REPO}"` and passed
# while six tests failed from anywhere else, so the test concealed exactly the gap it was
# written to expose. Pinning the cwd to the repository is the one thing this step must not
# do.
FAST_SELECTOR="$(cat "${REPO}/scripts/fast-selector.txt" 2>/dev/null || true)"
if [ -z "${FAST_SELECTOR}" ]; then
  echo "scripts/fast-selector.txt is missing or empty; pytest reads -m \"\" as no filter" >&2
  exit 1
fi
step "fast test selector from an unrelated cwd" env -C "${OUTSIDE}" "${VENV}/bin/python" -m pytest \
  -m "${FAST_SELECTOR}" -q "${REPO}/tests"

echo ""
echo "=================================="
echo "acceptance: exit ${STATUS}"
exit "${STATUS}"
