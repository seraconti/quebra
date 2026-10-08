#!/usr/bin/env bash
# Run `make check` against a clean dependency resolve, the way the workflow does.
#
# `make check` runs whatever is installed in the working environment. The workflow installs
# the newest version every ">=" admits, into a fresh non-editable environment with no
# caches. Those are different inputs, so a green `make check` is not evidence about CI.
# Two failure modes this script catches and `make check` structurally cannot:
#
#   - A tool release changes behaviour. A formatter that gains a file type, or a linter
#     that gains a rule, turns a green tree red with no commit in between.
#   - An editable install resolves the package through a path hook rather than a real
#     directory, which changes which modules mypy follows and therefore whose errors it
#     reports.
#
# The mypy cache lives in the sandbox, never the repository: a cache written by an earlier
# configuration can report clean where a fresh checkout reports errors.
#
# Exit codes are printed per step rather than only propagated, so a failure names its step
# without reading the whole log.
set -u

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SANDBOX="$(mktemp -d)"
# Remove the sandbox on exit, including on failure; a venv with the full dependency set is
# several hundred megabytes. Keep it with SANDBOX_KEEP=1 when a failure needs inspecting.
cleanup() {
  if [ -n "${SANDBOX_KEEP:-}" ]; then
    echo "sandbox kept at ${SANDBOX}"
  else
    rm -rf "${SANDBOX}"
  fi
}
trap cleanup EXIT
VENV="${SANDBOX}/venv"
CACHE="${SANDBOX}/mypy_cache"
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

echo "repo:    ${REPO}"
echo "sandbox: ${SANDBOX}"

python3 -m venv "${VENV}"
"${VENV}/bin/pip" install --quiet --upgrade pip

# Non-editable, matching the workflow. `pip install -e` is the wrong test here.
step "install .[dev], clean resolve" "${VENV}/bin/pip" install --quiet "${REPO}[dev]"
if [ "${STATUS}" -ne 0 ]; then
  echo "install failed; the remaining steps would report nothing"
  exit 1
fi

echo ""
echo "resolved:"
for tool in ruff mypy deptry; do
  printf "  %-14s %s\n" "${tool}" "$("${VENV}/bin/${tool}" --version 2>&1 | head -1)"
done
printf "  %-14s %s\n" numpy "$("${VENV}/bin/python" -c 'import numpy; print(numpy.__version__)')"
printf "  %-14s %s\n" python "$("${VENV}/bin/python" -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')"

# Every tool runs from the repository root with the same arguments the Makefile uses, so a
# disagreement between this and `make check` is the environment and nothing else.
cd "${REPO}"

step "ruff check"   "${VENV}/bin/ruff" check .
step "ruff format"  "${VENV}/bin/ruff" format --check .
step "mypy"         "${VENV}/bin/mypy" --cache-dir "${CACHE}" --no-incremental src/quebra/core
step "lint-imports" "${VENV}/bin/lint-imports" --no-cache
step "deptry"       "${VENV}/bin/deptry" .
FAST_SELECTOR="$(cat "${REPO}/scripts/fast-selector.txt" 2>/dev/null || true)"
if [ -z "${FAST_SELECTOR}" ]; then
  echo "scripts/fast-selector.txt is missing or empty; pytest reads -m \"\" as no filter" >&2
  exit 1
fi
step "pytest"       "${VENV}/bin/python" -m pytest -q \
  -m "${FAST_SELECTOR}"

echo ""
echo "=================================="
echo "check-ci: exit ${STATUS}"
exit "${STATUS}"
