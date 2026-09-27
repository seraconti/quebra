"""Session-wide setup for the test suite: import paths, hypothesis, and two shared probes.

`quebra` itself is NOT provided here. It comes from the installed distribution, which is
the entire point of dropping `pythonpath` when `pyproject.toml` landed: if the suite could
import the package from `src/`, `pip install` would never actually be exercised and the
claim this phase exists to make would go unverified.

What this file does provide is the two trees that are deliberately NOT in the wheel:

  jobs/   the researcher's analysis configuration. It stays at the repository root so that
          `output/` is never written inside site-packages and so a study's dependencies
          (joblib) never become dependencies of the toolkit.
  tests/  the suite itself, which imports its own fixtures package.

pytest's default `prepend` import mode puts the directory containing the rootdir conftest
on `sys.path`, so simply existing here is enough. It is spelled out rather than left
implicit because a reader coming from the deleted `pytest.ini` will look for the setting
that replaced `pythonpath` and should find this explanation instead of a bare file.

Beyond the import paths this file holds two things the suite shares:

  the collection record   every item as it was collected, before any selector deselected
                          it, so a guard that is itself selected sees the whole run rather
                          than only what its own `-m` expression kept.
  `requires_rscript`      the single Rscript probe every `r`-marked test uses.
"""

from pathlib import Path

import pytest
from hypothesis import HealthCheck, Verbosity, settings

REPO_ROOT = Path(__file__).resolve().parent

# Hypothesis is DERANDOMISED here, and that is a contract rather than a preference.
#
# This repository treats an unseeded generator as a defect: `checks/_permutation` raises on
# a `None` rng because defaulting to OS entropy made three consecutive runs on identical
# input return p = 0.3860 / 0.4040 / 0.3790 while the provenance record stayed unchanged.
# A property test that draws fresh examples each run is the same failure in a different
# place - a suite that passes today and fails tomorrow with no diff between them.
#
# `derandomize=True` makes each test's example set a deterministic function of its own
# source, and the example database is disabled so a local `.hypothesis/` cache cannot make
# one machine's run differ from another's.
settings.register_profile(
    "quebra",
    settings(
        derandomize=True,
        database=None,
        max_examples=100,
        deadline=None,  # a KM fit on a large draw is slower than the default 200 ms
        suppress_health_check=[HealthCheck.function_scoped_fixture],
        verbosity=Verbosity.normal,
    ),
)
settings.load_profile("quebra")


@pytest.fixture
def in_repo(monkeypatch):
    """Run the test with the repository as the working directory.

    Needed because `repo_root()` and `resolve_data_root()` are now defined BY the working
    directory: `repo_root()` walks up looking for a project marker, and `resolve_data_root()`
    walks up looking for `quebra.toml`. That was the fix for an installed wheel reporting
    its repo root as `<venv>/lib/python3.13`, and it is correct - but it means a test that
    asserts something about THIS repository has to say which directory it means.

    Without it the suite passed only when pytest happened to be invoked from the repository.
    Measured: six tests failed under `cd /tmp && pytest <repo>/tests`, which is precisely
    the "unrelated working directory" the installability phase claims to support. The
    acceptance script pinned the cwd itself, so the gate hid the gap instead of catching it.

    `Path(__file__).parent` and not `Path.cwd()`: this file's location is stable, which is
    the whole point.

    FUTURE: worth revisiting whether `repo_root()` and `resolve_data_root()` should take the
    starting directory as an argument, defaulted to the cwd, instead of reading it from
    global state. Call sites that need a specific root could then say so, and this fixture
    would not be necessary. That is a signature change across ~17 call sites, so it is not
    a small edit.
    """
    monkeypatch.chdir(REPO_ROOT)
    return REPO_ROOT


# Every collected item, recorded before anything can deselect it.
#
# Marker and keyword deselection both happen in `pytest_collection_modifyitems`, so
# `session.items` is the post-deselection list and a guard reading it is only as total as
# the selector that ran. `pytest_itemcollected` fires per item during collection itself,
# which is strictly earlier than any `modifyitems` implementation, so the record it builds
# does not depend on hook ordering between plugins.
COLLECTED_PRE_DESELECTION: pytest.StashKey[list[pytest.Item]] = pytest.StashKey()


def pytest_itemcollected(item: pytest.Item) -> None:
    """Append one collected item to the session's pre-deselection record."""
    item.session.stash.setdefault(COLLECTED_PRE_DESELECTION, []).append(item)


@pytest.fixture(scope="session")
def collected_pre_deselection(request) -> list[pytest.Item]:
    """Every item collected from under this conftest's tree, deselected ones included.

    Kept on the session stash rather than in a module global so two sessions in one
    process cannot pool their items.
    """
    return request.session.stash.setdefault(COLLECTED_PRE_DESELECTION, [])


@pytest.fixture(scope="session")
def requires_rscript() -> str:
    """Absolute path to `Rscript`, or a skip when the machine has none.

    The one probe every `r`-marked test shares, so the tier does not accumulate a private
    copy of it per file. It defers to the bridge's own `rscript_path()`, which is what
    stops the suite and the pipeline disagreeing about whether R is present.

    A missing interpreter skips. A present interpreter without the `copula` package does
    NOT skip: `c3_serial_copula.R` stops on a missing `copula` and the bridge raises, and
    that is the failure the `r` marker promises to report rather than hide. Never satisfy
    this fixture with a mocked path: a test that goes green without R is evidence about
    the mock.
    """
    from quebra.analyzers.checks.c3_serial_copula import rscript_path

    executable = rscript_path()
    if executable is None:
        pytest.skip("Rscript is not on PATH, and the `r` tier needs a real interpreter")
    return executable
