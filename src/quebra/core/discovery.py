"""Find jobs and read their declared identity and family WITHOUT importing them.

Two things that could be encoded in the filesystem are instead
declared in the file:

- **`JOB_ID`** — the logical name. `job.include("jobs/active/t2star_q1_070423.py")` made every
  reorganisation a breaking change to every composite that included it.
- **`JOB_FAMILY`** — the category. A directory would mean recategorising costs moving
  a file, which under logical IDs must stop being meaningful.

**Nothing here imports a job.** Importing a job file BUILDS its graph — that is the whole
design of the tool — so a discovery pass that imported would execute every job in the tree
just to read two strings, and one broken job would take the listing down with it. The
constants are read statically with `ast`, the way `tests/test_bench_isolation.py` already
walks sources without importing them.

**Not `pkgutil`**, which `quebraplan.md` 3.5 suggests. It enumerates *importable* packages,
which contradicts the no-import rule above, and `jobs/` is deliberately outside the wheel
with `sys.path` manipulation forbidden, so it could not reach the
tree from an installed copy anyway.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path

JOB_ID_CONSTANT = "JOB_ID"
JOB_FAMILY_CONSTANT = "JOB_FAMILY"
JOB_SWEEP_CONSTANT = "JOB_SWEEP"


@dataclass(frozen=True, slots=True)
class DiscoveredJob:
    """One job file's declared identity, read without executing it."""

    job_id: str
    family: str | None
    path: Path
    # Whether a bare `run --all` includes this job. Composites `include` other jobs, so
    # sweeping them re-runs every sub-job, and `independence_survey` re-runs the R bridge on
    # every eligible cell. The total has not been measured.
    # Declared, not encoded in the DIRECTORY (`jobs/active` swept, `jobs/composite` not).
    # Dropping the distinction entirely silently widens `--all` from 9 jobs to 13.
    sweep: bool = True


class DuplicateJobId(ValueError):
    """Two files claim the same `JOB_ID`.

    An error, not a last-one-wins: the ID is what `include` resolves against, so a silent
    winner would make a composite depend on whichever file the walk happened to reach first.
    """


def _bool_constant(tree: ast.Module, name: str) -> bool | None:
    """The value of a module-level `NAME = True/False`, or None if not declared."""
    for node in tree.body:
        if isinstance(node, ast.Assign):
            targets = [t.id for t in node.targets if isinstance(t, ast.Name)]
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            targets = [node.target.id]
        else:
            continue
        if name in targets:
            value = node.value
            if isinstance(value, ast.Constant) and isinstance(value.value, bool):
                return value.value
    return None


def _string_constant(tree: ast.Module, name: str) -> str | None:
    """The value of a module-level `NAME = "literal"`, or None.

    Only a plain string literal counts. A computed value would have to be executed to be
    read, which is exactly what this module exists to avoid.
    """
    for node in tree.body:
        if isinstance(node, ast.Assign):
            targets = [t.id for t in node.targets if isinstance(t, ast.Name)]
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            targets = [node.target.id]
        else:
            continue
        if name in targets:
            value = node.value
            if isinstance(value, ast.Constant) and isinstance(value.value, str):
                return value.value
    return None


class UnreadableJobFile(ValueError):
    """A file under `jobs/` could not be parsed.

    Raised rather than skipped. The old directory glob handed every file to the importer, so
    a syntax error surfaced as a failed job in the batch; swallowing it here would make a
    broken job silently VANISH from `--all` instead - a job that stops running and says
    nothing is the worst of the three outcomes.
    """


def read_declaration(path: Path) -> DiscoveredJob | None:
    """Read one job file's declarations. None if it declares no `JOB_ID`.

    A file that cannot be PARSED raises; a file that parses and declares no `JOB_ID` is not a
    job and is skipped. The two are different and only the second is silent.
    """
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except SyntaxError as exc:
        raise UnreadableJobFile(
            f"{path} could not be parsed and so cannot be discovered: {exc}. A job that "
            f"stops running without saying so is worse than one that fails."
        ) from exc
    except OSError as exc:
        raise UnreadableJobFile(f"{path} could not be read: {exc}") from exc
    job_id = _string_constant(tree, JOB_ID_CONSTANT)
    if job_id is None:
        return None
    sweep = _bool_constant(tree, JOB_SWEEP_CONSTANT)
    return DiscoveredJob(
        job_id=job_id,
        family=_string_constant(tree, JOB_FAMILY_CONSTANT),
        path=path,
        sweep=True if sweep is None else sweep,
    )


def discover(jobs_root: Path) -> dict[str, DiscoveredJob]:
    """Every declared job under `jobs_root`, keyed by `JOB_ID`.

    Files without a `JOB_ID` are skipped silently — `jobs/bench/` holds study modules that are
    not jobs. Files WITH one that collide raise.
    """
    found: dict[str, DiscoveredJob] = {}
    for path in sorted(jobs_root.rglob("*.py")):
        if path.name == "__init__.py":
            continue
        declared = read_declaration(path)
        if declared is None:
            continue
        if declared.job_id in found:
            first = found[declared.job_id].path
            raise DuplicateJobId(
                f"JOB_ID {declared.job_id!r} is claimed by two files:\n"
                f"  {first}\n  {declared.path}\n"
                f"An ID is what `include` resolves against; a silent winner would make a "
                f"composite depend on whichever file the walk reached first."
            )
        found[declared.job_id] = declared
    return found


def resolve(job_id: str, jobs_root: Path) -> Path:
    """The file declaring `job_id`, or a ValueError naming what was available."""
    found = discover(jobs_root)
    if job_id not in found:
        raise ValueError(
            f"no job declares JOB_ID {job_id!r} under {jobs_root}. "
            f"Declared: {', '.join(sorted(found)) or '<none>'}"
        )
    return found[job_id].path


def swept(jobs_root: Path) -> list[DiscoveredJob]:
    """Jobs a bare `run --all` should run: everything not opting out via `JOB_SWEEP = False`."""
    return sorted(
        (j for j in discover(jobs_root).values() if j.sweep), key=lambda j: j.job_id
    )


def by_family(jobs_root: Path, family: str) -> list[DiscoveredJob]:
    """Every job declaring `JOB_FAMILY == family`, sorted by ID."""
    return sorted(
        (j for j in discover(jobs_root).values() if j.family == family),
        key=lambda j: j.job_id,
    )
