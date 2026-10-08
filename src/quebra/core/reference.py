"""References: resolvable step inputs.

A step input is a `Reference` - something the runner turns into a concrete value
via `resolve(context)`. Two kinds share that one contract:

  - `LocalRef`  - another node's result in the SAME job's run (generalizes the
    old NodeHandle): resolve = read it from the run's results map.
  - `ArtifactRef` - a materialized artifact of an *included* job's node: resolve =
    locate (and if absent, produce) that artifact via an injected locator.

This is the self-similarity the composite layer is built on: a composite consumes
a sub-job exactly the way a step consumes an upstream result - both are References
in the same input list, resolved through the same call. The locator is injected on
the `ResolutionContext` (not baked into `ArtifactRef`), so the identity-keyed
dir-glob strategy can be swapped for another without touching this module.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Protocol

if TYPE_CHECKING:
    from quebra.core.job import Job, _IncludedJob


class Reference(Protocol):
    """A resolvable step input. `resolve(context)` returns the typed value."""

    def resolve(self, context: ResolutionContext) -> object: ...


@dataclass(slots=True)
class LocatedArtifact:
    """The outcome of locating one included artifact: its value plus the
    provenance facts the composite records about how it was obtained."""

    value: object
    artifact_hash: str
    mode: str  # "cached" | "fresh"
    prov_dir_rel: str


# A locator turns an ArtifactRef into a LocatedArtifact, given the run context.
# Injected via ResolutionContext.locate so the strategy is swappable.
ArtifactLocator = Callable[["ArtifactRef", "ResolutionContext"], LocatedArtifact]


@dataclass(slots=True)
class ResolutionContext:
    """Everything a Reference needs to resolve during one run.

    `results` backs LocalRef; the run dirs + flags + injected `locate` back
    ArtifactRef, whose located artifacts are memoized in `artifacts` (keyed by
    the ref object's id) so a referenced sub-job runs at most once per ref.
    """

    results: dict[str, object]
    locate: ArtifactLocator | None = None
    artifacts: dict[int, LocatedArtifact] = field(default_factory=dict)
    subjobs_dir: Path | None = None
    job_out_dir: Path | None = None
    # The single top-level output pool, threaded UNCHANGED through nested runs so a
    # locator at any depth searches (and reuses from) the same pool.
    pool_root: Path | None = None
    reuse_deps: bool = False
    data_root: Path | None = None
    # Resolved dataset root + this run's reuse-rule inputs (an included artifact is
    # reusable only if its identity AND commit match and the tree is clean).
    dataset_root: Path | None = None
    git_commit: str = "nogit"
    tree_clean: bool = False


@dataclass(slots=True)
class LocalRef:
    """Reference to another node's result within the same job (was NodeHandle).

    `fn_name`/`kwargs` are carried only for provenance and inspect labels; the
    value comes from the run's results map keyed by `node_id`.
    """

    node_id: str
    job_ref: Job
    fn_name: str
    kwargs: dict[str, object]

    def resolve(self, context: ResolutionContext) -> object:
        return context.results[self.node_id]


@dataclass(slots=True)
class ArtifactRef:
    """Reference to a materialized artifact of an included job's node.

    The structured (job locator, node id) pair: `included` carries the job
    locator (alias, sub-job, source path), `node_name` the node id within it.
    The locator strategy lives on the context, not here, so it can be swapped.
    """

    included: _IncludedJob
    node_name: str

    def resolve(self, context: ResolutionContext) -> object:
        key = id(self)
        located = context.artifacts.get(key)
        if located is None:
            if context.locate is None:
                raise RuntimeError(
                    "ArtifactRef resolved without a locator on the context"
                )
            located = context.locate(self, context)
            context.artifacts[key] = located
        return located.value
