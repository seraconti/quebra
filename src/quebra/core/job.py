from __future__ import annotations

import hashlib
import importlib.util
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Callable

if TYPE_CHECKING:
    # Annotation-only (`_FigureSink.plot_class`, `Job.figure`'s PlotClass). Importing it at
    # runtime pulled `plots.base`, which imports matplotlib AND plotly at module scope, into
    # every `import quebra.core.job` - 813 modules against 648. The graph runtime has no
    # business loading a rendering stack.
    from quebra.plots.base import BasePlot

import pandas as pd

from quebra.loaders.registry import load as load_dataframe
from quebra.core.closure import code_closure, parameter_row
from quebra.core import discovery
from quebra.core.dataset import Dataset
from quebra.core.identity import Identity, content_hash
from quebra.core.paths import repo_root, resolve_dataset_path, resolve_repo_path
from quebra.core.reference import ArtifactRef, LocalRef, Reference
from quebra.core.types import Norm
from quebra.provenance import hash_string
from quebra.schemas.ramsey_series import RamseySeriesSchema


@dataclass(slots=True)
class _DAGNode:
    node_id: str
    fn: Callable[..., object]
    fn_name: str
    inputs: list[Reference]
    kwargs: dict[str, object]


@dataclass(slots=True)
class _FigureSink:
    plot_class: type[BasePlot]
    input: LocalRef
    targets: list[str]
    name: str


@dataclass(slots=True)
class _MaterializeSink:
    node: LocalRef
    name: str


# In-progress include chain (resolved absolute paths). Re-entering a path means a
# cycle: composites include each other at module-import time, so a cycle would
# otherwise recurse into exec_module forever and die with a useless RecursionError.
_IMPORT_STACK: list[Path] = []


def _rel_to_repo(p: Path) -> str:
    try:
        return str(p.relative_to(repo_root()))
    except ValueError:
        return str(p)


def _import_job(path: Path) -> Job:
    """Import a job module by file path and return its top-level `job` (with
    `job_file` attached), mirroring main._module_from_path. Used by Job.include.

    Detects include cycles (A→B→A) and self-includes at import time: re-entering a
    path already on the in-progress stack raises with the full chain."""
    resolved = path.resolve()
    if resolved in _IMPORT_STACK:
        chain = " -> ".join(_rel_to_repo(p) for p in [*_IMPORT_STACK, resolved])
        raise ValueError(f"include cycle: {chain}")
    _IMPORT_STACK.append(resolved)
    try:
        module_name = (
            f"subjob_{re.sub(r'[^A-Za-z0-9_]', '_', resolved.stem)}_"
            f"{hashlib.sha256(str(resolved).encode('utf-8')).hexdigest()[:8]}"
        )
        spec = importlib.util.spec_from_file_location(module_name, resolved)
        if spec is None or spec.loader is None:
            raise ValueError(f"Cannot import job module from {resolved}")
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)
        job = getattr(module, "job", None)
        if job is None:
            raise ValueError(f"Job module {resolved} has no top-level 'job'")
        job.job_file = resolved
        return job
    finally:
        _IMPORT_STACK.pop()


@dataclass(slots=True)
class _IncludedJob:
    """A sub-job pulled into a composite via Job.include.

    Every figure()/materialize() sink always persists its input node as a .pkl
    (ref-able by node id), regardless of `figures`. `figures` only controls whether
    a fresh nested run ALSO renders the figure sinks' PDFs (`figures=True`) or not
    (`figures=False`, the default).
    """

    alias: str
    path: Path
    job: Job
    figures: bool
    composite: Job

    def _referenceable(self) -> set[str]:
        """Node names whose results the composite can `ref` (persisted artifacts)."""
        names: set[str] = set()
        for sink in self.job.sinks:
            if isinstance(sink, _MaterializeSink):
                names.add(sink.name)
            elif isinstance(sink, _FigureSink):
                names.add(sink.input.node_id)
        return names

    def ref(self, node_name: str) -> ArtifactRef:
        """A reference to this sub-job's `node_name` artifact, usable directly as
        a composite step input. It is a structured (job locator, node id) pair -
        no node registered in the composite, no name mangling; the runner resolves
        it via the locator on the resolution context. Ownership (that this ref is
        used only in the composite that created the include) is checked at step()
        time, where a foreign ref would actually leak in."""
        available = self._referenceable()
        if node_name not in available:
            raise ValueError(
                f"Sub-job '{self.job.name}' (alias '{self.alias}') does not persist a node "
                f"'{node_name}'. Referenceable outputs: {sorted(available) or '<none>'}. "
                f"(Figure inputs are always persisted under their node id; for other "
                f"nodes, add an explicit job.materialize(node, '{node_name}') in the sub-job.)"
            )
        return ArtifactRef(included=self, node_name=node_name)


# The runner and Job.build_identity classify DAG nodes by fn.__name__; these
# internal loader names are therefore reserved (Job.step rejects them).
_LOAD_NODE_FN_NAMES = frozenset({"_load_dataset", "_load_dataframe_raw"})


def _load_dataset(dataset: Dataset) -> Norm:
    """Read a file and hand it to a schema. Nothing here is specific to any experiment.

    Requiring `timestamp` and `frequency` columns and a resolvable run start here would make
    `job.load` the Ramsey path wearing a generic name, unusable by a job about anything else.
    That work lives in `quebra.schemas.ramsey_series.RamseySeriesSchema`, which is the DEFAULT
    and not a requirement.

    Two schema shapes, dispatched below:

    - `to_norm(frame, dataset) -> Norm` owns the conversion outright. This is the
      extension point: supply one and the default never runs.
    - `validate(frame, dataset=...) -> DataFrame` checks and cleans, then the default
      normaliser runs on the result. `track912Schema` is this shape.
    """
    frame = load_dataframe(dataset.path, meta=dict(dataset.loader_kwargs))
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("Dataset loader must return a pandas.DataFrame")

    frame = frame.copy()
    # Attach declared metadata as columns so a schema can key device-specific validation
    # off them. Only when declared: an unconditional `int(dataset.qubit)` raises
    # `int() argument must be...` on `qubit=None`, naming neither `qubit` nor `Dataset`,
    # before any schema can intervene - and would force even a `to_norm` schema that never
    # looks at a qubit to declare one.
    if dataset.qubit is not None and "qubit_id" not in frame.columns:
        frame["qubit_id"] = int(dataset.qubit)
    if dataset.device is not None and "device" not in frame.columns:
        frame["device"] = dataset.device

    schema = dataset.schema
    if schema is not None:
        if hasattr(schema, "to_norm"):
            return schema.to_norm(frame, dataset)
        if hasattr(schema, "validate"):
            frame = schema.validate(frame, dataset=dataset)
        else:
            raise TypeError(
                f"{schema!r} is not a usable schema: it needs a to_norm(frame, dataset) "
                f"or a validate(frame, dataset=...) method. See docs/WRITING_A_SCHEMA.md."
            )

    return RamseySeriesSchema.to_norm(frame, dataset)


def _load_dataframe_raw(dataset: Dataset) -> pd.DataFrame:
    df = load_dataframe(dataset.path, meta=dict(dataset.loader_kwargs))
    if not isinstance(df, pd.DataFrame):
        raise TypeError("Raw loader must return a pandas.DataFrame")
    return df.copy()


def _safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    cleaned = cleaned.strip("_")
    return cleaned or "node"


class Job:
    def __init__(self, name: str) -> None:
        self.name = name
        self.dag: dict[str, _DAGNode] = {}
        self.sinks: list[_FigureSink | _MaterializeSink] = []
        self._node_counts: dict[str, int] = {}
        # Composite jobs: sub-jobs pulled in via include(); empty for normal jobs.
        self.includes: list[_IncludedJob] = []
        # Set at import time (main._module_from_path / _import_job); required to
        # compute the code hash. Memoized identity caches (compute once per Job).
        self.job_file: Path | None = None
        self._job_code_hash: str | None = None
        self._identity: Identity | None = None
        self._identity_root: Path | None = None

    def job_code_hash(self) -> str:
        """The identity's `code` contribution. Named for the three things it folds:

        1. the text of this job file,
        2. the digest of every `quebra.*` module the job's step functions can reach, by
           dotted module name and hashed from file bytes,
        3. the arguments each step was called with.

        The job file alone would not be the code that produced the result: an analyzer could
        change completely while a file-only digest asserted nothing had. Item 3 is needed
        because two family members can share a definition file and differ only by a
        parameter, which would otherwise collide.

        `quebra.core.closure` explains why item 2 is an import closure rather than the whole
        package (blast radius) and rather than `sys.modules` (under `run --all`, job N would
        inherit the union of jobs 1..N).

        Memoized so a composite and its own run don't re-read/re-hash the same source."""
        if self._job_code_hash is None:
            if self.job_file is None:
                raise ValueError("Job.job_code_hash() requires job_file to be set")
            parts = [Path(self.job_file).read_text(encoding="utf-8")]
            # Sorted by module name inside `code_closure`, so the fold order is a property of
            # the graph and not of DAG iteration order.
            nodes = list(self.dag.values())
            for name, digest in code_closure([n.fn for n in nodes]).items():
                parts.append(f"{name}={digest}")
            parts.extend(parameter_row(nodes))
            self._job_code_hash = hash_string("\n".join(parts))
        return self._job_code_hash

    def build_identity(self, dataset_root: Path) -> Identity:
        """This job's content identity, folding code + data + children (memoized).

        `data` is the content hash of every dataset the job loads (keyed by
        dataset-root-relative path); `children` are the identities of included
        sub-jobs, so a composite's identity changes iff a child's does."""
        if self._identity is not None:
            if self._identity_root != dataset_root:
                raise ValueError(
                    "Job.build_identity() called with a different dataset_root than "
                    "the memoized one; identity would be silently stale"
                )
            return self._identity

        data: dict[str, str] = {}
        for node in self.dag.values():
            if node.fn.__name__ not in _LOAD_NODE_FN_NAMES:
                continue
            ds = node.kwargs.get("dataset")
            if ds is None:
                continue
            # `_DAGNode.kwargs` is `dict[str, object]` by design - the DAG stores arbitrary
            # step arguments - so this needs narrowing. It RAISES rather than skipping: a
            # skipped load node would drop that dataset's content hash from `data` silently,
            # and an identity that quietly stops covering an input is the one failure this
            # whole scheme exists to prevent. Before the type fix a non-Dataset reached
            # `ds.path` and raised AttributeError; this keeps that loudness and names the
            # cause. Matches `runner._dataset_of`.
            if not isinstance(ds, Dataset):
                raise TypeError(
                    f"load node {node.node_id!r} carries a non-Dataset 'dataset' kwarg: "
                    f"{type(ds).__name__}. Its content hash cannot enter the identity."
                )
            path = resolve_dataset_path(ds.path, dataset_root)
            try:
                key = str(path.relative_to(dataset_root))
            except ValueError:
                # dataset outside dataset_root: fall back to the absolute path.
                # Machine-specific - flagged for the reuse-key work in a later
                # increment; does not occur for the in-repo-parent datasets.
                key = str(path)
            data[key] = content_hash(path)
        children = tuple(
            inc.job.build_identity(dataset_root).digest for inc in self.includes
        )
        self._identity = Identity(
            code=self.job_code_hash(),
            data=tuple(sorted(data.items())),
            children=children,
        )
        self._identity_root = dataset_root
        return self._identity

    def include(
        self, path: str | Path, alias: str | None = None, figures: bool = False
    ) -> _IncludedJob:
        """Pull another job in as a dependency, to `ref` one of its outputs.

        The composite runs the included job through the normal run_job (so its
        datasets/provenance resolve exactly as a standalone run); every figure()
        sink's input is always persisted for `ref`, and `figures` only controls
        whether a fresh nested run also renders the figure sinks' PDFs.
        """
        # A logical JOB_ID first, a path only as a fallback. The ID names the
        # graph node and content hashes determine identity; where the file sits is neither.
        # Resolve BEFORE _import_job so the import works from any CWD, and store the absolute
        # path on the _IncludedJob (the runner's locator resolves reuse by the sub-job's
        # identity, not by re-reading this path).
        text = str(path)
        sub_path: Path | None = None
        if "/" not in text and not text.endswith(".py"):
            jobs_root = repo_root() / "jobs"
            if jobs_root.is_dir():
                sub_path = discovery.resolve(text, jobs_root)
        if sub_path is None:
            sub_path = resolve_repo_path(path)
        if not sub_path.exists():
            raise FileNotFoundError(
                f"include not found: '{path}' resolved to '{sub_path}'. "
                "Prefer a logical JOB_ID (e.g. 't2star_q1_070423'); a repo-root-relative "
                "path still works."
            )
        sub_job = _import_job(sub_path)  # raises on an include cycle / self-include
        included = _IncludedJob(
            alias=alias or sub_job.name,
            path=sub_path,
            job=sub_job,
            figures=figures,
            composite=self,
        )
        self.includes.append(included)
        return included

    def _allocate_node_id(self, base_name: str) -> str:
        slug = _safe_name(base_name)
        count = self._node_counts.get(slug, 0) + 1
        self._node_counts[slug] = count
        node_id = slug if count == 1 else f"{slug}_{count}"
        while node_id in self.dag:
            count += 1
            self._node_counts[slug] = count
            node_id = f"{slug}_{count}"
        return node_id

    def _register_node(
        self,
        fn: Callable[..., object],
        inputs: list[Reference],
        kwargs: dict[str, object],
        base_name: str,
    ) -> LocalRef:
        # Same-job discipline applies to LocalRefs (a node's result belongs to one
        # run). ArtifactRefs are cross-job by construction - that is exactly the
        # case the old blanket `job_ref is not self` lock wrongly rejected, forcing
        # the sentinel; they are allowed here.
        for input_ref in inputs:
            if isinstance(input_ref, LocalRef):
                if input_ref.job_ref is not self:
                    raise ValueError(
                        "A LocalRef input must belong to the same Job instance"
                    )
                # Membership, not just ownership. `LocalRef` is a public dataclass, so a
                # caller can hand one in naming a node that does not exist yet. Without
                # this, a FORWARD reference is accepted and registration order stops being
                # a topological order: `job.step(fn, LocalRef(node_id="g", ...), name="f")`
                # followed by `job.step(fn, f, name="g")` builds a real cycle through public
                # calls, which `runner._toposort` then reports one node of. A reference to
                # an id never allocated at all surfaced later as a bare `KeyError` from
                # inside the traversal rather than an error where the mistake was made.
                if input_ref.node_id not in self.dag:
                    raise ValueError(
                        f"LocalRef names node {input_ref.node_id!r}, which this job has "
                        f"not registered. Inputs must be references returned by an earlier "
                        f"step; a forward reference would make the graph cyclic. "
                        f"Registered: {sorted(self.dag)}"
                    )
            elif isinstance(input_ref, ArtifactRef):
                # An ArtifactRef may only be consumed by the composite that created
                # its include - this is where a foreign ref would actually leak in.
                if input_ref.included.composite is not self:
                    raise ValueError(
                        "ArtifactRef belongs to a different composite than this job"
                    )
            else:
                raise TypeError(
                    f"step inputs must be References (LocalRef/ArtifactRef), "
                    f"got {type(input_ref).__name__}"
                )
        node_id = self._allocate_node_id(base_name)
        self.dag[node_id] = _DAGNode(
            node_id=node_id, fn=fn, fn_name=base_name, inputs=inputs, kwargs=kwargs
        )
        return LocalRef(node_id=node_id, job_ref=self, fn_name=base_name, kwargs=kwargs)

    def load(self, dataset: Dataset) -> LocalRef:
        return self._register_node(
            fn=_load_dataset,
            inputs=[],
            kwargs={"dataset": dataset},
            base_name="load",
        )

    def load_df(self, dataset: Dataset) -> LocalRef:
        """Load a dataset as a raw pandas.DataFrame node (for side tables and auxiliary files).

        This registers a node that yields a DataFrame (not the normalized mapping).
        """
        return self._register_node(
            fn=_load_dataframe_raw,
            inputs=[],
            kwargs={"dataset": dataset},
            base_name="load_df",
        )

    def step(
        self,
        fn: Callable[..., object],
        *inputs: Reference,
        name: str | None = None,
        **kwargs: object,
    ) -> LocalRef:
        if getattr(fn, "__name__", None) in _LOAD_NODE_FN_NAMES:
            raise ValueError(
                f"Step function name '{fn.__name__}' is reserved for internal "
                f"dataset-load nodes: the runner and build_identity classify DAG "
                f"nodes by fn.__name__, so a user step with this name would be "
                f"silently misclassified. Rename the function."
            )
        step_name: str = name or str(getattr(fn, "__name__", fn.__class__.__name__))
        return self._register_node(
            fn=fn, inputs=list(inputs), kwargs=dict(kwargs), base_name=step_name
        )

    def figure(
        self,
        PlotClass: type[BasePlot],
        input: LocalRef,
        targets: list[str],
        title: str = "",
    ) -> None:
        self._require_own_node(input, "Figure input")
        sink_name = _safe_name(title) if title else f"fig_{input.node_id}"
        self.sinks.append(
            _FigureSink(
                plot_class=PlotClass, input=input, targets=list(targets), name=sink_name
            )
        )

    def _require_own_node(self, ref: object, what: str) -> None:
        """Ownership AND membership, for every site that accepts a node reference.

        `LocalRef` is a public dataclass, so a caller can build one naming a node this job
        never registered. Checking ownership alone let a ghost sink reference through: the
        run then created its timestamped output directory before dying as a `KeyError` from
        inside the traversal. Both SINK sites call this. `_register_node` still carries its
        own inlined pair of raises, because its messages name the offending input rather than
        the sink; the two paths are therefore consistent by review, not by construction.
        """
        if not isinstance(ref, LocalRef) or ref.job_ref is not self:
            raise ValueError(f"{what} must be a same-job node (LocalRef)")
        if ref.node_id not in self.dag:
            raise ValueError(
                f"{what} names node {ref.node_id!r}, which this job has not registered. "
                f"Registered: {sorted(self.dag)}"
            )

    def materialize(self, node: LocalRef, name: str) -> None:
        self._require_own_node(node, "Materialized node")
        self.sinks.append(_MaterializeSink(node=node, name=_safe_name(name)))
