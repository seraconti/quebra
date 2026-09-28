from __future__ import annotations

import dataclasses
import json
import pickle
import sys
from datetime import datetime
from pathlib import Path

from quebra.core.job import Job, _FigureSink, _DAGNode, _LOAD_NODE_FN_NAMES
from quebra.core.identity import content_hash
from quebra.core.dataset import Dataset
from quebra.core.paths import repo_root, resolve_data_root, resolve_dataset_path
from quebra.core.reference import (
    ArtifactRef,
    LocalRef,
    LocatedArtifact,
    ResolutionContext,
)
from quebra.provenance import (
    build_prov_record,
    get_git_commit,
    hash_file,
    is_tree_clean,
    save_prov,
)
from quebra.core._artifact_guard import load_artifact
from quebra.plots.targets import RENDER_TARGETS


def _require(value: Path | None, name: str) -> Path:
    """Narrow an optional context field that this code path requires.

    `ResolutionContext`'s run-directory fields are optional because a context can be built
    for LocalRef resolution alone. Every ArtifactRef path populates them, so a None here is a
    programming error and deserves to say which field was missing.
    """
    if value is None:
        raise ValueError(
            f"ResolutionContext.{name} is required to resolve an artifact reference but was "
            f"not set. This is a runner bug: run_job populates it before resolution."
        )
    return value


def _dataset_of(node: _DAGNode) -> Dataset:
    """The `dataset` kwarg of a load node, narrowed.

    `_DAGNode.kwargs` is `dict[str, object]` because the DAG stores arbitrary step arguments;
    the value model is deliberately untyped. A load node whose `dataset`
    is not a Dataset is a defect, so this raises rather than silently skipping.
    """
    dataset = node.kwargs.get("dataset")
    if not isinstance(dataset, Dataset):
        raise TypeError(
            f"load node carries a non-Dataset 'dataset' kwarg: {type(dataset).__name__}"
        )
    return dataset


def _job_file(job: Job) -> Path:
    # Prefer explicit job_file attached at import time (main._module_from_path sets this).
    explicit = getattr(job, "job_file", None)
    if explicit is not None:
        return Path(explicit).resolve()

    # Fallback: attempt to locate via job.__module__ (legacy behavior)
    module = sys.modules.get(getattr(job, "__module__", ""))
    module_file = getattr(module, "__file__", None) if module is not None else None
    if module_file is None:
        raise ValueError(
            "Cannot resolve job file: attach 'job.job_file' when importing job modules"
        )
    return Path(module_file).resolve()


def _format_step(node: _DAGNode) -> str:
    if node.kwargs:
        kwargs = ", ".join(f"{key}={value!r}" for key, value in node.kwargs.items())
        return f"{node.fn_name}({kwargs})"
    return node.fn_name


def _fmt_dataset_path(p: Path, dataset_root: Path) -> str:
    """Stringify a dataset path dataset_root-relative when possible, else absolute."""
    try:
        return str(p.relative_to(dataset_root))
    except Exception:
        return str(p)


def _dataset_load_nodes(job: Job) -> dict[str, _DAGNode]:
    """Load nodes (main + companions) that carry a Dataset kwarg, by node id."""
    return {
        node_id: node
        for node_id, node in job.dag.items()
        if node.fn.__name__ in _LOAD_NODE_FN_NAMES
        and node.kwargs.get("dataset") is not None
    }


def _toposort(job: Job, root_ids: list[str]) -> list[str]:
    ordered: list[str] = []
    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(node_id: str) -> None:
        if node_id in visited:
            return
        if node_id in visiting:
            raise ValueError(f"Cycle detected in job DAG at node '{node_id}'")
        node = job.dag[node_id]
        visiting.add(node_id)
        for input_ref in node.inputs:
            # Only LocalRefs point at in-job nodes; ArtifactRefs are external roots.
            if isinstance(input_ref, LocalRef):
                visit(input_ref.node_id)
        visiting.remove(node_id)
        visited.add(node_id)
        ordered.append(node_id)

    for node_id in root_ids:
        visit(node_id)
    return ordered


def _ancestors(job: Job, root_id: str) -> set[str]:
    reachable: set[str] = set()

    def visit(node_id: str) -> None:
        if node_id in reachable:
            return
        reachable.add(node_id)
        for input_ref in job.dag[node_id].inputs:
            if isinstance(input_ref, LocalRef):
                visit(input_ref.node_id)

    visit(root_id)
    return reachable


def _check_render_targets(job: Job, included_by: tuple[str, ...] = ()) -> None:
    """Refuse a figure sink that names a render target nobody registered, before any work.

    Looked up at render time instead, the name fails as a bare KeyError after every step has
    run. A target exists only once `plots.targets.register_target` has registered it. A
    composite checks every job it includes, at any depth, whether or not it renders their
    figures: a sub-job that declares a missing target is a broken job file either way, and
    finding out after the composite has written its run dir is too late.
    """
    for inc in job.includes:
        _check_render_targets(inc.job, (*included_by, f"{job.name!r} as {inc.alias!r}"))
    for sink in job.sinks:
        if not isinstance(sink, _FigureSink):
            continue
        unknown = sorted(set(sink.targets) - set(RENDER_TARGETS))
        if unknown:
            via = (
                f", included by {' via '.join(reversed(included_by))}"
                if included_by
                else ""
            )
            raise ValueError(
                f"figure {sink.name!r} in job {job.name!r}{via} asks for render target(s) "
                f"{unknown}, which are not registered; the registered targets are "
                f"{sorted(RENDER_TARGETS)}. A target exists only once "
                "plots.targets.register_target has registered it."
            )


def _check_sink_artifact_names(job: Job) -> None:
    """Reject two sinks that would write the same artifact name from different nodes.

    Both sink kinds are keyed on `name`, and they share one namespace: a materialize writes
    `{name}.pkl`, a figure renders `{name}_{target}.pdf`, and BOTH write
    `provenance/{name}.prov.json`. So a figure title and a materialize name that coincide
    overwrite one another's provenance record even though their artifacts do not collide.

    Keying a figure on its SOURCE NODE instead cannot detect this: with the basename and the
    comparison value both taken from `sink.input.node_id`, the mismatch test is `x != x` and
    the check is inert for every figure.

    `name` is already safe-named - `_safe_name` maps each run of non-`[A-Za-z0-9_.-]` to a
    single `_` - which is what makes the collision reachable: the titles "a b" and "a_b" read
    as two figures and land on one filename.

    Keyed on `(kind, source node)`, not on the source node alone, because the KIND changes the
    record. A figure's record carries `targets_rendered` and a figure label; a materialize's
    carries neither. So a figure and a materialize sharing a name overwrite each other even
    when they read the SAME node - the second declared wins, and a run that rendered PDFs ends
    up with a record asserting it rendered none.

    Rejected at run start, before any output dir exists. Two sinks of the same kind on the
    same node are genuinely harmless: that is one artifact requested twice.
    """
    seen: dict[str, tuple[str, str]] = {}
    for sink in job.sinks:
        if isinstance(sink, _FigureSink):
            entry = ("figure", sink.input.node_id)
        else:
            entry = ("materialize", sink.node.node_id)
        existing = seen.get(sink.name)
        if existing is not None and existing != entry:
            raise ValueError(
                f"sink artifact name collision: '{sink.name}' would be written by two "
                f"different sinks - {existing[0]} of '{existing[1]}' and {entry[0]} of "
                f"'{entry[1]}'. They share provenance/{sink.name}.prov.json, so the second "
                f"would overwrite the first's record. Rename one sink."
            )
        seen[sink.name] = entry


def _expected_sink_pkls(job: Job) -> set[str]:
    """The `.pkl` files a COMPLETE run of this job leaves in its output directory.

    Mode-independent, which is what makes it usable as a completeness test: a figure sink
    persists its input under the SOURCE NODE's id whether or not a PDF is rendered, and a
    materialize persists under its own name.

    The sink loop writes these one at a time, so a run that raised on its third sink still
    holds a perfectly readable provenance record from its first. Identity, commit and
    tree-clean all match on that record, so without this test the gate reads a half-written
    run as reusable and every later invocation skips the job - leaving the missing sinks
    permanently unproduced. `run --all` makes it likely, because it catches per job and
    carries on, so the partial directory survives the batch.
    """
    names: set[str] = set()
    for sink in job.sinks:
        if isinstance(sink, _FigureSink):
            names.add(f"{sink.input.node_id}.pkl")
        else:
            names.add(f"{sink.name}.pkl")
    return names


def _read_prov_reuse_fields(
    run_dir: Path,
) -> tuple[str | None, str | None, bool]:
    """Read (identity, git_commit, tree_clean_at_build) from any prov record in
    `run_dir` (all sinks of a run share them), or (None, None, False) if
    unreadable. Missing tree_clean (older records) reads as False → not reusable."""
    prov_files = sorted((run_dir / "provenance").glob("*.prov.json"))
    if not prov_files:
        return None, None, False
    try:
        record = json.loads(prov_files[0].read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None, None, False
    return (
        record.get("identity"),
        record.get("git_commit"),
        bool(record.get("tree_clean", False)),
    )


def _reuse_eligible_dir(
    candidates: list[Path], identity: str, git_commit: str, tree_clean: bool
) -> Path | None:
    """The newest candidate run reusable under the gate, else None (→ re-run fresh).

    Reuse requires ALL of: matching content identity, matching git commit, the
    consumer's tree clean NOW, and the artifact PRODUCED on a clean tree. Commit +
    clean-tree (both ends) stand in for a dependency hash on imported/shared code
    (which the identity's code component does not capture): they guarantee the
    artifact was built from - and is being reused under - exactly the committed
    code. Any mismatch re-runs - commit-safe, not content-safe against code edits.
    `candidates` is sorted oldest→newest.

    `"nogit"` never matches, even against itself. It is the sentinel for "no commit could be
    read", and two runs that both failed to read a commit are not two runs at the same commit.
    The pair `("nogit", tree_clean=True)` is reachable - a repository with no commits answers
    `status --porcelain` cleanly while `rev-parse HEAD` fails - and admitting it would reduce
    the gate to an identity match alone. That matters because commit-plus-clean-tree is the
    stand-in for the code the identity does not cover: a plot class reaches a run through
    `job.figure`, not through a step function, so it is absent from the closure and an edit to
    it changes every rendered figure without moving the digest.
    """
    if not tree_clean or git_commit == "nogit":
        return None
    for run_dir in reversed(candidates):
        rec_identity, rec_commit, rec_tree_clean = _read_prov_reuse_fields(run_dir)
        if rec_identity == identity and rec_commit == git_commit and rec_tree_clean:
            return run_dir
    return None


def _locate_artifact(ref: ArtifactRef, context: ResolutionContext) -> LocatedArtifact:
    """The identity-keyed dir-glob locator (the strategy injected on the context).

    Finds the included sub-job's `{node_name}.pkl` - under --reuse-deps from a
    prior standalone or nested run whose identity AND commit match and whose tree
    is clean, otherwise by fresh-running the sub-job in figures-as-materialize mode
    nested under this composite's subjobs_output/ - loads it, and returns the value
    plus provenance facts. The strategy is injected on the context, so it can be
    swapped without changing ArtifactRef.
    """
    inc = ref.included
    node_name = ref.node_name
    # These are `Path | None` on ResolutionContext because a context can be built for
    # LocalRef resolution alone, where no run directory exists. On the ArtifactRef path all
    # four are populated - `run_job` sets them together at its single construction site.
    # Asserting that here rather than suppressing the type error turns a latent
    # `NoneType has no attribute glob` deep in the search into a named failure at the
    # boundary. `job_out_dir` is guarded too: it feeds the `not in d.parents` self-read
    # check, and `None` is a legal operand there, so an unset value would silently DISABLE
    # that guard and permit a wrong reuse rather than crash - the worst of the four.
    subjobs_dir = _require(context.subjobs_dir, "subjobs_dir")
    job_out_dir = _require(context.job_out_dir, "job_out_dir")
    pool_root = _require(context.pool_root, "pool_root")
    dataset_root = _require(context.dataset_root, "dataset_root")

    # Identity-keyed: the sub-job's content identity (memoized on the job) names its
    # output dir, so reuse keys on what the run IS, not merely its source.
    sub_identity = inc.job.build_identity(dataset_root).digest
    identity_short = sub_identity[:6]
    artifact_rel = f"{node_name}.pkl"
    dir_glob = f"{inc.job.name}_{identity_short}_*"

    def _cached_runs() -> list[Path]:
        # A reusable run is one whose dir holds {node_name}.pkl - a standalone run
        # that materialized it, or ANY composite's nested subjobs_output at any
        # depth. Search the whole pool recursively (`**/` matches every depth,
        # including direct children of the pool root); the dir basename ends with
        # the run timestamp, so sorting by name => newest last.
        #
        # `job_out_dir not in d.parents` excludes THIS run's own in-progress
        # subtree so a second include of the same sub-job can't self-read the first
        # include's fresh write. Under nesting `job_out_dir` is the *innermost*
        # run's dir, so this is a per-nested-run guarantee: a sibling branch's
        # already-FINISHED artifact (written earlier in the same top invocation,
        # under an ancestor composite) is a valid candidate - resolution is
        # sequential (no partial-write race) and the identity+commit+clean-tree
        # gate still applies, so that is content-correct diamond dedup, not a leak.
        candidates = pool_root.glob(f"**/{dir_glob}")
        return sorted(
            (
                d
                for d in candidates
                if (d / artifact_rel).exists() and job_out_dir not in d.parents
            ),
            key=lambda d: d.name,
        )

    cached = _cached_runs()
    eligible = (
        _reuse_eligible_dir(
            cached, sub_identity, context.git_commit, context.tree_clean
        )
        if context.reuse_deps
        else None
    )
    if eligible is not None:
        produced_dir = eligible
        mode = "cached"
    else:
        if cached and context.reuse_deps:
            print(
                f"[NOTE] '{inc.job.name}' has a cached '{node_name}' but it is not "
                f"reuse-eligible (identity/commit/clean-tree); re-running.",
                flush=True,
            )
        elif cached:
            print(
                f"[NOTE] composite is re-running sub-job '{inc.job.name}' which has a "
                f"cached '{node_name}'; pass --reuse-deps to reuse a cache that "
                f"matches the current identity + commit on a clean tree.",
                flush=True,
            )
        run_job(
            inc.job,
            subjobs_dir,
            force=True,
            data_root=context.data_root,
            render_figures=inc.figures,
            reuse_deps=context.reuse_deps,  # so a nested locator can reuse too
            pool_root=pool_root,  # every depth searches the same pool
        )
        produced = sorted(
            (d for d in subjobs_dir.glob(dir_glob) if (d / artifact_rel).exists()),
            key=lambda d: d.name,
        )
        if not produced:
            raise FileNotFoundError(
                f"sub-job '{inc.job.name}' did not persist '{artifact_rel}' under {subjobs_dir}. "
                f"Reference a figure-input node or an explicit materialize sink."
            )
        produced_dir = produced[-1]
        mode = "fresh"

    artifact_path = produced_dir / artifact_rel
    with artifact_path.open("rb") as handle:
        obj = load_artifact(handle)
    artifact_hash = f"sha256:{hash_file(artifact_path)}"
    # Record the sub-job's provenance dir relative to the top output/ pool root:
    # stable and copy-pasteable at every nesting depth (consistent with how depth-1
    # already renders it), never a machine-specific absolute path.
    prov_path = produced_dir / "provenance"
    try:
        prov_rel = str(prov_path.relative_to(pool_root))
    except ValueError:
        prov_rel = str(prov_path)
    return LocatedArtifact(
        value=obj, artifact_hash=artifact_hash, mode=mode, prov_dir_rel=prov_rel
    )


def _collect_includes_prov(
    job: Job, ordered_ids: list[str], context: ResolutionContext
) -> list[dict[str, object]]:
    """One prov entry per distinct ArtifactRef consumed, in first-seen toposort
    order (deterministic, independent of resolution timing)."""
    includes_prov: list[dict[str, object]] = []
    seen: set[int] = set()
    for node_id in ordered_ids:
        for input_ref in job.dag[node_id].inputs:
            if not isinstance(input_ref, ArtifactRef) or id(input_ref) in seen:
                continue
            seen.add(id(input_ref))
            located = context.artifacts[id(input_ref)]
            includes_prov.append(
                {
                    "alias": input_ref.included.alias,
                    "job_name": input_ref.included.job.name,
                    "node_name": input_ref.node_name,
                    "artifact_hash": located.artifact_hash,
                    "mode": located.mode,
                    "subjob_prov_dir": located.prov_dir_rel,
                }
            )
    return includes_prov


def _emit_prov(
    *,
    job: Job,
    root_node_id: str,
    ordered_ids: list[str],
    job_file: Path,
    repo: Path,
    dataset_root: Path,
    resolved_datasets: dict[str, Path],
    job_code_hash_full: str,
    git_commit: str,
    identity: str,
    tree_clean: bool,
    includes_prov: list[dict[str, object]],
    job_out_dir: Path,
    prov_name: str,
    targets: list[str],
    figure_node_label: str | None = None,
) -> None:
    """Build and save one sink's provenance record.

    Datasets reachable from `root_node_id` are looked up in `resolved_datasets`
    (the run's single resolution point - the same paths the loader opened, so the
    recorded hash describes the file actually loaded) and hashed via `content_hash`
    (memoized per process). Shared by the figure and materialize sinks so their two
    prov passes cannot drift apart. Only `prov_name`, `targets`, and
    `figure_node_label` differ between them. Hash failures propagate - a
    placeholder hash is never recorded.
    """
    ancestors = _ancestors(job, root_node_id)
    dataset_paths: list[Path] = []
    dataset_hashes: list[str] = []
    for node_id in job.dag:
        if node_id not in ancestors or node_id not in resolved_datasets:
            continue
        p = resolved_datasets[node_id]
        dataset_paths.append(p)
        dataset_hashes.append(f"sha256:{content_hash(p)}")

    try:
        job_file_rendered: Path | str = job_file.relative_to(repo)
    except ValueError:
        job_file_rendered = job_file

    record = build_prov_record(
        job_file=job_file_rendered,
        job_code_hash=f"sha256:{job_code_hash_full}",
        dataset_paths=[_fmt_dataset_path(p, dataset_root) for p in dataset_paths],
        dataset_hashes=dataset_hashes,
        git_commit=git_commit,
        identity=identity,
        tree_clean=tree_clean,
        pipeline_steps=[
            _format_step(job.dag[node_id])
            for node_id in ordered_ids
            if node_id in ancestors and job.dag[node_id].fn.__name__ != "_load_dataset"
        ],
        targets=targets,
        node_name=prov_name,
        figure_node_label=figure_node_label,
        includes=includes_prov,
    )
    save_prov(record, job_out_dir, prov_name)


def run_job(
    job: Job,
    out_dir: Path,
    force: bool = False,
    data_root: Path | None = None,
    render_figures: bool = True,
    reuse_deps: bool = False,
    pool_root: Path | None = None,
) -> None:
    # The reuse pool root: `out_dir` at the top level (== output/), threaded
    # unchanged into nested runs so every depth searches one shared pool.
    pool_root = out_dir if pool_root is None else pool_root
    job_file = _job_file(job)
    job.job_file = job_file  # ensure set for job_code_hash()/build_identity()
    repo = repo_root()
    # Routed through `resolve_data_root` rather than resolved here. `Path(data_root).resolve()`
    # accepts a root that does not exist, and `resolve_dataset_path` then falls back to the
    # repo root for any relative path - which every job in `jobs/` declares - so the run would
    # complete against the repository tree and record ITS dataset hashes. Passing the value in
    # is what makes a named-but-missing root a demand that fails. `None` behaves as before:
    # `resolve_data_root` falls through to `quebra.toml` and the user data directory.
    dataset_root = resolve_data_root(data_root)
    job_code_hash_full = job.job_code_hash()
    git_commit = get_git_commit()
    tree_clean = is_tree_clean()
    is_composite = bool(job.includes)

    _check_sink_artifact_names(job)  # fail-fast before any output dir exists
    _check_render_targets(job)

    # Resolve every dataset path once (fail-fast on a missing file, before any
    # output dir exists) - this mapping feeds BOTH the execution loop (loader) and
    # _emit_prov (hashing), so the recorded hash describes the file actually loaded.
    resolved_datasets: dict[str, Path] = {
        node_id: resolve_dataset_path(_dataset_of(node).path, dataset_root)
        for node_id, node in _dataset_load_nodes(job).items()
    }

    # Content identity (code + dataset content + sub-job identities): names the
    # output dir and keys reuse. build_identity is a function of static job
    # structure only, so a missing dataset raises here, before mkdir.
    identity = job.build_identity(dataset_root).digest
    identity_short = identity[:6]

    # Reuse gate (standalone jobs; composites always run fresh). Skip only when a
    # prior run has the SAME identity AND the SAME commit AND the tree is clean -
    # any mismatch (edited shared code → commit differs; dirty tree) re-runs fresh.
    #
    # A candidate must also be COMPLETE: every sink's artifact present, not just a matching
    # directory name. The composite path filters the same way, in `_cached_runs`. Without it a
    # run that died mid-sink-loop stays eligible, because its first sink's provenance record
    # is readable and matches.
    if not is_composite and not force:
        expected_pkls = _expected_sink_pkls(job)
        candidates = sorted(
            (
                run_dir
                for run_dir in out_dir.glob(f"{job.name}_{identity_short}_*")
                if all((run_dir / name).exists() for name in expected_pkls)
            ),
            key=lambda d: d.name,
        )
        if _reuse_eligible_dir(candidates, identity, git_commit, tree_clean):
            print(
                f"Skipping {job.name}: reusable output exists "
                f"(identity + commit match, clean tree)"
            )
            return

    # Create timestamped folder: always unique when creating or force-rerunning
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    job_out_dir = out_dir / f"{job.name}_{identity_short}_{timestamp}"
    job_out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Running {job.name} -> {job_out_dir.relative_to(out_dir)}")

    sink_ids = [
        sink.input.node_id if isinstance(sink, _FigureSink) else sink.node.node_id
        for sink in job.sinks
    ]
    ordered_ids = _toposort(job, sink_ids)
    results: dict[str, object] = {}

    # One context resolves every input the same way: LocalRefs read `results`,
    # ArtifactRefs (composite → included sub-job) run/locate via the injected
    # locator. A referenced sub-job runs lazily, at most once per ref (memoized).
    context = ResolutionContext(
        results=results,
        locate=_locate_artifact,
        subjobs_dir=job_out_dir / "subjobs_output",
        job_out_dir=job_out_dir,
        pool_root=pool_root,
        reuse_deps=reuse_deps,
        data_root=data_root,
        dataset_root=dataset_root,
        git_commit=git_commit,
        tree_clean=tree_clean,
    )

    for node_id in ordered_ids:
        node = job.dag[node_id]
        inputs = [input_ref.resolve(context) for input_ref in node.inputs]
        kwargs = node.kwargs
        if node_id in resolved_datasets:
            # Hand the loader a COPY of the Dataset carrying the resolved absolute
            # path; node.kwargs stays pristine so relative paths (not machine-
            # specific absolutes) render into pipeline_steps/Mermaid labels.
            kwargs = {
                **kwargs,
                "dataset": dataclasses.replace(
                    _dataset_of(node), path=resolved_datasets[node_id]
                ),
            }
        results[node_id] = node.fn(*inputs, **kwargs)

    includes_prov = _collect_includes_prov(job, ordered_ids, context)

    for sink in job.sinks:
        if isinstance(sink, _FigureSink):
            result = results[sink.input.node_id]
            if render_figures:
                plot = sink.plot_class(name=sink.name)
                for target_name in sink.targets:
                    target = RENDER_TARGETS[target_name]
                    target(plot, result, job_out_dir)
                # Persist the figure's input too, so a standalone job's figure
                # output is reusable by a composite (not only figures-as-materialize runs).
                with (job_out_dir / f"{sink.input.node_id}.pkl").open("wb") as handle:
                    pickle.dump(result, handle)
                prov_name = sink.name
                rendered_targets = sink.targets
                fig_label = f"{sink.plot_class.__name__}\\ngit:{git_commit}"
            else:
                # Figures-as-materialize: persist the figure's input under its node
                # id so a composite can ref it by that name; emit no PDF.
                prov_name = sink.input.node_id
                with (job_out_dir / f"{prov_name}.pkl").open("wb") as handle:
                    pickle.dump(result, handle)
                rendered_targets = []
                fig_label = None
            _emit_prov(
                job=job,
                root_node_id=sink.input.node_id,
                ordered_ids=ordered_ids,
                job_file=job_file,
                repo=repo,
                dataset_root=dataset_root,
                resolved_datasets=resolved_datasets,
                job_code_hash_full=job_code_hash_full,
                git_commit=git_commit,
                identity=identity,
                tree_clean=tree_clean,
                includes_prov=includes_prov,
                job_out_dir=job_out_dir,
                prov_name=prov_name,
                targets=rendered_targets,
                figure_node_label=fig_label,
            )
            continue

        result = results[sink.node.node_id]
        with (job_out_dir / f"{sink.name}.pkl").open("wb") as handle:
            pickle.dump(result, handle)
        _emit_prov(
            job=job,
            root_node_id=sink.node.node_id,
            ordered_ids=ordered_ids,
            job_file=job_file,
            repo=repo,
            dataset_root=dataset_root,
            resolved_datasets=resolved_datasets,
            job_code_hash_full=job_code_hash_full,
            git_commit=git_commit,
            identity=identity,
            tree_clean=tree_clean,
            includes_prov=includes_prov,
            job_out_dir=job_out_dir,
            prov_name=sink.name,
            targets=[],
        )
