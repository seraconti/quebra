from __future__ import annotations

import hashlib
import json
import re
import subprocess
from pathlib import Path


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def hash_string(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


# Both git helpers below anchor on the working directory and share this timeout. They answer
# the commit and clean-tree conditions of the reuse rule - `_reuse_eligible_dir` requires a
# commit match AND a clean tree - so they must describe the SAME repository or those two
# conditions mean nothing. `Path.cwd()` is what makes that true: this module ships inside the
# wheel, so `__file__` is site-packages once installed, where one helper would read a
# repository the other never saw.
#
# The timeout bounds the damage rather than fixing a known hang: both commands are local and
# should answer in milliseconds, so a wait means something is wrong (a stale lock, a stalled
# filesystem) and a provenance helper is the wrong place to block a run indefinitely. A
# timeout reads as "could not answer", which the callers already handle.
_GIT_TIMEOUT_S = 30.0


def _git(*args: str) -> str | None:
    """stdout of a git command run in the cwd, or None if git could not answer.

    Every failure mode collapses to None on purpose, and callers must treat None as "do not
    know" rather than as an answer. `OSError` covers a missing or non-executable git and a
    deleted cwd; `SubprocessError` covers a non-zero exit and the timeout.
    """
    try:
        completed = subprocess.run(
            ["git", *args],
            cwd=Path.cwd(),
            check=True,
            capture_output=True,
            text=True,
            timeout=_GIT_TIMEOUT_S,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return completed.stdout


def get_git_commit() -> str:
    """The commit of the PROJECT being run, or "nogit"."""
    out = _git("rev-parse", "--short", "HEAD")
    return (out.strip() or "nogit") if out is not None else "nogit"


def is_tree_clean() -> bool:
    """True iff there are no uncommitted *tracked* changes.

    `--untracked-files=no` so a stray untracked scratch file does not disable artifact
    reuse; the claim is only "no uncommitted tracked changes". "Do not know" reports as
    dirty, so reuse never fires on a tree this cannot vouch for.
    """
    out = _git("status", "--porcelain", "--untracked-files=no")
    return out is not None and out.strip() == ""


def build_prov_record(
    job_file: str | Path,
    job_code_hash: str,
    dataset_paths: list[str],
    dataset_hashes: list[str],
    git_commit: str,
    pipeline_steps: list[str],
    targets: list[str],
    node_name: str,
    identity: str | None = None,
    tree_clean: bool = True,
    figure_node_label: str | None = None,
    includes: list[dict[str, object]] | None = None,
) -> dict[str, object]:
    # keep backwards-compatible top-level fields for the primary dataset
    primary_path = dataset_paths[0] if dataset_paths else ""
    primary_hash = dataset_hashes[0] if dataset_hashes else ""
    return {
        "node_name": node_name,
        "job_file": str(job_file),
        # NOT a digest of `job_file`: it folds that file's text, the digests of every
        # module the steps can reach, and the arguments each step was called with. The
        # name says so, because `sha256sum` on `job_file` will not reproduce it.
        "job_code_hash": job_code_hash,
        "dataset_path": str(primary_path),
        "dataset_hash": primary_hash,
        "dataset_paths": list(dataset_paths),
        "dataset_hashes": list(dataset_hashes),
        # content identity of the run (code + data + sub-job identities); folds
        # nothing time- or machine-specific. git_commit below is separate lineage.
        "identity": identity,
        "git_commit": git_commit,
        # whether the working tree was clean when this artifact was produced. An
        # artifact built on a dirty tree is never reuse-eligible: it does not
        # correspond to any commit's code (see _reuse_eligible_dir).
        "tree_clean": tree_clean,
        "pipeline_steps": list(pipeline_steps),
        "targets_rendered": list(targets),
        "figure_node_label": figure_node_label,
        # composite jobs: opaque upstream sub-jobs (see save_prov mermaid).
        "includes": list(includes) if includes else [],
    }


# Semantics colour for analyst-declared inputs, matching the paper figures. Kept here
# rather than imported from plots.theme: theme is the render layer, provenance is not.
_DECLARED_INPUT_FILL = "#FDF3D2"
_DECLARED_INPUT_STROKE = "#B8912A"


def _short_hash(value: str, length: int = 6) -> str:
    stripped = value.removeprefix("sha256:")
    return stripped[:length]


def _mermaid_label(text: str) -> str:
    # Mermaid has no backslash escape; a literal quote inside a quoted label is written
    # as the entity #quot;. Labels come from repr() of step kwargs, which switches to
    # double quotes as soon as a value contains an apostrophe - so this is reachable.
    escaped = text.replace("\n", "\\n").replace('"', "#quot;")
    return f'"{escaped}"'


def save_prov(record: dict[str, object], out_dir: Path, node_name: str) -> None:
    # All provenance lives under a provenance/ subdir to keep the job output
    # dir uncluttered (figures and .pkl artifacts stay at the top level).
    prov_dir = out_dir / "provenance"
    prov_dir.mkdir(parents=True, exist_ok=True)

    json_path = prov_dir / f"{node_name}.prov.json"
    md_path = prov_dir / f"{node_name}.prov.md"

    json_path.write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    md_path.write_text(_mermaid_graph(record, node_name), encoding="utf-8")


def _mermaid_id_factory():
    """Allocate canonical, collision-safe mermaid node ids from human names
    (replacing the old throwaway A/B/C positional letters)."""
    used: dict[str, int] = {}

    def make(base: str) -> str:
        slug = re.sub(r"[^0-9A-Za-z]+", "_", base).strip("_").lower() or "n"
        used[slug] = used.get(slug, 0) + 1
        return slug if used[slug] == 1 else f"{slug}_{used[slug]}"

    return make


def _mermaid_graph(record: dict[str, object], node_name: str) -> str:
    """One DAG renderer for every job (standalone or composite).

    Source nodes - loaded datasets and/or references to included sub-jobs - feed
    the pipeline steps, which feed the sink node. Node ids are canonical slugs of
    the real names (not positional A/B/C letters); labels are human, with machine
    hashes shown only as short tags. An included reference is a first-class node
    that click-throughs to the sub-job's own provenance graph.

    Node SHAPE carries meaning (the grammar the paper figures already use):
    `[/data/]`, `[process]`, `{decision}`, `[[analyst-declared input]]`,
    `[(sub-job reference)]`, and `-.->` for an optional edge. Only the categories
    that map to real record fields are emitted; decision and declared-input shapes
    are reserved, and appear in the legend so a reader can decode the whole grammar
    from one file. Step kwargs stay on the step's own label: they are part of that
    step's identity, not a separate input, and no node for them exists in job.dag.
    """
    dataset_paths = record.get("dataset_paths") or (
        [record.get("dataset_path")] if record.get("dataset_path") else []
    )
    dataset_hashes = record.get("dataset_hashes") or (
        [record.get("dataset_hash")] if record.get("dataset_hash") else []
    )
    includes = record.get("includes") or []
    steps = [str(step) for step in record.get("pipeline_steps", [])]
    git_commit = str(record["git_commit"])  # always set by build_prov_record
    figure_node_label = record.get("figure_node_label")

    nid = _mermaid_id_factory()
    lines = ["graph LR"]
    click_lines: list[str] = []
    source_ids: list[str] = []

    for p, h in zip(dataset_paths, dataset_hashes, strict=True):
        name = Path(str(p)).name
        node_id = nid(f"ds_{name}")
        # Built in two steps on purpose. Nesting this f-string inside the outer one put a
        # backslash escape inside an f-string EXPRESSION, which is a SyntaxError on Python
        # 3.11 and only legal from 3.12. pyproject declares requires-python >= 3.11, so on
        # the oldest supported interpreter the module would not even import. It went unseen
        # because this machine runs 3.13; `ruff --target-version py311` is what caught it.
        label = f"{name}\n{_short_hash(str(h))}"
        lines.append(f"  {node_id}[/{_mermaid_label(label)}/]")
        source_ids.append(node_id)

    for inc in includes:
        alias = str(inc.get("alias", "sub"))
        job_name = str(inc.get("job_name", ""))
        ref_node = str(inc.get("node_name", ""))
        short = _short_hash(str(inc.get("artifact_hash", "")))
        node_id = nid(f"ref_{alias}_{ref_node}")
        label = f"{alias}: {job_name}:{ref_node}\\n{short}"
        lines.append(f"  {node_id}[({_mermaid_label(label)})]")
        source_ids.append(node_id)
        prov_ref = inc.get("subjob_prov_dir")
        if prov_ref:
            # runner-generated relative path; strip any quote defensively so the
            # click line can't be broken by the URL
            safe_ref = str(prov_ref).replace('"', "")
            click_lines.append(f'  click {node_id} "{safe_ref}"')

    if not source_ids:
        # Absent data is still data: same shape, so the grammar has no exception.
        node_id = nid("source")
        lines.append(f"  {node_id}[/{_mermaid_label('(no inputs)')}/]")
        source_ids.append(node_id)

    previous_ids = source_ids
    for step in steps:
        fn_name = step.split("(", 1)[0]
        node_id = nid(f"step_{fn_name}")
        lines.append(f"  {node_id}[{_mermaid_label(step)}]")
        for src in previous_ids:
            lines.append(f"  {src} --> {node_id}")
        previous_ids = [node_id]

    if isinstance(figure_node_label, str) and figure_node_label.strip():
        final_label = figure_node_label
    else:
        final_label = f"{node_name}\\ngit:{git_commit}"
    sink_id = nid(f"sink_{node_name}")
    lines.append(f"  {sink_id}[{_mermaid_label(final_label)}]")
    for src in previous_ids:
        lines.append(f"  {src} --> {sink_id}")

    lines.extend(click_lines)
    lines.extend(_mermaid_legend())
    return "\n".join(lines) + "\n"


def _mermaid_legend() -> list[str]:
    """The grammar key, appended to every graph so one file decodes itself.

    Every legend node declared at the start of a line also appears on the left of a
    solid edge: the graph must stay fully wired, and tests/test_provenance_graph.py
    asserts exactly that over all declared nodes. Mermaid places a disconnected
    subgraph wherever its layout engine likes, so this is "at the end of the file",
    not "at the bottom of the picture".
    """
    return [
        f"  classDef declared fill:{_DECLARED_INPUT_FILL},stroke:{_DECLARED_INPUT_STROKE}",
        '  subgraph legend["grammar"]',
        "    direction LR",
        '    lg_data[/"data"/] --> lg_step["process"]',
        '    lg_step --> lg_decision{"decision"}',
        '    lg_declared[["analyst-declared input"]] --> lg_step',
        '    lg_ref[("sub-job reference")] --> lg_step',
        '    lg_decision -.-> lg_optional["optional edge"]',
        "  end",
        "  class lg_declared declared",
    ]
