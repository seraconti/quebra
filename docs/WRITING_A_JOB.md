# Writing a job

A job is a Python file that declares what to load, what to compute, and what to draw.
Nothing in it runs when the file is imported. The file builds a graph; `quebra run` walks
it. That separation is what lets `quebra inspect` show you the whole computation, with the
arguments that affect every result, without touching a dataset.

You do not need to work inside this repository to write one. Install QUEBRA from a clone
(`pip install .`; there is no PyPI release yet), put a `.py` file anywhere, and run it by
path.

## The smallest job that does something

Nothing here is ours: your own data file, your own schema, four lines of computation.

Say you have a `my_data.csv` with a `value` column. Put this job file beside it:

```python
"""One file in, one number out."""

from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.core.types import Norm
from quebra.schemas.null import NullSchema

job = Job("my_first_job")

readings = Dataset(path="my_data.csv", schema=NullSchema)
norm = job.load(readings)


def mean_value(n: Norm) -> float:
    return float(n["value"].mean())


average = job.step(mean_value, norm, name="mean_value")
job.materialize(average, name="mean_value")
```

Run it:

```bash
quebra inspect my_first_job.py            # see the graph
quebra run my_first_job.py --data-root .  # then compute it
```

`inspect` prints the two nodes and their arguments and never touches a file, so it needs no
data root. `run` does: `--data-root .` says the dataset paths are relative to the current
directory. `run` writes `output/my_first_job_<hash>_<timestamp>/mean_value.pkl` holding the
mean, with a provenance record beside it.

**CSV is only what this example happens to use.** Files are read by extension, and `.csv`,
`.pickle` / `.pkl`, `.yaml` / `.yml` and `.h5` / `.hdf5` ship with the tool — `quebra inspect`
prints the list under `Loaders:`. Any other format is fine, it just needs a loader, which is
a function returning a DataFrame and a decorator, and can live in the job file itself:

```python
from quebra.loaders.registry import register_loader


@register_loader(".tsv")
def _load_tsv(path, meta):
    return pd.read_csv(path, sep="\t", **meta)
```

If a shipped reader just needs options, skip the loader and pass them as
`Dataset(loader_kwargs=...)` — a semicolon-separated file is `loader_kwargs={"sep": ";"}`.
Keep run metadata out of there and in `extra`, which is a separate field for that purpose.

`NullSchema` copies every column into the `Norm` under its own name and asks nothing of your
data. It is the starting point, not the destination: a schema is what maps your columns onto
the names the analyzers read, and writing one is the normal way to point this tool at your
own records. That is [WRITING_A_SCHEMA.md](WRITING_A_SCHEMA.md).

`job.load_df(dataset)` skips the schema entirely and hands your step a plain
`pandas.DataFrame`. It suits a side table — the bench jobs read a CSV of results that
way — but a step taking a raw frame has to know that file's column names, so it does not
compose the way a `Norm` does.

## The same job against a real record

The example above supplies its own schema. The other way is to use one that ships with the
tool, which is what our own jobs do — `track912Schema` validates this group's Ramsey column
conventions, and the default normaliser then produces the `t_rel_s` / `delta_hz` arrays the
analyzers read:

```python
"""One dataset, one T2* fit."""

from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.schemas.track912 import track912Schema
import quebra.analyzers.t2star as t2star

job = Job("my_first_job")

reads = Dataset(
    path="data/real_private/6D2S/070423_6D2S_qubit1.pickle",
    schema=track912Schema,
    qubit=1,
    device="6D2S",
)

norm = job.load(reads)


def fit_t2star(n):
    return t2star.run(t2star.make_inputs_from_norm(n))


result = job.step(fit_t2star, norm, name="t2star")
job.materialize(result, name="t2star_result")
```

That path and that schema are ours, and the pickle is not distributed. To run the analyzers
on your own records, write a schema that maps your columns onto what they read — that is
[WRITING_A_SCHEMA.md](WRITING_A_SCHEMA.md), and it is a much smaller job than it looks.

## The four things a job file can say

**`job.load(dataset)`** reads a file through the loader registry and hands the frame to the
dataset's schema, which returns a `Norm` — the mapping of named arrays every analyzer takes.
With no schema you get the default normaliser, which expects a `timestamp`/`frequency`
series. The path is relative to the data root, not to your file. See "Where data comes from"
below.

**`job.load_df(dataset)`** does the same read but skips the schema and yields the raw
`pandas.DataFrame`.

**`job.step(fn, *inputs, name=..., **kwargs)`** adds a node. `fn` must be a pure function of
its inputs: the graph caches and re-runs steps by content hash, so a step that reads a clock,
a file, or a global cannot be cached correctly. Keyword arguments are part of the node's
identity and appear on the provenance graph, which is why anything that changes the answer
belongs in a kwarg rather than a closure.

**`job.materialize(node, name=...)`** writes a node's result to the run directory as a
pickle, with a provenance record beside it.

**`job.figure(PlotClass, node, targets=[...], title=...)`** renders a node. The plot class
draws and nothing else; every number it shows must already be a field of the artifact it is
given.

## Where data comes from

A `Dataset.path` is relative to the **data root**. Two of the four mechanisms are *demands*
and two are *candidates*.

A demand is a root you named for this run. If it does not exist, the run fails; it never
falls through to another mechanism:

1. `--data-root` on the command line
2. the `QUEBRA_DATA_ROOT` environment variable

A candidate is where to look when you named nothing. First existing hit wins:

3. `[tool.quebra] data_root` in a `quebra.toml`, in the working directory or any directory
   above it
4. a `platformdirs` user data directory

The distinction matters because a relative `Dataset.path` that misses the data root is then
tried against the **project root** (see below), and most in-repo paths exist there. So a
typo'd `--data-root` that fell through would not fail - it would quietly analyse the project
tree and record *its* dataset hashes as if they were the ones you asked for.

If none resolves, QUEBRA raises `DataRootNotFound` listing every location it tried. It never
guesses a relative path, because a guess would silently point at the wrong tree and the first
symptom would be a dataset hash that disagrees with a published record.

For a project of your own, the simplest answer is a `quebra.toml` beside your job files:

```toml
[tool.quebra]
data_root = "../data"
```

A relative `data_root` anchors on the file that declares it, not on where you happen to be
standing.

## What a run produces

`quebra run` writes to `./output/<job>_<identity>_<timestamp>/`, containing the figures, any
materialised artifacts, and `provenance/` with a record in both JSON and Markdown. Use
`--output-root` to write elsewhere.

The `<identity>` is a content hash of the job file's source, **the library code its steps
reach**, the resolved arguments those steps were called with, the datasets it loaded, and the
identities of any sub-jobs. Two runs with the same identity computed the same thing.

The library half matters: editing an analyzer your job calls moves the identity, because it
changed what produced the numbers. Editing a module your job does not reach does not — the
hash covers the static import closure of your steps, not the whole package.

### Publishing a figure

`output/` is gitignored — the artifacts run to hundreds of megabytes. The provenance beside
them is kilobytes and carries the identity, the dataset digests, the git commit and the
pipeline steps, which is everything needed to say what produced a figure without shipping any
of it. When a figure appears in a paper, commit that:

```bash
make promote RUN=output/t2star_q1_070423_<identity>_<timestamp> NOTE="thesis ch4 fig 3"
```

It writes `published/<job>_<identity>/` with the provenance records and a `PROMOTED.toml`
saying where the figure appears. **Only the provenance is copied** — never an artifact, or
the rule keeping `output/` out of the repository would mean nothing.

A run made against an uncommitted working tree is refused, because no commit reproduces it.
Commit first and re-run, or pass `PROMOTE_FLAGS=--allow-dirty` to promote it anyway with
`tree_clean = false` recorded in the file.

The `dataset_hash` in a promoted record is the same digest
`data/real_private/MANIFEST.toml` lists, so a reader can check which records a published
figure rests on without receiving any of them.

## Rules that are not style

**A step is a pure function.** No disk access, no `matplotlib`, no clock, no global state.
I/O belongs in loaders; drawing belongs in plot classes.

**Physical quantities carry unit suffixes**: `_hz`, `_rel_s`, `_unix_s`, `_utc_dt`. Unit and
clock mismatch is the expensive bug class in this domain, and a suffix is what makes it
visible at a call site.

**Errors are raised, not swallowed.** A silent fallback produces a wrong-but-plausible
number, which in a reproducibility tool is worse than a crash.

**Results are typed dataclasses, not dicts.** Define the type in the same file as the step
that returns it.

## Composing jobs

A job can include another and reference its results, naming it by its `JOB_ID`:

```python
sub = job.include("other_job", alias="upstream")
node = job.step(combine, sub.ref("some_node"), name="combined")
```

The included job runs through the normal runner, so its datasets and provenance resolve
exactly as they would standalone. Its identity is consumed by the parent, so the parent's
identity changes when the child's does.

Composites declare `JOB_SWEEP = False`, so a bare `quebra run --all` skips them — sweeping a
composite re-runs every sub-job. They live in `jobs/composite/` for readability, but the
directory decides nothing; the declaration does.

## What a job file declares

Three module-level constants, read statically — discovery never imports your job, because
importing one builds its graph:

```python
JOB_ID = "my_first_job"     # the logical name `job.include` resolves. Move the file freely.
JOB_FAMILY = "t2star"       # the subject. `quebra run --all --family t2star` selects it.
JOB_SWEEP = True            # default. False keeps it out of a bare `run --all`.
```

`JOB_ID` is what makes a job's location stop being part of its identity: a composite naming
`"t2star_q1_070423"` keeps working when that file moves. A job with no `JOB_ID` is not
discovered, so `run --all` will not find it.

`docs/JOBS.md` is the generated table of every declared job and the parameters that
distinguish it — regenerate with `python scripts/make_job_manifest.py`.
