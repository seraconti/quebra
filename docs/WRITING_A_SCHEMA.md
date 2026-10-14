# Writing a schema

A schema is the one piece of QUEBRA that is about *your* data. Everything else — the graph,
the caching, the provenance, the analyzers — works on a `Norm`, which is a mapping of named
numpy arrays plus a `meta` dict. A schema turns your file into one.

There is no registration step and no base class to inherit. `job.load` looks for a method by
name, so any class with the right method is a schema.

## Start with none

`NullSchema` copies every column into the `Norm` under its own name:

```python
from quebra.core.dataset import Dataset
from quebra.schemas.null import NullSchema

Dataset(path="my_data.csv", schema=NullSchema)
```

That is enough to run a job. You write a real schema when you want the next two things:
columns named what the analyzers expect, and the same checks applied to every file in a
family rather than remembered each time.

## The two shapes

`job.load` checks for these in order and stops at the first it finds.

**`to_norm(frame, dataset) -> Norm`** — you own the conversion outright:

```python
import pandas as pd

from quebra.core.dataset import Dataset
from quebra.core.types import Norm


class ReadingsSchema:
    @classmethod
    def to_norm(cls, frame: pd.DataFrame, dataset: Dataset) -> Norm:
        return Norm(
            {
                "t_rel_s": frame["elapsed"].to_numpy(dtype=float),
                "value": frame["reading"].to_numpy(dtype=float),
                "meta": {"dataset_id": dataset.path.stem},
            }
        )
```

Use this for anything that is not already a timestamp/frequency series — event logs, count
data, a table with your own column names.

**`validate(frame, dataset=...) -> DataFrame`** — you check and clean a frame and hand it
back, and `RamseySeriesSchema` runs on the result to produce the `Norm`. Use this when your
data *is* a timestamp/frequency series and you only need to rename or filter your way into
it. It is the `pandera.DataFrameModel` shape.

An object with neither method is rejected by name at load time.

## A worked example: `RamseySeriesSchema`

`src/quebra/schemas/ramsey_series.py` is the longest schema in the repo and a good model for
a `to_norm` of your own, because its problem is the ordinary one: a file of readings that has
to become time and signal arrays with enough metadata to be reproducible.

It expects a `timestamp` column and a `frequency` column, both numeric, and a resolvable run
start — either `Dataset.extra["run_start_unix_s"]` or a `DDMMYY_` filename prefix, which it
reads as local midnight and warns about. From those it produces `t_rel_s`, `delta_hz`
(frequency relative to the mean), `raw_frequency_hz` and `meta`, and passes through
`normalised chi-square`, `T2star`, `T2star error` and `frequency error` (as
`frequency_error_hz`) when they are present.

Three things in it are worth copying:

- **It refuses rather than guesses.** No run start means an error naming the two ways to
  supply one, not a fallback that would silently date the run to 1970.
- **Optional columns are optional.** Present, they are carried; absent, nothing breaks.
- **Everything the run depends on lands in `meta`.** That is what makes the artifact
  reproducible later.

It is also the schema `job.load` falls through to once a `validate` schema has run, which is
the route every Ramsey job in `jobs/active/` actually takes: they name `track912Schema`, it
renames this group's columns onto `timestamp`/`frequency`, and `RamseySeriesSchema` turns the
result into a `Norm`. It is the default for `Dataset.schema=None` too, but no shipped job
relies on that.

## What `Dataset` carries

```python
Dataset(
    path,                # the file, relative to the data root
    schema=None,         # your schema class
    qubit=None,          # optional; added to the frame as a qubit_id column
    device=None,         # optional; added to the frame as a device column
    duration_h=None,     # optional; run length - carried into meta by the shipped
                         #   schemas, though NullSchema does not copy it
    extra={},            # run metadata: reaches your schema and provenance
    loader_kwargs={},    # options for the file reader, e.g. {"sep": ";"}
)
```

`qubit` and `device` become columns only if you set them, so a schema that has no use for
them can ignore them entirely.

> **PS — `extra` and `loader_kwargs` are not interchangeable.** `extra` is metadata: it
> reaches your schema as `dataset.extra` and is conventionally copied into `Norm["meta"]`,
> so it ends up in provenance. `loader_kwargs` goes to the reader for that file extension
> (`pd.read_csv(path, **loader_kwargs)`) and is never read back. A run name is `extra`; a
> separator or an HDF5 key is `loader_kwargs`.

## What `Norm` has to contain

`Norm` is a `MutableMapping` (`src/quebra/core/types.py`). Nothing enforces a key set, so the
contract is with whatever consumes it:

- **your own steps** — any keys you like;
- **the shipped analyzers** — `t2star.make_inputs_from_norm` and friends expect
  `RamseySeriesSchema`'s keys, so a custom schema feeding a shipped analyzer must produce
  them;
- **`meta`** — always include it, with at least a `dataset_id`. Provenance and figure titles
  read it, and an absent one surfaces late and unhelpfully.

## Other examples in this repo

- `src/quebra/schemas/null.py` — the shortest possible `to_norm`.
- `src/quebra/schemas/track912.py` — a `validate` schema, `pandera` model, inherits
  `BaseQubitSchema`.
- `src/quebra/schemas/calibration_log.py` — a `to_norm` turning a calibration log into
  `CalibrationEvent` objects. The closest template for data that is not a time series.

## Testing yours

`tests/test_load_dataset_contract.py` is the pattern: write a small file under `tmp_path`,
load it through `_load_dataset`, assert on the `Norm`. It needs no private data and runs in
about a second.
