from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(slots=True)
class Dataset:
    """A file to read, plus what is needed to interpret it.

    `extra` and `loader_kwargs` are BOTH free-form dicts and it matters which one a key
    goes in. `extra` is run metadata: it survives into `Norm["meta"]` and reaches
    provenance. `loader_kwargs` is passed to the reader for this file's extension and is
    never read back.

    Merging them is a defect rather than a simplification: one field would be forwarded to
    the loader AND read back as metadata, so `extra={'run_start_unix_s': ...}` - the value
    `_load_dataset`'s own error message tells you to set - would reach
    `pd.read_csv(path, **extra)` and raised `unexpected keyword argument`. It only ever
    worked because every in-repo caller used `.pickle`, whose loader discards the dict.
    """

    path: Path | str
    schema: Any = None
    qubit: int | None = None
    device: str | None = None
    duration_h: float | None = None
    extra: dict[str, Any] = field(default_factory=dict)
    loader_kwargs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.path = Path(self.path)
        # keep path normalized; a side table is a separate Dataset
