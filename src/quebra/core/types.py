from __future__ import annotations

from collections.abc import Iterator, MutableMapping
from dataclasses import dataclass, field
from typing import Any


# Type aliases - semantic labels for array roles in the pipeline
TimeSeries = Any  # 1-D np.ndarray indexed by t_rel_s
Parameter = (
    Any  # 1-D np.ndarray broadcast-constant (scalar value repeated to match t_rel_s)
)

# Rows of the file that the load drops before any step runs, recorded in `norm["meta"]`
# and copied onto the results that describe the record: by a schema's `validate`, and by
# the default normaliser's finite timestamp and frequency mask.
LOAD_COUNT_KEYS = ("n_rows_dropped_by_schema", "n_rows_dropped_invalid")


@dataclass
class Measurement:
    """A single scalar measurement with its timestamp."""

    t_unix_s: float
    value: float


@dataclass
class CalibrationEvent:
    """One calibration event from a calibration log dataset."""

    t_event_unix_s: float
    payload: dict = field(default_factory=dict)
    pre: Measurement | None = None  # reserved for future companion work
    post: Measurement | None = None  # reserved for future companion work


class Norm(MutableMapping):
    """Normalized dataset mapping flowing between pipeline steps.

    Behaves exactly like a plain dict (MutableMapping). All canonical key names
    use unit suffixes (_hz, _rel_s, _unix_s). Use dict(norm) to get a plain dict.

    The ``events`` attribute holds named lists of typed events (e.g. CalibrationEvent)
    and is not accessible via the MutableMapping interface.
    """

    __slots__ = ("_data", "events")

    def __init__(
        self,
        data: dict[str, Any] | None = None,
        events: dict[str, list] | None = None,
    ) -> None:
        self._data: dict[str, Any] = dict(data) if data is not None else {}
        self.events: dict[str, list] = events if events is not None else {}

    def __getitem__(self, key: str) -> Any:
        return self._data[key]

    def __setitem__(self, key: str, value: Any) -> None:
        self._data[key] = value

    def __delitem__(self, key: str) -> None:
        del self._data[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._data)

    def __len__(self) -> int:
        return len(self._data)

    def __repr__(self) -> str:
        return f"Norm({list(self._data.keys())})"

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> Norm:
        return cls(d)
