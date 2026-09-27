"""Unpickle staleness guard for materialized dataclass artifacts.

A pickle written before a field was added to its dataclass lacks that field: unpickling restores ``__dict__`` directly, bypassing ``__init__`` defaults,
so a composite reusing such an artifact would crash mid-render (AttributeError on
a factory field) or silently draw a class-level default. Loading one must instead
fail loudly at the pickle boundary ("errors are raised, not swallowed"). This guard
covers the unpickle path; construction-time completeness is enforced separately by
each panel dataclass's __post_init__ (unpickling bypasses __post_init__, hence both).
"""

from __future__ import annotations

import dataclasses
from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:
    from _typeshed import DataclassInstance


class StaleArtifactGuard:
    """Mixin for non-slots dataclasses: validate state on unpickle.

    The required key set derives from ``dataclasses.fields()`` - never a
    hand-maintained list - so it tracks field additions automatically. A valid
    builder-produced pickle always carries every field in ``__dict__``
    (default_factory fields included); a pickle predating a field does not.
    """

    def __setstate__(self, state: dict[str, object]) -> None:
        cls = type(self)
        if not isinstance(state, dict):
            # e.g. the (dict, slots_dict) tuple a slots=True dataclass would emit -
            # fail with the guard's error, not an AttributeError on .keys().
            raise ValueError(
                f"stale {cls.__name__} artifact: unexpected pickle state of type "
                f"{type(state).__name__} - re-run the sub-job (run the composite "
                "without --reuse-deps)."
            )
        # `cast`: every concrete subclass IS a dataclass, but the mixin itself is not,
        # so mypy cannot see `cls` as one. Narrowing at runtime would be a lie - a
        # non-dataclass subclass is a programming error, not a case to handle.
        missing = {
            f.name for f in dataclasses.fields(cast("type[DataclassInstance]", cls))
        } - state.keys()
        if missing:
            raise ValueError(
                f"stale {cls.__name__} artifact: pickle lacks field(s) "
                f"{sorted(missing)} - it predates a field added to this artifact. "
                "Re-run the sub-job (run the composite without --reuse-deps)."
            )
        self.__dict__.update(state)
