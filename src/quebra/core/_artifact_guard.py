"""Unpickle compatibility for materialized dataclass artifacts: moved modules, and stale ones.

A pickle written before a field was added to its dataclass lacks that field: unpickling restores
the instance state directly, bypassing ``__init__`` defaults, so a consumer of such an artifact
would crash mid-render (AttributeError on a missing slot or factory field) or silently read a
class-level default. Loading one must instead fail loudly at the pickle boundary ("errors are
raised, not swallowed"). Two mechanisms cover it. `StaleArtifactGuard` is a mixin that checks its
own class in ``__setstate__``. `load_artifact` checks every dataclass reachable from the loaded
object through dicts, lists, tuples, sets and dataclass fields, whatever its base class; it does
not look inside arrays or frames, and it cannot judge an `init=False` field. Construction-time completeness is enforced separately by each
panel dataclass's __post_init__ (unpickling bypasses __post_init__, hence both).

A pickle also records the fully qualified module path of every class it holds, so moving a
module leaves every artifact written before the move unreadable, even when the class itself is
unchanged. `load_artifact` resolves those old paths through `MODULE_ALIASES`.

WHAT AN ALIAS MAY AND MAY NOT DO. It maps a MOVED package to its new home, and only when the new
home defines a class of the same name. It never renames a class and never maps one schema onto
another: a class that was renamed or restructured is a different object, and pretending
otherwise would load old field values into a new shape. Those stay unreadable, with an error
that says so. A class that was moved AND gained a field loads, and is then refused as stale.

The pipeline reads artifacts in one place, the runner's composite transport, and its reuse rule
admits only runs at the current commit on a clean tree, which are written under current paths.
So the alias serves artifacts opened by hand; the runner uses the same loader so that every read
gets the same completeness check.
"""

from __future__ import annotations

import dataclasses
import importlib
import pickle
import sys
from typing import IO, TYPE_CHECKING, Any, cast

import numpy as np

# Old top-level package -> where it lives now. The lookup key is only ever the FIRST dotted
# component of a pickled module path, so an entry covers one whole package that moved;
# renaming a single module inside a package cannot be expressed here. A missing entry fails
# loudly at load.
#
# These are the packages that sat at the repository root before `src/quebra/` existed.
MODULE_ALIASES = {
    "analyzers": "quebra.analyzers",
    "core": "quebra.core",
    "loaders": "quebra.loaders",
    "panels": "quebra.panels",
    "plots": "quebra.plots",
    "schemas": "quebra.schemas",
    "transforms": "quebra.transforms",
}


class _AliasingUnpickler(pickle.Unpickler):
    """Resolve a class under its moved package path, if and only if it is the same class name."""

    def find_class(self, module: str, name: str) -> Any:
        head, _, rest = module.partition(".")
        target = MODULE_ALIASES.get(head)
        if target is None:
            return super().find_class(module, name)
        # The package's new home must import. If it does not, the install is broken, which
        # is not a rename, and its own error says so.
        importlib.import_module(target)
        moved = f"{target}.{rest}" if rest else target
        try:
            importlib.import_module(moved)
        except ModuleNotFoundError as exc:
            # Only the moved module itself being absent means "renamed, not moved". A missing
            # dependency or a broken import inside it is a real error and propagates as is.
            if exc.name is None or not (
                moved == exc.name or moved.startswith(f"{exc.name}.")
            ):
                raise
            raise ModuleNotFoundError(
                f"No module named {module!r}: the artifact names {module}.{name}, and its "
                f"moved home {moved} does not exist either, so it was renamed or "
                "restructured, not moved. Re-run the job that wrote it.",
                name=module,
            ) from exc
        owner: object = sys.modules[moved]
        for part in name.split("."):
            if not hasattr(owner, part):
                raise AttributeError(
                    f"the artifact names {module}.{name}; {moved} exists but defines no "
                    f"{name}, so the class was renamed or restructured, not moved. Re-run "
                    "the job that wrote it."
                )
            owner = getattr(owner, part)
        return super().find_class(moved, name)


def _slot_names(cls: type) -> set[str]:
    names: set[str] = set()
    for klass in cls.__mro__:
        slots = klass.__dict__.get("__slots__", ())
        names.update((slots,) if isinstance(slots, str) else slots)
    return names


def _missing_fields(obj: object) -> list[str]:
    """Fields the instance's class declares that its restored state does not hold.

    A class-level default does not count as present: that default is exactly what a stale
    pickle would silently read. An `init=False` field is exempt, because `__init__` never
    writes one, so a fresh instance lacks it too and stale cannot be told from current.
    """
    state = getattr(obj, "__dict__", None)
    slots = _slot_names(type(obj))
    missing = []
    for f in dataclasses.fields(cast("DataclassInstance", obj)):
        if not f.init:
            continue
        if state is not None and f.name in state:
            continue
        if f.name in slots:
            try:
                object.__getattribute__(obj, f.name)
                continue
            except AttributeError:
                pass
        missing.append(f.name)
    return missing


_LEAVES = (str, bytes, int, float, complex, bool, type(None), np.ndarray, np.generic)


def _require_complete(root: object) -> None:
    """Raise if any dataclass reachable from `root` lacks a field its class declares now."""
    seen: set[int] = set()
    stack = [root]
    while stack:
        obj = stack.pop()
        if isinstance(obj, _LEAVES) or id(obj) in seen:
            continue
        seen.add(id(obj))
        if isinstance(obj, dict):
            stack.extend(obj.keys())
            stack.extend(obj.values())
        elif isinstance(obj, (list, tuple, set, frozenset)):
            stack.extend(obj)
        elif dataclasses.is_dataclass(obj) and not isinstance(obj, type):
            missing = _missing_fields(obj)
            if missing:
                raise ValueError(
                    f"stale {type(obj).__name__} artifact: pickle lacks field(s) "
                    f"{sorted(missing)} - it predates a field added to this class, and "
                    "reading it would return a class default or fail mid-render. Re-run the "
                    "job that wrote it (for a composite: without --reuse-deps)."
                )
            stack.extend(
                getattr(obj, f.name)
                for f in dataclasses.fields(obj)
                if hasattr(obj, f.name)
            )


def load_artifact(handle: IO[bytes]) -> Any:
    """Unpickle a materialized artifact: resolve moved packages, then refuse a stale one."""
    obj = _AliasingUnpickler(handle).load()
    _require_complete(obj)
    return obj


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
        # An `init=False` field is exempt, as in `load_artifact`'s walk: `__init__` never
        # writes one, so a fresh instance's state lacks it too.
        missing = {
            f.name
            for f in dataclasses.fields(cast("type[DataclassInstance]", cls))
            if f.init
        } - state.keys()
        if missing:
            raise ValueError(
                f"stale {cls.__name__} artifact: pickle lacks field(s) "
                f"{sorted(missing)} - it predates a field added to this artifact. "
                "Re-run the sub-job (run the composite without --reuse-deps)."
            )
        self.__dict__.update(state)
