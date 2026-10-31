"""Join a record's declared Rabi frequency onto its normalized reads.

The Rabi frequency of a record is a DECLARED value: one row of a table supplied with the
data, matched exactly on (device, run_date, qubit). There is no time column and no "latest
value before" search, so a record without a row of its own raises rather than borrowing a
neighbour's, and a key that names two rows raises rather than picking one.

The table carries `rabi_khz`, a cyclic frequency in kHz, as given by its source; the join
writes `rabi_hz`, the cyclic hertz the fidelity analyzer reads, broadcast over the reads, and
records where it came from in `meta["rabi_source"]` and `meta["rabi_declared_key"]`. A Rabi
column declared in an angular unit is refused, naming it: x = Delta / Omega would move by
2 pi. A column declares one when its name holds a `rad` or `angular` token beside a token
starting `rabi` or `omega`, tokens being the runs of letters and digits.
"""

from __future__ import annotations

import re
from collections.abc import Mapping

import numpy as np
import pandas as pd

KEY_COLUMNS = ("device", "run_date", "qubit")
VALUE_COLUMN = "rabi_khz"
SOURCE_LABEL = "declared"
_TOKEN = re.compile(r"[a-z0-9]+")
_ANGULAR_TOKENS = frozenset({"rad", "angular"})


def _is_angular_rabi(column: str) -> bool:
    """Whether a column's name declares a Rabi frequency in an angular unit."""
    tokens = _TOKEN.findall(column.lower())
    names_rabi = any(t.startswith(("rabi", "omega")) for t in tokens)
    return names_rabi and not _ANGULAR_TOKENS.isdisjoint(tokens)


def declared_key(device: str, run_date: str, qubit: int) -> str:
    """The key as one string, the form it takes in metadata and in error messages."""
    return f"{device}/{run_date}/qubit{qubit}"


def join_declared_rabi(
    norm: Mapping[str, object],
    table: pd.DataFrame,
    *,
    device: str,
    run_date: str,
    qubit: int,
) -> dict[str, object]:
    """Return `norm` with `rabi_hz` taken from the one table row matching the key."""
    angular = [str(c) for c in table.columns if _is_angular_rabi(str(c))]
    if angular:
        raise ValueError(
            f"the declared Rabi table carries {angular}, an angular unit; the join reads "
            f"the cyclic frequency in kHz from {VALUE_COLUMN!r} only"
        )
    missing = [c for c in (*KEY_COLUMNS, VALUE_COLUMN) if c not in table.columns]
    if missing:
        raise KeyError(f"the declared Rabi table has no column(s) {missing}")
    if "t_rel_s" not in norm:
        raise KeyError(
            "the declared Rabi join needs 't_rel_s' in the normalized mapping"
        )

    if "rabi_hz" in norm:
        raise KeyError(
            "the normalized mapping already carries 'rabi_hz'; a record's Rabi frequency "
            "has one source, and the declared join will not overwrite another"
        )

    key = declared_key(device, run_date, qubit)
    match = (
        (table["device"].astype(str) == str(device))
        & (table["run_date"].astype(str) == str(run_date))
        & (table["qubit"].astype(int) == int(qubit))
    )
    rows = table[match]
    if len(rows) == 0:
        raise KeyError(
            f"no declared Rabi frequency for {key}: every record needs a row of its own, "
            f"and none is borrowed from another date or qubit"
        )
    if len(rows) > 1:
        raise ValueError(
            f"{len(rows)} declared Rabi rows for {key}: the key must name exactly one row"
        )
    rabi_khz = float(rows[VALUE_COLUMN].iloc[0])
    if not np.isfinite(rabi_khz) or rabi_khz <= 0.0:
        raise ValueError(
            f"the declared Rabi frequency for {key} is {rabi_khz!r} kHz; "
            f"it must be finite and positive"
        )

    out = dict(norm)
    n_reads = len(np.asarray(norm["t_rel_s"]))
    out["rabi_hz"] = np.full(n_reads, rabi_khz * 1e3)
    meta = dict(norm.get("meta") or {})
    meta["rabi_source"] = SOURCE_LABEL
    meta["rabi_declared_key"] = key
    out["meta"] = meta
    return out
