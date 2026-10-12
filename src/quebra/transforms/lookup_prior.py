"""Validation of a unix-seconds scalar.

`check_unix_s` is what the default normaliser runs on an explicit
`Dataset.extra['run_start_unix_s']` (`schemas/ramsey_series.py`). Its public import path is
`quebra.transforms.lookup_prior.check_unix_s`, which callers and the import contract
(`pyproject.toml`) name, so the module keeps that name.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def check_unix_s(value: object, *, label: str = "unix_s") -> float:
    """Return `value` as float unix seconds, or raise.

    Raises TypeError when `value` is not a number and ValueError when it is not finite or
    falls outside 2000-01-01 to 2100-01-01, the range that catches a millisecond or
    nanosecond value passed as seconds.
    """

    try:
        unix_s = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{label} must be a number of seconds, got {value!r}") from exc

    if not np.isfinite(unix_s):
        raise ValueError(f"{label} must be finite, got {unix_s!r}")

    lower = pd.Timestamp("2000-01-01").value / 1e9
    upper = pd.Timestamp("2100-01-01").value / 1e9
    if not (lower <= unix_s <= upper):
        raise ValueError(
            f"{label}={unix_s:.6g} looks implausible; expected seconds since epoch in range "
            f"[{lower:.6g}, {upper:.6g}]."
        )
    return unix_s
