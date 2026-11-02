#!/usr/bin/env python3
"""Count, from check-ledger artifacts on disk, where C1 and C2 gave no answer and where the
trend group is CvM alone. Read-only: it opens the artifacts and writes nothing.

    python scripts/count_trend_answers.py <ledger run folder> [<ledger run folder> ...]

Each folder is one ledger job's run; every `*_check_ledger_data.pkl` in it is one record's
ledger. The counting and its definitions are `independence_survey.trend_group_counts`.

Artifacts load through `load_artifact`, which refuses a stale one. If any is refused the
script names each, prints no count and exits 2: a count over some of the ledgers would read
as a count over all of them.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from quebra.analyzers.independence_survey import trend_group_counts
from quebra.core._artifact_guard import load_artifact
from quebra.ledger_recipe import LEDGER_DATA_SUFFIX


def load_ledgers(folders: list[Path]) -> tuple[list[object], list[tuple[Path, str]]]:
    """Every record ledger in `folders`, and the files refused as stale with the reason."""
    ledgers: list[object] = []
    refused: list[tuple[Path, str]] = []
    for folder in folders:
        paths = sorted(folder.glob(f"*{LEDGER_DATA_SUFFIX}.pkl"))
        if not paths:
            raise FileNotFoundError(f"{folder} holds no *{LEDGER_DATA_SUFFIX}.pkl")
        for path in paths:
            try:
                with path.open("rb") as handle:
                    ledgers.append(load_artifact(handle))
            except ValueError as exc:
                if not str(exc).startswith("stale "):
                    raise
                refused.append((path, str(exc)))
    return ledgers, refused


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("folders", nargs="+", type=Path)
    args = parser.parse_args(argv)
    ledgers, refused = load_ledgers(args.folders)
    if refused:
        for path, reason in refused:
            print(f"refused {path}: {reason}")
        print(
            f"refused {len(refused)} of {len(refused) + len(ledgers)} ledgers as stale; "
            "no count printed"
        )
        return 2
    counts = trend_group_counts(ledgers)
    with pd.option_context("display.width", 200, "display.max_colwidth", 120):
        print(f"{len(ledgers)} ledgers from {len(args.folders)} folders\n")
        print("Cells where neither C1 nor C2 answered, by bag, clock and message:")
        print(counts.no_answer.to_string(index=False))
        print("\nThe trend group, per bag and clock:")
        print(
            counts.cvm_alone.to_string(index=False, float_format=lambda v: f"{v:.4f}")
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
