"""Tellability wiring: one tellability step per metric record, materialised for reading.

Apart from `recipes.py`, which never imports it, so the tellability modules stay out of the
identity closure of every job that does not ask for them, the check ledgers and the
Kaplan-Meier jobs among them. A job calls the configure function it already uses, which
returns its nodes, and then the matching `wire_*_tellability` below.

Every tellability step reads the same reads as the carve beside it, under the same gap rule,
and every setting is a step argument, so it reaches the provenance label.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

from quebra.analyzers import fidelity, tellability
from quebra.analyzers.fidelity import FidelityResult
from quebra.analyzers.t2star import T2StarResult
from quebra.analyzers.windows import WindowsResult
from quebra.core.job import Job

# One seed and one permutation count for every dependence panel. B + 1 = 1000 stays well
# above the panel's five members over alpha = 0.05 (`min_p`).
TELLABILITY_SEED = 20261014
TELLABILITY_PERMUTATIONS = 999

SIGMA_NOT_DERIVED = "not computed: sigma propagation not derived"


def _thresholds(
    rows: Sequence[Sequence[object]],
) -> tuple[tuple[str, float, bool], ...]:
    return tuple((str(r[0]), float(r[1]), bool(r[2])) for r in rows)  # type: ignore[arg-type]


def _t2star_tellability(
    result: T2StarResult,
    windows_result: WindowsResult,
    *,
    thresholds: Sequence[Sequence[object]],
    gap_mult: float,
    unsure_k: float,
    n_permutations: int,
    seed: int,
) -> tellability.TellabilityResult:
    frame = result.frame
    declared = "t2star_error_s" in frame.columns
    return tellability.run(
        tellability.TellabilityInputs(
            dataset_id=str(result.meta["dataset_id"]),
            metric_key="t2star",
            value_key="t2star_s",
            t_rel_s=frame["t_rel_s"].to_numpy(dtype=float),
            values=frame["t2star_s"].to_numpy(dtype=float),
            sigma=frame["t2star_error_s"].to_numpy(dtype=float) if declared else None,
            sigma_status=tellability.SIGMA_DECLARED
            if declared
            else tellability.NOT_DECLARED,
            thresholds=_thresholds(thresholds),
            windows=windows_result,
            unsure_k=unsure_k,
            gap_mult=gap_mult,
            clip_floor=None,
            n_permutations=n_permutations,
            seed=seed,
        )
    )


def _frequency_tellability(
    norm: Mapping[str, object],
    *,
    value_key: str,
    sigma_key: str,
    gap_mult: float,
    n_permutations: int,
    seed: int,
) -> tellability.TellabilityResult:
    meta = norm.get("meta", {})
    sigma = norm.get(sigma_key)
    return tellability.run(
        tellability.TellabilityInputs(
            dataset_id=str(meta.get("dataset_id", ""))
            if isinstance(meta, Mapping)
            else "",
            metric_key="frequency",
            value_key=value_key,
            t_rel_s=np.asarray(norm["t_rel_s"], dtype=float),
            values=np.asarray(norm[value_key], dtype=float),
            sigma=None if sigma is None else np.asarray(sigma, dtype=float),
            sigma_status=tellability.SIGMA_DECLARED
            if sigma is not None
            else tellability.NOT_DECLARED,
            thresholds=(),
            windows=None,
            unsure_k=1.0,
            gap_mult=gap_mult,
            clip_floor=None,
            n_permutations=n_permutations,
            seed=seed,
        )
    )


def _infidelity_tellability(
    result: FidelityResult,
    windows_result: WindowsResult,
    *,
    thresholds: Sequence[Sequence[object]],
    gap_mult: float,
    clip_floor: float,
    n_permutations: int,
    seed: int,
) -> tellability.TellabilityResult:
    return tellability.run(
        tellability.TellabilityInputs(
            dataset_id=str(result.meta["dataset_id"]),
            metric_key="infidelity",
            value_key="infidelity",
            t_rel_s=result.frame["t_rel_s"].to_numpy(dtype=float),
            values=fidelity.panel_series(result),
            sigma=None,
            sigma_status=SIGMA_NOT_DERIVED,
            thresholds=_thresholds(thresholds),
            windows=windows_result,
            unsure_k=float(windows_result.meta["k"]),
            gap_mult=gap_mult,
            clip_floor=clip_floor,
            n_permutations=n_permutations,
            seed=seed,
        )
    )


def wire_t2star_tellability(
    job: Job,
    refs: Mapping[str, object],
    *,
    prefix: str,
    n_permutations: int = TELLABILITY_PERMUTATIONS,
    seed: int = TELLABILITY_SEED,
) -> object:
    """Tellability of the T2* record, on the fit and the windows of `configure_t2star_job`."""
    node = job.step(
        _t2star_tellability,
        refs["t2star"],
        refs["windows"],
        name="tellability_t2star",
        thresholds=refs["thresholds"],
        gap_mult=refs["gap_mult"],
        unsure_k=refs["k"],
        n_permutations=n_permutations,
        seed=seed,
    )
    job.materialize(node, name=f"{prefix}_tellability")
    return node


def wire_frequency_tellability(
    job: Job,
    refs: Mapping[str, object],
    *,
    prefix: str,
    value_key: str = "raw_frequency_hz",
    sigma_key: str = "frequency_error_hz",
    n_permutations: int = TELLABILITY_PERMUTATIONS,
    seed: int = TELLABILITY_SEED,
) -> object:
    """Tellability of the fitted frequency of a Ramsey record, with no thresholds, on the
    filtered reads of `configure_ramsey_job`. `value_key` is recorded on the artifact."""
    node = job.step(
        _frequency_tellability,
        refs["final_filtered"],
        name="tellability_frequency",
        value_key=value_key,
        sigma_key=sigma_key,
        gap_mult=refs["gap_mult"],
        n_permutations=n_permutations,
        seed=seed,
    )
    job.materialize(node, name=f"{prefix}_tellability")
    return node


def wire_infidelity_tellability(
    job: Job,
    refs: Mapping[str, object],
    *,
    prefix: str,
    n_permutations: int = TELLABILITY_PERMUTATIONS,
    seed: int = TELLABILITY_SEED,
) -> object:
    """Tellability of the infidelity of a Ramsey record, on the fidelity result and the
    fidelity windows of `configure_ramsey_job`. The infidelity's sigma is not derived, so
    every row that reads sigma says so."""
    node = job.step(
        _infidelity_tellability,
        refs["fidelity_raw"],
        refs["fidelity_windows"],
        name="tellability_infidelity",
        thresholds=refs["fidelity_thresholds"],
        gap_mult=refs["gap_mult"],
        clip_floor=fidelity.INFIDELITY_FLOOR,
        n_permutations=n_permutations,
        seed=seed,
    )
    job.materialize(node, name=f"{prefix}_tellability")
    return node
