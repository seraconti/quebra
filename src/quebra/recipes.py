from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Callable

import pandas as pd

import quebra.analyzers.fidelity as fidelity
import quebra.analyzers.t2star as t2star
import quebra.analyzers.windows as windows
from quebra.analyzers.allan import run as run_allan
from quebra.analyzers.fidelity import FidelityResult
from quebra.analyzers.fidelity import make_inputs_from_norm as _fidelity_make_inputs
from quebra.analyzers.fidelity import run as _run_fidelity
from quebra.analyzers.t2star import T2STAR_DEFAULT_LADDER, T2StarResult
from quebra.analyzers.tlf import run as run_tlf
from quebra.analyzers.windows import DEFAULT_GAP_MULT, WindowsResult
from quebra.core.dataset import Dataset
from quebra.core.job import Job
from quebra.panels.within_calibration import (
    WithinCalibrationPanel,
    WithinCalibrationPanelData,
)
from quebra.plots.allan_plot import AllanPlot
from quebra.plots.tlf_plot import TLFPlot
from quebra.transforms.filter import run as run_filter
from quebra.transforms.interpolate import run as run_interpolate
from quebra.transforms.lookup_prior import check_unix_s

RAMSEY_CONFIG: dict[str, object] = {
    "filter": {
        "apply_chi_squared": True,
        "chi_squared_threshold": 0.8,
        "apply_sigma": True,
        "sigma_factor": 3.5,
        "apply_frequency_window": False,
        "frequency_window_hz": {"min": 1.0e6, "max": 3.0e6},
    },
    "allan": {"modes": ["raw"], "taus_mode": "all", "min_points_for_allan": 8},
    "fidelity": {"use_angular_frequency": False},
    "dataset_profile": "overnight",
}


def _copy_config(profile: str) -> dict[str, object]:
    config = dict(RAMSEY_CONFIG)
    config["dataset_profile"] = profile
    config["filter"] = dict(RAMSEY_CONFIG["filter"])
    config["filter"]["frequency_window_hz"] = dict(
        RAMSEY_CONFIG["filter"]["frequency_window_hz"]
    )
    config["allan"] = dict(RAMSEY_CONFIG["allan"])
    config["fidelity"] = dict(RAMSEY_CONFIG["fidelity"])
    return config


def run_start_unix_s_from_hdf5(path: str | Path) -> float:
    """Derive run_start_unix_s from core-tools HDF5 measurement_time.

    The HDF5 attribute is a naive local timestamp; we keep it naive-as-UTC
    and validate with check_unix_s.
    """
    import h5py

    file_path = Path(path)
    with h5py.File(file_path, "r") as handle:
        measurement_time = handle.attrs.get("measurement_time")

    if measurement_time is None:
        raise KeyError(f"HDF5 file {file_path} is missing measurement_time attribute")

    run_start_local_dt = pd.Timestamp(str(measurement_time))
    return check_unix_s(float(run_start_local_dt.value / 1e9), label="run_start_unix_s")


def _filter_step(
    config: Mapping[str, object],
) -> Callable[[dict[str, object]], dict[str, object]]:
    def step(norm: dict[str, object]) -> dict[str, object]:
        return run_filter(norm, config, None)

    return step


def _interpolate_step(
    config: Mapping[str, object],
) -> Callable[[dict[str, object]], dict[str, object]]:
    def step(norm: dict[str, object]) -> dict[str, object]:
        return run_interpolate(norm, config)

    return step


def _allan_step(config: Mapping[str, object]) -> Callable[[dict[str, object]], object]:
    def step(
        norm: dict[str, object],
        fractional: bool = False,
        carrier_col: str = "frequency",
    ) -> object:
        return run_allan(norm, config, fractional=fractional, carrier_col=carrier_col)

    return step


def _fidelity_step(
    config: Mapping[str, object],
) -> Callable[[dict[str, object]], object]:
    def step(norm: dict[str, object]) -> object:
        return _run_fidelity(_fidelity_make_inputs(norm, config))

    return step


def _tlf_step() -> Callable[..., dict[str, object]]:
    def step(norm: dict[str, object], seed: int) -> dict[str, object]:
        # TLF analysis on filtered (but NOT interpolated) data to preserve noise metrics.
        # Input norm is the final filtered stage from filter step, with uninterpolated timestamps.
        if "raw_frequency_hz" in norm:
            values_hz = norm["raw_frequency_hz"]
        elif "delta_hz" in norm:
            values_hz = norm["delta_hz"]
        else:
            raise KeyError(
                "TLF analysis requires 'raw_frequency_hz' or 'delta_hz' in normalized mapping"
            )

        # timestamps (seconds, relative to start) required for dynamics computation
        if "t_rel_s" not in norm:
            raise KeyError(
                "TLF analysis requires 't_rel_s' (relative seconds) in normalized mapping for dynamics"
            )
        timestamps = norm["t_rel_s"]

        result = run_tlf(values_hz, timestamps, seed=seed)
        return {
            "result": result,
            "values_hz": values_hz,
            "meta": dict(norm.get("meta", {})),
        }

    return step


def _final_stage(bundle: object) -> dict[str, object]:
    from quebra.transforms.filter import FilterResult

    if isinstance(bundle, FilterResult):
        return dict(bundle.final_norm)
    # legacy dict path (unmigrated callers)
    return dict(bundle["stages"][bundle["final_stage"]])


def _fidelity_windows(result: FidelityResult, gap_mult: float) -> WindowsResult:
    """Carve in-spec windows on the infidelity series.

    The ladder is data-derived, so this calls the SAME pure `panel_thresholds` on the
    SAME clipped series the panel adapter uses - the two nodes cannot disagree.
    """
    series = fidelity.panel_series(result)
    return windows.run(
        windows.WindowsInputs(
            t_rel_s=result.frame["t_rel_s"].to_numpy(dtype=float),
            values=series,
            thresholds=fidelity.panel_thresholds(series),
            dataset_id=str(result.meta.get("dataset_id", "")),
            gap_mult=gap_mult,
        )
    )


# ONE definition, imported by the jobs. The literal was in three files after the seed fix,
# which is a drift risk of exactly the kind this repo's provenance rules exist to prevent:
# two jobs quoting different seeds while both labels claim to be reproducible.
#
# NOTE, carried as open work: one seed serves every window. `paired_permutation_test` keys
# `default_rng(seed)` on the seed alone, so two windows of equal n draw an IDENTICAL
# permutation sequence and their Monte Carlo error does not average down when
# `for_windows` takes the median. Per-window seeds are the fix; it is a design change, not
# a rename.
XI_SEED = 20260813

# GaussianMixture initialises by k-means, so an unseeded fit makes `is_bimodal` and every
# dwell statistic a fresh random variable per call. Passed as a step kwarg, so it reaches the
# parameter row and the provenance label.
TLF_SEED = 20260902


def _fidelity_panel_data(
    result: FidelityResult, window_result: WindowsResult, xi_seed: int
) -> WithinCalibrationPanelData:
    return fidelity.make_panel_data(
        result,
        windows=window_result.windows,
        reads=window_result.reads,
        gap_spans_s=window_result.diagnostics["gap_spans_s"],
        # `windows.run` always records this, so a missing key is a broken artifact rather
        # than an old one; defaulting it to False would silently redraw the panel in the
        # other mode. This repo raises instead of falling back.
        use_uncertainty=bool(window_result.meta["use_uncertainty"]),
        dataset_id=str(result.meta.get("dataset_id", "")),
        xi_seed=xi_seed,
    )


def configure_ramsey_job(
    job: Job,
    dataset: object,
    *,
    profile: str,
    include_fidelity: bool,
    include_tlf: bool = False,
    allan_fractional: bool = False,
    allan_carrier_col: str = "qubit_frequency_hz",
    xi_seed: int = XI_SEED,
    tlf_seed: int = TLF_SEED,
    figure_prefix: str | None = None,
) -> None:
    config = _copy_config(profile)
    # dataset may be a Dataset (to load) or an already-registered Reference
    # (LocalRef/ArtifactRef - a node whose result is loaded/enriched upstream).
    if hasattr(dataset, "resolve"):
        raw = dataset
    else:
        raw = job.load(dataset)

    # Decide fidelity inclusion based on device family when caller didn't explicitly
    # disable it. Device-aware branching keeps the job semantics clear per dataset.
    dev = None
    # if raw was created from a Dataset node, that original Dataset lives in kwargs
    if hasattr(raw, "kwargs") and isinstance(raw.kwargs.get("dataset"), Dataset):
        dev = raw.kwargs.get("dataset").device
    elif hasattr(dataset, "device"):
        dev = dataset.device
    if dev is not None and isinstance(dev, str) and dev.lower().startswith("2x2"):
        include_fidelity = False
        include_tlf = False
        allan_fractional = False
    filtered = job.step(_filter_step(config), raw, name="filter")
    final_filtered = job.step(_final_stage, filtered, name="final_filter_stage")
    interpolated = job.step(
        _interpolate_step(config), final_filtered, name="interpolate"
    )
    allan = job.step(
        _allan_step(config),
        interpolated,
        name="allan",
        fractional=allan_fractional,
        carrier_col=allan_carrier_col,
    )

    prefix = figure_prefix or job.name
    job.figure(
        AllanPlot, allan, targets=["static", "academic"], title=f"{prefix} Allan"
    )

    if include_fidelity:
        fidelity_raw = job.step(
            _fidelity_step(config), final_filtered, name="fidelity_raw"
        )
        # Fidelity is carved and drawn from the FILTERED reads, never the interpolated
        # ones. Interpolation puts points on a uniform grid, which erases exactly the
        # read gaps the carve exists to find: on 100423 the interpolated fidelity
        # reported gaps=0 across a real 30-minute gap, because 99.8% of its points
        # were manufactured. Allan still uses `interpolated` - it needs uniform
        # sampling and is not downstream of a carve.
        #
        # Panel data is a STEP, not a draw-time side effect: the carve has to be a node
        # so its gap policy reaches provenance and its tables can be materialized.
        fidelity_windows = job.step(
            _fidelity_windows,
            fidelity_raw,
            name="fidelity_windows",
            gap_mult=DEFAULT_GAP_MULT,
        )
        fidelity_panel = job.step(
            _fidelity_panel_data,
            fidelity_raw,
            fidelity_windows,
            name="fidelity_panel_data",
            # Declared, not defaulted: see the note in jobs/active/t2star_q1_070423.py.
            xi_seed=xi_seed,
        )
        job.figure(
            WithinCalibrationPanel,
            fidelity_panel,
            targets=["static", "academic"],
            title=f"{prefix} Fidelity",
        )
        job.materialize(fidelity_raw, name=f"{prefix}_fidelity_raw")

    if include_tlf:
        tlf = job.step(_tlf_step(), final_filtered, name="tlf", seed=tlf_seed)
        job.figure(TLFPlot, tlf, targets=["static", "academic"], title=f"{prefix} TLF")
        job.materialize(tlf, name=f"{prefix}_tlf")


# ---------------------------------------------------------------------------- the T2* family

# The T2* family's threshold ladder, in SI seconds.
#
# DERIVED from the analyzer's own default rather than retyped. A first version of this
# declared its own comprehension under the comment "one ladder, one place", which was false:
# `analyzers.t2star.T2STAR_DEFAULT_LADDER` already held the same ten rungs, so three copies
# became two and nothing pinned them equal. They coincide today, and a silent divergence
# would put the survey's columns out of correspondence with the panel's.
#
# A family that needs a different ladder passes `thresholds=` - which is what the analyzer's
# own docstring asks production jobs to do, and what keeps this a default rather than a rule.
T2STAR_THRESHOLDS: list[tuple[str, float, bool]] = list(T2STAR_DEFAULT_LADDER)


def configure_t2star_job(
    job: "Job",
    *,
    dataset: "Dataset",
    prefix: str,
    thresholds: list[tuple[str, float, bool]] | None = None,
    gap_mult: float = 10.0,
    k: float = 1.0,
    use_uncertainty: bool = True,
    shape_min_reads: int = 5,
    xi_seed: int = XI_SEED,
) -> None:
    """Wire the whole T2* within-calibration graph onto `job`.

    Collapses the family: two job files that were
    byte-identical once the date and the run duration were normalised. What is left in each
    job file is its parameter row, which is the part a reader should actually diff.

    Deliberately does NOT call `configure_ramsey_job` and adds NO interpolate node. Windows
    are carved from the filtered reads: the gap policy is meaningless on a uniform grid, and
    an interpolated point is not an observation, so a window must never be built from one.
    `tests/test_windows_not_interpolated.py` pins that as a property of the DAG.

    `gap_mult`, `k`, `use_uncertainty`, `shape_min_reads` and `xi_seed` are passed as step
    KWARGS rather than captured, so they reach the provenance label - the `allan` pattern,
    not the `filter` pattern. `runner` builds that label from `node.kwargs`, so an argument
    left to its default would be invisible to provenance; they are therefore passed
    explicitly below even where the value equals the default.
    """
    ladder = T2STAR_THRESHOLDS if thresholds is None else thresholds

    def _t2star_run(norm: object) -> T2StarResult:
        return t2star.run(t2star.make_inputs_from_norm(norm))  # type: ignore[arg-type]

    def _windows_run(
        result: T2StarResult,
        gap_mult: float,
        k: float,
        use_uncertainty: bool,
        thresholds: list[tuple[str, float, bool]],
    ) -> WindowsResult:
        return windows.run(
            windows.make_inputs_from_frame(
                result.frame,
                time_col="t_rel_s",
                value_col="t2star_s",
                sigma_col="t2star_error_s" if use_uncertainty else None,
                thresholds=thresholds,
                dataset_id=str(result.meta.get("dataset_id", "")),
                gap_mult=gap_mult,
                k=k,
                use_uncertainty=use_uncertainty,
            )
        )

    def _t2star_panel_data(
        result: T2StarResult,
        window_result: WindowsResult,
        shape_min_reads: int,
        use_uncertainty: bool,
        xi_seed: int,
        thresholds: list[tuple[str, float, bool]],
    ) -> WithinCalibrationPanelData:
        return t2star.make_panel_data(
            result,
            windows=window_result.windows,
            reads=window_result.reads,
            gap_spans_s=window_result.diagnostics["gap_spans_s"],
            thresholds=thresholds,
            shape_min_reads=shape_min_reads,
            use_uncertainty=use_uncertainty,
            xi_seed=xi_seed,
        )

    main_node = job.load(dataset)
    filtered = job.step(_filter_step(RAMSEY_CONFIG), main_node, name="t2star_filter")
    final = job.step(_final_stage, filtered, name="t2star_final_filter_stage")
    result = job.step(_t2star_run, final, name="t2star")
    window_node = job.step(
        _windows_run,
        result,
        name="windows",
        gap_mult=gap_mult,
        k=k,
        use_uncertainty=use_uncertainty,
        # A kwarg, not the closure capture it was: the ladder defines every rung of the
        # panel, and this function's own docstring says an argument left out of node.kwargs
        # is invisible to the provenance label. It was the one parameter breaking that rule.
        thresholds=ladder,
    )
    panel = job.step(
        _t2star_panel_data,
        result,
        window_node,
        name="t2star_panel_data",
        shape_min_reads=shape_min_reads,
        use_uncertainty=use_uncertainty,
        xi_seed=xi_seed,
        thresholds=ladder,
    )
    job.materialize(window_node, name=f"{prefix}_windows")
    job.figure(
        WithinCalibrationPanel,
        panel,
        targets=["static", "academic"],
        title=f"{prefix}_t2star",
    )
