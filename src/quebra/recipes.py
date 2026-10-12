from __future__ import annotations

from collections.abc import Mapping
import re
from pathlib import Path
from typing import Callable

import pandas as pd

import quebra.analyzers.event_table as event_table
import quebra.analyzers.fidelity as fidelity
import quebra.analyzers.kaplan_meier as kaplan_meier
import quebra.analyzers.survival as survival
import quebra.analyzers.t2star as t2star
import quebra.analyzers.windows as windows
from quebra.analyzers.allan import run as run_allan
from quebra.analyzers.event_table import EventTables
from quebra.analyzers.fidelity import FidelityResult
from quebra.analyzers.fidelity import make_inputs_from_norm as _fidelity_make_inputs
from quebra.analyzers.fidelity import run as _run_fidelity
from quebra.analyzers.kaplan_meier import CONF_LEVEL, KaplanMeierSet
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
from quebra.transforms.declared_rabi import join_declared_rabi
from quebra.transforms.filter import run as run_filter
from quebra.transforms.lookup_prior import check_unix_s

# Unless a comment cites a source, a value here was ported unchanged from the predecessor
# pipeline's "overnight" profile, and no source for it is located.
RAMSEY_CONFIG: dict[str, object] = {
    "filter": {
        # Ported; the guard it switches is cited in `transforms/filter.py`.
        "apply_chi_squared": True,
        # Capannelli et al., Supplementary notes, section III.C, p. 9: drift datasets scoring
        # a normalised chi-squared of 0.8 or more are filtered out. The released records hold
        # the points left after that filter (Supplementary Figs. 10 and 11, pp. 10-11).
        "chi_squared_threshold": 0.8,
    },
    # Ported; no source located. Allan computes the raw detuning series only.
    "allan": {"modes": ["raw"], "taus_mode": "all", "min_points_for_allan": 8},
    # Ported. It has no effect on the result: the gate-fidelity formula depends on detuning
    # and Rabi frequency only through their ratio, which a common factor of 2*pi leaves alone.
    "fidelity": {"use_angular_frequency": False},
    # Ported. Overwritten per job by `configure_ramsey_job(profile=...)`.
    "dataset_profile": "overnight",
}


def _copy_config(profile: str) -> dict[str, object]:
    config = dict(RAMSEY_CONFIG)
    config["dataset_profile"] = profile
    config["filter"] = dict(RAMSEY_CONFIG["filter"])
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


def _fidelity_windows(
    result: FidelityResult,
    gap_mult: float,
    thresholds: list[tuple[str, float, bool]],
) -> WindowsResult:
    """Carve windows on the infidelity series, against the declared ladder.

    The panel-data step receives the same `thresholds` kwarg, so the two nodes carve and
    render one ladder; both read the series through `fidelity.panel_series`.
    """
    series = fidelity.panel_series(result)
    return windows.run(
        windows.WindowsInputs(
            t_rel_s=result.frame["t_rel_s"].to_numpy(dtype=float),
            values=series,
            thresholds=list(thresholds),
            dataset_id=str(result.meta["dataset_id"]),
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


def _event_tables(window_result: WindowsResult) -> EventTables:
    return event_table.event_tables_from_carve(window_result)


def _kaplan_meier_set(tables: EventTables, conf_level: float) -> KaplanMeierSet:
    return kaplan_meier.kaplan_meier_set(tables, conf_level=conf_level)


def wire_kaplan_meier(
    job: "Job", window_node: object, *, conf_level: float, node_suffix: str = ""
) -> tuple[object, object]:
    """Wire the event tables and the ONE Kaplan-Meier node off a carve; return both refs.

    The panel and every other consumer read this node, so one set of windows cannot ship
    two bands. `conf_level` is a step kwarg, so it reaches the provenance label.
    """
    tables = job.step(_event_tables, window_node, name=f"event_tables{node_suffix}")
    km = job.step(
        _kaplan_meier_set,
        tables,
        name=f"kaplan_meier{node_suffix}",
        conf_level=conf_level,
    )
    return tables, km


def wire_survival(
    job: "Job",
    tables_node: object,
    *,
    conf_level: float,
    quantiles: tuple[float, ...] = survival.DEFAULT_QUANTILES,
    rmst_tau_s: float | None = None,
    turnbull_tol: float = survival.TURNBULL_TOL,
    turnbull_max_iter: int = survival.TURNBULL_MAX_ITER,
    node_suffix: str = "",
) -> dict[str, object]:
    """Wire the survival estimators off an event-table node; return their refs by name.

    Nelson-Aalen, the curve summaries, Turnbull and the placement bracket, each one node
    over every threshold and side. Nothing here adds a sink: a job that wants these
    results materializes the refs it uses. Every choice that changes a number is a step
    kwarg, so it reaches the identity and the provenance label; `rmst_tau_s` is ONE kwarg,
    shared by the summaries and the bracket. `conf_level` has no default: pass the value
    given to `wire_kaplan_meier`, since the summaries rebuild the same product-limit band.
    """
    return {
        "nelson_aalen": job.step(
            survival.nelson_aalen_set,
            tables_node,
            name=f"nelson_aalen{node_suffix}",
            conf_level=conf_level,
        ),
        "curve_summaries": job.step(
            survival.curve_summaries_set,
            tables_node,
            name=f"curve_summaries{node_suffix}",
            quantiles=tuple(quantiles),
            rmst_tau_s=rmst_tau_s,
            conf_level=conf_level,
        ),
        "turnbull": job.step(
            survival.turnbull_set,
            tables_node,
            name=f"turnbull{node_suffix}",
            tol=turnbull_tol,
            max_iter=turnbull_max_iter,
        ),
        "placement_bracket": job.step(
            survival.placement_bracket_set,
            tables_node,
            name=f"placement_bracket{node_suffix}",
            rmst_tau_s=rmst_tau_s,
        ),
    }


def _fidelity_panel_data(
    result: FidelityResult,
    window_result: WindowsResult,
    km: KaplanMeierSet,
    xi_seed: int,
    thresholds: list[tuple[str, float, bool]],
) -> WithinCalibrationPanelData:
    return fidelity.make_panel_data(
        result,
        windows=window_result.windows_in_spec,
        reads=window_result.reads,
        gap_spans_s=window_result.diagnostics["gap_spans_s"],
        kaplan_meier=km,
        thresholds=thresholds,
        # `windows.run` always records this, so a missing key is a broken artifact rather
        # than an old one; defaulting it to False would silently redraw the panel in the
        # other mode. This repo raises instead of falling back.
        use_uncertainty=bool(window_result.meta["use_uncertainty"]),
        xi_seed=xi_seed,
    )


# `DDMMYY_<device>_qubit<N>[_<condition>]`: the record-name convention `ramsey_series` and
# `scripts/make_data_manifest.py` also read.
_RECORD_NAME = re.compile(r"^(\d{2})(\d{2})(\d{2})_(\w+?)_qubit(\d+)(?:_(\w+))?$")


def _declared_rabi_key(dataset: object, rabi_table: object) -> tuple[str, str, int]:
    """(device, run_date, qubit) for the declared Rabi join, or raise at graph build.

    The key comes from the record's own `Dataset` declaration, never from a loaded node, and
    each of its three parts must agree with the record's file name, so a key copied from
    another job cannot silently join another record's Rabi frequency.
    """
    if not isinstance(dataset, Dataset):
        raise TypeError(
            "include_fidelity=True needs the record's Dataset, not a loaded node: its "
            "device, qubit and extra['run_date'] key the declared Rabi table"
        )
    if rabi_table is None:
        raise ValueError(
            "include_fidelity=True needs rabi_table=Dataset(...), the table that declares "
            "this record's Rabi frequency"
        )
    if not isinstance(rabi_table, Dataset):
        raise TypeError(
            f"rabi_table must be a Dataset, got {type(rabi_table).__name__}"
        )
    run_date = dataset.extra.get("run_date")
    missing = [
        name
        for name, value in (
            ("device", dataset.device),
            ("qubit", dataset.qubit),
            ("extra['run_date']", run_date),
        )
        if value is None
    ]
    if missing:
        raise ValueError(
            f"the declared Rabi join needs Dataset {', '.join(missing)} for {dataset.path}"
        )
    stem = Path(dataset.path).stem
    match = _RECORD_NAME.match(stem)
    if match is None:
        raise ValueError(
            f"cannot check the Rabi key against {stem!r}: the file name is not "
            f"DDMMYY_<device>_qubit<N>"
        )
    day, month, year, named_device, named_qubit, _condition = match.groups()
    declared = (str(dataset.device), str(run_date), int(dataset.qubit))
    named = (named_device, f"20{year}-{month}-{day}", int(named_qubit))
    for part, said, found in zip(("device", "run_date", "qubit"), declared, named):
        if said != found:
            raise ValueError(
                f"the Dataset declares {part}={said!r} but the file {stem!r} names "
                f"{found!r}; the declared Rabi row would belong to another record"
            )
    return declared


def configure_ramsey_job(
    job: Job,
    dataset: object,
    *,
    profile: str,
    include_fidelity: bool,
    include_tlf: bool = False,
    allan_fractional: bool = False,
    allan_carrier_col: str = "qubit_frequency_hz",
    rabi_table: Dataset | None = None,
    xi_seed: int = XI_SEED,
    tlf_seed: int = TLF_SEED,
    conf_level: float = CONF_LEVEL,
    figure_prefix: str | None = None,
) -> None:
    config = _copy_config(profile)
    # Checked before any node is registered, so a mis-declared job fails at build time.
    if rabi_table is not None and not include_fidelity:
        raise ValueError(
            "rabi_table was given but include_fidelity=False; nothing would read it"
        )
    rabi_key = _declared_rabi_key(dataset, rabi_table) if include_fidelity else None
    # dataset may be a Dataset (to load) or an already-registered Reference
    # (LocalRef/ArtifactRef - a node whose result is loaded/enriched upstream).
    if hasattr(dataset, "resolve"):
        raw = dataset
    else:
        raw = job.load(dataset)

    # `include_fidelity`, `include_tlf` and `allan_fractional` mean what the caller says;
    # nothing here overrides them from the dataset's device. Every analysis reads the
    # filtered reads themselves, never a resampled series.
    filtered = job.step(_filter_step(config), raw, name="filter")
    final_filtered = job.step(_final_stage, filtered, name="final_filter_stage")
    allan = job.step(
        _allan_step(config),
        final_filtered,
        name="allan",
        fractional=allan_fractional,
        carrier_col=allan_carrier_col,
    )

    prefix = figure_prefix or job.name
    job.figure(
        AllanPlot, allan, targets=["static", "academic"], title=f"{prefix} Allan"
    )

    if include_fidelity:
        # Validated by `_declared_rabi_key` before any node was registered.
        device, run_date, qubit = rabi_key  # type: ignore[misc]
        # The join sits on the filtered reads and feeds fidelity only: Allan and TLF never
        # depend on the Rabi table.
        with_rabi = job.step(
            join_declared_rabi,
            final_filtered,
            job.load_df(rabi_table),
            name="declared_rabi",
            device=device,
            run_date=run_date,
            qubit=qubit,
        )
        fidelity_raw = job.step(_fidelity_step(config), with_rabi, name="fidelity_raw")
        ladder = list(fidelity.FIDELITY_LADDER)
        # Fidelity is carved and drawn from the FILTERED reads. A resampled series puts
        # points on a uniform grid, which erases exactly the read gaps the carve exists to
        # find.
        #
        # Panel data is a STEP, not a draw-time side effect: the carve has to be a node
        # so its gap policy reaches provenance and its tables can be materialized.
        fidelity_windows = job.step(
            _fidelity_windows,
            fidelity_raw,
            name="fidelity_windows",
            gap_mult=DEFAULT_GAP_MULT,
            thresholds=ladder,
        )
        _fidelity_tables, fidelity_km = wire_kaplan_meier(
            job, fidelity_windows, conf_level=conf_level, node_suffix="_fidelity"
        )
        fidelity_panel = job.step(
            _fidelity_panel_data,
            fidelity_raw,
            fidelity_windows,
            fidelity_km,
            name="fidelity_panel_data",
            # Declared, not defaulted: see the note in jobs/active/t2star_q1_070423.py.
            xi_seed=xi_seed,
            thresholds=ladder,
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

# The T2* carve parameters, the defaults of both `wire_t2star_carve` and
# `configure_t2star_job`, so a T2* panel and a check ledger cannot default to different
# windows.
T2STAR_GAP_MULT = 10.0
T2STAR_K = 1.0
T2STAR_USE_UNCERTAINTY = True


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
            dataset_id=str(result.meta["dataset_id"]),
            gap_mult=gap_mult,
            k=k,
            use_uncertainty=use_uncertainty,
        )
    )


def _t2star_panel_data(
    result: T2StarResult,
    window_result: WindowsResult,
    km: KaplanMeierSet,
    shape_min_reads: int,
    use_uncertainty: bool,
    xi_seed: int,
    thresholds: list[tuple[str, float, bool]],
) -> WithinCalibrationPanelData:
    return t2star.make_panel_data(
        result,
        windows=window_result.windows_in_spec,
        reads=window_result.reads,
        gap_spans_s=window_result.diagnostics["gap_spans_s"],
        kaplan_meier=km,
        thresholds=thresholds,
        shape_min_reads=shape_min_reads,
        use_uncertainty=use_uncertainty,
        xi_seed=xi_seed,
    )


def wire_t2star_carve(
    job: "Job",
    *,
    dataset: "Dataset",
    thresholds: list[tuple[str, float, bool]] | None = None,
    gap_mult: float = T2STAR_GAP_MULT,
    k: float = T2STAR_K,
    use_uncertainty: bool = T2STAR_USE_UNCERTAINTY,
    node_suffix: str = "",
) -> tuple[object, object]:
    """Wire load, filter, T2* fit and window carve onto `job`; return (fit, windows) nodes.

    The one definition of the T2* carve. `configure_t2star_job` draws a panel on top of it
    and `ledger_recipe.configure_check_ledger_job` scores its windows, so the ledger and the
    panel cannot carve differently. `node_suffix` is appended to every step name (the load
    node is numbered by the job instead), so one job can carve several records with readable
    node ids.

    No interpolate node, deliberately. Windows are carved from the filtered reads: the gap
    policy is meaningless on a uniform grid, and an interpolated point is not an observation,
    so a window must never be built from one. `tests/test_windows_not_interpolated.py` pins
    that as a property of the DAG.

    Every carve parameter is passed as a step KWARG, so it reaches the provenance label even
    where it equals the default.
    """
    ladder = T2STAR_THRESHOLDS if thresholds is None else thresholds
    main_node = job.load(dataset)
    filtered = job.step(
        _filter_step(RAMSEY_CONFIG), main_node, name=f"t2star_filter{node_suffix}"
    )
    final = job.step(
        _final_stage, filtered, name=f"t2star_final_filter_stage{node_suffix}"
    )
    result = job.step(_t2star_run, final, name=f"t2star{node_suffix}")
    window_node = job.step(
        _windows_run,
        result,
        name=f"windows{node_suffix}",
        gap_mult=gap_mult,
        k=k,
        use_uncertainty=use_uncertainty,
        thresholds=ladder,
    )
    return result, window_node


def configure_t2star_job(
    job: "Job",
    *,
    dataset: "Dataset",
    prefix: str,
    thresholds: list[tuple[str, float, bool]] | None = None,
    gap_mult: float = T2STAR_GAP_MULT,
    k: float = T2STAR_K,
    use_uncertainty: bool = T2STAR_USE_UNCERTAINTY,
    shape_min_reads: int = 5,
    xi_seed: int = XI_SEED,
    conf_level: float = CONF_LEVEL,
) -> None:
    """Wire the whole T2* within-calibration graph onto `job`: the carve, then the panel.

    Collapses the family: two job files that were
    byte-identical once the date and the run duration were normalised. What is left in each
    job file is its parameter row, which is the part a reader should actually diff.

    Deliberately does NOT call `configure_ramsey_job`: the carve comes from
    `wire_t2star_carve`, which adds no interpolate node.

    `shape_min_reads` and `xi_seed` are passed as step KWARGS rather than captured, so they
    reach the provenance label - the `allan` pattern, not the `filter` pattern. `runner`
    builds that label from `node.kwargs`, so an argument left to its default would be
    invisible to provenance; they are therefore passed explicitly below even where the value
    equals the default.
    """
    ladder = T2STAR_THRESHOLDS if thresholds is None else thresholds
    result, window_node = wire_t2star_carve(
        job,
        dataset=dataset,
        thresholds=ladder,
        gap_mult=gap_mult,
        k=k,
        use_uncertainty=use_uncertainty,
    )
    _tables, km = wire_kaplan_meier(job, window_node, conf_level=conf_level)
    panel = job.step(
        _t2star_panel_data,
        result,
        window_node,
        km,
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
