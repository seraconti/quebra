from __future__ import annotations

import pandera.pandas as pa
from pandera.typing import Series

from quebra.schemas.base import BaseQubitSchema


class track912Schema(BaseQubitSchema):
    # amp: fit amplitude, unitless.
    amp: Series[float] = pa.Field(alias="amp")
    # amp_error: amplitude standard error, unitless.
    amp_error: Series[float] = pa.Field(alias="amp error")
    # offset: fit offset, unitless.
    offset: Series[float] = pa.Field(alias="offset")
    # offset_error: offset standard error, unitless.
    offset_error: Series[float] = pa.Field(alias="offset error")
    # T2star_s: Ramsey coherence time in seconds.
    T2star_s: Series[float] = pa.Field(alias="T2star")
    # T2star_error_s: uncertainty on Ramsey coherence time in seconds.
    T2star_error_s: Series[float] = pa.Field(alias="T2star error")
    # frequency_hz: fitted Ramsey fringe frequency in hertz, of order MHz; not the
    # absolute qubit frequency.
    frequency_hz: Series[float] = pa.Field(alias="frequency")
    # frequency_error_hz: uncertainty on fitted frequency in hertz.
    frequency_error_hz: Series[float] = pa.Field(alias="frequency error")
    # phase_rad: fitted phase in radians.
    phase_rad: Series[float] = pa.Field(alias="phase")
    # phase_error_rad: uncertainty on fitted phase in radians.
    phase_error_rad: Series[float] = pa.Field(alias="phase_error")
    # normalised_chi_square: reduced chi-square of the fit.
    normalised_chi_square: Series[float] = pa.Field(alias="normalised chi-square")

    @classmethod
    def validate(cls, df, *args, **kwargs):
        """Perform device-aware cleaning before running the strict pandera validation.

        This method will drop any rows that lack required non-null fit fields.
        For devices that legitimately lack a Ramsey `frequency` (e.g. device families
        without a Ramsey frequency), the presence/absence should be handled by the
        downstream pipeline (job configuration). We do not relax schema constraints;
        instead we remove malformed rows so the strict schema can be enforced.
        """
        required = [
            "timestamp",
            "amp",
            "amp error",
            "offset",
            "offset error",
            "T2star",
            "T2star error",
            "phase",
            "phase_error",
            "normalised chi-square",
        ]

        df = df.copy()

        # Allow caller to provide dataset context (either a Dataset object or a dict)
        dataset = kwargs.pop("dataset", None)
        if dataset is not None:
            # Extract qubit and device from provided dataset context
            qubit = (
                getattr(dataset, "qubit", None)
                if not isinstance(dataset, dict)
                else dataset.get("qubit")
            )
            device = (
                getattr(dataset, "device", None)
                if not isinstance(dataset, dict)
                else dataset.get("device")
            )
            if qubit is not None and "qubit_id" not in df.columns:
                df["qubit_id"] = int(qubit)
            if device is not None and "device" not in df.columns:
                df["device"] = device

        # Always drop rows missing the core fit fields listed above.
        present_required = [c for c in required if c in df.columns]
        if present_required:
            df = df.dropna(subset=present_required)

        # After cleaning, defer to the base pandera validation to ensure types/aliases.
        return super().validate(df, *args, **kwargs)
