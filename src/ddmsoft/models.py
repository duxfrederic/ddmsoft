"""Typed, GUI-independent records shared by the backend layers."""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral, Real
from pathlib import Path

import numpy as np


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or value < 0:
        raise ValueError(f"{name} must be a non-negative integer")
    return int(value)


@dataclass(frozen=True)
class VideoMetadata:
    """Acquisition values needed to interpret a video."""

    path: Path
    frame_rate: float
    pixel_size: float
    temperature: float | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "path", Path(self.path))
        for value, name in ((self.frame_rate, "frame_rate"), (self.pixel_size, "pixel_size")):
            if not isinstance(value, Real) or not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be a finite positive number")
        if self.temperature is not None and (
            not isinstance(self.temperature, Real) or not np.isfinite(self.temperature)
        ):
            raise ValueError("temperature must be finite when provided")


@dataclass(frozen=True)
class DDMData:
    """A DDM matrix and its physical axes."""

    matrix: np.ndarray
    lag_times: np.ndarray
    q_values: np.ndarray
    source: Path | None = None

    def __post_init__(self) -> None:
        matrix = np.array(self.matrix, copy=True)
        lag_times = np.array(self.lag_times, dtype=float, copy=True)
        q_values = np.array(self.q_values, dtype=float, copy=True)
        if matrix.ndim != 2:
            raise ValueError("matrix must be two-dimensional")
        # Legacy directional matrices use NaN for bins with no pixels.
        if np.any(np.isinf(matrix)):
            raise ValueError("matrix must not contain infinite values")
        if lag_times.ndim != 1 or lag_times.size != matrix.shape[0]:
            raise ValueError("lag_times must have one value per matrix row")
        if q_values.ndim != 1 or q_values.size != matrix.shape[1]:
            raise ValueError("q_values must have one value per matrix column")
        if not np.all(np.isfinite(lag_times)) or not np.all(np.isfinite(q_values)):
            raise ValueError("lag_times and q_values must contain finite values")
        object.__setattr__(self, "matrix", matrix)
        object.__setattr__(self, "lag_times", lag_times)
        object.__setattr__(self, "q_values", q_values)
        if self.source is not None:
            object.__setattr__(self, "source", Path(self.source))

    @property
    def delta_times(self) -> np.ndarray:
        """Alias using the terminology in the legacy file suffix."""
        return self.lag_times

    @property
    def dts(self) -> np.ndarray:
        return self.lag_times

    @property
    def qs(self) -> np.ndarray:
        return self.q_values


@dataclass(frozen=True)
class FitRange:
    """Inclusive q and lag-time index bounds."""

    q_min: int
    q_max: int
    time_min: int
    time_max: int

    def __post_init__(self) -> None:
        values = (
            (self.q_min, "q_min"),
            (self.q_max, "q_max"),
            (self.time_min, "time_min"),
            (self.time_max, "time_max"),
        )
        for value, name in values:
            _index(value, name)
        if self.q_min > self.q_max:
            raise ValueError("q_min must not exceed q_max")
        if self.time_min > self.time_max:
            raise ValueError("time_min must not exceed time_max")

    def to_slices(self) -> tuple[slice, slice]:
        """Return inclusive bounds as `(q_slice, time_slice)` NumPy slices."""
        return slice(self.q_min, self.q_max + 1), slice(self.time_min, self.time_max + 1)

    @property
    def qmin(self) -> int:
        return self.q_min

    @property
    def qmax(self) -> int:
        return self.q_max

    @property
    def timemin(self) -> int:
        return self.time_min

    @property
    def timemax(self) -> int:
        return self.time_max


@dataclass(frozen=True)
class FitRequest:
    """Inputs for one model fit."""

    model_id: str
    initial_values: tuple[float | None, ...]
    fixed_flags: tuple[bool, ...]
    fit_range: FitRange

    def __post_init__(self) -> None:
        if not self.model_id or not isinstance(self.model_id, str):
            raise ValueError("model_id must be a non-empty string")
        values = tuple(self.initial_values)
        fixed = tuple(bool(value) for value in self.fixed_flags)
        if len(values) != len(fixed):
            raise ValueError("initial_values and fixed_flags must have equal lengths")
        for value in values:
            if value is not None and (not isinstance(value, Real) or not np.isfinite(value)):
                raise ValueError("initial values must be finite numbers or None")
        if not isinstance(self.fit_range, FitRange):
            raise TypeError("fit_range must be a FitRange")
        object.__setattr__(self, "initial_values", values)
        object.__setattr__(self, "fixed_flags", fixed)

    @property
    def range(self) -> FitRange:
        return self.fit_range


@dataclass(frozen=True)
class FitResult:
    """Structured output from fitting one matrix."""

    model_id: str
    q_values: np.ndarray
    amplitude: np.ndarray
    noise: np.ndarray
    model_parameters: tuple[np.ndarray, ...]
    correlation: np.ndarray
    fitted_matrix: np.ndarray
    convergence_status: tuple[bool, ...]
    messages: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        q_values = np.array(self.q_values, dtype=float, copy=True)
        amplitude = np.array(self.amplitude, dtype=float, copy=True)
        noise = np.array(self.noise, dtype=float, copy=True)
        correlation = np.array(self.correlation, dtype=float, copy=True)
        fitted_matrix = np.array(self.fitted_matrix, dtype=float, copy=True)
        status = tuple(bool(value) for value in self.convergence_status)
        messages = tuple(str(value) for value in self.messages)
        if q_values.ndim != 1:
            raise ValueError("q_values must be one-dimensional")
        q_count = q_values.size
        if any(array.ndim != 1 or array.size != q_count for array in (amplitude, noise)):
            raise ValueError("amplitude and noise must have one value per q value")
        if correlation.ndim != 2 or fitted_matrix.shape != correlation.shape:
            raise ValueError("correlation and fitted_matrix must be equal-shaped matrices")
        if correlation.shape[1] != q_count:
            raise ValueError("fit matrices must have one column per q value")
        parameters = tuple(np.array(parameter, dtype=float, copy=True) for parameter in self.model_parameters)
        if any(parameter.ndim != 1 or parameter.size != q_count for parameter in parameters):
            raise ValueError("every model parameter must have one value per q value")
        if len(status) != q_count:
            raise ValueError("convergence_status must have one value per q value")
        if messages and len(messages) != q_count:
            raise ValueError("messages must have one value per q value when provided")
        object.__setattr__(self, "q_values", q_values)
        object.__setattr__(self, "amplitude", amplitude)
        object.__setattr__(self, "noise", noise)
        object.__setattr__(self, "model_parameters", parameters)
        object.__setattr__(self, "correlation", correlation)
        object.__setattr__(self, "fitted_matrix", fitted_matrix)
        object.__setattr__(self, "convergence_status", status)
        object.__setattr__(self, "messages", messages)

    @property
    def fitted_q_values(self) -> np.ndarray:
        return self.q_values

    @property
    def converged(self) -> tuple[bool, ...]:
        return self.convergence_status
