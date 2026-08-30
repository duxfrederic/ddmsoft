"""Structured, GUI-independent fitting for DDM matrices."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from math import factorial

import numpy as np
from scipy.optimize import minimize

from .models import DDMData, FitRange, FitRequest, FitResult

ModelFunction = Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray]
ProgressCallback = Callable[[int, int], None]


@dataclass(frozen=True)
class ParameterDefinition:
    """Stable parameter identity and its UI-independent default value."""

    identifier: str
    display_name: str
    default: float | None


@dataclass(frozen=True)
class ModelDefinition:
    """Complete definition of one supported relaxation model."""

    identifier: str
    display_name: str
    function: ModelFunction
    parameters: tuple[ParameterDefinition, ...]

    @property
    def defaults(self) -> tuple[float | None, ...]:
        return tuple(parameter.default for parameter in self.parameters)

    @property
    def parameter_names(self) -> tuple[str, ...]:
        return tuple(parameter.identifier for parameter in self.parameters)

    @property
    def physical_parameters(self) -> tuple[ParameterDefinition, ...]:
        return self.parameters[:-2]


class FitCancelled(ValueError):
    """A cooperative cancellation request stopped fitting before completion."""


def _tau(q_values: np.ndarray | float, times: np.ndarray) -> np.ndarray:
    q = np.asarray(q_values)
    times = np.asarray(times)
    if q.ndim == 1 and times.ndim == 1:
        return times[:, None] * q[None, :] ** 2
    return times * q**2


def _q_times(q_values: np.ndarray | float, times: np.ndarray) -> np.ndarray:
    q = np.asarray(q_values)
    times = np.asarray(times)
    if q.ndim == 1 and times.ndim == 1:
        return times[:, None] * q[None, :]
    return times * q


def _cumulant(
    order: int, parameters: np.ndarray, q_values: np.ndarray, times: np.ndarray
) -> np.ndarray:
    tau = _tau(q_values, times)
    result = np.ones_like(tau, dtype=float)
    for index, cumulant in enumerate(parameters[1:order], start=2):
        result += (-1) ** index * cumulant * tau**index / factorial(index)
    return result * np.exp(-parameters[0] * tau)


def _single_exponential(
    parameters: np.ndarray, q_values: np.ndarray, times: np.ndarray
) -> np.ndarray:
    return np.exp(-parameters[0] * _tau(q_values, times))


def _stretch(parameters: np.ndarray, q_values: np.ndarray, times: np.ndarray) -> np.ndarray:
    return np.exp(-(parameters[0] * _tau(q_values, times)) ** parameters[1])


def _double_exponential(
    parameters: np.ndarray, q_values: np.ndarray, times: np.ndarray
) -> np.ndarray:
    tau = _tau(q_values, times)
    diffusion_1, diffusion_2, beta, weight = parameters
    return weight * np.exp(-diffusion_1 * tau) + (1 - weight) * np.exp(
        -(diffusion_2 * tau) ** beta
    )


def _exponential_flow(
    parameters: np.ndarray, q_values: np.ndarray, times: np.ndarray
) -> np.ndarray:
    diffusion, flow = parameters
    return np.exp(-diffusion * _tau(q_values, times)) * np.cos(
        _q_times(q_values, times) * flow
    )


def _stretched_flow(
    parameters: np.ndarray, q_values: np.ndarray, times: np.ndarray
) -> np.ndarray:
    diffusion, flow, beta = parameters
    return np.exp(-(diffusion * _tau(q_values, times)) ** beta) * np.cos(
        _q_times(q_values, times) * flow
    )


def _double_exponential_flow(
    parameters: np.ndarray, q_values: np.ndarray, times: np.ndarray
) -> np.ndarray:
    flow = parameters[4]
    return _double_exponential(parameters[:4], q_values, times) * np.cos(
        _q_times(q_values, times) * flow
    )


def _parameter(
    identifier: str, display_name: str, default: float | None
) -> ParameterDefinition:
    return ParameterDefinition(identifier, display_name, default)


_AMPLITUDE = _parameter("amplitude", "A", None)
_BACKGROUND = _parameter("background", "B", None)


def _definition(
    identifier: str,
    display_name: str,
    function: ModelFunction,
    physical: Sequence[ParameterDefinition],
) -> ModelDefinition:
    return ModelDefinition(identifier, display_name, function, (*physical, _AMPLITUDE, _BACKGROUND))


MODEL_REGISTRY: dict[str, ModelDefinition] = {
    "cumulant_1": _definition(
        "cumulant_1",
        "Single exponential decay",
        _single_exponential,
        (_parameter("diffusion", "Diffusion coefficient", 2e-12),),
    ),
    "cumulant_2": _definition(
        "cumulant_2",
        "Cumulants up to second order",
        lambda parameters, q, times: _cumulant(2, parameters, q, times),
        (
            _parameter("diffusion", "Diffusion coefficient", 2e-12),
            _parameter("cumulant_2", "2nd cumulant", 0.0),
        ),
    ),
    "cumulant_3": _definition(
        "cumulant_3",
        "Cumulants up to third order",
        lambda parameters, q, times: _cumulant(3, parameters, q, times),
        (
            _parameter("diffusion", "Diffusion coefficient", 2e-12),
            _parameter("cumulant_2", "2nd cumulant", 0.0),
            _parameter("cumulant_3", "3rd cumulant", 0.0),
        ),
    ),
    "stretch": _definition(
        "stretch",
        "Stretched exponential",
        _stretch,
        (
            _parameter("diffusion", "Diffusion coefficient", 2e-12),
            _parameter("stretch", "Stretch parameter", 1.0),
        ),
    ),
    "dblexp_2ndstretched": _definition(
        "dblexp_2ndstretched",
        "Double exponential, second stretched",
        _double_exponential,
        (
            _parameter("diffusion_1", "Diffusion coefficient 1", 2e-12),
            _parameter("diffusion_2", "Diffusion coefficient 2", 2e-13),
            _parameter("stretch", "Stretch coefficient (2)", 1.0),
            _parameter("weight", "Weighting parameter", 1.0),
        ),
    ),
    "expcos": _definition(
        "expcos",
        "Exponential and flow",
        _exponential_flow,
        (
            _parameter("diffusion", "Diffusion coefficient", 2e-12),
            _parameter("flow", "Effective flow speed", 0.0),
        ),
    ),
    "expcosstretch": _definition(
        "expcosstretch",
        "Stretched exponential and flow",
        _stretched_flow,
        (
            _parameter("diffusion", "Diffusion coefficient", 2e-12),
            _parameter("flow", "Effective flow speed", 0.0),
            _parameter("stretch", "Stretch parameter", 1.0),
        ),
    ),
    "dblexpcosstretch": _definition(
        "dblexpcosstretch",
        "Double exponential with flow",
        _double_exponential_flow,
        (
            _parameter("diffusion_1", "Diffusion coefficient 1", 2e-12),
            _parameter("diffusion_2", "Diffusion coefficient 2", 2e-13),
            _parameter("stretch", "Stretch coefficient (2)", 1.0),
            _parameter("weight", "Weighting parameter", 1.0),
            _parameter("flow", "Effective flow speed", 0.0),
        ),
    ),
}


def get_model(model_id: str) -> ModelDefinition:
    """Return a registered model by stable identifier."""
    try:
        return MODEL_REGISTRY[model_id]
    except KeyError as error:
        raise ValueError(f"unknown fit model: {model_id}") from error


def default_fit_request(model_id: str, fit_range: FitRange) -> FitRequest:
    """Build a request from registry defaults without sharing mutable state."""
    model = get_model(model_id)
    return FitRequest(model_id, model.defaults, (False,) * len(model.parameters), fit_range)


def estimate_amplitude_background(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Estimate A and B for every q column without changing the matrix."""
    values = np.asarray(matrix, dtype=float)
    if values.ndim != 2 or values.shape[0] < 3:
        raise ValueError("at least three time rows are required to estimate A and B")
    background = values[0].copy()
    amplitude = values[-3:].mean(axis=0) - background
    return amplitude, background


def evaluate_model(
    model_id: str,
    parameters: Sequence[float],
    q_values: Iterable[float] | float,
    lag_times: Iterable[float],
) -> np.ndarray:
    """Evaluate a registered correlation model from physical parameters."""
    model = get_model(model_id)
    values = np.asarray(tuple(parameters), dtype=float)
    expected = len(model.physical_parameters)
    if values.size != expected:
        raise ValueError(f"{model_id} expects {expected} physical parameters")
    q = np.asarray(q_values, dtype=float)
    times = np.asarray(tuple(lag_times), dtype=float)
    return np.asarray(model.function(values, q, times), dtype=float)


def _resolved_initial_values(model: ModelDefinition, request: FitRequest) -> np.ndarray:
    if len(request.initial_values) != len(model.parameters):
        raise ValueError(
            f"{request.model_id} expects {len(model.parameters)} initial values, "
            f"got {len(request.initial_values)}"
        )
    resolved: list[float | None] = []
    for index, (definition, value) in enumerate(zip(model.parameters, request.initial_values)):
        if value is None:
            value = definition.default if index < len(model.parameters) - 2 else 0.0
        resolved.append(value)
    values = np.asarray(resolved, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("fit initial values must resolve to finite numbers")
    return values


def _fit_one_q(
    model: ModelDefinition,
    q_value: float,
    times: np.ndarray,
    observations: np.ndarray,
    initial: np.ndarray,
    fixed: np.ndarray,
    scales: np.ndarray,
) -> tuple[np.ndarray, bool, str]:
    denominator = np.median(times) + times
    denominator = np.where(denominator > 0, denominator, 1.0)

    def objective(candidate: np.ndarray) -> float:
        values = np.asarray(candidate, dtype=float) * scales
        values[fixed] = initial[fixed]
        with np.errstate(all="ignore"):
            correlation = model.function(values[:-2], np.asarray(q_value), times)
        predicted = values[-2] * (1.0 - correlation) + values[-1]
        residual = (observations - predicted) ** 2 / denominator
        if not np.all(np.isfinite(residual)):
            return 1e300
        return float(np.sum(residual))

    try:
        optimized = minimize(
            objective,
            initial / scales,
            method="Nelder-Mead",
            options={"maxiter": 10000, "maxfev": 10000},
        )
    except (ArithmeticError, FloatingPointError, ValueError) as error:
        return initial.copy(), False, str(error)
    values = np.asarray(optimized.x, dtype=float) * scales
    if values.shape != initial.shape or not np.all(np.isfinite(values)):
        return initial.copy(), False, "optimizer returned non-finite parameters"
    values[fixed] = initial[fixed]
    message = "" if optimized.success else str(optimized.message)
    return values, bool(optimized.success), message


def _parameter_scales(model: ModelDefinition) -> np.ndarray:
    scales = []
    for parameter in model.parameters:
        if parameter.identifier.startswith("diffusion"):
            scales.append(1e-12)
        elif parameter.identifier.startswith("cumulant"):
            scales.append(1e-26)
        elif parameter.identifier == "flow":
            scales.append(1e-7)
        else:
            scales.append(1.0)
    return np.asarray(scales)


def fit_ddm(
    data: DDMData,
    request: FitRequest,
    *,
    progress: ProgressCallback | None = None,
    cancel: Callable[[], bool] | object | None = None,
) -> FitResult:
    """Fit one inclusive q/time selection and return per-q diagnostics."""
    if not isinstance(data, DDMData):
        raise TypeError("data must be a DDMData instance")
    if not isinstance(request, FitRequest):
        raise TypeError("request must be a FitRequest instance")
    model = get_model(request.model_id)
    scales = _parameter_scales(model)
    if len(request.fixed_flags) != len(model.parameters):
        raise ValueError(
            f"{request.model_id} expects {len(model.parameters)} fixed flags, "
            f"got {len(request.fixed_flags)}"
        )
    q_slice, time_slice = request.fit_range.to_slices()
    if request.fit_range.q_max >= data.q_values.size:
        raise ValueError("fit q range exceeds available q values")
    if request.fit_range.time_max >= data.lag_times.size:
        raise ValueError("fit time range exceeds available lag times")
    if any(not isinstance(value, bool) for value in request.fixed_flags):
        raise ValueError("fixed_flags must contain booleans")

    q_values = data.q_values[q_slice]
    times = data.lag_times[time_slice]
    observations = data.matrix[time_slice, q_slice]
    amplitude_estimate, background_estimate = estimate_amplitude_background(data.matrix)
    initial_values = _resolved_initial_values(model, request)
    fixed = np.asarray(request.fixed_flags, dtype=bool)
    amplitudes: list[float] = []
    backgrounds: list[float] = []
    physical_parameters = [[] for _ in model.physical_parameters]
    correlations: list[np.ndarray] = []
    fitted_matrices: list[np.ndarray] = []
    statuses: list[bool] = []
    messages: list[str] = []

    for index, q_value in enumerate(q_values):
        if _cancel_requested(cancel):
            raise FitCancelled("fit cancelled before completing all q values")
        initial = initial_values.copy()
        if request.initial_values[-2] is None:
            initial[-2] = amplitude_estimate[request.fit_range.q_min + index]
        if request.initial_values[-1] is None:
            initial[-1] = background_estimate[request.fit_range.q_min + index]
        fitted, success, message = _fit_one_q(
            model, q_value, times, observations[:, index], initial, fixed, scales
        )
        with np.errstate(all="ignore"):
            correlation = np.asarray(model.function(fitted[:-2], q_value, times), dtype=float)
        fitted_matrix = fitted[-2] * (1.0 - correlation) + fitted[-1]
        amplitudes.append(float(fitted[-2]))
        backgrounds.append(float(fitted[-1]))
        for parameter, values in zip(fitted[:-2], physical_parameters):
            values.append(float(parameter))
        correlations.append(correlation)
        fitted_matrices.append(fitted_matrix)
        statuses.append(success)
        messages.append(message)
        if progress is not None:
            progress(index + 1, q_values.size)

    return FitResult(
        model_id=model.identifier,
        q_values=q_values,
        amplitude=np.asarray(amplitudes),
        noise=np.asarray(backgrounds),
        model_parameters=tuple(np.asarray(values) for values in physical_parameters),
        correlation=np.column_stack(correlations),
        fitted_matrix=np.column_stack(fitted_matrices),
        convergence_status=tuple(statuses),
        messages=tuple(messages),
    )


def _cancel_requested(cancel: Callable[[], bool] | object | None) -> bool:
    if cancel is None:
        return False
    if callable(cancel):
        return bool(cancel())
    is_set = getattr(cancel, "is_set", None)
    if callable(is_set):
        return bool(is_set())
    raise TypeError("cancel must be callable or provide an is_set() method")


__all__ = [
    "MODEL_REGISTRY",
    "FitCancelled",
    "ModelDefinition",
    "ParameterDefinition",
    "default_fit_request",
    "estimate_amplitude_background",
    "evaluate_model",
    "fit_ddm",
    "get_model",
]
