"""Qt-free CONTIN-style regularized distribution fitting."""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from .science import stokes_einstein_radius


class CONTINError(ValueError):
    """Base class for invalid CONTIN inputs."""


class CONTINCancelled(CONTINError):
    """A cooperative cancellation request stopped a CONTIN run."""


@dataclass(frozen=True)
class CONTINResult:
    """All candidate alpha solutions and the selected minimum-residual result."""

    gamma_range: np.ndarray
    alphas: np.ndarray
    tau: np.ndarray
    ddmdata: np.ndarray
    alpha_residuals: np.ndarray
    alpha_g: np.ndarray
    alpha_ddmfit: np.ndarray
    alpha_amplitude: np.ndarray
    alpha_noise: np.ndarray
    selected_index: int
    sizes: np.ndarray | None = None

    def __post_init__(self) -> None:
        gamma = np.array(self.gamma_range, dtype=float, copy=True)
        alphas = np.array(self.alphas, dtype=float, copy=True)
        tau = np.array(self.tau, dtype=float, copy=True)
        data = np.array(self.ddmdata, dtype=float, copy=True)
        residuals = np.array(self.alpha_residuals, dtype=float, copy=True)
        distributions = np.array(self.alpha_g, dtype=float, copy=True)
        fits = np.array(self.alpha_ddmfit, dtype=float, copy=True)
        amplitudes = np.array(self.alpha_amplitude, dtype=float, copy=True)
        noises = np.array(self.alpha_noise, dtype=float, copy=True)
        if (
            gamma.ndim != 1
            or gamma.size < 3
            or not np.all(np.isfinite(gamma))
            or np.any(gamma <= 0)
            or np.any(np.diff(gamma) <= 0)
        ):
            raise CONTINError("gamma_range must contain at least three increasing positive values")
        if (
            alphas.ndim != 1
            or alphas.size < 1
            or np.any(alphas < 0)
            or not np.all(np.isfinite(alphas))
            or (alphas.size > 1 and np.any(np.diff(alphas) <= 0))
        ):
            raise CONTINError(
                "alphas must be a strictly increasing finite non-negative sequence"
            )
        if (
            tau.ndim != 1
            or tau.size < 2
            or np.any(tau < 0)
            or not np.all(np.isfinite(tau))
        ):
            raise CONTINError("tau must contain at least two finite non-negative values")
        if data.shape != tau.shape:
            raise CONTINError("ddmdata must have one value per tau")
        alpha_count = alphas.size
        if (
            residuals.shape != (alpha_count,)
            or amplitudes.shape != (alpha_count,)
            or noises.shape != (alpha_count,)
        ):
            raise CONTINError("CONTIN scalar candidate arrays must have one value per alpha")
        if distributions.shape != (alpha_count, gamma.size):
            raise CONTINError("alpha_g must have shape (alpha_count, gamma_count)")
        if fits.shape != (alpha_count, tau.size):
            raise CONTINError("alpha_ddmfit must have shape (alpha_count, tau_count)")
        if not all(
            np.all(np.isfinite(array))
            for array in (residuals, distributions, fits, amplitudes, noises)
        ):
            raise CONTINError("CONTIN candidate results must be finite")
        if isinstance(self.selected_index, bool) or not isinstance(
            self.selected_index, (int, np.integer)
        ):
            raise CONTINError("selected_index must be an alpha index")
        selected_index = int(self.selected_index)
        if not 0 <= selected_index < alpha_count:
            raise CONTINError("selected_index is outside the alpha candidates")
        sizes = None if self.sizes is None else np.array(self.sizes, dtype=float, copy=True)
        if sizes is not None and (
            sizes.shape != gamma.shape or not np.all(np.isfinite(sizes))
        ):
            raise CONTINError("sizes must have one finite value per gamma")
        for name, value in (
            ("gamma_range", gamma),
            ("alphas", alphas),
            ("tau", tau),
            ("ddmdata", data),
        ):
            object.__setattr__(self, name, value)
        for name, value in (
            ("alpha_residuals", residuals),
            ("alpha_g", distributions),
            ("alpha_ddmfit", fits),
            ("alpha_amplitude", amplitudes),
            ("alpha_noise", noises),
        ):
            object.__setattr__(self, name, value)
        object.__setattr__(self, "selected_index", selected_index)
        object.__setattr__(self, "sizes", sizes)

    @property
    def selection_method(self) -> str:
        return "minimum residual"

    @property
    def chosen_alpha(self) -> float:
        return float(self.alphas[self.selected_index])

    @property
    def g(self) -> np.ndarray:
        return self.alpha_g[self.selected_index]

    @property
    def ddmfit(self) -> np.ndarray:
        return self.alpha_ddmfit[self.selected_index]

    @property
    def amplitude(self) -> float:
        return float(self.alpha_amplitude[self.selected_index])

    @property
    def noise(self) -> float:
        return float(self.alpha_noise[self.selected_index])

    def with_particle_sizes(
        self, temperature_kelvin: float, viscosity_pa_s: float
    ) -> CONTINResult:
        """Return a copy with each diffusion rate converted to a radius in metres."""
        sizes = gamma_to_radius(self.gamma_range, temperature_kelvin, viscosity_pa_s)
        return replace(self, sizes=sizes)


def regularized_square_sum(
    candidate: np.ndarray,
    data: np.ndarray,
    convolution: np.ndarray,
    alpha: float,
    weights: np.ndarray,
) -> float:
    """Return the legacy CONTIN weighted residual plus curvature penalties."""
    amplitude, noise = candidate[:2]
    distribution = np.abs(candidate[2:])
    scale = amplitude + noise
    edge_penalty = 0.1 * scale * (abs(distribution[0]) + abs(distribution[-1]))
    fitted = amplitude * (1.0 - convolution @ distribution) + noise
    residual = np.sum(weights * (data - fitted) ** 2)
    curvature = scale * alpha**2 * np.sum(np.diff(distribution, 2) ** 2)
    score = residual + curvature + edge_penalty
    return float(score) if np.isfinite(score) else 1e300


def _as_alphas(alpha: float | Iterable[float]) -> np.ndarray:
    if isinstance(alpha, (str, bytes)):
        raise CONTINError("alpha must be a number or iterable of numbers")
    if np.isscalar(alpha):
        values = np.asarray([alpha], dtype=float)
    else:
        values = np.asarray(tuple(alpha), dtype=float)
    if (
        values.ndim != 1
        or values.size == 0
        or not np.all(np.isfinite(values))
        or np.any(values < 0)
        or (values.size > 1 and np.any(np.diff(values) <= 0))
    ):
        raise CONTINError(
            "alpha values must be a strictly increasing finite non-negative sequence"
        )
    return values


def _validate_ranges(
    gamma_range: Iterable[float], alpha: float | Iterable[float]
) -> tuple[np.ndarray, np.ndarray]:
    gamma = np.asarray(tuple(gamma_range), dtype=float)
    if gamma.ndim != 1 or gamma.size < 3 or not np.all(np.isfinite(gamma)):
        raise CONTINError("gamma_range must contain at least three finite values")
    if np.any(gamma <= 0) or np.any(np.diff(gamma) <= 0):
        raise CONTINError("gamma_range must be strictly increasing and positive")
    return gamma, _as_alphas(alpha)


def contin_ranges(
    gamma_min: float,
    gamma_max: float,
    gamma_count: int,
    alpha_min: float,
    alpha_max: float,
    alpha_count: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Build and validate the linear ranges used by the legacy CONTIN dialog."""
    try:
        gamma_min = float(gamma_min)
        gamma_max = float(gamma_max)
        alpha_min = float(alpha_min)
    except (TypeError, ValueError) as error:
        raise CONTINError("CONTIN range endpoints must be numeric") from error
    try:
        alpha_max = float(alpha_max)
    except (TypeError, ValueError) as error:
        raise CONTINError("CONTIN range endpoints must be numeric") from error
    if (
        not all(np.isfinite(value) for value in (gamma_min, gamma_max, alpha_min, alpha_max))
        or gamma_min <= 0
        or gamma_max <= gamma_min
        or alpha_min < 0
        or alpha_max < alpha_min
    ):
        raise CONTINError("CONTIN ranges must have finite ordered endpoints")
    if (
        isinstance(gamma_count, bool)
        or not isinstance(gamma_count, (int, np.integer))
        or gamma_count < 3
        or isinstance(alpha_count, bool)
        or not isinstance(alpha_count, (int, np.integer))
        or alpha_count < 1
    ):
        raise CONTINError("gamma_count must be at least 3 and alpha_count must be positive")
    if alpha_count > 1 and alpha_max == alpha_min:
        raise CONTINError("alpha range must have distinct endpoints for multiple candidates")
    return (
        np.linspace(gamma_min, gamma_max, int(gamma_count)),
        np.linspace(alpha_min, alpha_max, int(alpha_count)),
    )


def run_contin(
    tau: Iterable[float],
    ddmdata: Iterable[float],
    gamma_range: Iterable[float],
    *,
    alpha: float | Iterable[float] = 0.1,
    g0: Iterable[float] | None = None,
    weights: Iterable[float] | None = None,
    maxiter: int = 10,
    tol: float = 1e-3,
    progress: Callable[[int, int], None] | None = None,
    cancel: Callable[[], bool] | object | None = None,
) -> CONTINResult:
    """Run a regularized alpha scan and return all candidates explicitly."""
    tau_values = np.asarray(tuple(tau), dtype=float)
    data = np.asarray(tuple(ddmdata), dtype=float)
    if tau_values.ndim != 1 or data.shape != tau_values.shape or tau_values.size < 2:
        raise CONTINError("tau and ddmdata must be equal-length one-dimensional arrays")
    if (
        np.any(tau_values < 0)
        or not np.all(np.isfinite(tau_values))
        or not np.all(np.isfinite(data))
    ):
        raise CONTINError(
            "tau and ddmdata must contain finite non-negative tau values"
        )
    gamma, alphas = _validate_ranges(gamma_range, alpha)
    if (
        isinstance(maxiter, bool)
        or not isinstance(maxiter, (int, np.integer))
        or maxiter < 1
    ):
        raise CONTINError("maxiter must be a positive integer")
    if not np.isfinite(tol) or tol < 0:
        raise CONTINError("tol must be a finite non-negative number")
    if weights is None:
        weight_values = 1.0 / np.sqrt(np.arange(1, data.size + 1, dtype=float))
    else:
        weight_values = np.asarray(tuple(weights), dtype=float)
        if (
            weight_values.shape != data.shape
            or not np.all(np.isfinite(weight_values))
            or np.any(weight_values <= 0)
        ):
            raise CONTINError(
                "weights must be finite and positive with one value per data point"
            )
    if data.size < 3:
        raise CONTINError("at least three data points are required")
    if g0 is None:
        initial_distribution = np.ones(gamma.size, dtype=float)
        initial_distribution[[0, -1]] = 0.0
    else:
        initial_distribution = np.asarray(tuple(g0), dtype=float)
        if (
            initial_distribution.shape != gamma.shape
            or not np.all(np.isfinite(initial_distribution))
        ):
            raise CONTINError("g0 must contain one finite value per gamma")
    if np.sum(np.abs(initial_distribution)) == 0:
        raise CONTINError("g0 must not be all zero")
    initial_distribution = np.abs(initial_distribution)
    initial_distribution /= np.sum(initial_distribution)
    convolution = np.exp(-tau_values[:, None] * gamma[None, :])
    amplitude = np.mean(data[-2:]) - data[0]
    noise = data[0]
    residuals: list[float] = []
    distributions: list[np.ndarray] = []
    fits: list[np.ndarray] = []
    amplitudes: list[float] = []
    noises: list[float] = []

    for index, alpha_value in enumerate(alphas):
        if _cancel_requested(cancel):
            raise CONTINCancelled("CONTIN cancelled before the next alpha")
        candidate = np.concatenate(([amplitude, noise], initial_distribution))
        previous_residual: float | None = None
        residual = np.inf
        for _ in range(int(maxiter)):
            if _cancel_requested(cancel):
                raise CONTINCancelled("CONTIN cancelled during optimizer iterations")

            def check_cancel(candidate_values: np.ndarray) -> None:
                if _cancel_requested(cancel):
                    raise CONTINCancelled("CONTIN cancelled during an optimizer iteration")

            solution = minimize(
                regularized_square_sum,
                candidate,
                args=(data, convolution, float(alpha_value), weight_values),
                method="Nelder-Mead",
                callback=check_cancel,
            )
            candidate = np.abs(np.asarray(solution.x, dtype=float))
            residual = float(solution.fun)
            if (
                previous_residual is not None
                and abs(previous_residual - residual) <= abs(residual) * tol
            ):
                break
            previous_residual = residual
        if not np.all(np.isfinite(candidate)) or not np.isfinite(residual):
            raise CONTINError(f"non-finite CONTIN result for alpha {alpha_value}")
        residuals.append(residual)
        distributions.append(candidate[2:].copy())
        fits.append(convolution @ candidate[2:])
        amplitudes.append(float(candidate[0]))
        noises.append(float(candidate[1]))
        if progress is not None:
            progress(index + 1, alphas.size)

    selected_index = int(np.argmin(residuals))
    return CONTINResult(
        gamma_range=gamma,
        alphas=alphas,
        tau=tau_values,
        ddmdata=data,
        alpha_residuals=np.asarray(residuals),
        alpha_g=np.asarray(distributions),
        alpha_ddmfit=np.asarray(fits),
        alpha_amplitude=np.asarray(amplitudes),
        alpha_noise=np.asarray(noises),
        selected_index=selected_index,
    )


def run_contin_scan(
    tau: Iterable[float],
    ddmdata: Iterable[float],
    gamma_min: float,
    gamma_max: float,
    gamma_count: int,
    alpha_min: float,
    alpha_max: float,
    alpha_count: int,
    **kwargs: object,
) -> CONTINResult:
    """Run CONTIN from validated min/max/count controls."""
    gamma, alphas = contin_ranges(
        gamma_min,
        gamma_max,
        gamma_count,
        alpha_min,
        alpha_max,
        alpha_count,
    )
    return run_contin(tau, ddmdata, gamma, alpha=alphas, **kwargs)


def gamma_to_radius(
    gamma: Iterable[float], temperature_kelvin: float, viscosity_pa_s: float
) -> np.ndarray:
    """Convert CONTIN diffusion rates to hydrodynamic radii in metres."""
    values = np.asarray(tuple(gamma), dtype=float)
    if values.ndim != 1 or not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise CONTINError("gamma must contain finite positive diffusion rates")
    try:
        return np.asarray(
            [
                stokes_einstein_radius(value, temperature_kelvin, viscosity_pa_s)
                for value in values
            ]
        )
    except ValueError as error:
        raise CONTINError(str(error)) from error


def export_contin(
    path: str | Path,
    result: CONTINResult,
    *,
    video: str | Path | None = None,
    q: float | None = None,
    all_alphas: bool = True,
) -> Path:
    """Export selected or all alpha candidates without any GUI dependency."""
    output = Path(path)
    indices = range(result.alphas.size) if all_alphas else (result.selected_index,)
    with output.open("w", encoding="utf-8", newline="\n") as save_file:
        save_file.write(f"Matrix:\t{video if video is not None else ''}\n")
        if q is not None:
            save_file.write(f"Wavenumber [1/m]:\t{q:.03e}\n")
        save_file.write(f"selection method:\t{result.selection_method}\n")
        save_file.write(f"selected alpha:\t{result.chosen_alpha:.03e}\n")
        for index in indices:
            save_file.write("###############################\n")
            save_file.write(
                f"alpha: {result.alphas[index]:.03e} with residual "
                f"{result.alpha_residuals[index]:.05e}\n"
            )
            save_file.write(f"amplitude:\t{result.alpha_amplitude[index]:.03e}\n")
            save_file.write(f"noise:\t{result.alpha_noise[index]:.03e}\n\n")
            save_file.write("Distribution of decay rates\n")
            if result.sizes is None:
                save_file.write("Gamma [m^2/s]\tintensity []\n")
                for gamma, value in zip(result.gamma_range, result.alpha_g[index]):
                    save_file.write(f"{gamma:.03e}\t{value:.03e}\n")
            else:
                save_file.write("Gamma [m^2/s]\tsize [nm]\tintensity []\n")
                for gamma, size, value in zip(
                    result.gamma_range, result.sizes, result.alpha_g[index]
                ):
                    save_file.write(f"{gamma:.03e}\t{size * 1e9:.03e}\t{value:.03e}\n")
            save_file.write("\ntau q^2 [s/m^2]\texp data\tCONTIN fit\n")
            intensity = (
                result.alpha_amplitude[index] * (1.0 - result.alpha_ddmfit[index])
                + result.alpha_noise[index]
            )
            for tau, measured, fitted in zip(result.tau, result.ddmdata, intensity):
                save_file.write(f"{tau:.03e}\t{measured:.03e}\t{fitted:.03e}\n")
    return output


def _cancel_requested(cancel: Callable[[], bool] | object | None) -> bool:
    if cancel is None:
        return False
    if callable(cancel):
        return bool(cancel())
    is_set = getattr(cancel, "is_set", None)
    if callable(is_set):
        return bool(is_set())
    raise TypeError("cancel must be callable or provide an is_set() method")


# The explicit function is the public replacement for the legacy generator.
CONTIN = run_contin


__all__ = [
    "CONTIN",
    "CONTINCancelled",
    "CONTINError",
    "CONTINResult",
    "contin_ranges",
    "export_contin",
    "gamma_to_radius",
    "regularized_square_sum",
    "run_contin",
    "run_contin_scan",
]
