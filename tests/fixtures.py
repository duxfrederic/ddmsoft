"""Independent deterministic fixtures for backend and scientific tests."""

from __future__ import annotations

import math
from collections.abc import Iterable
from pathlib import Path

import numpy as np

from ddmsoft.models import DDMData

FIT_MODEL_IDS = (
    "cumulant_1",
    "cumulant_2",
    "cumulant_3",
    "stretch",
    "dblexp_2ndstretched",
    "expcos",
    "expcosstretch",
    "dblexpcosstretch",
)

CUMULANT_REFERENCE_Q = 2.0e6
CUMULANT_PDI = 0.03
CUMULANT_THIRD_NORMALIZED = -0.01

LEGACY_SUFFIXES = ("_DDM_matrix.npy", "_deltaTs.npy", "_QS.npy")


def constant_frames(count: int = 5, shape: tuple[int, int] = (8, 8), value: float = 1.0) -> np.ndarray:
    """Return identical frames."""
    return np.full((count, *shape), value, dtype=float)


def changing_intensity_frames(
    count: int = 5, shape: tuple[int, int] = (8, 8)
) -> np.ndarray:
    """Return uniform frames with a deterministic linear intensity change."""
    levels = np.arange(count, dtype=float)[:, None, None]
    return np.broadcast_to(levels, (count, *shape)).copy()


def translated_sinusoidal_grating(
    count: int = 6, shape: tuple[int, int] = (16, 16), cycles: int = 3
) -> np.ndarray:
    """Return a grating translated by one pixel per frame."""
    _, x = np.indices(shape)
    frames = [np.sin(2.0 * np.pi * cycles * (x - shift) / shape[1]) for shift in range(count)]
    return np.asarray(frames, dtype=float)


def orthogonal_gratings(count: int = 6, shape: tuple[int, int] = (16, 16)) -> np.ndarray:
    """Return two orthogonal gratings whose phases change independently."""
    y, x = np.indices(shape)
    frames = [
        np.sin(2.0 * np.pi * (x - index) / shape[1])
        + 0.7 * np.sin(2.0 * np.pi * (y - 2 * index) / shape[0])
        for index in range(count)
    ]
    return np.asarray(frames, dtype=float)


def random_frames(
    count: int = 5, shape: tuple[int, int] = (8, 8), seed: int = 1234
) -> np.ndarray:
    """Return seeded standard-normal frames."""
    return np.random.default_rng(seed).standard_normal((count, *shape))


def direct_ddm_reference(frames: Iterable[np.ndarray], lags: Iterable[int] | None = None) -> np.ndarray:
    """Calculate full-plane DDM power independently of production averaging code."""
    stack = np.asarray(list(frames), dtype=float)
    if stack.ndim != 3 or stack.shape[0] < 2:
        raise ValueError("frames must contain at least two equally shaped 2D frames")
    if lags is None:
        lag_values = np.arange(1, stack.shape[0], dtype=int)
    else:
        lag_values = np.asarray(list(lags), dtype=int)
    if lag_values.ndim != 1 or np.any(lag_values < 1) or np.any(lag_values >= stack.shape[0]):
        raise ValueError("lags must be in [1, frame_count - 1]")
    transformed = np.fft.fft2(stack, axes=(1, 2))
    return np.asarray(
        [
            np.mean(np.abs(transformed[lag:] - transformed[:-lag]) ** 2, axis=0)
            for lag in lag_values
        ]
    )


def model_parameters(model_id: str) -> tuple[float, ...]:
    """Return stable, non-degenerate parameters in each model's physical units.

    The cumulant values use the reference decay rate at q=2e6 m^-1, with
    PDI=0.03 and normalized mu3=-0.01.
    """
    diffusion = 2.0e-12
    gamma = diffusion * CUMULANT_REFERENCE_Q**2
    values = {
        "cumulant_1": (diffusion,),
        "cumulant_2": (diffusion, CUMULANT_PDI * gamma**2),
        "cumulant_3": (
            diffusion,
            CUMULANT_PDI * gamma**2,
            CUMULANT_THIRD_NORMALIZED * gamma**3,
        ),
        "stretch": (2.0e-12, 0.82),
        "dblexp_2ndstretched": (2.0e-12, 6.0e-13, 0.78, 0.65),
        "expcos": (2.0e-12, 1.0e-7),
        "expcosstretch": (2.0e-12, 1.0e-7, 0.82),
        "dblexpcosstretch": (2.0e-12, 6.0e-13, 0.78, 0.65, 1.0e-7),
    }
    try:
        return values[model_id]
    except KeyError as error:
        raise ValueError(f"unknown fit model: {model_id}") from error


def _correlation(model_id: str, parameters: tuple[float, ...], q: np.ndarray, times: np.ndarray) -> np.ndarray:
    tau = times[:, None] * q[None, :] ** 2
    if model_id.startswith("cumulant_"):
        gamma = parameters[0] * q[None, :] ** 2
        exponent = -gamma * times[:, None]
        for order, cumulant in enumerate(parameters[1:], start=2):
            exponent += (-1) ** order * cumulant * times[:, None] ** order / math.factorial(order)
        return np.exp(exponent)
    if model_id == "stretch":
        return np.exp(-(parameters[0] * tau) ** parameters[1])
    if model_id in {"dblexp_2ndstretched", "dblexpcosstretch"}:
        diffusion_1, diffusion_2, beta, weight = parameters[:4]
        result = weight * np.exp(-diffusion_1 * tau) + (1 - weight) * np.exp(
            -(diffusion_2 * tau) ** beta
        )
        if model_id == "dblexpcosstretch":
            result *= np.cos(q[None, :] * parameters[4] * times[:, None])
        return result
    if model_id in {"expcos", "expcosstretch"}:
        diffusion, flow = parameters[:2]
        decay = np.exp(-diffusion * tau)
        if model_id == "expcosstretch":
            decay = np.exp(-(diffusion * tau) ** parameters[2])
        return decay * np.cos(q[None, :] * flow * times[:, None])
    raise ValueError(f"unknown fit model: {model_id}")


def _q_scaled_cumulant_correlation(
    model_id: str, q: np.ndarray, times: np.ndarray
) -> np.ndarray:
    """Generate cumulants with fixed PDI and normalized third cumulant."""
    diffusion = 2.0e-12
    gamma = diffusion * q[None, :] ** 2
    second = CUMULANT_PDI * gamma**2
    exponent = -gamma * times[:, None] + second * times[:, None] ** 2 / 2.0
    if model_id == "cumulant_3":
        third = CUMULANT_THIRD_NORMALIZED * gamma**3
        exponent -= third * times[:, None] ** 3 / 6.0
    return np.exp(exponent)


def generate_model_data(
    model_id: str,
    q_values: np.ndarray | None = None,
    lag_times: np.ndarray | None = None,
    noise: float = 0.0,
    seed: int = 5678,
) -> DDMData:
    """Generate a deterministic DDM matrix for each legacy fit model."""
    parameters = model_parameters(model_id)
    q = np.asarray(q_values if q_values is not None else np.linspace(1.0e6, 4.0e6, 4), dtype=float)
    times = np.asarray(
        lag_times if lag_times is not None else np.linspace(0.01, 0.20, 7), dtype=float
    )
    amplitude = np.linspace(1.0, 1.3, q.size)
    background = np.linspace(0.02, 0.05, q.size)
    correlation = (
        _q_scaled_cumulant_correlation(model_id, q, times)
        if model_id in {"cumulant_2", "cumulant_3"}
        else _correlation(model_id, parameters, q, times)
    )
    matrix = amplitude[None, :] * (1.0 - correlation) + background[None, :]
    if noise:
        matrix = matrix + np.random.default_rng(seed).normal(0.0, noise, matrix.shape)
    return DDMData(matrix=matrix, lag_times=times, q_values=q)


def write_legacy_matrix(directory: Path, stem: str, data: DDMData) -> tuple[Path, Path, Path]:
    """Write one legacy three-file matrix set and return its paths."""
    target = Path(directory) / "ddm_matrices"
    target.mkdir(parents=True, exist_ok=True)
    arrays = (data.matrix, data.lag_times, data.q_values)
    paths = tuple(target / f"{stem}{suffix}" for suffix in LEGACY_SUFFIXES)
    for path, array in zip(paths, arrays):
        np.save(path, array)
    return paths


def read_legacy_matrix(directory: Path, stem: str) -> DDMData:
    """Read one legacy three-file matrix set without a production I/O helper."""
    paths = tuple(Path(directory) / "ddm_matrices" / f"{stem}{suffix}" for suffix in LEGACY_SUFFIXES)
    if not all(path.is_file() for path in paths):
        missing = [str(path) for path in paths if not path.is_file()]
        raise FileNotFoundError("incomplete legacy matrix set: " + ", ".join(missing))
    return DDMData(matrix=np.load(paths[0]), lag_times=np.load(paths[1]), q_values=np.load(paths[2]))
