"""OpenCV video input and compatibility DDM computation.

The engine accepts any iterable yielding two-dimensional frames, so file input
and deterministic tests use the same computation path.  It intentionally keeps
the first implementation straightforward rather than optimizing FFT storage.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator
from numbers import Real
from pathlib import Path
from typing import Protocol

import cv2
import numpy as np
from scipy.signal.windows import tukey

from .io import save_matrix_set
from .models import DDMData


class ProgressCallback(Protocol):
    def __call__(self, stage: str, completed: int, total: int) -> None: ...


class FrameIterator(Protocol):
    def __iter__(self) -> Iterator[np.ndarray]: ...


class DDMEngineError(ValueError):
    """Base class for invalid DDM input or computation errors."""


class VideoReadError(DDMEngineError):
    """A video cannot be opened or decoded into consistent frames."""


class ComputationCancelled(DDMEngineError):
    """A cooperative cancellation request stopped computation before output."""


def _positive_number(value: float, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite positive number")
    result = float(value)
    if not np.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be a finite positive number")
    return result


def _positive_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _cancel_requested(cancel: Callable[[], bool] | object | None) -> bool:
    if cancel is None:
        return False
    if callable(cancel):
        return bool(cancel())
    is_set = getattr(cancel, "is_set", None)
    if callable(is_set):
        return bool(is_set())
    raise TypeError("cancel must be callable or provide an is_set() method")


def _report(progress: ProgressCallback | None, stage: str, completed: int, total: int) -> None:
    if progress is not None:
        progress(stage, completed, total)


def _as_grayscale(frame: np.ndarray, *, path: Path | None = None) -> np.ndarray:
    array = np.asarray(frame)
    if array.ndim == 2:
        gray = array
    elif array.ndim == 3 and array.shape[2] == 1:
        gray = array[:, :, 0]
    elif array.ndim == 3 and array.shape[2] == 3:
        gray = cv2.cvtColor(array, cv2.COLOR_BGR2GRAY)
    elif array.ndim == 3 and array.shape[2] == 4:
        gray = cv2.cvtColor(array, cv2.COLOR_BGRA2GRAY)
    else:
        location = f" in {path}" if path is not None else ""
        raise VideoReadError(f"decoded frame must be a 2D grayscale or 3/4-channel image{location}")
    gray = np.asarray(gray)
    try:
        valid = np.issubdtype(gray.dtype, np.number) and np.all(np.isfinite(gray))
    except TypeError:
        valid = False
    if not valid:
        location = f" in {path}" if path is not None else ""
        raise VideoReadError(f"decoded frame contains invalid values{location}")
    return gray


class OpenCVFrameReader:
    """Iterable OpenCV reader with validation and guaranteed capture release."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)

    def __iter__(self) -> Iterator[np.ndarray]:
        capture = cv2.VideoCapture(str(self.path))
        count = 0
        shape: tuple[int, int] | None = None
        try:
            if not capture.isOpened():
                raise VideoReadError(f"could not open video: {self.path}")
            while True:
                ok, frame = capture.read()
                if not ok:
                    break
                gray = _as_grayscale(frame, path=self.path)
                if shape is None:
                    shape = gray.shape
                elif gray.shape != shape:
                    raise VideoReadError(
                        f"video frames have inconsistent shapes in {self.path}: "
                        f"expected {shape}, got {gray.shape}"
                    )
                count += 1
                yield gray
            if count == 0:
                raise VideoReadError(f"video contains no decodable frames: {self.path}")
        except cv2.error as error:
            raise VideoReadError(f"could not decode video {self.path}: {error}") from error
        finally:
            capture.release()


def read_video_frames(path: str | Path) -> OpenCVFrameReader:
    """Return a validated-on-consumption frame iterator for an AVI/video path."""
    return OpenCVFrameReader(path)


def log_spaced_lags(frame_count: int, points_per_decade: float = 20) -> np.ndarray:
    """Return unique increasing integer lags in ``[1, frame_count - 1]``."""
    if isinstance(frame_count, bool) or not isinstance(frame_count, (int, np.integer)):
        raise TypeError("frame_count must be an integer")
    frame_count = int(frame_count)
    if frame_count < 2:
        raise ValueError("at least two frames are required")
    points = _positive_number(points_per_decade, "points_per_decade")
    if not points.is_integer():
        raise ValueError("points_per_decade must be an integer")
    number = max(1, int(np.ceil(np.log10(frame_count) * points)))
    raw = np.logspace(0.0, np.log10(frame_count), num=number + 1, endpoint=False)
    lags = np.unique(raw.astype(int))
    lags = lags[(lags >= 1) & (lags < frame_count)]
    if lags.size == 0 or lags[0] != 1:
        lags = np.insert(lags, 0, 1)
    return lags.astype(int, copy=False)


def tukey_two_dimensional(shape: tuple[int, int], alpha: float) -> np.ndarray:
    """Create a separable Tukey window for an explicitly rectangular shape."""
    if len(shape) != 2 or any(size <= 0 for size in shape):
        raise ValueError("shape must contain two positive dimensions")
    if not np.isfinite(alpha) or not 0 <= alpha <= 1:
        raise ValueError("Tukey alpha must be between 0 and 1")
    return np.outer(tukey(shape[0], alpha), tukey(shape[1], alpha))


class RadialAverager:
    """Legacy-compatible radial or opposite-direction sector averaging."""

    def __init__(self, shape: tuple[int, int], sectors: int = 1) -> None:
        if len(shape) != 2 or shape[0] != shape[1] or shape[0] < 2:
            raise ValueError("DDM currently requires square frames with at least 2 pixels per side")
        self.shape = tuple(int(size) for size in shape)
        self.sectors = _positive_integer(sectors, "sectors")
        size = self.shape[0]
        self.distances = np.sqrt(
            np.fft.fftfreq(size)[:, None] ** 2 + np.fft.fftfreq(size)[None, :] ** 2
        )
        self.bins = np.arange(size // 2 + 1, dtype=float) / size
        if self.sectors == 1:
            self._masks = (self.distances <= 0.5,)
        else:
            with np.errstate(divide="ignore", invalid="ignore"):
                args = np.arctan(
                    np.fft.fftfreq(size)[None, :] / np.fft.fftfreq(size)[:, None]
                ) + np.pi / 2
            edges = np.arange(-0.5, self.sectors + 0.01) / self.sectors * np.pi
            args[args > (self.sectors - 0.5) / self.sectors * np.pi] -= np.pi
            edges[0] -= 1e-10
            self._masks = tuple(
                ((args > edges[index]) & (args <= edges[index + 1]) & (self.distances <= 0.5))
                | np.isnan(args)
                for index in range(self.sectors)
            )
        self._counts = tuple(
            np.histogram(self.distances[mask], self.bins)[0] for mask in self._masks
        )

    @property
    def q_bin_centers(self) -> np.ndarray:
        return (self.bins[:-1] + self.bins[1:]) / 2.0

    def __call__(self, spectrum: np.ndarray) -> tuple[np.ndarray, ...]:
        values = np.asarray(spectrum)
        if values.shape != self.shape:
            raise ValueError(f"spectrum shape must be {self.shape}, got {values.shape}")
        curves = []
        for mask, counts in zip(self._masks, self._counts):
            weighted = np.histogram(self.distances[mask], self.bins, weights=values[mask])[0]
            curves.append(
                np.divide(
                    weighted,
                    counts,
                    out=np.zeros_like(weighted, dtype=float),
                    where=counts != 0,
                )
            )
        return tuple(curves)


def _validate_frames(
    frames: Iterable[np.ndarray],
    cancel: Callable[[], bool] | object | None,
    progress: ProgressCallback | None,
) -> list[np.ndarray]:
    validated: list[np.ndarray] = []
    shape: tuple[int, int] | None = None
    for index, frame in enumerate(frames, start=1):
        if _cancel_requested(cancel):
            raise ComputationCancelled("DDM computation cancelled while reading frames")
        gray = _as_grayscale(frame)
        if shape is None:
            shape = gray.shape
        elif gray.shape != shape:
            raise ValueError(f"frames have inconsistent shapes: expected {shape}, got {gray.shape}")
        if shape[0] != shape[1] or shape[0] < 2:
            raise ValueError("DDM currently requires square frames with at least 2 pixels per side")
        validated.append(gray.astype(float, copy=False))
        _report(progress, "frame_read", index, 0)
    if len(validated) < 2:
        raise ValueError("at least two frames are required")
    return validated


def compute_ddm(
    frames: Iterable[np.ndarray],
    frame_rate: float,
    pixel_size: float,
    *,
    max_couples: int = 300,
    points_per_decade: float = 20,
    sectors: int = 1,
    progress: ProgressCallback | None = None,
    cancel: Callable[[], bool] | object | None = None,
    tukey_alpha: float | None = None,
) -> DDMData | tuple[DDMData, ...]:
    """Compute isotropic or directional DDM data from an in-memory iterator."""
    frame_rate = _positive_number(frame_rate, "frame_rate")
    pixel_size = _positive_number(pixel_size, "pixel_size")
    if (
        isinstance(max_couples, bool)
        or not isinstance(max_couples, (int, np.integer))
        or max_couples < 0
    ):
        raise ValueError("max_couples must be a non-negative integer")
    sectors = _positive_integer(sectors, "sectors")
    source_frames = _validate_frames(frames, cancel, progress)
    frame_count = len(source_frames)
    lags = log_spaced_lags(frame_count, points_per_decade)
    window = (
        None
        if tukey_alpha is None
        else tukey_two_dimensional(source_frames[0].shape, tukey_alpha)
    )
    transformed = []
    for index, frame in enumerate(source_frames, start=1):
        if _cancel_requested(cancel):
            raise ComputationCancelled("DDM computation cancelled during FFT")
        transformed.append(np.fft.fft2(frame if window is None else frame * window))
        _report(progress, "frame_fft", index, frame_count)

    averager = RadialAverager(source_frames[0].shape, sectors)
    matrices = [
        np.empty((lags.size, averager.q_bin_centers.size), dtype=float)
        for _ in range(sectors)
    ]
    for lag_index, lag in enumerate(lags):
        if _cancel_requested(cancel):
            raise ComputationCancelled("DDM computation cancelled during lag averaging")
        increment = (
            1
            if max_couples == 0
            else max(1, int(np.ceil((frame_count - int(lag)) / max_couples)))
        )
        starts = range(0, frame_count - int(lag), increment)
        accumulated = np.zeros(averager.shape, dtype=float)
        used = 0
        for start in starts:
            if _cancel_requested(cancel):
                raise ComputationCancelled("DDM computation cancelled during lag averaging")
            difference = transformed[start + int(lag)] - transformed[start]
            accumulated += np.abs(difference) ** 2
            used += 1
        if used == 0:
            raise DDMEngineError(f"no frame pairs available for lag {lag}")
        for sector, curve in enumerate(averager(accumulated / used)):
            matrices[sector][lag_index] = curve
        _report(progress, "lag_average", lag_index + 1, lags.size)

    q_values = 2.0 * np.pi * averager.q_bin_centers / pixel_size
    results = tuple(
        DDMData(matrix=matrix, lag_times=lags / frame_rate, q_values=q_values)
        for matrix in matrices
    )
    return results[0] if sectors == 1 else results


def compute_video_ddm(
    path: str | Path,
    frame_rate: float,
    pixel_size: float,
    **kwargs: object,
) -> DDMData | tuple[DDMData, ...]:
    """Compute DDM data from a video without trusting its reported frame count."""
    return compute_ddm(read_video_frames(path), frame_rate, pixel_size, **kwargs)


def compute_and_save_video_ddm(
    path: str | Path,
    output_prefix: str | Path,
    frame_rate: float,
    pixel_size: float,
    **kwargs: object,
) -> tuple[Path, ...]:
    """Compute fully before saving one or more legacy matrix sets."""
    result = compute_video_ddm(path, frame_rate, pixel_size, **kwargs)
    datasets = (result,) if isinstance(result, DDMData) else result
    prefix = Path(output_prefix)
    if len(datasets) == 1:
        return save_matrix_set(prefix, datasets[0])
    prefix_name = prefix.name
    for suffix in ("_DDM_matrix.npy", "_deltaTs.npy", "_QS.npy"):
        if prefix_name.endswith(suffix):
            prefix_name = prefix_name[: -len(suffix)]
            break
    paths: list[Path] = []
    for index, dataset in enumerate(datasets):
        angle = index * 180.0 / len(datasets)
        paths.extend(save_matrix_set(prefix.with_name(f"{prefix_name}_{angle:.1f}_"), dataset))
    return tuple(paths)


__all__ = [
    "ComputationCancelled",
    "DDMEngineError",
    "FrameIterator",
    "OpenCVFrameReader",
    "ProgressCallback",
    "RadialAverager",
    "VideoReadError",
    "compute_and_save_video_ddm",
    "compute_ddm",
    "compute_video_ddm",
    "log_spaced_lags",
    "read_video_frames",
    "tukey_two_dimensional",
]
