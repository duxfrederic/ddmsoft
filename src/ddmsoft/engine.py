"""OpenCV video input and compatibility DDM computation.

The engine accepts any iterable yielding two-dimensional frames, so file input
and deterministic tests use the same computation path.  It intentionally keeps
the first implementation straightforward rather than optimizing FFT storage.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass
from numbers import Real
from pathlib import Path
from struct import Struct
from struct import error as StructError
from tempfile import TemporaryDirectory
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


@dataclass(frozen=True)
class _RawAVIInfo:
    width: int
    height: int
    bits_per_pixel: int
    stream_id: bytes
    movi_start: int
    movi_end: int


def _scan_avi_chunks(
    stream, start: int, end: int, state: dict[str, object], stream_id: bytes | None = None
) -> None:
    stream.seek(start)
    while stream.tell() + 8 <= end:
        chunk_id = stream.read(4)
        chunk_size_data = stream.read(4)
        if len(chunk_id) != 4 or len(chunk_size_data) != 4:
            raise ValueError("truncated AVI chunk header")
        chunk_size = Struct("<I").unpack(chunk_size_data)[0]
        data_start = stream.tell()
        data_end = data_start + chunk_size
        if data_end > end:
            raise ValueError("AVI chunk extends beyond its containing list")
        if chunk_id == b"LIST":
            if chunk_size < 4:
                raise ValueError("AVI LIST chunk has no type")
            list_type = stream.read(4)
            content_start = stream.tell()
            if list_type == b"movi":
                state["movi"] = (content_start, data_end)
                return
            child_stream_id = stream_id
            if list_type == b"strl":
                child_stream_id = f"{state['stream_count']:02d}".encode("ascii")
                state["stream_count"] = int(state["stream_count"]) + 1
            _scan_avi_chunks(stream, content_start, data_end, state, child_stream_id)
        elif chunk_id == b"strh" and chunk_size >= 8:
            header = stream.read(min(chunk_size, 56))
            if header[:4] == b"vids":
                state["video_stream_id"] = stream_id or b"00"
        elif chunk_id == b"strf" and state.get("video_stream_id") == stream_id:
            format_header = stream.read(min(chunk_size, 40))
            if len(format_header) >= 20:
                width, height, _, bits_per_pixel, compression = Struct("<iiHHI").unpack_from(
                    format_header, 4
                )
                state["format"] = (width, height, bits_per_pixel, compression)
        stream.seek(data_end + (chunk_size & 1))


def _raw_avi_info(path: Path) -> _RawAVIInfo | None:
    """Return format details for uncompressed DIB AVI files, if recognized."""
    try:
        with path.open("rb") as stream:
            header = stream.read(12)
            if len(header) != 12 or header[:4] != b"RIFF" or header[8:] != b"AVI ":
                return None
            riff_size = Struct("<I").unpack_from(header, 4)[0]
            state: dict[str, object] = {"stream_count": 0}
            _scan_avi_chunks(stream, 12, min(8 + riff_size, path.stat().st_size), state)
    except (OSError, StructError, ValueError):
        return None
    format_values = state.get("format")
    movi_values = state.get("movi")
    stream_id = state.get("video_stream_id")
    if (
        not isinstance(format_values, tuple)
        or len(format_values) != 4
        or not isinstance(movi_values, tuple)
        or len(movi_values) != 2
        or not isinstance(stream_id, bytes)
    ):
        return None
    width, height, bits_per_pixel, compression = format_values
    if (
        not isinstance(width, int)
        or not isinstance(height, int)
        or not isinstance(bits_per_pixel, int)
        or not isinstance(compression, int)
        or width <= 0
        or height == 0
        or compression != 0
        or bits_per_pixel not in (8, 24, 32)
    ):
        return None
    movi_start, movi_end = movi_values
    if not isinstance(movi_start, int) or not isinstance(movi_end, int):
        return None
    return _RawAVIInfo(width, height, bits_per_pixel, stream_id, movi_start, movi_end)


def _raw_avi_payloads(stream, start: int, end: int, stream_id: bytes) -> Iterator[bytes]:
    stream.seek(start)
    while stream.tell() + 8 <= end:
        chunk_id = stream.read(4)
        chunk_size_data = stream.read(4)
        if len(chunk_id) != 4 or len(chunk_size_data) != 4:
            raise ValueError("truncated AVI frame chunk header")
        chunk_size = Struct("<I").unpack(chunk_size_data)[0]
        data_start = stream.tell()
        data_end = data_start + chunk_size
        if data_end > end:
            raise ValueError("AVI frame chunk extends beyond the movie list")
        if chunk_id == b"LIST":
            if chunk_size < 4:
                raise ValueError("AVI LIST chunk has no type")
            yield from _raw_avi_payloads(stream, data_start + 4, data_end, stream_id)
        elif (
            len(stream_id) == 2
            and chunk_id[:2] == stream_id
            and chunk_id[2:] in (b"db", b"dc")
        ):
            payload = stream.read(chunk_size)
            if len(payload) != chunk_size:
                raise ValueError("truncated AVI frame payload")
            yield payload
        stream.seek(data_end + (chunk_size & 1))


class RawAVIFrameReader:
    """Reader for uncompressed DIB AVI files that avoids a crashing OpenCV path."""

    def __init__(self, path: str | Path, info: _RawAVIInfo) -> None:
        self.path = Path(path)
        self.info = info

    def __iter__(self) -> Iterator[np.ndarray]:
        height = abs(self.info.height)
        channels = self.info.bits_per_pixel // 8
        row_bytes = self.info.width * channels
        row_stride = (row_bytes + 3) & ~3
        expected_size = row_stride * height
        count = 0
        shape: tuple[int, int] | None = None
        try:
            with self.path.open("rb") as stream:
                for payload in _raw_avi_payloads(
                    stream, self.info.movi_start, self.info.movi_end, self.info.stream_id
                ):
                    if len(payload) < expected_size:
                        raise ValueError("raw AVI frame is shorter than its declared dimensions")
                    rows = np.frombuffer(payload[:expected_size], dtype=np.uint8).reshape(
                        height, row_stride
                    )[:, :row_bytes]
                    frame = rows.reshape(height, self.info.width, channels)
                    if self.info.height > 0:
                        frame = frame[::-1]
                    gray = _as_grayscale(np.array(frame, copy=True), path=self.path)
                    if shape is None:
                        shape = gray.shape
                    elif gray.shape != shape:
                        raise VideoReadError(
                            f"video frames have inconsistent shapes in {self.path}: "
                            f"expected {shape}, got {gray.shape}"
                        )
                    count += 1
                    yield gray
        except (OSError, StructError, ValueError) as error:
            raise VideoReadError(f"could not decode video {self.path}: {error}") from error
        if count == 0:
            raise VideoReadError(f"video contains no decodable frames: {self.path}")


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


def read_video_frames(path: str | Path) -> FrameIterator:
    """Return a validated-on-consumption frame iterator for an AVI/video path."""
    video_path = Path(path)
    raw_info = _raw_avi_info(video_path)
    if raw_info is not None:
        return RawAVIFrameReader(video_path, raw_info)
    return OpenCVFrameReader(video_path)


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


def _cache_frames(
    frames: Iterable[np.ndarray],
    cancel: Callable[[], bool] | object | None,
    progress: ProgressCallback | None,
    path: Path,
) -> tuple[int, tuple[int, int]]:
    count = 0
    shape: tuple[int, int] | None = None
    with path.open("wb") as cache:
        for index, frame in enumerate(frames, start=1):
            if _cancel_requested(cancel):
                raise ComputationCancelled("DDM computation cancelled while reading frames")
            gray = _as_grayscale(frame)
            if shape is None:
                shape = gray.shape
            elif gray.shape != shape:
                raise ValueError(
                    f"frames have inconsistent shapes: expected {shape}, got {gray.shape}"
                )
            if shape[0] != shape[1] or shape[0] < 2:
                raise ValueError(
                    "DDM currently requires square frames with at least 2 pixels per side"
                )
            np.asarray(gray, dtype=np.float32).tofile(cache)
            count = index
            _report(progress, "frame_read", index, 0)
    if count < 2 or shape is None:
        raise ValueError("at least two frames are required")
    return count, shape


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
    cache_directory: str | Path | None = None,
) -> DDMData | tuple[DDMData, ...]:
    """Compute isotropic or directional DDM data using a disk-backed cache."""
    frame_rate = _positive_number(frame_rate, "frame_rate")
    pixel_size = _positive_number(pixel_size, "pixel_size")
    if (
        isinstance(max_couples, bool)
        or not isinstance(max_couples, (int, np.integer))
        or max_couples < 0
    ):
        raise ValueError("max_couples must be a non-negative integer")
    sectors = _positive_integer(sectors, "sectors")
    with TemporaryDirectory(prefix=".ddmsoft-cache-", dir=cache_directory) as temporary:
        temporary_directory = Path(temporary)
        frame_path = temporary_directory / "frames.dat"
        frame_count, frame_shape = _cache_frames(frames, cancel, progress, frame_path)
        lags = log_spaced_lags(frame_count, points_per_decade)
        window = (
            None
            if tukey_alpha is None
            else tukey_two_dimensional(frame_shape, tukey_alpha)
        )
        source_frames = np.memmap(
            frame_path,
            dtype=np.float32,
            mode="r",
            shape=(frame_count, *frame_shape),
        )
        fft_path = temporary_directory / "transformed.dat"
        transformed: np.memmap | None = None
        try:
            transformed = np.memmap(
                fft_path,
                dtype=np.complex64,
                mode="w+",
                shape=(frame_count, *frame_shape),
            )
            for index, frame in enumerate(source_frames, start=1):
                if _cancel_requested(cancel):
                    raise ComputationCancelled("DDM computation cancelled during FFT")
                transformed[index - 1] = np.asarray(
                    np.fft.fft2(frame if window is None else frame * window),
                    dtype=np.complex64,
                )
                _report(progress, "frame_fft", index, frame_count)
        finally:
            del source_frames
            frame_path.unlink(missing_ok=True)

        try:
            averager = RadialAverager(frame_shape, sectors)
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
        finally:
            if transformed is not None:
                transformed.flush()
                del transformed


def compute_video_ddm(
    path: str | Path,
    frame_rate: float,
    pixel_size: float,
    **kwargs: object,
) -> DDMData | tuple[DDMData, ...]:
    """Compute DDM data from a video without trusting its reported frame count."""
    video_path = Path(path)
    options = dict(kwargs)
    options.setdefault("cache_directory", video_path.parent)
    return compute_ddm(read_video_frames(video_path), frame_rate, pixel_size, **options)


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
    "RawAVIFrameReader",
    "VideoReadError",
    "compute_and_save_video_ddm",
    "compute_ddm",
    "compute_video_ddm",
    "log_spaced_lags",
    "read_video_frames",
    "tukey_two_dimensional",
]
