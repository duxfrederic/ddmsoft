"""File loading and export helpers for the legacy DDMSoft data format.

This module deliberately contains no GUI code.  The three NumPy files used by
the original application remain the persistence format:
``*_DDM_matrix.npy``, ``*_deltaTs.npy``, and ``*_QS.npy``.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .models import DDMData, VideoMetadata
from .science import stokes_einstein_radius

LEGACY_SUFFIXES = ("_DDM_matrix.npy", "_deltaTs.npy", "_QS.npy")
CSV_SUFFIXES = ("_DDM_matrix.csv", "_deltaTs.csv", "_QS.csv")
DDM_MATRICES_DIRECTORY = "ddm_matrices"


class DDMIOError(Exception):
    """Base class for errors raised by the modern I/O layer."""

    def __init__(self, message: str, *, path: Path | None = None) -> None:
        super().__init__(message)
        self.path = path


class MetadataError(DDMIOError):
    """Base class for acquisition metadata errors."""


class MetadataDirectoryError(MetadataError):
    """The requested acquisition directory does not exist or is not a directory."""


class MetadataFileError(MetadataError):
    """A metadata file cannot be read or contains malformed content."""


class MissingMetadataError(MetadataError):
    """A required metadata file or value is absent."""


class InvalidMetadataError(MetadataError):
    """A required metadata value is not a finite positive number."""


class AmbiguousMetadataError(MetadataError):
    """Metadata files cannot be unambiguously assigned to videos."""


class MatrixError(DDMIOError):
    """Base class for legacy matrix loading errors."""


class IncompleteMatrixError(MatrixError):
    """One or more files in a legacy matrix set are missing."""

    def __init__(self, base: Path, missing: Iterable[Path]) -> None:
        missing_paths = tuple(Path(path) for path in missing)
        text = ", ".join(str(path) for path in missing_paths)
        super().__init__(f"incomplete legacy matrix set for {base}: missing {text}", path=base)
        self.base = base
        self.missing = missing_paths


class MatrixLoadError(MatrixError):
    """A complete legacy matrix set contains invalid NumPy data."""


@dataclass(frozen=True)
class MatrixFileSet:
    """The three files making up one legacy matrix dataset."""

    base: Path
    matrix: Path
    delta_times: Path
    q_values: Path

    @property
    def display_name(self) -> str:
        return self.matrix.name

    def load(self) -> DDMData:
        """Load this set and validate its array dimensions through ``DDMData``."""
        try:
            matrix = np.load(self.matrix, allow_pickle=False)
            delta_times = np.load(self.delta_times, allow_pickle=False)
            q_values = np.load(self.q_values, allow_pickle=False)
            return DDMData(
                matrix=matrix,
                lag_times=delta_times,
                q_values=q_values,
                source=self.matrix,
            )
        except (OSError, ValueError, TypeError) as error:
            raise MatrixLoadError(
                f"could not load legacy matrix set {self.base}: {error}", path=self.base
            ) from error


def parse_metadata_file(path: str | Path) -> dict[str, str]:
    """Parse an acquisition text file, ignoring blank and comment lines.

    Only the first colon is structural.  This preserves values such as URLs,
    timestamps, or Windows paths that contain additional colons.
    """
    metadata_path = Path(path)
    if not metadata_path.is_file():
        raise MissingMetadataError(
            f"metadata file does not exist: {metadata_path}", path=metadata_path
        )

    values: dict[str, str] = {}
    try:
        with metadata_path.open("r", encoding="utf-8") as metadata_file:
            for line_number, line in enumerate(metadata_file, start=1):
                content = line.strip()
                if not content or content.startswith("#"):
                    continue
                if ":" not in content:
                    raise MetadataFileError(
                        f"malformed metadata at {metadata_path}:{line_number}: expected key: value",
                        path=metadata_path,
                    )
                key, value = content.split(":", 1)
                key = key.strip().lower()
                value = value.strip()
                if not key:
                    raise MetadataFileError(
                        f"malformed metadata at {metadata_path}:{line_number}: empty key",
                        path=metadata_path,
                    )
                values[key] = value
    except UnicodeError as error:
        raise MetadataFileError(
            f"metadata file is not valid UTF-8: {metadata_path}", path=metadata_path
        ) from error
    return values


def _metadata_number(values: Mapping[str, str], key: str, path: Path) -> float:
    value = values.get(key)
    if value is None or not value:
        raise MissingMetadataError(f"missing required metadata value '{key}' in {path}", path=path)
    try:
        number = float(value)
    except ValueError as error:
        raise InvalidMetadataError(
            f"metadata value '{key}' is not numeric in {path}: {value!r}", path=path
        ) from error
    if not np.isfinite(number) or number <= 0:
        raise InvalidMetadataError(
            f"metadata value '{key}' must be finite and positive in {path}: {value!r}", path=path
        )
    return number


def metadata_for_video(video: str | Path, metadata_path: str | Path) -> VideoMetadata:
    """Read and validate metadata for one video."""
    video_path = Path(video)
    path = Path(metadata_path)
    values = parse_metadata_file(path)
    frame_rate = _metadata_number(values, "framerate", path)
    pixel_size = _metadata_number(values, "pixelsize", path)
    temperature = None
    if values.get("temperature"):
        try:
            temperature = float(values["temperature"])
        except ValueError as error:
            raise InvalidMetadataError(
                f"metadata value 'temperature' is not numeric in {path}: {values['temperature']!r}",
                path=path,
            ) from error
        if not np.isfinite(temperature):
            raise InvalidMetadataError(f"temperature must be finite in {path}", path=path)
    return VideoMetadata(video_path, frame_rate, pixel_size, temperature)


def _files_by_suffix(directory: Path, suffix: str) -> list[Path]:
    return sorted(
        (path for path in directory.iterdir() if path.is_file() and path.suffix.lower() == suffix),
        key=lambda path: path.name.casefold(),
    )


def load_directory(directory: str | Path) -> dict[Path, VideoMetadata]:
    """Load all AVI acquisition metadata using legacy assignment rules.

    One text file is shared by every video.  When there is one text file per
    video, names must match by stem.  Any other combination is ambiguous and
    raises ``AmbiguousMetadataError`` rather than silently pairing sorted files.
    """
    root = Path(directory)
    if not root.is_dir():
        raise MetadataDirectoryError(f"acquisition directory does not exist: {root}", path=root)
    videos = sorted(
        (path for path in root.iterdir() if path.is_file() and path.suffix.lower() == ".avi"),
        key=lambda path: path.name.casefold(),
    )
    if not videos:
        raise MissingMetadataError(f"no AVI videos found in {root}", path=root)
    metadata_files = _files_by_suffix(root, ".txt")
    if not metadata_files:
        raise MissingMetadataError(f"no acquisition metadata file found in {root}", path=root)

    assignments: dict[Path, Path]
    if len(metadata_files) == 1:
        assignments = {video: metadata_files[0] for video in videos}
    elif len(metadata_files) == len(videos) and len(videos) > 1:
        by_stem = {path.stem.casefold(): path for path in metadata_files}
        if len(by_stem) != len(metadata_files) or any(
            video.stem.casefold() not in by_stem for video in videos
        ):
            raise AmbiguousMetadataError(
                f"metadata files do not match video names in {root}; expected one file per video", path=root
            )
        assignments = {video: by_stem[video.stem.casefold()] for video in videos}
    else:
        raise AmbiguousMetadataError(
            f"found {len(metadata_files)} metadata files for {len(videos)} videos in {root}", path=root
        )
    return {
        video: metadata_for_video(video, metadata_path)
        for video, metadata_path in assignments.items()
    }


def discover_matrix_sets(
    directory: str | Path, *, strict: bool = True
) -> tuple[MatrixFileSet, ...]:
    """Discover legacy matrix sets below ``directory/ddm_matrices``.

    With ``strict=False`` incomplete candidates are explicitly skipped, matching
    the legacy catalog behavior.  Strict mode is useful for validation and
    raises ``IncompleteMatrixError`` on the first incomplete candidate.
    """
    root = Path(directory) / DDM_MATRICES_DIRECTORY
    if not root.is_dir():
        return ()
    sets: list[MatrixFileSet] = []
    for matrix_path in sorted(
        root.glob(f"*{LEGACY_SUFFIXES[0]}"), key=lambda path: path.name.casefold()
    ):
        base = Path(str(matrix_path)[: -len(LEGACY_SUFFIXES[0])])
        paths = tuple(base.parent / f"{base.name}{suffix}" for suffix in LEGACY_SUFFIXES)
        missing = [path for path in paths if not path.is_file()]
        if missing:
            if strict:
                raise IncompleteMatrixError(base, missing)
            continue
        sets.append(MatrixFileSet(base, paths[0], paths[1], paths[2]))
    return tuple(sets)


def load_matrices(directory: str | Path, *, strict: bool = False) -> dict[Path, DDMData]:
    """Load discovered matrices keyed by their full matrix-file path."""
    return {
        matrix_set.matrix: matrix_set.load()
        for matrix_set in discover_matrix_sets(directory, strict=strict)
    }


def display_names(matrices: Iterable[MatrixFileSet]) -> dict[str, Path]:
    """Map human-readable names to full matrix paths without collisions."""
    sets = tuple(matrices)
    result: dict[str, Path] = {}
    for matrix_set in sets:
        name = matrix_set.display_name
        if name in result:
            name = f"{matrix_set.matrix.parent.parent.name}/{name}"
        result[name] = matrix_set.matrix
    return result


def _as_ddm_data(data: DDMData | Sequence[np.ndarray]) -> DDMData:
    if isinstance(data, DDMData):
        return data
    try:
        matrix, delta_times, q_values = data
    except (TypeError, ValueError) as error:
        raise TypeError("matrix data must be DDMData or a three-item sequence") from error
    return DDMData(matrix=matrix, lag_times=delta_times, q_values=q_values)


def _prefix(path: str | Path, suffixes: Iterable[str]) -> Path:
    result = Path(path)
    name = result.name
    for suffix in suffixes:
        if name.endswith(suffix):
            return result.with_name(name[: -len(suffix)])
    return result


def save_matrix_csv(
    output_prefix: str | Path, data: DDMData | Sequence[np.ndarray]
) -> tuple[Path, Path, Path]:
    """Write the three matrix exports and return their exact paths."""
    ddm_data = _as_ddm_data(data)
    prefix = _prefix(output_prefix, CSV_SUFFIXES)
    paths = tuple(prefix.with_name(prefix.name + suffix) for suffix in CSV_SUFFIXES)
    for path, array in zip(paths, (ddm_data.matrix, ddm_data.lag_times, ddm_data.q_values)):
        np.savetxt(path, array, fmt="%.6e", delimiter="\t")
    return paths


def save_matrix_set(
    output_prefix: str | Path, data: DDMData | Sequence[np.ndarray]
) -> tuple[Path, Path, Path]:
    """Atomically write a legacy three-file NumPy matrix set.

    Temporary files are created beside the destinations and are only replaced
    after all three arrays have been serialized successfully.  This prevents a
    cancelled or failed computation from looking like a complete dataset.
    """
    ddm_data = _as_ddm_data(data)
    prefix = _prefix(output_prefix, LEGACY_SUFFIXES)
    paths = tuple(prefix.with_name(prefix.name + suffix) for suffix in LEGACY_SUFFIXES)
    temporary_paths: list[Path] = []
    try:
        for index, array in enumerate((ddm_data.matrix, ddm_data.lag_times, ddm_data.q_values)):
            temporary = paths[index].with_name(f".{paths[index].name}.tmp")
            temporary_paths.append(temporary)
            with temporary.open("wb") as temporary_file:
                np.save(temporary_file, array, allow_pickle=False)
        for temporary, destination in zip(temporary_paths, paths):
            temporary.replace(destination)
    except (OSError, ValueError, TypeError):
        for temporary in temporary_paths:
            temporary.unlink(missing_ok=True)
        raise
    return paths


def save_partitioned_matrix_sets(
    output_prefix: str | Path,
    results: Iterable[object],
) -> tuple[Path, ...]:
    """Save time-dependent results using their actual starting-frame names."""
    paths: list[Path] = []
    prefix = Path(output_prefix)
    for result in results:
        start_frame = result.start_frame
        data = result.data
        datasets = (data,) if isinstance(data, DDMData) else tuple(data)
        partition_prefix = prefix.with_name(f"{prefix.name}__i={start_frame}__")
        if len(datasets) == 1:
            paths.extend(save_matrix_set(partition_prefix, datasets[0]))
        else:
            for index, dataset in enumerate(datasets):
                angle = index * 180.0 / len(datasets)
                paths.extend(
                    save_matrix_set(
                        partition_prefix.with_name(f"{partition_prefix.name}{angle:.1f}_"),
                        dataset,
                    )
                )
    return tuple(paths)


def save_autocorrelation_csv(
    output_prefix: str | Path,
    data: DDMData | Sequence[np.ndarray],
    amplitude: np.ndarray,
    noise: np.ndarray,
    q_values: Iterable[float],
) -> tuple[Path, Path, Path]:
    """Export refined correlation data for the requested physical q range."""
    ddm_data = _as_ddm_data(data)
    amplitude = np.asarray(amplitude, dtype=float)
    noise = np.asarray(noise, dtype=float)
    selected_q = np.asarray(tuple(q_values), dtype=float)
    if (
        amplitude.ndim != 1
        or noise.ndim != 1
        or amplitude.size != ddm_data.q_values.size
        or noise.size != ddm_data.q_values.size
    ):
        raise ValueError("amplitude and noise must contain one value per matrix q value")
    if selected_q.ndim != 1 or selected_q.size == 0:
        raise ValueError("q_values must contain at least one q value")
    q_indices = np.asarray([int(np.argmin(np.abs(ddm_data.q_values - q))) for q in selected_q])
    matrix = 1.0 - (ddm_data.matrix[:, q_indices] - noise[q_indices]) / amplitude[q_indices]
    prefix = _prefix(output_prefix, ("_autocorrelationmatrix.csv", "_qs.csv", "_dts.csv"))
    paths = (
        prefix.with_name(prefix.name + "_autocorrelationmatrix.csv"),
        prefix.with_name(prefix.name + "_qs.csv"),
        prefix.with_name(prefix.name + "_dts.csv"),
    )
    np.savetxt(paths[0], matrix, fmt="%.6e", delimiter="\t")
    np.savetxt(paths[1], ddm_data.q_values[q_indices], fmt="%.6e", delimiter="\t")
    np.savetxt(paths[2], ddm_data.lag_times, fmt="%.6e", delimiter="\t")
    return paths


def save_fit_text(
    path: str | Path,
    q_values: Iterable[float],
    amplitude: Iterable[float],
    noise: Iterable[float],
    parameters: Sequence[Iterable[float]],
    parameter_names: Sequence[str],
    *,
    viscosity: float | None = None,
    temperature: float | None = None,
) -> Path:
    """Write the legacy tab-delimited fit export using an explicit suffix."""
    qs = np.asarray(tuple(q_values), dtype=float)
    amplitudes = np.asarray(tuple(amplitude), dtype=float)
    noises = np.asarray(tuple(noise), dtype=float)
    parameter_arrays = tuple(np.asarray(tuple(parameter), dtype=float) for parameter in parameters)
    names = tuple(parameter_names)
    count = qs.size
    if amplitudes.size != count or noises.size != count or len(parameter_arrays) != len(names):
        raise ValueError("fit export columns have inconsistent lengths")
    if any(parameter.size != count for parameter in parameter_arrays):
        raise ValueError("every fit parameter must contain one value per q value")
    radius: np.ndarray | None = None
    if (viscosity is None) != (temperature is None):
        raise ValueError("viscosity and temperature must be provided together")
    if viscosity is not None and temperature is not None:
        temperature_kelvin = temperature + 273.15
        radius = np.asarray(
            [
                stokes_einstein_radius(value, temperature_kelvin, viscosity) * 1e9
                for value in parameter_arrays[0]
            ]
        )
    output = Path(path)
    if output.suffix.lower() not in {".txt", ".csv"}:
        output = output.with_name(output.name + ".txt")
    header = ["q [m^-1]", "A", "B", *names]
    if radius is not None:
        header.append("R_H eff (nm)")
    with output.open("w", encoding="utf-8", newline="\n") as save_file:
        save_file.write("\t".join(header) + "\n")
        for index in range(count):
            values = [qs[index], amplitudes[index], noises[index]]
            values.extend(parameter[index] for parameter in parameter_arrays)
            if radius is not None:
                values.append(radius[index])
            save_file.write("\t".join(f"{value:.3e}" for value in values) + "\n")
    return output


# Names used by the first package migration are intentionally simple aliases.
parse_metadata = parse_metadata_file
load_video_metadata = load_directory
save_matrix = save_matrix_csv
save_fit = save_fit_text


__all__ = [
    "CSV_SUFFIXES",
    "DDM_MATRICES_DIRECTORY",
    "LEGACY_SUFFIXES",
    "AmbiguousMetadataError",
    "DDMIOError",
    "IncompleteMatrixError",
    "InvalidMetadataError",
    "MatrixError",
    "MatrixFileSet",
    "MatrixLoadError",
    "MetadataDirectoryError",
    "MetadataError",
    "MetadataFileError",
    "MissingMetadataError",
    "discover_matrix_sets",
    "display_names",
    "load_directory",
    "load_matrices",
    "load_video_metadata",
    "metadata_for_video",
    "parse_metadata",
    "parse_metadata_file",
    "save_autocorrelation_csv",
    "save_fit",
    "save_fit_text",
    "save_matrix",
    "save_matrix_csv",
    "save_matrix_set",
    "save_partitioned_matrix_sets",
]
