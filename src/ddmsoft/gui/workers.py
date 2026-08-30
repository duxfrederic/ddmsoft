"""Background computation workers for the Qt application."""

from __future__ import annotations

import traceback
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Event

from PySide6.QtCore import QObject, Signal, Slot

from ..contin import CONTINCancelled, CONTINResult, run_contin_scan
from ..engine import ComputationCancelled, compute_video_ddm, read_video_frames
from ..fitting import FitCancelled, fit_ddm, get_model
from ..io import (
    LEGACY_SUFFIXES,
    save_fit_text,
    save_matrix_set,
    save_partitioned_matrix_sets,
)
from ..media import FFmpegCancelled, concatenate_videos
from ..models import DDMData, FitRange, FitRequest, FitResult, VideoMetadata
from ..time_dependent import compute_time_dependent_ddm

ProgressCallback = Callable[[str, int, int], None]
CancelCallback = Callable[[], bool]
WorkFunction = Callable[[ProgressCallback, CancelCallback], object]


@dataclass(frozen=True)
class VideoComputationRequest:
    """Plain inputs for one sequential video-computation job."""

    videos: tuple[VideoMetadata, ...]
    max_couples: int
    points_per_decade: float
    sectors: int
    recompute: bool


@dataclass(frozen=True)
class ComputationResult:
    """Paths produced by a completed computation job."""

    paths: tuple[Path, ...]
    processed_videos: tuple[Path, ...]
    kept_videos: tuple[Path, ...]


@dataclass(frozen=True)
class FitComputationRequest:
    """Plain inputs for one matrix-fitting job."""

    matrix_path: Path
    data: DDMData
    fit_request: FitRequest


@dataclass(frozen=True)
class FitComputationResult:
    """A fit result retaining the matrix identity and request that produced it."""

    matrix_path: Path
    fit_request: FitRequest
    fit: FitResult


@dataclass(frozen=True)
class BatchFitRequest:
    """Inputs for fitting and exporting a selected set of matrices."""

    matrices: tuple[tuple[Path, DDMData], ...]
    fit_request: FitRequest
    output_prefix: Path
    continue_on_failure: bool = True
    overwrite: bool = False
    viscosity: float | None = None
    temperature: float | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "matrices",
            tuple((Path(path), data) for path, data in self.matrices),
        )
        object.__setattr__(self, "output_prefix", Path(self.output_prefix))


@dataclass(frozen=True)
class BatchFitResult:
    """Completed batch fits and per-matrix failures."""

    fits: tuple[FitComputationResult, ...]
    output_paths: tuple[Path, ...]
    failures: tuple[tuple[Path, str], ...]


@dataclass(frozen=True)
class TimeDependentComputationRequest:
    """Inputs for sequential partitioned DDM computation."""

    videos: tuple[VideoMetadata, ...]
    partitions: int
    max_couples: int
    points_per_decade: float
    sectors: int
    recompute: bool = False


@dataclass(frozen=True)
class TimeDependentComputationResult:
    """Partitioned matrix paths produced for each source video."""

    paths: tuple[Path, ...]
    processed_videos: tuple[Path, ...]


@dataclass(frozen=True)
class CONTINComputationRequest:
    """Inputs for one q-index CONTIN scan."""

    matrix_path: Path
    data: DDMData
    q_index: int
    fit_range: FitRange
    gamma_min: float
    gamma_max: float
    gamma_count: int
    alpha_min: float
    alpha_max: float
    alpha_count: int
    maxiter: int


@dataclass(frozen=True)
class CONTINComputationResult:
    """A CONTIN result retaining its matrix and q identities."""

    matrix_path: Path
    q_index: int
    fit_range: FitRange
    result: CONTINResult


@dataclass(frozen=True)
class VideoConcatenationRequest:
    """Inputs for one ffmpeg concat job."""

    videos: tuple[Path, ...]
    output: Path
    overwrite: bool = False


@dataclass(frozen=True)
class VideoConcatenationResult:
    """Output path from a completed video concat job."""

    output: Path


class ComputationWorker(QObject):
    """Run one plain callable on a ``QThread`` with cooperative cancellation."""

    progress = Signal(str, int, int)
    status = Signal(str)
    result = Signal(object)
    failure = Signal(str, str)
    cancelled = Signal()
    finished = Signal()

    def __init__(self, work: WorkFunction) -> None:
        super().__init__()
        self._work = work
        self._cancel_event = Event()

    def request_cancel(self) -> None:
        """Request cancellation without requiring a queued worker call."""
        self._cancel_event.set()

    def is_cancelled(self) -> bool:
        return self._cancel_event.is_set()

    @Slot()
    def run(self) -> None:
        try:
            value = self._work(self._report, self.is_cancelled)
            self.result.emit(value)
        except (ComputationCancelled, CONTINCancelled, FFmpegCancelled, FitCancelled):
            self.cancelled.emit()
        except Exception as error:  # noqa: BLE001 - workers must relay all failures
            message = str(error) or error.__class__.__name__
            self.failure.emit(message, traceback.format_exc())
        finally:
            self.finished.emit()

    def _report(self, stage: str, completed: int, total: int) -> None:
        self.status.emit(_stage_label(stage))
        self.progress.emit(stage, completed, total)


def run_video_computation(
    request: VideoComputationRequest,
    progress: ProgressCallback,
    cancel: CancelCallback,
) -> ComputationResult:
    """Compute and transactionally save videos sequentially."""
    if not request.videos:
        raise ValueError("at least one video is required")
    if request.sectors < 1:
        raise ValueError("sectors must be positive")

    total = len(request.videos) * 1_000
    paths: list[Path] = []
    processed: list[Path] = []
    kept: list[Path] = []
    for video_index, metadata in enumerate(request.videos):
        base = video_index * 1_000
        prefixes = _output_prefixes(metadata.path, request.sectors)
        targets = tuple(path for prefix in prefixes for path in _matrix_paths(prefix))
        existing = tuple(path for path in targets if path.exists())
        if not request.recompute and existing:
            if len(existing) != len(targets):
                raise FileExistsError(
                    f"incomplete output exists for {metadata.path}; select recompute to replace it"
                )
            progress("keeping", base + 1_000, total)
            paths.extend(targets)
            kept.append(metadata.path)
            continue

        _check_cancel(cancel)
        progress("loading", base + 10, total)

        def engine_progress(
            stage: str, completed: int, stage_total: int, video_base: int = base
        ) -> None:
            local = _stage_progress(stage, completed, stage_total)
            progress(stage, video_base + local, total)

        result = compute_video_ddm(
            metadata.path,
            metadata.frame_rate,
            metadata.pixel_size,
            max_couples=request.max_couples,
            points_per_decade=request.points_per_decade,
            sectors=request.sectors,
            progress=engine_progress,
            cancel=cancel,
        )
        _check_cancel(cancel)
        datasets = (result,) if isinstance(result, DDMData) else tuple(result)
        if len(datasets) != len(prefixes):
            raise ValueError(
                f"computation returned {len(datasets)} dataset(s), expected {len(prefixes)}"
            )

        output_directory = metadata.path.parent / "ddm_matrices"
        output_directory.mkdir(parents=True, exist_ok=True)
        progress("saving", base + 900, total)
        _check_cancel(cancel)
        with TemporaryDirectory(prefix=".ddmsoft-stage-", dir=output_directory) as staging:
            staged_paths: list[tuple[Path, Path]] = []
            for prefix, dataset in zip(prefixes, datasets):
                staged_prefix = Path(staging) / prefix.name
                staged = _matrix_paths(save_matrix_set(staged_prefix, dataset))
                staged_paths.extend(zip(staged, _matrix_paths(prefix)))
            _check_cancel(cancel)
            _commit_staged(staged_paths, recompute=request.recompute)
        paths.extend(targets)
        processed.append(metadata.path)
        progress("saving", base + 990, total)

    progress("complete", total, total)
    return ComputationResult(tuple(paths), tuple(processed), tuple(kept))


def run_fit(
    request: FitComputationRequest,
    progress: ProgressCallback,
    cancel: CancelCallback,
) -> FitComputationResult:
    """Fit a selected matrix while retaining its stable catalog identity."""
    if not request.matrix_path:
        raise ValueError("matrix_path is required")
    result = fit_ddm(
        request.data,
        request.fit_request,
        progress=lambda completed, total: progress("fitting", completed, total),
        cancel=cancel,
    )
    return FitComputationResult(request.matrix_path, request.fit_request, result)


def run_batch_fit(
    request: BatchFitRequest,
    progress: ProgressCallback,
    cancel: CancelCallback,
) -> BatchFitResult:
    """Fit and export matrices sequentially with explicit failure policy."""
    if not request.matrices:
        raise ValueError("at least one matrix is required")
    model = get_model(request.fit_request.model_id)
    targets = batch_fit_target_paths(
        request.output_prefix, (path for path, _ in request.matrices)
    )
    existing = tuple(path for path in targets if path.exists())
    if existing and not request.overwrite:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"batch fit output already exists: {names}")
    request.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    results: list[FitComputationResult] = []
    output_paths: list[Path] = []
    failures: list[tuple[Path, str]] = []
    total = len(request.matrices)
    for index, (matrix_path, data) in enumerate(request.matrices):
        if cancel():
            raise FitCancelled("batch fit cancelled")
        progress(f"batch fitting {matrix_path.name}", index, total)
        try:
            fit = fit_ddm(data, request.fit_request, cancel=cancel)
            target = targets[index]
            save_fit_text(
                target,
                fit.q_values,
                fit.amplitude,
                fit.noise,
                fit.model_parameters,
                model.parameter_names[:-2],
                viscosity=request.viscosity,
                temperature=request.temperature,
            )
        except FitCancelled:
            raise
        except Exception as error:  # noqa: BLE001 - batch policy reports each matrix failure
            failures.append((matrix_path, str(error) or error.__class__.__name__))
            progress(f"batch failed {matrix_path.name}", index + 1, total)
            if not request.continue_on_failure:
                break
            continue
        results.append(FitComputationResult(matrix_path, request.fit_request, fit))
        output_paths.append(target)
        progress(f"batch complete {matrix_path.name}", index + 1, total)
    progress("complete", total, total)
    return BatchFitResult(tuple(results), tuple(output_paths), tuple(failures))


def run_time_dependent_computation(
    request: TimeDependentComputationRequest,
    progress: ProgressCallback,
    cancel: CancelCallback,
) -> TimeDependentComputationResult:
    """Compute and save every requested source-frame partition sequentially."""
    if not request.videos:
        raise ValueError("at least one video is required")
    if request.partitions < 1:
        raise ValueError("partitions must be positive")
    total = len(request.videos) * request.partitions
    paths: list[Path] = []
    processed: list[Path] = []
    for video_index, metadata in enumerate(request.videos):
        offset = video_index * request.partitions

        def partition_progress(
            stage: str,
            completed: int,
            partition_total: int,
            video_name: str = metadata.path.name,
            video_offset: int = offset,
        ) -> None:
            progress(f"time dependent {video_name}", video_offset + completed, total)

        datasets = compute_time_dependent_ddm(
            read_video_frames(metadata.path),
            metadata.frame_rate,
            metadata.pixel_size,
            request.partitions,
            max_couples=request.max_couples,
            points_per_decade=request.points_per_decade,
            sectors=request.sectors,
            progress=partition_progress,
            cancel=cancel,
        )
        if cancel():
            raise ComputationCancelled("time-dependent DDM cancelled before saving")
        output_prefix = metadata.path.parent / "ddm_matrices" / metadata.path.stem
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(prefix=".ddmsoft-stage-", dir=output_prefix.parent) as staging:
            staged_paths = save_partitioned_matrix_sets(Path(staging) / output_prefix.name, datasets)
            if cancel():
                raise ComputationCancelled("time-dependent DDM cancelled before committing")
            staged_targets = tuple(
                (staged_path, output_prefix.parent / staged_path.name)
                for staged_path in staged_paths
            )
            existing = tuple(target for _, target in staged_targets if target.exists())
            if existing and not request.recompute:
                names = ", ".join(str(path) for path in existing)
                raise FileExistsError(
                    f"time-dependent output already exists: {names}; select recompute to replace it"
                )
            current_targets = {target for _, target in staged_targets}
            stale_paths = (
                tuple(
                    path
                    for path in output_prefix.parent.glob("*.npy")
                    if path.name.startswith(f"{output_prefix.name}__i=")
                    and path not in current_targets
                )
                if request.recompute
                else ()
            )
            _commit_staged(
                staged_targets,
                recompute=request.recompute,
                stale_paths=stale_paths,
            )
            committed_paths = [target for _, target in staged_targets]
        paths.extend(committed_paths)
        processed.append(metadata.path)
        progress(f"time dependent {metadata.path.name}", offset + request.partitions, total)
    progress("complete", total, total)
    return TimeDependentComputationResult(tuple(paths), tuple(processed))


def run_contin_computation(
    request: CONTINComputationRequest,
    progress: ProgressCallback,
    cancel: CancelCallback,
) -> CONTINComputationResult:
    """Run CONTIN for one inclusive time range and q index."""
    if request.q_index >= request.data.q_values.size:
        raise ValueError("CONTIN q index exceeds available q values")
    _, time_slice = request.fit_range.to_slices()
    if request.q_index < 0:
        raise ValueError("CONTIN q index must be non-negative")
    q_value = request.data.q_values[request.q_index]
    tau = request.data.lag_times[time_slice] * q_value**2
    ddmdata = request.data.matrix[time_slice, request.q_index]
    result = run_contin_scan(
        tau,
        ddmdata,
        request.gamma_min,
        request.gamma_max,
        request.gamma_count,
        request.alpha_min,
        request.alpha_max,
        request.alpha_count,
        maxiter=request.maxiter,
        progress=lambda completed, total: progress("contin", completed, total),
        cancel=cancel,
    )
    return CONTINComputationResult(request.matrix_path, request.q_index, request.fit_range, result)


def run_video_concatenation(
    request: VideoConcatenationRequest,
    progress: ProgressCallback,
    cancel: CancelCallback,
) -> VideoConcatenationResult:
    """Run safe ffmpeg concatenation in the shared worker abstraction."""
    if cancel():
        raise ComputationCancelled("video concatenation cancelled")
    progress("concatenating", 0, 1)
    output = concatenate_videos(
        request.videos,
        request.output,
        overwrite=request.overwrite,
        cancel=cancel,
    )
    progress("concatenating", 1, 1)
    return VideoConcatenationResult(output)


def batch_fit_target_paths(
    output_prefix: str | Path, matrix_paths: Iterable[str | Path]
) -> tuple[Path, ...]:
    """Return deterministic, collision-checked text targets for matrix paths."""
    prefix = Path(output_prefix)
    if prefix.name.endswith(".txt"):
        prefix = prefix.with_name(prefix.name[:-4])
    paths = tuple(
        prefix.with_name(f"{prefix.name}_{_matrix_stem(path)}.txt") for path in matrix_paths
    )
    if len(set(paths)) != len(paths):
        raise ValueError("batch matrix names produce duplicate fit output paths")
    return paths


def _matrix_stem(path: str | Path) -> str:
    stem = Path(path).stem
    return stem.removesuffix("_DDM_matrix")


def _output_prefixes(path: Path, sectors: int) -> tuple[Path, ...]:
    prefix = path.parent / "ddm_matrices" / path.stem
    if sectors == 1:
        return (prefix,)
    return tuple(
        prefix.with_name(f"{prefix.name}_{index * 180.0 / sectors:.1f}_")
        for index in range(sectors)
    )


def _matrix_paths(prefix: Path | Sequence[Path]) -> tuple[Path, ...]:
    if isinstance(prefix, Path):
        return tuple(prefix.with_name(prefix.name + suffix) for suffix in LEGACY_SUFFIXES)
    return tuple(Path(path) for path in prefix)


def _commit_staged(
    staged_paths: Sequence[tuple[Path, Path]],
    *,
    recompute: bool,
    stale_paths: Sequence[Path] = (),
) -> None:
    committed: list[Path] = []
    backups: list[tuple[Path, Path]] = []
    try:
        for index, target in enumerate(stale_paths):
            backup = staged_paths[0][0].parent / f".stale-{index}-{target.name}"
            target.replace(backup)
            backups.append((backup, target))
        for index, (staged, target) in enumerate(staged_paths):
            backup = staged.parent / f".backup-{index}-{target.name}"
            if target.exists():
                if not recompute:
                    raise FileExistsError(f"output already exists: {target}")
                target.replace(backup)
                backups.append((backup, target))
            staged.replace(target)
            committed.append(target)
    except Exception:
        for target in committed:
            target.unlink(missing_ok=True)
        for backup, target in reversed(backups):
            if backup.exists():
                backup.replace(target)
        raise


def _check_cancel(cancel: CancelCallback) -> None:
    if cancel():
        raise ComputationCancelled("DDM computation cancelled")


def _stage_progress(stage: str, completed: int, total: int) -> int:
    if stage == "frame_read":
        return 50
    if stage == "frame_fft":
        return 100 + _scaled_progress(completed, total, 300)
    if stage == "lag_average":
        return 400 + _scaled_progress(completed, total, 500)
    return 100


def _scaled_progress(completed: int, total: int, span: int) -> int:
    if total <= 0:
        return 0
    return max(0, min(span, int(span * completed / total)))


def _stage_label(stage: str) -> str:
    if stage.startswith("time dependent "):
        return f"Time-dependent DDM: {stage.removeprefix('time dependent ')}"
    if stage.startswith("batch fitting "):
        return f"Fitting {stage.removeprefix('batch fitting ')}"
    if stage.startswith("batch complete "):
        return f"Saved fit for {stage.removeprefix('batch complete ')}"
    if stage.startswith("batch failed "):
        return f"Fit failed for {stage.removeprefix('batch failed ')}"
    return {
        "loading": "Loading frames",
        "frame_read": "Loading frames",
        "frame_fft": "Computing FFT",
        "lag_average": "Temporal averaging",
        "saving": "Saving DDM matrices",
        "keeping": "Keeping existing matrices",
        "fitting": "Fitting correlation curves",
        "contin": "Computing CONTIN candidates",
        "concatenating": "Concatenating videos",
        "complete": "Computation complete",
    }.get(stage, stage.replace("_", " ").capitalize())


__all__ = [
    "BatchFitRequest",
    "BatchFitResult",
    "CONTINComputationRequest",
    "CONTINComputationResult",
    "ComputationResult",
    "ComputationWorker",
    "FitComputationRequest",
    "FitComputationResult",
    "TimeDependentComputationRequest",
    "TimeDependentComputationResult",
    "VideoComputationRequest",
    "VideoConcatenationRequest",
    "VideoConcatenationResult",
    "batch_fit_target_paths",
    "run_batch_fit",
    "run_contin_computation",
    "run_fit",
    "run_time_dependent_computation",
    "run_video_computation",
    "run_video_concatenation",
]
