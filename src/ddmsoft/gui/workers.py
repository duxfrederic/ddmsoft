"""Background computation workers for the Qt application."""

from __future__ import annotations

import traceback
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Event

from PySide6.QtCore import QObject, Signal, Slot

from ..engine import ComputationCancelled, compute_video_ddm
from ..io import LEGACY_SUFFIXES, save_matrix_set
from ..models import DDMData, VideoMetadata

ProgressCallback = Callable[[str, int, int], None]
CancelCallback = Callable[[], bool]
WorkFunction = Callable[[ProgressCallback, CancelCallback], object]


@dataclass(frozen=True)
class VideoComputationRequest:
    """Plain inputs for one sequential video-computation job."""

    videos: tuple[VideoMetadata, ...]
    max_couples: int
    points_per_decade: int | float
    sectors: int
    recompute: bool


@dataclass(frozen=True)
class ComputationResult:
    """Paths produced by a completed computation job."""

    paths: tuple[Path, ...]
    processed_videos: tuple[Path, ...]
    kept_videos: tuple[Path, ...]


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
            if self.is_cancelled():
                self.cancelled.emit()
            else:
                self.result.emit(value)
        except ComputationCancelled:
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


def _commit_staged(staged_paths: Sequence[tuple[Path, Path]], *, recompute: bool) -> None:
    committed: list[Path] = []
    backups: list[tuple[Path, Path]] = []
    try:
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
    return {
        "loading": "Loading frames",
        "frame_read": "Loading frames",
        "frame_fft": "Computing FFT",
        "lag_average": "Temporal averaging",
        "saving": "Saving DDM matrices",
        "keeping": "Keeping existing matrices",
        "complete": "Computation complete",
    }.get(stage, stage.replace("_", " ").capitalize())


__all__ = [
    "ComputationResult",
    "ComputationWorker",
    "VideoComputationRequest",
    "run_video_computation",
]
