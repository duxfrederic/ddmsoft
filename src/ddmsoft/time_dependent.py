"""Time-dependent DDM partitioning and sequential computation."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .engine import ComputationCancelled, compute_ddm
from .models import DDMData, VideoMetadata


@dataclass(frozen=True)
class PartitionedDDM:
    """One time-dependent result and its actual source-frame interval."""

    start_frame: int
    stop_frame: int
    data: DDMData | tuple[DDMData, ...]

    @property
    def name(self) -> str:
        return f"__i={self.start_frame}__"


def partition_frame_ranges(frame_count: int, partitions: int) -> tuple[tuple[int, int], ...]:
    """Return array-split-equivalent half-open ranges covering every frame."""
    if isinstance(frame_count, bool) or not isinstance(frame_count, (int, np.integer)):
        raise ValueError("frame_count must be an integer")
    if isinstance(partitions, bool) or not isinstance(partitions, (int, np.integer)):
        raise ValueError("partitions must be an integer")
    frame_count = int(frame_count)
    partitions = int(partitions)
    if frame_count < 1:
        raise ValueError("frame_count must be positive")
    if partitions < 1:
        raise ValueError("partitions must be positive")
    if partitions > frame_count:
        raise ValueError("partitions cannot exceed frame_count")
    base_size, remainder = divmod(frame_count, partitions)
    ranges = []
    start = 0
    for index in range(partitions):
        stop = start + base_size + (index < remainder)
        ranges.append((start, stop))
        start = stop
    return tuple(ranges)


def compute_time_dependent_ddm(
    frames: Iterable[np.ndarray],
    frame_rate: float,
    pixel_size: float,
    partitions: int,
    *,
    max_couples: int = 300,
    points_per_decade: int | float = 20,
    sectors: int = 1,
    progress: Callable[[str, int, int], None] | None = None,
    cancel: Callable[[], bool] | object | None = None,
) -> tuple[PartitionedDDM, ...]:
    """Compute each partition sequentially, preserving actual frame starts."""
    source = list(frames)
    ranges = partition_frame_ranges(len(source), partitions)
    results: list[PartitionedDDM] = []
    for index, (start, stop) in enumerate(ranges):
        if _cancel_requested(cancel):
            raise ComputationCancelled("time-dependent DDM cancelled before a partition")
        subset = source[start:stop]
        if len(subset) < 2:
            raise ValueError(
                f"partition {index} contains {len(subset)} frame; at least two are required per DDM"
            )
        result = compute_ddm(
            subset,
            frame_rate,
            pixel_size,
            max_couples=max_couples,
            points_per_decade=points_per_decade,
            sectors=sectors,
            cancel=cancel,
        )
        results.append(PartitionedDDM(start, stop, result))
        if progress is not None:
            progress("partition", index + 1, len(ranges))
    return tuple(results)


def compute_time_dependent_videos(
    videos: Sequence[VideoMetadata],
    reader: Callable[[Path], Iterable[np.ndarray]],
    partitions: int,
    **kwargs: object,
) -> dict[Path, tuple[PartitionedDDM, ...]]:
    """Process each video with its own metadata rather than shared parameters."""
    return {
        metadata.path: compute_time_dependent_ddm(
            reader(metadata.path),
            metadata.frame_rate,
            metadata.pixel_size,
            partitions,
            **kwargs,
        )
        for metadata in videos
    }


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
    "PartitionedDDM",
    "compute_time_dependent_ddm",
    "compute_time_dependent_videos",
    "partition_frame_ranges",
]
