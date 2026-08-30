from __future__ import annotations

import numpy as np
import pytest

from ddmsoft.engine import ComputationCancelled
from ddmsoft.models import DDMData
from ddmsoft.time_dependent import (
    PartitionedDDM,
    compute_time_dependent_ddm,
    partition_frame_ranges,
)


def test_partition_frame_ranges_cover_every_source_frame_without_empty_ranges():
    ranges = partition_frame_ranges(10, 3)

    assert ranges == ((0, 4), (4, 7), (7, 10))
    assert [frame for start, stop in ranges for frame in range(start, stop)] == list(range(10))


@pytest.mark.parametrize(
    "frame_count, partitions, message",
    [(0, 1, "positive"), (3, 0, "positive"), (3, 4, "cannot exceed")],
)
def test_partition_frame_ranges_reject_invalid_counts(frame_count, partitions, message):
    with pytest.raises(ValueError, match=message):
        partition_frame_ranges(frame_count, partitions)


def test_time_dependent_ddm_preserves_actual_partition_starts_and_reports_progress(monkeypatch):
    frames = [np.full((2, 2), value, dtype=float) for value in range(10)]
    calls: list[tuple[int, ...]] = []
    progress: list[tuple[str, int, int]] = []

    def fake_compute(subset, *args, **kwargs):
        calls.append(tuple(int(frame[0, 0]) for frame in subset))
        return DDMData(np.ones((2, 1)), np.array([0.1, 0.2]), np.array([1.0]))

    monkeypatch.setattr("ddmsoft.time_dependent.compute_ddm", fake_compute)
    results = compute_time_dependent_ddm(
        frames,
        30.0,
        1e-6,
        3,
        progress=lambda *event: progress.append(event),
    )

    assert calls == [(0, 1, 2, 3), (4, 5, 6), (7, 8, 9)]
    assert [(item.start_frame, item.stop_frame) for item in results] == [
        (0, 4),
        (4, 7),
        (7, 10),
    ]
    assert progress == [("partition", 1, 3), ("partition", 2, 3), ("partition", 3, 3)]


def test_time_dependent_ddm_cancels_before_next_partition(monkeypatch):
    calls = 0

    def fake_compute(*args, **kwargs):
        nonlocal calls
        calls += 1
        return DDMData(np.ones((2, 1)), np.array([0.1, 0.2]), np.array([1.0]))

    monkeypatch.setattr("ddmsoft.time_dependent.compute_ddm", fake_compute)
    with pytest.raises(ComputationCancelled):
        compute_time_dependent_ddm(
            [np.zeros((2, 2)) for _ in range(6)],
            30.0,
            1e-6,
            3,
            cancel=lambda: True,
        )
    assert calls == 0


def test_time_dependent_ddm_cancels_while_materializing_video():
    yielded = 0

    def frames():
        nonlocal yielded
        for _ in range(10):
            yielded += 1
            yield np.zeros((2, 2))

    with pytest.raises(ComputationCancelled, match="reading frames"):
        compute_time_dependent_ddm(
            frames(),
            30.0,
            1e-6,
            2,
            cancel=lambda: yielded >= 2,
        )
    assert yielded == 2


def test_time_dependent_ddm_rejects_short_partitions_before_computation(monkeypatch):
    monkeypatch.setattr(
        "ddmsoft.time_dependent.compute_ddm",
        lambda *args, **kwargs: pytest.fail("no partition should be computed"),
    )

    with pytest.raises(ValueError, match="at least two frames"):
        compute_time_dependent_ddm(
            [np.zeros((2, 2)) for _ in range(5)],
            30.0,
            1e-6,
            3,
        )


def test_partitioned_ddm_name_contains_source_frame_start():
    result = PartitionedDDM(
        7,
        10,
        DDMData(np.ones((2, 1)), np.array([0.1, 0.2]), np.array([1.0])),
    )

    assert result.name == "__i=7__"
