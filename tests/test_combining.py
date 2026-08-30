from __future__ import annotations

import numpy as np
import pytest

from ddmsoft.combining import CombinationError, average_ddm, average_groups, merge_ddm
from ddmsoft.io import save_partitioned_matrix_sets
from ddmsoft.models import DDMData, VideoMetadata
from ddmsoft.time_dependent import (
    compute_time_dependent_ddm,
    compute_time_dependent_videos,
    partition_frame_ranges,
)
from .fixtures import constant_frames


def _data(times, values, q=(1.0, 2.0)):
    matrix = np.asarray(values, dtype=float)
    return DDMData(matrix, np.asarray(times, dtype=float), np.asarray(q, dtype=float))


def test_merge_and_average_do_not_mutate_inputs():
    fast = _data([1.0, 2.0, 3.0, 4.0], [[1, 2], [2, 3], [3, 4], [4, 5]])
    slow = _data([3.0, 5.0, 7.0], [[3, 6], [4, 7], [5, 8]])
    fast_before = fast.matrix.copy()
    slow_before = slow.matrix.copy()
    merged = merge_ddm([fast, slow])
    averaged = average_ddm([fast, fast])
    assert np.array_equal(fast.matrix, fast_before)
    assert np.array_equal(slow.matrix, slow_before)
    assert merged.matrix.shape[1] == 2
    assert np.array_equal(averaged.matrix, fast.matrix)


def test_average_groups_by_lag_grid_and_names_results():
    first = _data([1.0, 2.0], [[1, 2], [3, 4]])
    second = _data([1.0, 2.0], [[3, 4], [5, 6]])
    third = _data([1.0, 3.0], [[5, 6], [7, 8]])
    groups = average_groups([first, third, second])
    assert [group.name for group in groups] == ["average_0", "average_1"]
    assert np.array_equal(groups[0].data.matrix, np.array([[2, 3], [4, 5]]))
    assert np.array_equal(groups[1].data.matrix, third.matrix)


def test_combination_rejects_incompatible_axes_and_zero_scaling():
    first = _data([1.0, 2.0], [[1, 2], [3, 4]])
    with pytest.raises(CombinationError, match="q grid"):
        average_ddm([first, _data([1.0, 2.0], [[1, 2], [3, 4]], q=(1.0, 3.0))])
    with pytest.raises(CombinationError, match="strictly increasing"):
        average_ddm([first, _data([1.0, 1.0], [[1, 2], [3, 4]])])
    with pytest.raises(CombinationError, match="zero"):
        merge_ddm(
            [
                _data([1.0, 2.0, 3.0], [[1, 1], [2, 2], [3, 3]]),
                _data([2.0, 4.0], [[0, 1], [2, 3]]),
            ]
        )


def test_partition_ranges_account_for_every_frame():
    ranges = partition_frame_ranges(10, 3)
    assert ranges == ((0, 4), (4, 7), (7, 10))
    assert [frame for start, stop in ranges for frame in range(start, stop)] == list(range(10))
    with pytest.raises(ValueError, match="exceed"):
        partition_frame_ranges(2, 3)


def test_time_dependent_results_use_actual_starts():
    frames = constant_frames(10, (4, 4))
    results = compute_time_dependent_ddm(frames, 10, 1e-6, 3, points_per_decade=1)
    assert [(result.start_frame, result.stop_frame) for result in results] == [
        (0, 4),
        (4, 7),
        (7, 10),
    ]
    assert [result.name for result in results] == ["__i=0__", "__i=4__", "__i=7__"]


def test_time_dependent_video_processing_uses_each_video_metadata(tmp_path):
    first_path = tmp_path / "first.avi"
    second_path = tmp_path / "second.avi"
    source = {
        first_path: constant_frames(6, (4, 4)),
        second_path: constant_frames(6, (4, 4)),
    }
    metadata = (
        VideoMetadata(first_path, 10.0, 1e-6),
        VideoMetadata(second_path, 20.0, 1e-6),
    )
    results = compute_time_dependent_videos(
        metadata, source.__getitem__, 2, points_per_decade=1
    )
    assert set(results) == {first_path, second_path}
    assert results[first_path][0].data.lag_times[0] == 0.1
    assert results[second_path][0].data.lag_times[0] == 0.05


def test_partitioned_results_save_with_actual_start_tokens(tmp_path):
    results = compute_time_dependent_ddm(
        constant_frames(6, (4, 4)), 10, 1e-6, 2, points_per_decade=1
    )
    paths = save_partitioned_matrix_sets(tmp_path / "sample", results)
    assert paths[0].name == "sample__i=0___DDM_matrix.npy"
    assert paths[3].name == "sample__i=3___DDM_matrix.npy"
