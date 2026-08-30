import numpy as np
import pytest

from ddmsoft.models import DDMData, FitRange, FitRequest, FitResult, VideoMetadata


def test_ddm_data_owns_arrays_and_validates_axes(tmp_path):
    matrix = np.ones((2, 3))
    data = DDMData(matrix, [0.1, 0.2], [1.0, 2.0, 3.0], tmp_path / "sample.avi")
    matrix[0, 0] = 99
    assert data.matrix[0, 0] == 1
    assert data.source == tmp_path / "sample.avi"
    assert np.array_equal(data.dts, data.lag_times)
    with pytest.raises(ValueError, match="one value per matrix row"):
        DDMData(np.ones((2, 3)), [0.1], [1.0, 2.0, 3.0])


def test_contracts_are_qt_free_and_fit_range_is_inclusive():
    fit_range = FitRange(q_min=1, q_max=3, time_min=0, time_max=4)
    assert fit_range.to_slices() == (slice(1, 4), slice(0, 5))
    request = FitRequest("stretch", (1.0, None), (False, True), fit_range)
    assert request.range is fit_range
    assert request.initial_values == (1.0, None)
    assert request.fixed_flags == (False, True)
    with pytest.raises(ValueError, match="q_min"):
        FitRange(3, 1, 0, 1)


def test_fit_result_has_one_status_per_q():
    q_values = np.array([1.0, 2.0])
    result = FitResult(
        "stretch",
        q_values,
        np.ones(2),
        np.zeros(2),
        (np.ones(2), np.ones(2)),
        np.ones((3, 2)),
        np.ones((3, 2)),
        (True, False),
        ("ok", "failed"),
    )
    assert result.fitted_q_values is not q_values
    assert result.converged == (True, False)
    with pytest.raises(ValueError, match="convergence_status"):
        FitResult(
            "stretch",
            [1.0, 2.0],
            [1.0, 1.0],
            [0.0, 0.0],
            (),
            np.ones((1, 2)),
            np.ones((1, 2)),
            (True,),
        )


def test_video_metadata_normalizes_path(tmp_path):
    metadata = VideoMetadata(tmp_path / "movie.avi", 100.0, 1.0e-7)
    assert metadata.path == tmp_path / "movie.avi"
    with pytest.raises(ValueError, match="frame_rate"):
        VideoMetadata(tmp_path / "movie.avi", 0, 1.0e-7)
