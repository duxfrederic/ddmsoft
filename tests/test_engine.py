from __future__ import annotations

import numpy as np
import pytest

from ddmsoft.engine import (
    ComputationCancelled,
    VideoReadError,
    compute_ddm,
    log_spaced_lags,
)
from ddmsoft.io import save_matrix_set
from ddmsoft.models import DDMData
from .fixtures import constant_frames, direct_ddm_reference, random_frames


def _radial_reference(spectrum: np.ndarray) -> np.ndarray:
    size = spectrum.shape[0]
    frequencies = np.fft.fftfreq(size)
    distances = np.sqrt(frequencies[:, None] ** 2 + frequencies[None, :] ** 2)
    bins = np.arange(size // 2 + 1, dtype=float) / size
    return np.histogram(distances, bins, weights=spectrum)[0] / np.histogram(distances, bins)[0]


def test_short_videos_always_have_valid_lag_one():
    for count in range(2, 10):
        lags = log_spaced_lags(count)
        assert np.array_equal(lags, np.arange(1, count))


def test_invalid_lag_inputs_are_rejected():
    with pytest.raises(ValueError, match="at least two"):
        log_spaced_lags(1)
    with pytest.raises(ValueError, match="points_per_decade"):
        log_spaced_lags(4, 0)
    with pytest.raises(ValueError, match="integer"):
        log_spaced_lags(4, 1.5)


def test_constant_frames_produce_zero_ddm():
    result = compute_ddm(constant_frames(5, (8, 8)), 10, 1e-6)
    assert isinstance(result, DDMData)
    assert np.allclose(result.matrix, 0)
    assert np.array_equal(result.lag_times, log_spaced_lags(5) / 10)


def test_isotropic_output_matches_independent_full_fft_reference():
    frames = random_frames(6, (8, 8), seed=42)
    result = compute_ddm(frames, 20, 2e-6, points_per_decade=1, max_couples=0)
    assert isinstance(result, DDMData)
    full_power = direct_ddm_reference(frames, result.lag_times * 20)
    expected = np.vstack([_radial_reference(power) for power in full_power])[:, :4]
    assert np.allclose(result.matrix, expected)
    assert np.all(result.q_values > 0)


def test_rectangular_frames_are_rejected_before_fft():
    with pytest.raises(ValueError, match="square"):
        compute_ddm(random_frames(3, (4, 6)), 10, 1e-6)


def test_progress_and_cancellation_stop_without_returning_data():
    events: list[tuple[str, int, int]] = []
    cancelled = False

    def cancel() -> bool:
        return cancelled

    result = compute_ddm(
        random_frames(3, (4, 4)),
        10,
        1e-6,
        progress=lambda stage, completed, total: events.append((stage, completed, total)),
        cancel=cancel,
    )
    assert isinstance(result, DDMData)
    assert {event[0] for event in events} == {"frame_read", "frame_fft", "lag_average"}

    should_cancel = False

    def stop_after_fft(stage: str, completed: int, total: int) -> None:
        nonlocal should_cancel
        events.append((stage, completed, total))
        if stage == "frame_fft":
            should_cancel = True

    with pytest.raises(ComputationCancelled):
        compute_ddm(
            random_frames(3, (4, 4)),
            10,
            1e-6,
            progress=stop_after_fft,
            cancel=lambda: should_cancel,
        )


def test_directional_output_has_one_dataset_per_sector():
    frames = random_frames(4, (8, 8), seed=8)
    result = compute_ddm(frames, 10, 1e-6, sectors=2)
    assert isinstance(result, tuple)
    assert len(result) == 2
    assert all(isinstance(dataset, DDMData) for dataset in result)
    assert all(dataset.matrix.shape == (result[0].matrix.shape[0], 4) for dataset in result)


def test_atomic_matrix_set_save_uses_legacy_names(tmp_path):
    data = DDMData(np.ones((2, 2)), np.array([0.1, 0.2]), np.array([1.0, 2.0]))
    paths = save_matrix_set(tmp_path / "run", data)
    assert [path.name for path in paths] == ["run_DDM_matrix.npy", "run_deltaTs.npy", "run_QS.npy"]
    assert all(path.is_file() for path in paths)


def test_opencv_reader_releases_capture_on_failure(monkeypatch, tmp_path):
    import ddmsoft.engine as engine

    class FakeCapture:
        released = False

        def __init__(self, _path):
            pass

        def isOpened(self):
            return True

        def read(self):
            return True, np.zeros((3, 3), dtype=np.uint8)

        def release(self):
            self.released = True

    capture = FakeCapture("unused")
    monkeypatch.setattr(engine.cv2, "VideoCapture", lambda _path: capture)
    with pytest.raises(VideoReadError, match="inconsistent"):
        iterator = iter(engine.read_video_frames(tmp_path / "video.avi"))
        next(iterator)
        capture.read = lambda: (True, np.zeros((2, 3), dtype=np.uint8))
        next(iterator)
    assert capture.released
