from __future__ import annotations

import time
from pathlib import Path
from threading import Event
from unittest.mock import ANY

import numpy as np
import pytest
from PySide6.QtCore import QThread, QTimer
from PySide6.QtWidgets import QApplication

from ddmsoft.engine import ComputationCancelled
from ddmsoft.fitting import default_fit_request
from ddmsoft.gui.workers import (
    BatchFitRequest,
    ComputationWorker,
    CONTINComputationRequest,
    FitComputationRequest,
    TimeDependentComputationRequest,
    VideoComputationRequest,
    VideoConcatenationRequest,
    batch_fit_target_paths,
    run_batch_fit,
    run_contin_computation,
    run_fit,
    run_time_dependent_computation,
    run_video_computation,
    run_video_concatenation,
)
from ddmsoft.models import DDMData, FitRange, VideoMetadata
from ddmsoft.time_dependent import PartitionedDDM

from .fixtures import generate_model_data


@pytest.fixture
def qapp():
    instance = QApplication.instance() or QApplication([])
    yield instance
    for widget in instance.topLevelWidgets():
        widget.close()


def test_worker_keeps_qt_event_loop_responsive_and_reports_monotonic_progress(qapp):
    progress = []
    statuses = []
    results = []
    heartbeats = []

    def slow_work(report, cancel):
        for index in range(8):
            time.sleep(0.015)
            report("synthetic", index + 1, 8)
        return "finished"

    thread = QThread()
    worker = ComputationWorker(slow_work)
    worker.moveToThread(thread)
    worker.progress.connect(lambda *value: progress.append(value))
    worker.status.connect(statuses.append)
    worker.result.connect(results.append)
    done = False

    def mark_done():
        nonlocal done
        done = True

    thread.started.connect(worker.run)
    worker.finished.connect(thread.quit)
    worker.finished.connect(worker.deleteLater)
    thread.finished.connect(mark_done)
    heartbeat = QTimer()
    heartbeat.timeout.connect(lambda: heartbeats.append(time.monotonic()))
    heartbeat.start(2)
    thread.start()
    _wait_until(qapp, lambda: done)
    heartbeat.stop()
    thread.wait()

    assert len(heartbeats) >= 3
    assert [value[1] for value in progress] == sorted(value[1] for value in progress)
    assert progress[-1] == ("synthetic", 8, 8)
    assert statuses[-1] == "Synthetic"
    assert results == ["finished"]


def test_worker_relays_complete_traceback_on_failure(qapp):
    failures = []

    def fail(report, cancel):
        raise RuntimeError("synthetic worker failure")

    _run_worker(qapp, fail, failure=failures)

    assert failures[0][0] == "synthetic worker failure"
    assert "RuntimeError: synthetic worker failure" in failures[0][1]
    assert "Traceback" in failures[0][1]


def test_worker_cancellation_emits_cancelled_without_result(qapp):
    cancelled = []
    results = []

    def cancellable(report, cancel):
        for index in range(100):
            if cancel():
                raise ComputationCancelled("synthetic cancellation")
            time.sleep(0.005)
            report("synthetic", index + 1, 100)
        return "should not finish"

    _run_worker(qapp, cancellable, cancel_after=25, cancelled=cancelled, results=results)

    assert cancelled == [True]
    assert results == []


def test_video_job_keeps_existing_sets_and_requires_recompute_to_replace(tmp_path, monkeypatch):
    video = tmp_path / "sample.avi"
    video.touch()
    metadata = VideoMetadata(video, 30.0, 1e-6)
    request = VideoComputationRequest((metadata,), 300, 20, 1, False)
    calls = []

    def fake_compute(path, frame_rate, pixel_size, **kwargs):
        calls.append(path)
        return DDMData(np.ones((2, 2)), np.array([0.01, 0.1]), np.array([1e6, 2e6]))

    monkeypatch.setattr("ddmsoft.gui.workers.compute_video_ddm", fake_compute)
    first_progress = []
    first = run_video_computation(
        request, lambda *value: first_progress.append(value), lambda: False
    )
    second_progress = []
    second = run_video_computation(
        request, lambda *value: second_progress.append(value), lambda: False
    )
    replaced_progress = []
    replaced = run_video_computation(
        VideoComputationRequest((metadata,), 300, 20, 1, True),
        lambda *value: replaced_progress.append(value),
        lambda: False,
    )

    assert len(calls) == 2
    assert first.paths == second.paths == replaced.paths
    assert second.kept_videos == (video,)
    assert replaced.processed_videos == (video,)
    assert all(path.exists() for path in first.paths)
    for progress in (first_progress, second_progress, replaced_progress):
        assert [value[1] for value in progress] == sorted(value[1] for value in progress)


def test_fit_job_retains_matrix_identity_and_inclusive_request(tmp_path):
    data = generate_model_data("stretch")
    matrix_path = tmp_path / "first_DDM_matrix.npy"
    fit_request = FitComputationRequest(
        matrix_path,
        data,
        default_fit_request("stretch", FitRange(1, 3, 2, 6)),
    )
    progress = []

    result = run_fit(fit_request, lambda *value: progress.append(value), lambda: False)

    assert result.matrix_path == matrix_path
    assert result.fit_request.fit_range == FitRange(1, 3, 2, 6)
    assert np.array_equal(result.fit.q_values, data.q_values[1:4])
    assert result.fit.correlation.shape == (5, 3)
    assert progress[-1] == ("fitting", 3, 3)


def test_batch_fit_exports_each_matrix_and_prevents_silent_overwrite(tmp_path):
    first = generate_model_data("stretch")
    second = generate_model_data("stretch", noise=0.001)
    request = BatchFitRequest(
        ((tmp_path / "first_DDM_matrix.npy", first), (tmp_path / "second_DDM_matrix.npy", second)),
        default_fit_request("stretch", FitRange(0, 3, 0, 6)),
        tmp_path / "batch",
    )
    progress = []

    result = run_batch_fit(request, lambda *value: progress.append(value), lambda: False)

    assert len(result.fits) == 2
    assert result.failures == ()
    assert result.output_paths == batch_fit_target_paths(
        tmp_path / "batch", (path for path, _ in request.matrices)
    )
    assert all(path.is_file() for path in result.output_paths)
    assert progress[-1] == ("complete", 2, 2)
    with pytest.raises(FileExistsError, match="batch fit output already exists"):
        run_batch_fit(request, lambda *value: None, lambda: False)


def test_batch_fit_can_continue_after_one_matrix_failure(tmp_path, monkeypatch):
    first = generate_model_data("stretch")
    second = generate_model_data("stretch")
    first_path = tmp_path / "bad_DDM_matrix.npy"
    second_path = tmp_path / "good_DDM_matrix.npy"
    request = BatchFitRequest(
        ((first_path, first), (second_path, second)),
        default_fit_request("stretch", FitRange(0, 3, 0, 6)),
        tmp_path / "batch",
    )
    from ddmsoft.gui import workers

    original = workers.fit_ddm

    def fail_first(data, fit_request, **kwargs):
        if data is first:
            raise ValueError("synthetic bad matrix")
        return original(data, fit_request, **kwargs)

    monkeypatch.setattr(workers, "fit_ddm", fail_first)
    result = run_batch_fit(request, lambda *value: None, lambda: False)

    assert result.failures == ((first_path, "synthetic bad matrix"),)
    assert [item.matrix_path for item in result.fits] == [second_path]
    assert result.output_paths == (tmp_path / "batch_good.txt",)


def test_video_job_cancellation_never_commits_partial_outputs(tmp_path, monkeypatch):
    video = tmp_path / "sample.avi"
    video.touch()
    metadata = VideoMetadata(video, 30.0, 1e-6)
    request = VideoComputationRequest((metadata,), 300, 20, 2, True)

    def fake_compute(path, frame_rate, pixel_size, **kwargs):
        return (
            DDMData(np.ones((2, 2)), np.array([0.01, 0.1]), np.array([1e6, 2e6])),
            DDMData(np.ones((2, 2)), np.array([0.01, 0.1]), np.array([1e6, 2e6])),
        )

    monkeypatch.setattr("ddmsoft.gui.workers.compute_video_ddm", fake_compute)
    cancel = Event()

    def report(stage, completed, total):
        if stage == "saving":
            cancel.set()

    with pytest.raises(ComputationCancelled):
        run_video_computation(request, report, cancel.is_set)
    output_directory = tmp_path / "ddm_matrices"
    assert not tuple(output_directory.glob("*.npy"))


def test_time_dependent_worker_uses_each_video_metadata_and_frame_partition_names(
    tmp_path, monkeypatch
):
    first = tmp_path / "first.avi"
    second = tmp_path / "second.avi"
    first.touch()
    second.touch()
    videos = (
        VideoMetadata(first, 30.0, 1e-6),
        VideoMetadata(second, 60.0, 2e-6),
    )
    datasets = (
        PartitionedDDM(
            0,
            2,
            DDMData(np.ones((2, 1)), np.array([0.1, 0.2]), np.array([1.0])),
        ),
    )
    calls = []
    saved = []

    def fake_read(path):
        return [np.zeros((2, 2)) for _ in range(2)]

    def fake_compute(frames, frame_rate, pixel_size, partitions, **kwargs):
        calls.append((len(frames), frame_rate, pixel_size, partitions))
        return datasets

    def fake_save(prefix, results):
        saved.append((prefix, tuple(results)))
        path = prefix.with_name(f"{prefix.name}__i=0___DDM_matrix.npy")
        path.touch()
        return (path,)

    monkeypatch.setattr("ddmsoft.gui.workers.read_video_frames", fake_read)
    monkeypatch.setattr("ddmsoft.gui.workers.compute_time_dependent_ddm", fake_compute)
    monkeypatch.setattr("ddmsoft.gui.workers.save_partitioned_matrix_sets", fake_save)
    progress = []
    result = run_time_dependent_computation(
        TimeDependentComputationRequest(videos, 1, 300, 20, 1),
        lambda *event: progress.append(event),
        lambda: False,
    )

    assert calls == [(2, 30.0, 1e-6, 1), (2, 60.0, 2e-6, 1)]
    assert [prefix.name for prefix, _ in saved] == ["first", "second"]
    assert result.processed_videos == (first, second)
    assert len(result.paths) == 2
    assert progress[-1] == ("complete", 2, 2)


def test_time_dependent_worker_requires_recompute_and_rolls_back_failed_commit(
    tmp_path, monkeypatch
):
    video = tmp_path / "sample[1].avi"
    video.touch()
    metadata = VideoMetadata(video, 30.0, 1e-6)
    data = DDMData(np.ones((2, 1)), np.array([0.1, 0.2]), np.array([1.0]))
    datasets = (PartitionedDDM(0, 2, data),)
    monkeypatch.setattr("ddmsoft.gui.workers.read_video_frames", lambda path: [0, 1])
    monkeypatch.setattr(
        "ddmsoft.gui.workers.compute_time_dependent_ddm", lambda *args, **kwargs: datasets
    )
    output_directory = tmp_path / "ddm_matrices"
    output_directory.mkdir()
    existing = output_directory / "sample[1]__i=0___DDM_matrix.npy"
    existing.write_bytes(b"original")
    stale = output_directory / "sample[1]__i=9___DDM_matrix.npy"
    stale.write_bytes(b"stale")
    unrelated = output_directory / "sample1__i=9___DDM_matrix.npy"
    unrelated.write_bytes(b"unrelated")
    request = TimeDependentComputationRequest((metadata,), 1, 300, 20, 1)

    with pytest.raises(FileExistsError, match="select recompute"):
        run_time_dependent_computation(request, lambda *event: None, lambda: False)
    assert existing.read_bytes() == b"original"
    assert stale.read_bytes() == b"stale"

    original_replace = Path.replace

    def fail_second_replace(path, target):
        if Path(target).name.endswith("_deltaTs.npy") and not Path(path).name.startswith(".backup"):
            raise OSError("synthetic commit failure")
        return original_replace(path, target)

    monkeypatch.setattr(Path, "replace", fail_second_replace)
    with pytest.raises(OSError, match="synthetic commit failure"):
        run_time_dependent_computation(
            TimeDependentComputationRequest((metadata,), 1, 300, 20, 1, True),
            lambda *event: None,
            lambda: False,
        )
    assert existing.read_bytes() == b"original"
    assert stale.read_bytes() == b"stale"
    assert not (output_directory / "sample[1]__i=0___deltaTs.npy").exists()

    monkeypatch.setattr(Path, "replace", original_replace)
    run_time_dependent_computation(
        TimeDependentComputationRequest((metadata,), 1, 300, 20, 1, True),
        lambda *event: None,
        lambda: False,
    )
    assert not stale.exists()
    assert unrelated.read_bytes() == b"unrelated"


def test_contin_worker_returns_q_identity_and_scans_all_requested_candidates():
    tau = np.linspace(0.01, 0.2, 8)
    gamma = np.linspace(1.0, 5.0, 5)
    distribution = np.array([0.0, 0.2, 0.6, 0.2, 0.0])
    correlation = np.exp(-tau[:, None] * gamma[None, :]) @ distribution
    matrix = (1.4 * (1.0 - correlation) + 0.08)[:, None]
    path = Path("sample_DDM_matrix.npy")
    request = CONTINComputationRequest(
        path,
        DDMData(matrix, tau, np.array([1.0])),
        0,
        FitRange(0, 0, 0, tau.size - 1),
        1.0,
        5.0,
        5,
        0.01,
        1.0,
        3,
        1,
    )
    progress = []

    result = run_contin_computation(request, lambda *event: progress.append(event), lambda: False)

    assert result.matrix_path == path
    assert result.q_index == 0
    assert result.result.alphas.size == 3
    assert result.result.gamma_range.size == 5
    assert progress[-1] == ("contin", 3, 3)


def test_video_concatenation_worker_reports_output_and_honors_preflight_cancel(
    tmp_path, monkeypatch
):
    first = tmp_path / "first.avi"
    second = tmp_path / "second.avi"
    output = tmp_path / "combined.avi"
    calls = []

    def fake_concat(videos, destination, **kwargs):
        calls.append((videos, destination, kwargs))
        destination.touch()
        return destination

    monkeypatch.setattr("ddmsoft.gui.workers.concatenate_videos", fake_concat)
    progress = []
    result = run_video_concatenation(
        VideoConcatenationRequest((first, second), output),
        lambda *event: progress.append(event),
        lambda: False,
    )

    assert result.output == output
    assert calls == [
        ((first, second), output, {"overwrite": False, "cancel": ANY})
    ]
    assert progress == [("concatenating", 0, 1), ("concatenating", 1, 1)]

    with pytest.raises(ComputationCancelled):
        run_video_concatenation(
            VideoConcatenationRequest((first, second), output),
            lambda *event: None,
            lambda: True,
        )
    assert len(calls) == 1


def _run_worker(qapp, work, *, cancel_after=None, failure=None, cancelled=None, results=None):
    thread = QThread()
    worker = ComputationWorker(work)
    worker.moveToThread(thread)
    done = False

    def mark_done():
        nonlocal done
        done = True

    if failure is not None:
        worker.failure.connect(lambda *value: failure.append(value))
    if cancelled is not None:
        worker.cancelled.connect(lambda: cancelled.append(True))
    if results is not None:
        worker.result.connect(results.append)
    thread.started.connect(worker.run)
    worker.finished.connect(thread.quit)
    worker.finished.connect(worker.deleteLater)
    thread.finished.connect(mark_done)
    thread.start()
    if cancel_after is not None:
        QTimer.singleShot(cancel_after, lambda: worker.request_cancel())
    _wait_until(qapp, lambda: done)
    thread.wait()


def _wait_until(qapp, predicate, timeout=3.0):
    deadline = time.monotonic() + timeout
    while not predicate():
        qapp.processEvents()
        if time.monotonic() >= deadline:
            raise AssertionError("condition did not become true before timeout")
        time.sleep(0.001)
