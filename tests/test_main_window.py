from __future__ import annotations

import time

import numpy as np
import pytest
from PySide6.QtCore import QSettings, QTimer
from PySide6.QtWidgets import QApplication

from ddmsoft.engine import ComputationCancelled
from ddmsoft.gui.main_window import DDMMainWindow
from ddmsoft.gui.workers import ComputationResult
from ddmsoft.io import save_matrix_set
from ddmsoft.models import DDMData, VideoMetadata


def test_main_window_maps_legacy_workflow_controls(qapp):
    window = DDMMainWindow()

    assert [action.menu().title() for action in window.menuBar().actions()] == [
        "File",
        "Tools",
        "Plotting",
        "Batch, export",
        "More Fitting",
        "Help",
    ]
    assert window.computation_group.title() == "DDM matrix computation"
    assert window.fitting_group.title() == "DDM matrix fitting"
    assert window.plotting_group.title() == "Plotting"
    assert window.video_table.columnCount() == 4
    assert window.video_table.horizontalHeaderItem(0).text() == "Path / name"
    assert window.video_table.editTriggers() != 0
    assert window.max_couples_spin.value() == 300
    assert window.points_per_decade_spin.value() == 20
    assert window.direction_spin.value() == 1
    assert window.keep_existing_radio.isChecked()
    assert window.model_selector.count() == 8
    assert window.matrix_selector.currentText() == "No matrix loaded"
    assert window.q_min_value.text() == "index: 0; q: unavailable"
    assert window.time_max_value.text() == "index: 100; time: unavailable"
    window.close()


def test_main_window_range_controls_show_physical_values_and_preserve_order(qapp):
    window = DDMMainWindow()
    window.set_axis_values((1e6, 2e6, 3e6), (0.01, 0.1, 1.0, 10.0))

    assert window.q_max_slider.maximum() == 2
    assert window.time_max_slider.maximum() == 3
    assert window.q_max_value.text() == "index: 2; q: 3e+06"
    assert window.time_max_value.text() == "index: 3; time: 10"

    window.q_min_slider.setValue(window.q_max_slider.value())
    assert window.q_min_slider.value() < window.q_max_slider.value()
    window.time_max_slider.setValue(window.time_min_slider.value())
    assert window.time_min_slider.value() < window.time_max_slider.value()
    window.close()


def test_main_window_persists_geometry_and_last_directory(qapp, tmp_path):
    settings = QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    first = DDMMainWindow(settings=settings)
    first.directory_edit.setText(str(tmp_path))
    first.close()
    settings.sync()

    second = DDMMainWindow(settings=settings)
    assert second.directory_edit.text() == str(tmp_path)
    assert second.saveGeometry().size() > 0
    second.close()


def test_load_directory_populates_metadata_and_matrix_catalog(qapp, tmp_path, monkeypatch):
    root = _legacy_directory(tmp_path, "sample", q_count=3, time_count=4)
    settings = QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    window = DDMMainWindow(settings=settings)

    monkeypatch.setattr(
        "ddmsoft.gui.main_window.QFileDialog.getExistingDirectory", lambda *args: str(root)
    )
    assert window.choose_directory() == str(root)
    assert window.video_table.rowCount() == 1
    assert window.video_table.item(0, 0).text() == str(root / "sample.avi")
    assert window.video_table.item(0, 3).text() == "Valid"
    assert window.video_metadata() == (
        VideoMetadata(root / "sample.avi", 30.0, 1e-6),
    )
    assert window.process_button.isEnabled()
    assert window.matrix_selector.count() == 2
    matrix_path = root / "ddm_matrices" / "sample_DDM_matrix.npy"
    assert window.matrix_selector.itemData(1) == matrix_path
    assert window.selected_matrix_path == matrix_path
    assert window.q_max_slider.maximum() == 2
    assert window.time_max_slider.maximum() == 3
    assert window.q_max_value.text() == "index: 2; q: 3e+06"
    assert window.time_max_value.text() == "index: 3; time: 1"
    assert window.plot_matrix_button.isEnabled()
    window.close()


def test_invalid_metadata_edits_are_retained_and_block_processing(qapp, tmp_path):
    root = _legacy_directory(tmp_path, "sample")
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)

    frame_rate = window.video_table.item(0, 1)
    frame_rate.setText("not-a-number")
    assert frame_rate.text() == "not-a-number"
    assert "Invalid: frame rate" in window.video_table.item(0, 3).text()
    assert "row 1, frame rate" in window.status_label.text()
    assert not window.process_button.isEnabled()
    with pytest.raises(ValueError, match="row 1, frame rate"):
        window.video_metadata()

    frame_rate.setText("60")
    assert window.video_table.item(0, 3).text() == "Valid"
    assert window.process_button.isEnabled()
    window.close()


def test_invalid_loaded_metadata_identifies_the_error_without_clearing_rows(qapp, tmp_path):
    root = tmp_path
    (root / "sample.avi").touch()
    (root / "sample.txt").write_text("framerate: 0\npixelsize: 1e-6\n", encoding="utf-8")
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )

    assert not window.load_directory_data(root)
    assert window.video_table.rowCount() == 1
    assert "Invalid metadata" in window.video_table.item(0, 3).text()
    assert "framerate" in window.video_table.item(0, 3).toolTip()
    assert not window.process_button.isEnabled()
    window.close()


def test_matrix_selection_clamps_ranges_and_disables_without_selection(qapp, tmp_path):
    root = _legacy_directory(tmp_path, "a_matrix", q_count=5, time_count=5)
    _write_dataset(root, "b_matrix", q_count=2, time_count=3)
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)

    window.q_min_slider.setValue(3)
    window.time_min_slider.setValue(2)
    window.matrix_selector.setCurrentIndex(2)
    assert window.selected_matrix_path == root / "ddm_matrices" / "b_matrix_DDM_matrix.npy"
    assert window.q_min_slider.maximum() == 1
    assert window.q_min_slider.value() == 0
    assert window.q_max_slider.value() == 1
    assert window.time_min_slider.maximum() == 2
    assert window.time_min_slider.value() == 1
    assert window.time_max_slider.value() == 2
    assert window.q_max_value.text() == "index: 1; q: 2e+06"

    window.matrix_selector.setCurrentIndex(0)
    assert window.selected_matrix is None
    assert not window.fit_button.isEnabled()
    assert not window.plot_matrix_button.isEnabled()
    window.close()


def test_processing_uses_worker_and_refreshes_catalog_without_blocking(qapp, tmp_path, monkeypatch):
    root = _metadata_only_directory(tmp_path, "sample")
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)
    heartbeats = []

    def fake_job(request, report, cancel):
        for index in range(8):
            time.sleep(0.01)
            report("synthetic", index + 1, 8)
        output_directory = request.videos[0].path.parent / "ddm_matrices"
        output_directory.mkdir()
        paths = save_matrix_set(
            output_directory / request.videos[0].path.stem,
            DDMData(np.ones((2, 2)), np.array([0.01, 0.1]), np.array([1e6, 2e6])),
        )
        return ComputationResult(paths, tuple(video.path for video in request.videos), ())

    monkeypatch.setattr("ddmsoft.gui.main_window.run_video_computation", fake_job)
    heartbeat = QTimer()
    heartbeat.timeout.connect(lambda: heartbeats.append(time.monotonic()))
    heartbeat.start(2)
    window.start_processing()
    assert window.cancel_button.isEnabled()
    assert not window.load_button.isEnabled()
    _wait_for(qapp, lambda: not window._job_active)
    heartbeat.stop()

    assert len(heartbeats) >= 3
    assert window.progress_bar.value() == 100
    assert "Completed 1 video" in window.status_label.text()
    assert window.process_button.isEnabled()
    assert window.selected_matrix_path == root / "ddm_matrices" / "sample_DDM_matrix.npy"
    window.close()


def test_processing_failure_preserves_traceback_in_details_view(qapp, tmp_path, monkeypatch):
    root = _metadata_only_directory(tmp_path, "sample")
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)

    def fail_job(request, report, cancel):
        raise RuntimeError("synthetic GUI failure")

    monkeypatch.setattr("ddmsoft.gui.main_window.run_video_computation", fail_job)
    window.start_processing()
    _wait_for(qapp, lambda: not window._job_active)

    assert window.status_label.text() == "Processing failed: synthetic GUI failure"
    assert not window.error_details.isHidden()
    assert "RuntimeError: synthetic GUI failure" in window.error_details.toPlainText()
    window.close()


def test_processing_cancel_button_requests_cooperative_cancellation(qapp, tmp_path, monkeypatch):
    root = _metadata_only_directory(tmp_path, "sample")
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)

    def cancellable_job(request, report, cancel):
        for index in range(100):
            if cancel():
                raise ComputationCancelled("cancelled by GUI test")
            time.sleep(0.005)
            report("synthetic", index + 1, 100)
        raise AssertionError("job was not cancelled")

    monkeypatch.setattr("ddmsoft.gui.main_window.run_video_computation", cancellable_job)
    window.start_processing()
    QTimer.singleShot(25, window.cancel_processing)
    _wait_for(qapp, lambda: not window._job_active)

    assert window.status_label.text() == "Processing cancelled"
    assert window.progress_bar.value() < 100
    assert window.cancel_button.isEnabled() is False
    window.close()


def test_closing_window_waits_for_worker_cancellation(qapp, tmp_path, monkeypatch):
    root = _metadata_only_directory(tmp_path, "sample")
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)

    def cancellable_job(request, report, cancel):
        while not cancel():
            time.sleep(0.005)
        raise ComputationCancelled("cancelled while closing")

    monkeypatch.setattr("ddmsoft.gui.main_window.run_video_computation", cancellable_job)
    window.start_processing()
    window.close()

    assert window._job_active is False
    assert window._thread is None


def _legacy_directory(tmp_path, stem, *, q_count=3, time_count=3):
    root = tmp_path
    (root / f"{stem}.avi").touch()
    (root / f"{stem}.txt").write_text(
        "framerate: 30\npixelsize: 1e-6\n", encoding="utf-8"
    )
    _write_dataset(root, stem, q_count=q_count, time_count=time_count)
    return root


def _write_dataset(root, stem, *, q_count, time_count):
    matrix_directory = root / "ddm_matrices"
    matrix_directory.mkdir(exist_ok=True)
    matrix = np.arange(q_count * time_count, dtype=float).reshape(time_count, q_count)
    prefix = matrix_directory / stem
    np.save(f"{prefix}_DDM_matrix.npy", matrix)
    np.save(f"{prefix}_deltaTs.npy", np.linspace(0.01, 1.0, time_count))
    np.save(f"{prefix}_QS.npy", np.arange(1, q_count + 1, dtype=float) * 1e6)


def _metadata_only_directory(tmp_path, stem):
    root = tmp_path
    (root / f"{stem}.avi").touch()
    (root / f"{stem}.txt").write_text(
        "framerate: 30\npixelsize: 1e-6\n", encoding="utf-8"
    )
    return root


def _wait_for(qapp, predicate, timeout=3.0):
    deadline = time.monotonic() + timeout
    while not predicate():
        qapp.processEvents()
        if time.monotonic() >= deadline:
            raise AssertionError("condition did not become true before timeout")
        time.sleep(0.001)


@pytest.fixture
def qapp():
    instance = QApplication.instance() or QApplication([])
    yield instance
    for widget in instance.topLevelWidgets():
        widget.close()
