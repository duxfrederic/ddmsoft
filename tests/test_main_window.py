from __future__ import annotations

import time

import numpy as np
import pytest
from PySide6.QtCore import QSettings, QTimer
from PySide6.QtWidgets import QApplication, QDialog

from ddmsoft.contin import run_contin
from ddmsoft.engine import ComputationCancelled
from ddmsoft.fitting import default_fit_request, fit_ddm
from ddmsoft.gui.main_window import DDMMainWindow
from ddmsoft.gui.workers import (
    ComputationResult,
    CONTINComputationResult,
    FitComputationResult,
    TimeDependentComputationResult,
    VideoConcatenationResult,
    batch_fit_target_paths,
)
from ddmsoft.io import save_matrix_set
from ddmsoft.models import DDMData, FitRange, VideoMetadata
from ddmsoft.plotting import (
    AmplitudeNoiseDiffusionPlotController,
    CONTINPlotController,
    CorrelationPlotController,
    FitParameterPlotController,
    MatrixPlotController,
)

from .fixtures import generate_model_data, write_legacy_matrix


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


def test_load_matrix_only_directory_populates_matrix_catalog(qapp, tmp_path):
    root = tmp_path / "matrix_archive"
    root.mkdir()
    data = generate_model_data("stretch")
    write_legacy_matrix(root, "archived", data)
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )

    assert window.load_directory_data(root)
    matrix_path = root / "ddm_matrices" / "archived_DDM_matrix.npy"
    assert window.video_table.rowCount() == 0
    assert window.matrix_selector.count() == 2
    assert window.matrix_selector.itemData(1) == matrix_path
    assert window.selected_matrix_path == matrix_path
    assert "1 matrix/matrices" in window.status_label.text()
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


def test_directional_processing_discovers_selects_plots_and_fits_each_sector(
    qapp, tmp_path, monkeypatch
):
    root = _metadata_only_directory(tmp_path, "sample")
    columns = np.arange(16, dtype=float)[None, :]
    frames = [
        np.repeat(np.sin(2 * np.pi * (columns - shift / 4) / 8), 16, axis=0)
        for shift in range(12)
    ]
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)
    window.direction_spin.setValue(2)
    monkeypatch.setattr("ddmsoft.engine.read_video_frames", lambda path: iter(frames))
    window.start_processing()
    _wait_for(qapp, lambda: not window._job_active)

    assert window.matrix_selector.count() == 3
    first_path = window.matrix_selector.itemData(1)
    second_path = window.matrix_selector.itemData(2)
    assert not np.allclose(
        window._matrices[first_path].matrix,
        window._matrices[second_path].matrix,
    )
    window.matrix_selector.setCurrentIndex(1)
    assert window.selected_matrix_path == first_path
    assert isinstance(window.plot_selected_matrix(), MatrixPlotController)
    window.model_selector.setCurrentIndex(window.model_selector.findData("stretch"))
    window.start_fitting()
    _wait_for(qapp, lambda: not window._job_active)
    assert window.selected_fit is not None
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


def test_fit_workflow_uses_worker_range_and_keeps_modeless_plots_independent(qapp, tmp_path):
    root = _metadata_only_directory(tmp_path, "sample")
    write_legacy_matrix(root, "sample", generate_model_data("stretch"))
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)
    window.model_selector.setCurrentIndex(window.model_selector.findData("stretch"))
    window.q_min_slider.setValue(1)
    window.time_min_slider.setValue(2)

    window.start_fitting()
    _wait_for(qapp, lambda: not window._job_active)

    fit = window.selected_fit
    assert fit is not None
    assert window.selected_fit_range is not None
    assert window.selected_fit_range.q_min == 1
    assert window.selected_fit_range.q_max == 3
    assert window.selected_fit_range.time_min == 2
    assert window.selected_fit_range.time_max == 6
    assert np.array_equal(fit.q_values, window.selected_matrix.q_values[1:4])
    assert fit.correlation.shape == (5, 3)

    matrix_plot = window.plot_selected_matrix()
    correlation_plot = window.plot_selected_correlation()
    parameter_plot = window.plot_fitted_parameters()
    amplitude_plot = window.plot_amplitude_noise_diffusion()
    assert isinstance(matrix_plot, MatrixPlotController)
    assert isinstance(correlation_plot, CorrelationPlotController)
    assert isinstance(parameter_plot, FitParameterPlotController)
    assert isinstance(amplitude_plot, AmplitudeNoiseDiffusionPlotController)
    assert matrix_plot.fit_image is not None
    assert correlation_plot.fit_range == window.selected_fit_range
    assert len(parameter_plot.parameter_lines) == 4
    assert len(amplitude_plot.lines) == 3

    window.temperature_edit.setText("25")
    radius_plot = window.plot_amplitude_noise_diffusion()
    assert isinstance(radius_plot, AmplitudeNoiseDiffusionPlotController)
    assert radius_plot.radius is not None
    assert radius_plot.axes[2].get_title() == "Hydrodynamic radius"

    old_value = window.q_min_slider.value()
    window.q_min_slider.setValue(0)
    assert window.q_min_slider.value() != old_value
    assert window.q_min_slider.isEnabled()
    window.close()


def test_fit_result_for_previous_matrix_is_ignored(qapp, tmp_path):
    root = _legacy_directory(tmp_path, "a", q_count=3, time_count=4)
    _write_dataset(root, "b", q_count=3, time_count=4)
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)
    first_path = window.matrix_selector.itemData(1)
    first_data = window._matrices[first_path]
    first_request = default_fit_request("stretch", FitRange(0, 2, 0, 3))
    first_fit = fit_ddm(first_data, first_request)
    window.matrix_selector.setCurrentIndex(2)

    window._fit_worker_result(FitComputationResult(first_path, first_request, first_fit))

    assert window.selected_matrix_path != first_path
    assert window.selected_fit is None
    window.close()


def test_matrix_selection_cancellation_has_no_combine_side_effect(qapp, tmp_path, monkeypatch):
    root = _legacy_directory(tmp_path, "a", q_count=3, time_count=4)
    _write_dataset(root, "b", q_count=3, time_count=4)
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)
    monkeypatch.setattr(window, "_select_matrix_paths", lambda **kwargs: None)

    assert not window.average_selected_matrices()
    assert not tuple((root / "ddm_matrices").glob("average_result*.npy"))
    window.close()


def test_average_refreshes_matrix_catalog_and_rejects_incompatible_inputs(
    qapp, tmp_path, monkeypatch
):
    root = _legacy_directory(tmp_path, "a", q_count=3, time_count=4)
    _write_dataset(root, "b", q_count=3, time_count=4)
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)
    selected = tuple(window._matrices)
    output_prefix = root / "ddm_matrices" / "average_result"
    monkeypatch.setattr(window, "_select_matrix_paths", lambda **kwargs: selected)
    monkeypatch.setattr(window, "_choose_output_prefix", lambda *args: output_prefix)

    assert window.average_selected_matrices()
    output_matrix = output_prefix.with_name(output_prefix.name + "_DDM_matrix.npy")
    assert output_matrix.is_file()
    assert window.selected_matrix_path == output_matrix

    incompatible_root = tmp_path / "incompatible"
    incompatible_root.mkdir()
    _legacy_directory(incompatible_root, "first", q_count=3, time_count=4)
    _write_dataset(incompatible_root, "second", q_count=2, time_count=4)
    assert window.load_directory_data(incompatible_root)
    incompatible = tuple(window._matrices)
    monkeypatch.setattr(window, "_select_matrix_paths", lambda **kwargs: incompatible)
    monkeypatch.setattr(
        window,
        "_choose_output_prefix",
        lambda *args: pytest.fail("output must not be requested after validation failure"),
    )
    assert not window.average_selected_matrices()
    assert "incompatible matrix width" in window.status_label.text()
    window.close()


def test_average_writes_one_matrix_set_per_compatible_lag_group(qapp, tmp_path, monkeypatch):
    root = _legacy_directory(tmp_path, "first", q_count=3, time_count=3)
    _write_dataset(root, "second", q_count=3, time_count=4)
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)
    selected = tuple(window._matrices)
    output_prefix = root / "ddm_matrices" / "average_result"
    confirmed: list[tuple] = []
    monkeypatch.setattr(window, "_select_matrix_paths", lambda **kwargs: selected)
    monkeypatch.setattr(window, "_choose_output_prefix", lambda *args: output_prefix)
    monkeypatch.setattr(
        window,
        "_confirm_output_targets",
        lambda targets: confirmed.append(tuple(targets)) or True,
    )

    assert window.average_selected_matrices()
    first_matrix = output_prefix.with_name(
        f"{output_prefix.name}_average_0_DDM_matrix.npy"
    )
    second_matrix = output_prefix.with_name(
        f"{output_prefix.name}_average_1_DDM_matrix.npy"
    )
    assert first_matrix.is_file()
    assert second_matrix.is_file()
    assert confirmed
    assert window.selected_matrix_path == first_matrix
    window.close()


def test_merge_refreshes_matrix_catalog_after_successful_write(qapp, tmp_path, monkeypatch):
    root = _metadata_only_directory(tmp_path, "fast")
    fast = generate_model_data("stretch")
    slow = DDMData(
        fast.matrix + 1.0,
        fast.lag_times * 2.0,
        fast.q_values,
    )
    write_legacy_matrix(root, "fast", fast)
    write_legacy_matrix(root, "slow", slow)
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)
    selected = tuple(window._matrices)
    output_prefix = root / "ddm_matrices" / "merged_result"
    monkeypatch.setattr(window, "_select_matrix_paths", lambda **kwargs: selected)
    monkeypatch.setattr(window, "_choose_output_prefix", lambda *args: output_prefix)

    assert window.merge_selected_matrices()
    output_matrix = output_prefix.with_name(output_prefix.name + "_DDM_matrix.npy")
    assert output_matrix.is_file()
    assert window.selected_matrix_path == output_matrix
    assert window.selected_matrix.matrix.shape[1] == fast.matrix.shape[1]
    window.close()


def test_exports_write_expected_dimensions_and_headers(qapp, tmp_path, monkeypatch):
    root = _metadata_only_directory(tmp_path, "sample")
    data = generate_model_data("stretch")
    write_legacy_matrix(root, "sample", data)
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)
    window.model_selector.setCurrentIndex(window.model_selector.findData("stretch"))
    request = default_fit_request("stretch", FitRange(0, 3, 0, 6))
    path = window.selected_matrix_path
    window._fit_results[path] = fit_ddm(data, request)
    window._fit_requests[path] = request
    window._update_matrix_action_state()
    export_root = tmp_path / "exports"

    def output_for(default, suffixes):
        if suffixes == (".txt",):
            return export_root / "fit"
        if suffixes == ("_autocorrelationmatrix.csv", "_qs.csv", "_dts.csv"):
            return export_root / "correlation"
        return export_root / "matrix"

    monkeypatch.setattr(window, "_choose_output_prefix", output_for)
    assert window.export_selected_matrix()
    assert window.export_selected_fit()
    assert window.export_selected_correlation()
    matrix_csv = np.loadtxt(export_root / "matrix_DDM_matrix.csv", delimiter="\t")
    correlation_csv = np.loadtxt(
        export_root / "correlation_autocorrelationmatrix.csv", delimiter="\t"
    )
    fit_lines = (export_root / "fit.txt").read_text(encoding="utf-8").splitlines()
    assert matrix_csv.shape == data.matrix.shape
    assert correlation_csv.shape == data.matrix.shape
    assert fit_lines[0] == "q [m^-1]\tA\tB\tdiffusion\tstretch"
    window.close()


def test_unwritable_export_is_reported_without_crashing(qapp, tmp_path, monkeypatch):
    root = _legacy_directory(tmp_path, "sample", q_count=3, time_count=4)
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)
    monkeypatch.setattr(window, "_choose_output_prefix", lambda *args: tmp_path / "blocked")

    def permission_denied(*args, **kwargs):
        raise PermissionError("synthetic unwritable output")

    monkeypatch.setattr("ddmsoft.gui.main_window.save_matrix_csv", permission_denied)

    assert not window.export_selected_matrix()
    assert "synthetic unwritable output" in window.status_label.text()
    window.close()


def test_about_action_shows_native_credits_dialog(qapp, monkeypatch):
    window = DDMMainWindow()
    calls = []
    monkeypatch.setattr(
        "ddmsoft.gui.main_window.QMessageBox.about",
        lambda *args: calls.append(args),
    )

    window.about_action.trigger()

    assert calls
    assert calls[0][1] == "About DDMSoft"
    window.close()


def test_batch_fit_main_window_stores_results_by_full_path(qapp, tmp_path, monkeypatch):
    root = _metadata_only_directory(tmp_path, "a")
    first = generate_model_data("stretch")
    second = generate_model_data("stretch", noise=0.001)
    write_legacy_matrix(root, "a", first)
    write_legacy_matrix(root, "b", second)
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)
    window.model_selector.setCurrentIndex(window.model_selector.findData("stretch"))
    selected = tuple(window._matrices)
    batch_prefix = tmp_path / "batch"

    class AcceptedBatchDialog:
        DialogCode = QDialog.DialogCode
        continue_on_failure = True

        def __init__(self, *args, **kwargs):
            pass

        def exec(self):
            return QDialog.DialogCode.Accepted

    AcceptedBatchDialog.selected_paths = selected
    AcceptedBatchDialog.output_prefix = batch_prefix
    monkeypatch.setattr("ddmsoft.gui.main_window.BatchFitDialog", AcceptedBatchDialog)
    monkeypatch.setattr(window, "_confirm_output_targets", lambda targets: True)
    assert window.start_batch_fitting()
    _wait_for(qapp, lambda: not window._job_active)

    assert set(window._fit_results) == set(selected)
    assert all(path.is_file() for path in batch_fit_target_paths(batch_prefix, selected))
    assert "2 succeeded" in window.status_label.text()
    window.close()


def test_time_dependent_main_window_runs_selected_videos_and_refreshes_catalog(
    qapp, tmp_path, monkeypatch
):
    root = _metadata_only_directory(tmp_path, "sample")
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)
    video = VideoMetadata(root / "sample.avi", 30.0, 1e-6)

    class AcceptedTimeDependentDialog:
        DialogCode = QDialog.DialogCode
        selected_videos = (video,)
        partitions = 2

        def __init__(self, *args, **kwargs):
            pass

        def exec(self):
            return QDialog.DialogCode.Accepted

    def fake_job(request, report, cancel):
        assert request.videos == (video,)
        output_directory = root / "ddm_matrices"
        output_directory.mkdir(exist_ok=True)
        paths = save_matrix_set(
            output_directory / "sample__i=0__",
            DDMData(np.ones((2, 2)), np.array([0.01, 0.1]), np.array([1e6, 2e6])),
        )
        return TimeDependentComputationResult(paths, (video.path,))

    monkeypatch.setattr("ddmsoft.gui.main_window.TimeDependentDialog", AcceptedTimeDependentDialog)
    monkeypatch.setattr("ddmsoft.gui.main_window.run_time_dependent_computation", fake_job)
    assert window.start_time_dependent_processing()
    _wait_for(qapp, lambda: not window._job_active)

    assert window.selected_matrix_path == root / "ddm_matrices" / "sample__i=0___DDM_matrix.npy"
    assert "1 video(s)" in window.status_label.text()
    window.close()


def test_video_concatenation_selects_arbitrary_avi_directory_without_metadata(
    qapp, tmp_path, monkeypatch
):
    first = tmp_path / "first.avi"
    second = tmp_path / "second.avi"
    first.touch()
    second.touch()
    output_prefix = tmp_path / "joined"
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )

    class AcceptedVideoDialog:
        DialogCode = QDialog.DialogCode
        selected_paths = (first, second)

        def __init__(self, videos, **kwargs):
            assert set(videos) == {first, second}

        def exec(self):
            return QDialog.DialogCode.Accepted

    monkeypatch.setattr(
        "ddmsoft.gui.main_window.QFileDialog.getExistingDirectory",
        lambda *args: str(tmp_path),
    )
    monkeypatch.setattr("ddmsoft.gui.main_window.VideoSelectionDialog", AcceptedVideoDialog)
    monkeypatch.setattr(window, "_choose_output_prefix", lambda *args: output_prefix)
    monkeypatch.setattr(
        "ddmsoft.gui.main_window.run_video_concatenation",
        lambda request, progress, cancel: VideoConcatenationResult(request.output),
    )

    assert window.concatenate_selected_videos()
    _wait_for(qapp, lambda: not window._job_active)
    assert window.status_label.text() == "Concatenated video: joined.avi"
    window.close()


def test_contin_main_window_keeps_independent_result_and_exports_candidates(
    qapp, tmp_path, monkeypatch
):
    root = _metadata_only_directory(tmp_path, "sample")
    data = generate_model_data("stretch")
    write_legacy_matrix(root, "sample", data)
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)
    path = window.selected_matrix_path
    result = run_contin(
        np.linspace(0.01, 0.2, data.matrix.shape[0]),
        np.linspace(0.1, 0.8, data.matrix.shape[0]),
        np.linspace(1.0, 5.0, 5),
        alpha=(0.01, 0.1),
        maxiter=1,
    )
    value = CONTINComputationResult(path, 1, FitRange(0, 3, 0, 6), result)
    window._contin_save_all[(path, 1)] = True
    window._contin_worker_result(value)

    assert window.selected_contin is result
    controller = window._plot_controllers[-1]
    assert isinstance(controller, CONTINPlotController)
    controller.set_alpha_index(1)
    output_prefix = tmp_path / "contin_export"
    monkeypatch.setattr(window, "_choose_output_prefix", lambda *args: output_prefix)
    assert window.export_selected_contin()
    output = output_prefix.with_suffix(".txt")
    assert output.is_file()
    assert output.read_text(encoding="utf-8").count("\nalpha:") == 2
    window.close()


def test_contin_selected_alpha_is_used_for_selected_export(qapp, tmp_path, monkeypatch):
    root = _metadata_only_directory(tmp_path, "sample")
    data = generate_model_data("stretch")
    write_legacy_matrix(root, "sample", data)
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)
    path = window.selected_matrix_path
    result = run_contin(
        np.linspace(0.01, 0.2, data.matrix.shape[0]),
        np.linspace(0.1, 0.8, data.matrix.shape[0]),
        np.linspace(1.0, 5.0, 5),
        alpha=(0.01, 0.1),
        maxiter=1,
    )
    window._contin_save_all[(path, 1)] = False
    window._contin_worker_result(
        CONTINComputationResult(path, 1, FitRange(0, 3, 0, 6), result)
    )
    controller = window._plot_controllers[-1]
    assert isinstance(controller, CONTINPlotController)
    controller.set_alpha_index(1)
    output_prefix = tmp_path / "selected_contin_export"
    monkeypatch.setattr(window, "_choose_output_prefix", lambda *args: output_prefix)

    assert window.export_selected_contin()
    text = output_prefix.with_suffix(".txt").read_text(encoding="utf-8")
    assert text.count("\nalpha:") == 1
    assert f"selected alpha:\t{result.alphas[1]:.03e}" in text
    window.close()


def test_contin_result_selection_follows_matrix_selection(qapp, tmp_path):
    root = _legacy_directory(tmp_path, "first", q_count=3, time_count=4)
    _write_dataset(root, "second", q_count=3, time_count=4)
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)
    first_path = window.matrix_selector.itemData(1)
    window._contin_results[(first_path, 0)] = object()
    window.matrix_selector.setCurrentIndex(2)
    assert window._selected_contin_key is None
    window.matrix_selector.setCurrentIndex(1)
    assert window._selected_contin_key == (first_path, 0)
    window.close()


def test_closing_window_ignores_queued_contin_result(qapp, tmp_path):
    root = _legacy_directory(tmp_path, "sample", q_count=3, time_count=4)
    window = DDMMainWindow(
        settings=QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    )
    assert window.load_directory_data(root)
    window._closing = True

    window._contin_worker_result(object())

    assert window._plot_controllers == []
    window._closing = False
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
