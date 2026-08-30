from __future__ import annotations

import pytest
from PySide6.QtCore import QSettings
from PySide6.QtWidgets import QApplication

from ddmsoft.gui.main_window import DDMMainWindow


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


@pytest.fixture
def qapp():
    instance = QApplication.instance() or QApplication([])
    yield instance
    for widget in instance.topLevelWidgets():
        widget.close()
