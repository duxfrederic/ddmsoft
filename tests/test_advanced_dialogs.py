from __future__ import annotations

import pytest
from PySide6.QtWidgets import QApplication, QDialog

from ddmsoft.gui.advanced_dialogs import CONTINDialog, TimeDependentDialog, VideoSelectionWidget
from ddmsoft.models import VideoMetadata


@pytest.fixture
def qapp():
    instance = QApplication.instance() or QApplication([])
    yield instance
    for widget in instance.topLevelWidgets():
        widget.close()


def test_video_selection_widget_filters_by_name_or_full_path(qapp, tmp_path):
    first = tmp_path / "first clip.avi"
    second = tmp_path / "second.avi"
    widget = VideoSelectionWidget({first: first.name, second: second.name})

    widget.search_edit.setText("second")

    assert widget.video_list.item(0).isHidden()
    assert not widget.video_list.item(1).isHidden()
    widget.select_all_visible()
    assert widget.selected_paths == (second,)
    widget.clear_selection()
    assert widget.selected_paths == ()


def test_time_dependent_dialog_returns_selected_per_video_metadata(qapp, tmp_path):
    first = tmp_path / "first.avi"
    second = tmp_path / "second.avi"
    dialog = TimeDependentDialog(
        {
            first: VideoMetadata(first, 30.0, 1e-6),
            second: VideoMetadata(second, 60.0, 2e-6),
        },
        partitions=4,
    )

    assert len(dialog.video_table.selectionModel().selectedRows()) == 2
    dialog.video_table.item(0, 1).setText("24")
    dialog.partition_spin.setValue(3)
    dialog._accept_time_dependent()

    assert dialog.result() == QDialog.DialogCode.Accepted
    assert dialog.partitions == 3
    assert dialog.selected_videos[0].frame_rate == 24.0
    assert dialog.selected_videos[1].pixel_size == 2e-6


def test_time_dependent_dialog_rejects_invalid_metadata(qapp, tmp_path):
    path = tmp_path / "sample.avi"
    dialog = TimeDependentDialog({path: VideoMetadata(path, 30.0, 1e-6)})
    dialog.video_table.item(0, 2).setText("0")

    dialog._accept_time_dependent()

    assert dialog.result() == QDialog.DialogCode.Rejected
    assert "positive" in dialog.error_label.text()


def test_contin_dialog_validates_ranges_and_exposes_counts(qapp):
    dialog = CONTINDialog(4, q_index=2)
    dialog.gamma_min_edit.setText("5")
    dialog.gamma_max_edit.setText("1")
    dialog._accept_contin()

    assert dialog.result() == QDialog.DialogCode.Rejected
    assert "ordered" in dialog.error_label.text()

    dialog.gamma_min_edit.setText("1")
    dialog.gamma_max_edit.setText("5")
    dialog.gamma_count_spin.setValue(5)
    dialog.alpha_count_spin.setValue(3)
    dialog.save_all_check.setChecked(False)
    dialog._accept_contin()

    assert dialog.result() == QDialog.DialogCode.Accepted
    assert dialog.q_index == 2
    assert dialog.gamma_count == 5
    assert dialog.alpha_count == 3
    assert not dialog.save_all
