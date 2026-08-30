from __future__ import annotations

import pytest
from PySide6.QtWidgets import QApplication, QDialog

from ddmsoft.gui.selection import (
    MatrixSelectionDialog,
    MatrixSelectionWidget,
    output_target_paths,
)


@pytest.fixture
def qapp():
    instance = QApplication.instance() or QApplication([])
    yield instance
    for widget in instance.topLevelWidgets():
        widget.close()


def _matrices(tmp_path):
    return {
        tmp_path / "a_DDM_matrix.npy": "a_DDM_matrix.npy",
        tmp_path / "b_DDM_matrix.npy": "b_DDM_matrix.npy",
        tmp_path / "third_DDM_matrix.npy": "third_DDM_matrix.npy",
    }


def test_matrix_selection_filters_and_selects_full_paths(qapp, tmp_path):
    matrices = _matrices(tmp_path)
    widget = MatrixSelectionWidget(matrices, selected_paths=(next(iter(matrices)),))

    assert widget.selected_paths == (next(iter(matrices)),)
    widget.search_edit.setText("third")
    widget.select_all_visible()
    assert widget.selected_paths == (tmp_path / "third_DDM_matrix.npy",)
    widget.clear_button.click()
    assert widget.selected_paths == ()
    widget.search_edit.clear()
    widget.select_all_button.click()
    assert set(widget.selected_paths) == set(matrices)


def test_selection_cancel_has_no_selected_result(qapp, tmp_path):
    dialog = MatrixSelectionDialog(_matrices(tmp_path), selected_paths=tuple(_matrices(tmp_path)))
    dialog.reject()

    assert dialog.result() == QDialog.DialogCode.Rejected
    assert dialog.selected_paths == ()


def test_selection_requires_at_least_one_matrix(qapp, tmp_path):
    dialog = MatrixSelectionDialog(_matrices(tmp_path))
    dialog._accept_selection()

    assert dialog.result() == QDialog.DialogCode.Rejected
    assert "at least one" in dialog.error_label.text()


def test_output_target_paths_show_exact_suffixes(tmp_path):
    prefix = tmp_path / "result"

    assert output_target_paths(
        prefix,
        ("_DDM_matrix.npy", "_deltaTs.npy", "_QS.npy"),
    ) == (
        tmp_path / "result_DDM_matrix.npy",
        tmp_path / "result_deltaTs.npy",
        tmp_path / "result_QS.npy",
    )
    assert output_target_paths(prefix, (".txt",)) == (tmp_path / "result.txt",)
