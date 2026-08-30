"""Reusable matrix selection and output-target dialogs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from .workers import batch_fit_target_paths


class MatrixSelectionWidget(QWidget):
    """Searchable multi-selection list whose item data is a full matrix path."""

    def __init__(
        self,
        matrices: Mapping[Path, str],
        *,
        selected_paths: Sequence[Path] = (),
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        selected = {Path(path) for path in selected_paths}
        layout = QVBoxLayout(self)
        self.search_edit = QLineEdit()
        self.search_edit.setObjectName("matrixSearchEdit")
        self.search_edit.setPlaceholderText("Search matrix names or paths")
        self.search_edit.setToolTip("Filter the matrix list without changing selection")
        layout.addWidget(self.search_edit)

        self.matrix_list = QListWidget()
        self.matrix_list.setObjectName("matrixSelectionList")
        self.matrix_list.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)
        for path, name in sorted(matrices.items(), key=lambda item: (item[1].casefold(), str(item[0]))):
            item = QListWidgetItem(name)
            item.setData(Qt.ItemDataRole.UserRole, Path(path))
            item.setToolTip(str(path))
            self.matrix_list.addItem(item)
            item.setSelected(Path(path) in selected)
        layout.addWidget(self.matrix_list, 1)

        actions = QHBoxLayout()
        self.select_all_button = QPushButton("Select All")
        self.select_all_button.setObjectName("selectAllMatricesButton")
        self.select_all_button.setToolTip("Select every currently visible matrix")
        self.clear_button = QPushButton("Clear")
        self.clear_button.setObjectName("clearMatrixSelectionButton")
        self.clear_button.setToolTip("Clear the matrix selection")
        actions.addWidget(self.select_all_button)
        actions.addWidget(self.clear_button)
        actions.addStretch(1)
        layout.addLayout(actions)
        self.search_edit.textChanged.connect(self._filter_items)
        self.select_all_button.clicked.connect(self.select_all_visible)
        self.clear_button.clicked.connect(self.clear_selection)

    @property
    def selected_paths(self) -> tuple[Path, ...]:
        return tuple(
            Path(item.data(Qt.ItemDataRole.UserRole))
            for item in self.matrix_list.selectedItems()
        )

    def select_all_visible(self) -> None:
        for index in range(self.matrix_list.count()):
            item = self.matrix_list.item(index)
            item.setSelected(not item.isHidden())

    def clear_selection(self) -> None:
        self.matrix_list.clearSelection()

    def _filter_items(self, text: str) -> None:
        query = text.casefold().strip()
        for index in range(self.matrix_list.count()):
            item = self.matrix_list.item(index)
            path = Path(item.data(Qt.ItemDataRole.UserRole))
            item.setHidden(bool(query) and query not in f"{item.text()} {path}".casefold())


class MatrixSelectionDialog(QDialog):
    """Modal matrix selector with cancellation that has no side effects."""

    def __init__(
        self,
        matrices: Mapping[Path, str],
        *,
        selected_paths: Sequence[Path] = (),
        title: str = "Select matrices",
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setMinimumSize(620, 420)
        self.selection = MatrixSelectionWidget(matrices, selected_paths=selected_paths)
        self.selected_paths: tuple[Path, ...] = ()
        layout = QVBoxLayout(self)
        layout.addWidget(self.selection, 1)
        self.error_label = QLabel()
        self.error_label.setObjectName("matrixSelectionError")
        self.error_label.setStyleSheet("color: #b00020")
        layout.addWidget(self.error_label)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._accept_selection)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _accept_selection(self) -> None:
        selected = self.selection.selected_paths
        if not selected:
            self.error_label.setText("Select at least one matrix.")
            return
        self.selected_paths = selected
        self.accept()


class OutputPathDialog(QDialog):
    """Collect an output prefix and show every exact target before writing."""

    def __init__(
        self,
        default_prefix: str | Path,
        suffixes: Sequence[str],
        *,
        title: str = "Output files",
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setMinimumWidth(650)
        self._suffixes = tuple(suffixes)
        self.output_prefix = Path(default_prefix)
        self.target_paths: tuple[Path, ...] = ()
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Output name or prefix"))
        self.prefix_edit = QLineEdit(str(default_prefix))
        self.prefix_edit.setObjectName("outputPrefixEdit")
        layout.addWidget(self.prefix_edit)
        layout.addWidget(QLabel("Exact files that will be written:"))
        self.preview_list = QListWidget()
        self.preview_list.setObjectName("outputTargetList")
        self.preview_list.setSelectionMode(QListWidget.SelectionMode.NoSelection)
        layout.addWidget(self.preview_list)
        self.error_label = QLabel()
        self.error_label.setObjectName("outputPathError")
        self.error_label.setStyleSheet("color: #b00020")
        layout.addWidget(self.error_label)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._accept_output)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.prefix_edit.textChanged.connect(self._update_preview)
        self._update_preview(self.prefix_edit.text())

    def _update_preview(self, text: str) -> None:
        self.error_label.clear()
        prefix = Path(text.strip()) if text.strip() else None
        self.target_paths = output_target_paths(prefix, self._suffixes) if prefix else ()
        self.preview_list.clear()
        for path in self.target_paths:
            self.preview_list.addItem(str(path))

    def _accept_output(self) -> None:
        text = self.prefix_edit.text().strip()
        if not text:
            self.error_label.setText("An output name is required.")
            return
        self.output_prefix = Path(text).expanduser()
        self._update_preview(str(self.output_prefix))
        if not self.target_paths:
            self.error_label.setText("The output name did not produce any target files.")
            return
        self.accept()


class BatchFitDialog(QDialog):
    """Select batch matrices, output prefix, and failure policy together."""

    def __init__(
        self,
        matrices: Mapping[Path, str],
        default_prefix: str | Path,
        *,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Fit and save matrices")
        self.setMinimumSize(700, 560)
        self.selection = MatrixSelectionWidget(matrices, selected_paths=tuple(matrices))
        self.selected_paths: tuple[Path, ...] = ()
        self.output_prefix = Path(default_prefix)
        self.continue_on_failure = True
        self.target_paths: tuple[Path, ...] = ()
        layout = QVBoxLayout(self)
        layout.addWidget(self.selection, 1)
        output_row = QHBoxLayout()
        output_row.addWidget(QLabel("Output prefix"))
        self.prefix_edit = QLineEdit(str(default_prefix))
        self.prefix_edit.setObjectName("batchOutputPrefixEdit")
        output_row.addWidget(self.prefix_edit, 1)
        layout.addLayout(output_row)
        self.continue_check = QCheckBox("Continue after a matrix fails")
        self.continue_check.setObjectName("batchContinueOnFailureCheck")
        self.continue_check.setChecked(True)
        self.continue_check.setToolTip("Keep fitting later matrices when one matrix cannot be fitted")
        layout.addWidget(self.continue_check)
        layout.addWidget(QLabel("Exact fit files that will be written:"))
        self.preview_list = QListWidget()
        self.preview_list.setObjectName("batchOutputTargetList")
        self.preview_list.setSelectionMode(QListWidget.SelectionMode.NoSelection)
        layout.addWidget(self.preview_list)
        self.error_label = QLabel()
        self.error_label.setObjectName("batchFitError")
        self.error_label.setStyleSheet("color: #b00020")
        layout.addWidget(self.error_label)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._accept_batch)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.prefix_edit.textChanged.connect(self._update_preview)
        self.selection.matrix_list.itemSelectionChanged.connect(
            lambda: self._update_preview(self.prefix_edit.text())
        )
        self._update_preview(self.prefix_edit.text())

    def _update_preview(self, text: str) -> None:
        self.error_label.clear()
        prefix = Path(text.strip()) if text.strip() else None
        try:
            self.target_paths = (
                batch_fit_target_paths(prefix, self.selection.selected_paths) if prefix else ()
            )
        except ValueError as error:
            self.target_paths = ()
            self.error_label.setText(str(error))
            return
        self.preview_list.clear()
        for path in self.target_paths:
            self.preview_list.addItem(str(path))

    def _accept_batch(self) -> None:
        selected = self.selection.selected_paths
        text = self.prefix_edit.text().strip()
        if not selected:
            self.error_label.setText("Select at least one matrix.")
            return
        if not text:
            self.error_label.setText("An output prefix is required.")
            return
        self.selected_paths = selected
        self.output_prefix = Path(text).expanduser()
        self.continue_on_failure = self.continue_check.isChecked()
        self._update_preview(str(self.output_prefix))
        if not self.target_paths:
            self.error_label.setText("The output prefix did not produce target files.")
            return
        self.accept()


def _target_paths(prefix: Path | None, suffixes: Sequence[str]) -> tuple[Path, ...]:
    if prefix is None:
        return ()
    paths = []
    for suffix in suffixes:
        base = prefix
        if base.name.endswith(suffix):
            base = base.with_name(base.name[: -len(suffix)])
        paths.append(base.with_name(base.name + suffix))
    return tuple(paths)


def output_target_paths(
    prefix: str | Path | None, suffixes: Sequence[str]
) -> tuple[Path, ...]:
    """Return exact files for an output prefix without touching the filesystem."""
    return _target_paths(Path(prefix) if prefix is not None else None, suffixes)


__all__ = [
    "BatchFitDialog",
    "MatrixSelectionDialog",
    "MatrixSelectionWidget",
    "OutputPathDialog",
    "output_target_paths",
]
