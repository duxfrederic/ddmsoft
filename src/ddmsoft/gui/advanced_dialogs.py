"""Validated dialogs for the advanced DDM workflows."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from math import isfinite
from pathlib import Path

from PySide6.QtCore import QItemSelectionModel, Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from ..contin import CONTINError, contin_ranges
from ..models import VideoMetadata


class VideoSelectionWidget(QWidget):
    """Searchable multi-selection list for video paths."""

    def __init__(
        self,
        videos: Mapping[Path, str],
        *,
        selected_paths: Sequence[Path] = (),
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        selected = {Path(path) for path in selected_paths}
        layout = QVBoxLayout(self)
        self.search_edit = QLineEdit()
        self.search_edit.setObjectName("videoSearchEdit")
        self.search_edit.setPlaceholderText("Search video names or paths")
        layout.addWidget(self.search_edit)
        self.video_list = QListWidget()
        self.video_list.setObjectName("videoSelectionList")
        self.video_list.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)
        for path, name in sorted(videos.items(), key=lambda item: (item[1].casefold(), str(item[0]))):
            item = QListWidgetItem(name)
            item.setData(Qt.ItemDataRole.UserRole, Path(path))
            item.setToolTip(str(path))
            self.video_list.addItem(item)
            item.setSelected(Path(path) in selected)
        layout.addWidget(self.video_list, 1)
        actions = QHBoxLayout()
        self.select_all_button = QPushButton("Select All")
        self.select_all_button.setObjectName("selectAllVideosButton")
        self.clear_button = QPushButton("Clear")
        self.clear_button.setObjectName("clearVideoSelectionButton")
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
            Path(item.data(Qt.ItemDataRole.UserRole)) for item in self.video_list.selectedItems()
        )

    def select_all_visible(self) -> None:
        for index in range(self.video_list.count()):
            item = self.video_list.item(index)
            item.setSelected(not item.isHidden())

    def clear_selection(self) -> None:
        self.video_list.clearSelection()

    def _filter_items(self, text: str) -> None:
        query = text.casefold().strip()
        for index in range(self.video_list.count()):
            item = self.video_list.item(index)
            path = Path(item.data(Qt.ItemDataRole.UserRole))
            item.setHidden(bool(query) and query not in f"{item.text()} {path}".casefold())


class VideoSelectionDialog(QDialog):
    """Select videos without changing application state until acceptance."""

    def __init__(
        self,
        videos: Mapping[Path, str],
        *,
        minimum: int = 1,
        title: str = "Select videos",
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setMinimumSize(620, 420)
        self.minimum = minimum
        self.selection = VideoSelectionWidget(videos)
        self.selected_paths: tuple[Path, ...] = ()
        layout = QVBoxLayout(self)
        layout.addWidget(self.selection, 1)
        self.error_label = QLabel()
        self.error_label.setObjectName("videoSelectionError")
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
        if len(selected) < self.minimum:
            self.error_label.setText(f"Select at least {self.minimum} video(s).")
            return
        self.selected_paths = selected
        self.accept()


class TimeDependentDialog(QDialog):
    """Select videos and validate their per-video metadata and partitions."""

    def __init__(
        self,
        videos: Mapping[Path, VideoMetadata],
        *,
        partitions: int = 10,
        title: str = "Time-dependent DDM",
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setMinimumSize(760, 540)
        self.selected_videos: tuple[VideoMetadata, ...] = ()
        self.partitions = partitions
        self._metadata = {Path(path): metadata for path, metadata in videos.items()}
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Select videos and verify each frame rate and pixel size."))
        self.search_edit = QLineEdit()
        self.search_edit.setObjectName("timeDependentSearchEdit")
        self.search_edit.setPlaceholderText("Search video names or paths")
        layout.addWidget(self.search_edit)
        self.video_table = QTableWidget(0, 3)
        self.video_table.setObjectName("timeDependentVideoTable")
        self.video_table.setHorizontalHeaderLabels(("Video", "Frame rate (fps)", "Pixel size (m)"))
        self.video_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.video_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        for row, (path, metadata) in enumerate(
            sorted(videos.items(), key=lambda item: str(item[0]).casefold())
        ):
            self.video_table.insertRow(row)
            path_item = QTableWidgetItem(str(path))
            path_item.setData(Qt.ItemDataRole.UserRole, Path(path))
            path_item.setFlags(path_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.video_table.setItem(row, 0, path_item)
            self.video_table.setItem(row, 1, QTableWidgetItem(f"{metadata.frame_rate:.12g}"))
            self.video_table.setItem(row, 2, QTableWidgetItem(f"{metadata.pixel_size:.12g}"))
        self.video_table.horizontalHeader().setStretchLastSection(True)
        self.video_table.horizontalHeader().setSectionResizeMode(
            0, self.video_table.horizontalHeader().ResizeMode.Stretch
        )
        layout.addWidget(self.video_table, 1)
        partition_form = QFormLayout()
        self.partition_spin = QSpinBox()
        self.partition_spin.setObjectName("partitionCountSpin")
        self.partition_spin.setRange(1, 10_000)
        self.partition_spin.setValue(partitions)
        self.partition_spin.setToolTip("Each partition must contain at least two source frames")
        partition_form.addRow("Number of portions", self.partition_spin)
        layout.addLayout(partition_form)
        actions = QHBoxLayout()
        self.select_all_button = QPushButton("Select All")
        self.clear_button = QPushButton("Clear")
        actions.addWidget(self.select_all_button)
        actions.addWidget(self.clear_button)
        actions.addStretch(1)
        layout.addLayout(actions)
        self.error_label = QLabel()
        self.error_label.setObjectName("timeDependentError")
        self.error_label.setStyleSheet("color: #b00020")
        layout.addWidget(self.error_label)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._accept_time_dependent)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.search_edit.textChanged.connect(self._filter_rows)
        self.select_all_button.clicked.connect(self._select_all_visible)
        self.clear_button.clicked.connect(self.video_table.clearSelection)
        self._select_all_visible()

    def _accept_time_dependent(self) -> None:
        rows = sorted({index.row() for index in self.video_table.selectionModel().selectedRows()})
        if not rows:
            self.error_label.setText("Select at least one video.")
            return
        records: list[VideoMetadata] = []
        try:
            for row in rows:
                path = Path(self.video_table.item(row, 0).data(Qt.ItemDataRole.UserRole))
                frame_rate = _positive_number(self.video_table.item(row, 1).text(), "frame rate")
                pixel_size = _positive_number(self.video_table.item(row, 2).text(), "pixel size")
                records.append(
                    VideoMetadata(
                        path,
                        frame_rate,
                        pixel_size,
                        self._metadata[path].temperature,
                    )
                )
        except ValueError as error:
            self.error_label.setText(str(error))
            return
        self.selected_videos = tuple(records)
        self.partitions = self.partition_spin.value()
        self.accept()

    def _select_all_visible(self) -> None:
        self.video_table.clearSelection()
        selection_model = self.video_table.selectionModel()
        for row in range(self.video_table.rowCount()):
            if not self.video_table.isRowHidden(row):
                selection_model.select(
                    self.video_table.model().index(row, 0),
                    QItemSelectionModel.SelectionFlag.Select
                    | QItemSelectionModel.SelectionFlag.Rows,
                )

    def _filter_rows(self, text: str) -> None:
        query = text.casefold().strip()
        for row in range(self.video_table.rowCount()):
            path = self.video_table.item(row, 0).text()
            self.video_table.setRowHidden(row, bool(query) and query not in path.casefold())


class CONTINDialog(QDialog):
    """Validated controls for one q-index CONTIN scan."""

    def __init__(
        self,
        q_count: int,
        *,
        q_index: int = 0,
        title: str = "CONTIN",
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.q_index = q_index
        self.gamma_min = 1e-13
        self.gamma_max = 1e-11
        self.gamma_count = 10
        self.alpha_min = 0.01
        self.alpha_max = 1.0
        self.alpha_count = 3
        self.maxiter = 10
        self.save_all = True
        layout = QVBoxLayout(self)
        form = QFormLayout()
        self.q_index_spin = QSpinBox()
        self.q_index_spin.setObjectName("continQIndexSpin")
        self.q_index_spin.setRange(0, max(q_count - 1, 0))
        self.q_index_spin.setValue(max(0, min(q_index, max(q_count - 1, 0))))
        form.addRow("q index", self.q_index_spin)
        self.gamma_min_edit = _line_edit("1e-13", "continGammaMinEdit")
        self.gamma_max_edit = _line_edit("1e-11", "continGammaMaxEdit")
        self.gamma_count_spin = _count_spin(10, 3, "continGammaCountSpin")
        self.alpha_min_edit = _line_edit("0.01", "continAlphaMinEdit")
        self.alpha_max_edit = _line_edit("1", "continAlphaMaxEdit")
        self.alpha_count_spin = _count_spin(3, 1, "continAlphaCountSpin")
        self.maxiter_spin = _count_spin(10, 1, "continMaxiterSpin")
        form.addRow("Smallest decay rate", self.gamma_min_edit)
        form.addRow("Biggest decay rate", self.gamma_max_edit)
        form.addRow("Number of decay rates", self.gamma_count_spin)
        form.addRow("Smallest alpha", self.alpha_min_edit)
        form.addRow("Biggest alpha", self.alpha_max_edit)
        form.addRow("Number of alpha candidates", self.alpha_count_spin)
        form.addRow("Maximum iterations", self.maxiter_spin)
        layout.addLayout(form)
        layout.addWidget(QLabel("Gamma and alpha use legacy linear spacing."))
        self.save_all_check = QCheckBox("Export all alpha candidates")
        self.save_all_check.setObjectName("continSaveAllCheck")
        self.save_all_check.setChecked(True)
        layout.addWidget(self.save_all_check)
        self.error_label = QLabel()
        self.error_label.setObjectName("continError")
        self.error_label.setStyleSheet("color: #b00020")
        layout.addWidget(self.error_label)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._accept_contin)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _accept_contin(self) -> None:
        try:
            gamma_min = _finite_number(self.gamma_min_edit.text(), "gamma minimum")
            gamma_max = _finite_number(self.gamma_max_edit.text(), "gamma maximum")
            alpha_min = _finite_number(self.alpha_min_edit.text(), "alpha minimum")
            alpha_max = _finite_number(self.alpha_max_edit.text(), "alpha maximum")
            gamma, alphas = contin_ranges(
                gamma_min,
                gamma_max,
                self.gamma_count_spin.value(),
                alpha_min,
                alpha_max,
                self.alpha_count_spin.value(),
            )
        except (CONTINError, ValueError) as error:
            self.error_label.setText(str(error))
            return
        self.q_index = self.q_index_spin.value()
        self.gamma_min, self.gamma_max = float(gamma[0]), float(gamma[-1])
        self.gamma_count = int(gamma.size)
        self.alpha_min, self.alpha_max = float(alphas[0]), float(alphas[-1])
        self.alpha_count = int(alphas.size)
        self.maxiter = self.maxiter_spin.value()
        self.save_all = self.save_all_check.isChecked()
        self.accept()


def _line_edit(value: str, object_name: str) -> QLineEdit:
    edit = QLineEdit(value)
    edit.setObjectName(object_name)
    return edit


def _count_spin(value: int, minimum: int, object_name: str) -> QSpinBox:
    spin = QSpinBox()
    spin.setObjectName(object_name)
    spin.setRange(minimum, 10_000)
    spin.setValue(value)
    return spin


def _finite_number(text: str, name: str) -> float:
    try:
        value = float(text.strip())
    except ValueError as error:
        raise ValueError(f"{name} must be numeric") from error
    if not isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _positive_number(text: str, name: str) -> float:
    value = _finite_number(text, name)
    if value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


__all__ = [
    "CONTINDialog",
    "TimeDependentDialog",
    "VideoSelectionDialog",
    "VideoSelectionWidget",
]
