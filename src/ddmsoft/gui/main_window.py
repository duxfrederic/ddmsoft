"""Responsive Qt shell for the DDMSoft laboratory workflow."""

from __future__ import annotations

from collections.abc import Sequence
from itertools import pairwise
from pathlib import Path

from PySide6.QtCore import QSettings, Qt
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QProgressBar,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSizePolicy,
    QSlider,
    QSpinBox,
    QTableWidget,
    QVBoxLayout,
    QWidget,
)

from ..fitting import MODEL_REGISTRY


class DDMMainWindow(QMainWindow):
    """Native, calculation-free shell for the DDMSoft workflow."""

    def __init__(self, *, settings: QSettings | None = None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("ddmMainWindow")
        self.setWindowTitle("DDMSoft")
        self.setMinimumSize(760, 620)
        self._settings = settings if settings is not None else QSettings("DDMSoft", "DDMSoft")
        self._q_values: tuple[float, ...] = ()
        self._lag_times: tuple[float, ...] = ()

        self._build_menus()
        content = self._build_content()
        self._install_responsive_central_widget(content)
        self._set_tab_order()
        self._restore_window_state()

    def _build_menus(self) -> None:
        file_menu = self.menuBar().addMenu("File")
        self.open_directory_action = QAction("Open directory", self)
        self.open_directory_action.setObjectName("openDirectoryAction")
        self.open_directory_action.setShortcut(QKeySequence.StandardKey.Open)
        self.open_directory_action.setToolTip("Choose a directory containing DDM videos")
        file_menu.addAction(self.open_directory_action)
        file_menu.addSeparator()
        self.exit_action = QAction("Exit", self)
        self.exit_action.setObjectName("exitAction")
        self.exit_action.setShortcut(QKeySequence.StandardKey.Quit)
        self.exit_action.triggered.connect(self.close)
        file_menu.addAction(self.exit_action)

        tools_menu = self.menuBar().addMenu("Tools")
        self.concatenate_action = QAction("Concatenate videos", self)
        self.concatenate_action.setObjectName("concatenateVideosAction")
        self.split_action = QAction("Split a video in N matrices (time dependant DDM)", self)
        self.split_action.setObjectName("splitVideoAction")
        tools_menu.addActions((self.concatenate_action, self.split_action))

        plotting_menu = self.menuBar().addMenu("Plotting")
        self.show_matrix_action = QAction("Show the DDM matrix (image)", self)
        self.show_matrix_action.setObjectName("showMatrixAction")
        self.plot_correlation_action = QAction("Plot some autocorrelation functions", self)
        self.plot_correlation_action.setObjectName("plotCorrelationAction")
        plotting_menu.addActions((self.show_matrix_action, self.plot_correlation_action))

        export_menu = self.menuBar().addMenu("Batch, export")
        self.save_matrix_action = QAction("Save the DDM matrix as a text file", self)
        self.save_fit_action = QAction("Save the current fit parameters", self)
        self.fit_all_action = QAction("Fit and save all the matrices", self)
        self.save_correlation_action = QAction("Save the correlation functions", self)
        export_menu.addActions(
            (
                self.save_matrix_action,
                self.save_fit_action,
                self.fit_all_action,
                self.save_correlation_action,
            )
        )

        fitting_menu = self.menuBar().addMenu("More Fitting")
        self.contin_action = QAction("CONTIN", self)
        self.contin_action.setObjectName("continAction")
        fitting_menu.addAction(self.contin_action)

        help_menu = self.menuBar().addMenu("Help")
        self.about_action = QAction("About...", self)
        self.about_action.setObjectName("aboutAction")
        help_menu.addAction(self.about_action)

    def _build_content(self) -> QWidget:
        content = QWidget()
        content.setObjectName("mainContent")
        content.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        layout = QVBoxLayout(content)
        layout.setContentsMargins(12, 10, 12, 12)
        layout.setSpacing(10)

        layout.addWidget(self._build_computation_group())
        layout.addWidget(self._build_fitting_group())
        layout.addWidget(self._build_plotting_group())

        status_layout = QHBoxLayout()
        self.status_label = QLabel("Idle")
        self.status_label.setObjectName("statusLabel")
        self.status_label.setToolTip("What the backend is currently doing")
        self.status_label.setWordWrap(True)
        self.status_label.setMinimumWidth(180)
        self.status_label.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred
        )
        self.progress_bar = QProgressBar()
        self.progress_bar.setObjectName("progressBar")
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setAccessibleName("Processing progress")
        status_layout.addWidget(self.status_label)
        status_layout.addWidget(self.progress_bar, 1)
        layout.addLayout(status_layout)
        return content

    def _build_computation_group(self) -> QGroupBox:
        group = QGroupBox("DDM matrix computation")
        group.setObjectName("computationGroup")
        self.computation_group = group
        layout = QVBoxLayout(group)

        directory_layout = QHBoxLayout()
        directory_label = QLabel("Directory containing the videos")
        self.directory_edit = QLineEdit(self._last_directory())
        self.directory_edit.setObjectName("directoryEdit")
        directory_label.setBuddy(self.directory_edit)
        self.directory_edit.setToolTip("Directory containing videos and existing DDM matrices")
        self.browse_button = QPushButton("Browse...")
        self.browse_button.setObjectName("browseButton")
        self.browse_button.setToolTip("Choose the directory containing the videos")
        directory_layout.addWidget(directory_label)
        directory_layout.addWidget(self.directory_edit, 1)
        directory_layout.addWidget(self.browse_button)
        layout.addLayout(directory_layout)

        buttons_layout = QHBoxLayout()
        self.load_button = QPushButton("Load")
        self.load_button.setObjectName("loadButton")
        self.load_button.setToolTip("Load video metadata and existing matrices")
        self.quit_button = QPushButton("Quit")
        self.quit_button.setObjectName("quitButton")
        self.quit_button.setToolTip("Close DDMSoft")
        self.quit_button.clicked.connect(self.close)
        buttons_layout.addWidget(self.load_button)
        buttons_layout.addWidget(self.quit_button)
        buttons_layout.addStretch(1)
        layout.addLayout(buttons_layout)

        recompute_layout = QHBoxLayout()
        self.keep_existing_radio = QRadioButton("Keep existing matrices")
        self.keep_existing_radio.setObjectName("keepExistingRadio")
        self.keep_existing_radio.setChecked(True)
        self.keep_existing_radio.setToolTip("Reuse matrix files that are already present")
        self.recompute_radio = QRadioButton("Re-compute, overwrite existing")
        self.recompute_radio.setObjectName("recomputeRadio")
        self.recompute_radio.setToolTip("Recalculate and replace existing matrix files")
        recompute_layout.addWidget(self.keep_existing_radio)
        recompute_layout.addWidget(self.recompute_radio)
        recompute_layout.addStretch(1)
        layout.addLayout(recompute_layout)

        parameters = QFormLayout()
        self.max_couples_spin = QSpinBox()
        self.max_couples_spin.setObjectName("maxCouplesSpin")
        self.max_couples_spin.setRange(0, 1_000_000)
        self.max_couples_spin.setValue(300)
        self.max_couples_spin.setSpecialValueText("all available")
        self.max_couples_spin.setToolTip(
            "Maximum frame couples used during averaging; zero uses all available couples"
        )
        parameters.addRow("Maximum frame couples", self.max_couples_spin)

        self.points_per_decade_spin = QSpinBox()
        self.points_per_decade_spin.setObjectName("pointsPerDecadeSpin")
        self.points_per_decade_spin.setRange(1, 1_000)
        self.points_per_decade_spin.setValue(20)
        self.points_per_decade_spin.setToolTip("Number of lag times considered per decade")
        parameters.addRow("Lag times per decade", self.points_per_decade_spin)
        layout.addLayout(parameters)

        table_label = QLabel("Video metadata (editable)")
        table_label.setToolTip("Correct metadata values before processing")
        layout.addWidget(table_label)
        self.video_table = QTableWidget(0, 4)
        self.video_table.setObjectName("videoTable")
        self.video_table.setHorizontalHeaderLabels(
            ("Path / name", "Frame rate (fps)", "Pixel size (m)", "Validation")
        )
        self.video_table.setEditTriggers(QAbstractItemView.EditTrigger.AllEditTriggers)
        self.video_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.video_table.setAlternatingRowColors(True)
        self.video_table.setMinimumHeight(110)
        self.video_table.setToolTip("Editable path, frame-rate, and pixel-size metadata")
        header = self.video_table.horizontalHeader()
        header.setStretchLastSection(True)
        header.setSectionResizeMode(0, header.ResizeMode.Stretch)
        layout.addWidget(self.video_table)

        process_layout = QHBoxLayout()
        self.process_button = QPushButton("Process")
        self.process_button.setObjectName("processButton")
        self.process_button.setToolTip("Compute DDM matrices for the listed videos")
        self.direction_label = QLabel("Direction-dependent dynamics: split into")
        self.direction_spin = QSpinBox()
        self.direction_spin.setObjectName("directionSpin")
        self.direction_spin.setRange(1, 360)
        self.direction_spin.setValue(1)
        self.direction_spin.setToolTip(
            "Number of directional sectors; one uses the normal isotropic calculation"
        )
        self.direction_suffix = QLabel("parts")
        process_layout.addWidget(self.process_button)
        process_layout.addWidget(self.direction_label)
        process_layout.addWidget(self.direction_spin)
        process_layout.addWidget(self.direction_suffix)
        process_layout.addStretch(1)
        layout.addLayout(process_layout)
        return group

    def _build_fitting_group(self) -> QGroupBox:
        group = QGroupBox("DDM matrix fitting")
        group.setObjectName("fittingGroup")
        self.fitting_group = group
        layout = QVBoxLayout(group)

        merge_layout = QHBoxLayout()
        self.merge_button = QPushButton("Merge matrices")
        self.merge_button.setObjectName("mergeButton")
        self.merge_button.setToolTip("Merge selected matrices into one matrix")
        self.average_button = QPushButton("Average matrices")
        self.average_button.setObjectName("averageButton")
        self.average_button.setToolTip("Average selected compatible matrices")
        merge_layout.addWidget(self.merge_button)
        merge_layout.addWidget(self.average_button)
        merge_layout.addStretch(1)
        layout.addLayout(merge_layout)

        selectors_layout = QHBoxLayout()
        matrix_form = QFormLayout()
        matrix_label = QLabel("Matrices to analyse")
        self.matrix_selector = QComboBox()
        self.matrix_selector.setObjectName("matrixSelector")
        self.matrix_selector.addItem("No matrix loaded", None)
        self.matrix_selector.setToolTip("Select a discovered DDM matrix")
        matrix_form.addRow(matrix_label, self.matrix_selector)
        model_form = QFormLayout()
        model_label = QLabel("Model for the autocorrelation function")
        self.model_selector = QComboBox()
        self.model_selector.setObjectName("modelSelector")
        for identifier, definition in MODEL_REGISTRY.items():
            self.model_selector.addItem(definition.display_name, identifier)
        self.model_selector.setToolTip("Select the registered autocorrelation model")
        model_form.addRow(model_label, self.model_selector)
        selectors_layout.addLayout(matrix_form, 1)
        selectors_layout.addLayout(model_form, 1)
        layout.addLayout(selectors_layout)

        q_layout, self.q_min_slider, self.q_min_value, self.q_max_slider, self.q_max_value = (
            self._build_range_row("q", "q minimum", "q maximum", 0, 127, 0, 127)
        )
        layout.addLayout(q_layout)
        (
            time_layout,
            self.time_min_slider,
            self.time_min_value,
            self.time_max_slider,
            self.time_max_value,
        ) = self._build_range_row("time", "time minimum", "time maximum", 0, 100, 0, 100)
        layout.addLayout(time_layout)

        actions_layout = QHBoxLayout()
        self.initial_guess_button = QPushButton("Initial guess for the fit")
        self.initial_guess_button.setObjectName("initialGuessButton")
        self.initial_guess_button.setToolTip(
            "Set initial values and fixed parameters for the model"
        )
        self.fit_button = QPushButton("Fit the selected matrix")
        self.fit_button.setObjectName("fitButton")
        self.fit_button.setToolTip("Fit the selected matrix over the selected q and time ranges")
        self.fitted_parameters_button = QPushButton("Show the fitted parameters")
        self.fitted_parameters_button.setObjectName("fittedParametersButton")
        self.fitted_parameters_button.setToolTip("Inspect parameters from the current fit")
        actions_layout.addWidget(self.initial_guess_button)
        actions_layout.addWidget(self.fit_button)
        actions_layout.addWidget(self.fitted_parameters_button)
        layout.addLayout(actions_layout)
        self._update_range_labels()
        return group

    def _build_range_row(
        self,
        axis_name: str,
        minimum_name: str,
        maximum_name: str,
        minimum: int,
        maximum: int,
        minimum_value: int,
        maximum_value: int,
    ) -> tuple[QHBoxLayout, QSlider, QLabel, QSlider, QLabel]:
        layout = QHBoxLayout()
        lower_label = QLabel(f"{axis_name} min")
        lower_slider = QSlider(Qt.Orientation.Horizontal)
        lower_slider.setObjectName(f"{axis_name}MinSlider")
        lower_slider.setRange(minimum, maximum)
        lower_slider.setValue(minimum_value)
        lower_slider.setToolTip(f"Select the inclusive {minimum_name} index")
        lower_slider.setAccessibleName(minimum_name)
        lower_value = QLabel()
        lower_value.setObjectName(f"{axis_name}MinValue")
        lower_value.setMinimumWidth(100)
        upper_label = QLabel(f"{axis_name} max")
        upper_slider = QSlider(Qt.Orientation.Horizontal)
        upper_slider.setObjectName(f"{axis_name}MaxSlider")
        upper_slider.setRange(minimum, maximum)
        upper_slider.setValue(maximum_value)
        upper_slider.setToolTip(f"Select the inclusive {maximum_name} index")
        upper_slider.setAccessibleName(maximum_name)
        upper_value = QLabel()
        upper_value.setObjectName(f"{axis_name}MaxValue")
        upper_value.setMinimumWidth(100)
        layout.addWidget(lower_label)
        layout.addWidget(lower_slider, 1)
        layout.addWidget(lower_value)
        layout.addSpacing(8)
        layout.addWidget(upper_label)
        layout.addWidget(upper_slider, 1)
        layout.addWidget(upper_value)
        if axis_name == "q":
            lower_slider.valueChanged.connect(self._q_min_changed)
            upper_slider.valueChanged.connect(self._q_max_changed)
        else:
            lower_slider.valueChanged.connect(self._time_min_changed)
            upper_slider.valueChanged.connect(self._time_max_changed)
        return layout, lower_slider, lower_value, upper_slider, upper_value

    def _build_plotting_group(self) -> QGroupBox:
        group = QGroupBox("Plotting")
        group.setObjectName("plottingGroup")
        self.plotting_group = group
        layout = QVBoxLayout(group)

        fields = QFormLayout()
        temperature_label = QLabel("(Optional) temperature (deg C)")
        self.temperature_edit = QLineEdit()
        self.temperature_edit.setObjectName("temperatureEdit")
        self.temperature_edit.setToolTip(
            "Optional temperature used when converting diffusion to hydrodynamic radius"
        )
        temperature_label.setBuddy(self.temperature_edit)
        fields.addRow(temperature_label, self.temperature_edit)
        viscosity_label = QLabel("(Optional) viscosity")
        self.viscosity_edit = QLineEdit("water")
        self.viscosity_edit.setObjectName("viscosityEdit")
        self.viscosity_edit.setToolTip("Use water or enter the dynamic viscosity in SI units")
        viscosity_label.setBuddy(self.viscosity_edit)
        fields.addRow(viscosity_label, self.viscosity_edit)
        layout.addLayout(fields)

        buttons = QHBoxLayout()
        self.plot_matrix_button = QPushButton("Plot the matrix and the fit")
        self.plot_matrix_button.setObjectName("plotMatrixButton")
        self.plot_matrix_button.setToolTip("Open independent measured and fitted matrix plots")
        self.plot_amplitude_button = QPushButton(
            "Plot the amplitude, the noise, the diffusion"
        )
        self.plot_amplitude_button.setObjectName("plotAmplitudeButton")
        self.plot_amplitude_button.setToolTip("Open amplitude, noise, and diffusion plots")
        buttons.addWidget(self.plot_matrix_button)
        buttons.addWidget(self.plot_amplitude_button)
        buttons.addStretch(1)
        layout.addLayout(buttons)
        return group

    def _install_responsive_central_widget(self, content: QWidget) -> None:
        screen = self.screen()
        available_height = screen.availableGeometry().height() if screen is not None else 900
        available_width = screen.availableGeometry().width() if screen is not None else 1200
        if available_height < 820 or available_width < 1_000:
            scroll = QScrollArea()
            scroll.setObjectName("mainScrollArea")
            scroll.setFrameShape(QScrollArea.Shape.NoFrame)
            scroll.setWidgetResizable(True)
            scroll.setWidget(content)
            self.setCentralWidget(scroll)
        else:
            self.setCentralWidget(content)
        width = min(max(self.minimumWidth(), 1_020), max(available_width - 40, 760))
        height = min(900, max(available_height - 40, 620))
        self.resize(width, height)

    def _set_tab_order(self) -> None:
        widgets = (
            self.directory_edit,
            self.browse_button,
            self.load_button,
            self.quit_button,
            self.keep_existing_radio,
            self.recompute_radio,
            self.max_couples_spin,
            self.points_per_decade_spin,
            self.video_table,
            self.process_button,
            self.direction_spin,
            self.merge_button,
            self.average_button,
            self.matrix_selector,
            self.model_selector,
            self.q_min_slider,
            self.q_max_slider,
            self.time_min_slider,
            self.time_max_slider,
            self.initial_guess_button,
            self.fit_button,
            self.fitted_parameters_button,
            self.temperature_edit,
            self.viscosity_edit,
            self.plot_matrix_button,
            self.plot_amplitude_button,
        )
        for current, following in pairwise(widgets):
            self.setTabOrder(current, following)

    def _last_directory(self) -> str:
        return str(self._settings.value("last_directory", str(Path.home())))

    def _restore_window_state(self) -> None:
        geometry = self._settings.value("geometry")
        if geometry is not None:
            self.restoreGeometry(geometry)

    def _q_min_changed(self, value: int) -> None:
        if value >= self.q_max_slider.value():
            self._set_slider_value(self.q_min_slider, min(value, self.q_max_slider.value() - 1))
        self._update_range_labels()

    def _q_max_changed(self, value: int) -> None:
        if value <= self.q_min_slider.value():
            self._set_slider_value(self.q_max_slider, max(value, self.q_min_slider.value() + 1))
        self._update_range_labels()

    def _time_min_changed(self, value: int) -> None:
        if value >= self.time_max_slider.value():
            self._set_slider_value(
                self.time_min_slider, min(value, self.time_max_slider.value() - 1)
            )
        self._update_range_labels()

    def _time_max_changed(self, value: int) -> None:
        if value <= self.time_min_slider.value():
            self._set_slider_value(
                self.time_max_slider, max(value, self.time_min_slider.value() + 1)
            )
        self._update_range_labels()

    @staticmethod
    def _set_slider_value(slider: QSlider, value: int) -> None:
        value = max(slider.minimum(), min(slider.maximum(), value))
        slider.blockSignals(True)
        slider.setValue(value)
        slider.blockSignals(False)

    def _update_range_labels(self) -> None:
        self.q_min_value.setText(self._range_text(self.q_min_slider.value(), self._q_values, "q"))
        self.q_max_value.setText(self._range_text(self.q_max_slider.value(), self._q_values, "q"))
        self.time_min_value.setText(
            self._range_text(self.time_min_slider.value(), self._lag_times, "time")
        )
        self.time_max_value.setText(
            self._range_text(self.time_max_slider.value(), self._lag_times, "time")
        )

    @staticmethod
    def _range_text(index: int, values: tuple[float, ...], name: str) -> str:
        if values and index < len(values):
            return f"index: {index}; {name}: {values[index]:.6g}"
        return f"index: {index}; {name}: unavailable"

    def set_axis_values(self, q_values: Sequence[float], lag_times: Sequence[float]) -> None:
        """Update slider bounds and physical-value labels for the selected matrix."""
        self._q_values = tuple(float(value) for value in q_values)
        self._lag_times = tuple(float(value) for value in lag_times)
        self._set_axis_range(self.q_min_slider, self.q_max_slider, len(self._q_values))
        self._set_axis_range(self.time_min_slider, self.time_max_slider, len(self._lag_times))
        self._update_range_labels()

    def _set_axis_range(self, lower: QSlider, upper: QSlider, count: int) -> None:
        if count < 2:
            lower.setRange(0, max(count - 1, 0))
            upper.setRange(0, max(count - 1, 0))
            lower.setValue(0)
            upper.setValue(0)
            return
        maximum = count - 1
        lower_value = min(lower.value(), maximum - 1)
        upper_value = min(max(upper.value(), lower_value + 1), maximum)
        lower.blockSignals(True)
        upper.blockSignals(True)
        lower.setRange(0, maximum)
        upper.setRange(0, maximum)
        lower.setValue(lower_value)
        upper.setValue(upper_value)
        lower.blockSignals(False)
        upper.blockSignals(False)

    def closeEvent(self, event: object) -> None:
        self._settings.setValue("geometry", self.saveGeometry())
        self._settings.setValue("last_directory", self.directory_edit.text())
        self._settings.sync()
        super().closeEvent(event)


def create_main_window(*, settings: QSettings | None = None) -> DDMMainWindow:
    """Create a shell window without starting the Qt event loop."""
    return DDMMainWindow(settings=settings)


__all__ = ["DDMMainWindow", "create_main_window"]
