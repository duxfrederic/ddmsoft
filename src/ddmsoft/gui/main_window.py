"""Responsive Qt shell for the DDMSoft laboratory workflow."""

from __future__ import annotations

from collections.abc import Sequence
from itertools import pairwise
from math import isfinite
from pathlib import Path

from PySide6.QtCore import QSettings, Qt, QThread
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSizePolicy,
    QSlider,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from ..fitting import MODEL_REGISTRY, default_fit_request, get_model
from ..io import (
    DDMIOError,
    MatrixFileSet,
    discover_matrix_sets,
    display_names,
    load_directory,
    load_matrices,
)
from ..models import DDMData, FitRange, FitRequest, FitResult, VideoMetadata
from ..plotting import (
    AmplitudeNoiseDiffusionPlotController,
    CorrelationPlotController,
    FitParameterPlotController,
    MatrixPlotController,
)
from .fit_dialog import InitialGuessDialog
from .workers import (
    ComputationResult,
    ComputationWorker,
    FitComputationRequest,
    FitComputationResult,
    VideoComputationRequest,
    run_fit,
    run_video_computation,
)


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
        self._matrices: dict[Path, DDMData] = {}
        self._fit_results: dict[Path, FitResult] = {}
        self._fit_requests: dict[Path, FitRequest] = {}
        self._fit_preferences: dict[
            str, tuple[tuple[float | None, ...], tuple[bool, ...]]
        ] = {}
        self._plot_controllers: list[object] = []
        self._updating_video_table = False
        self._selected_matrix_path: Path | None = None
        self._thread: QThread | None = None
        self._worker: ComputationWorker | None = None
        self._job_active = False
        self._job_kind: str | None = None
        self._progress_value = 0

        self._build_menus()
        content = self._build_content()
        self._install_responsive_central_widget(content)
        self._set_tab_order()
        self._connect_workflow_signals()
        self._update_processing_state()
        self._update_matrix_action_state()
        self._restore_window_state()

    def _connect_workflow_signals(self) -> None:
        self.open_directory_action.triggered.connect(
            lambda: self.choose_directory(load_after_select=True)
        )
        self.browse_button.clicked.connect(
            lambda: self.choose_directory(load_after_select=False)
        )
        self.load_button.clicked.connect(self.load_current_directory)
        self.process_button.clicked.connect(self.start_processing)
        self.cancel_button.clicked.connect(self.cancel_processing)
        self.video_table.cellChanged.connect(self._video_cell_changed)
        self.matrix_selector.currentIndexChanged.connect(self._matrix_selection_changed)
        self.initial_guess_button.clicked.connect(self.edit_initial_guess)
        self.fit_button.clicked.connect(self.start_fitting)
        self.fitted_parameters_button.clicked.connect(self.plot_fitted_parameters)
        self.plot_matrix_button.clicked.connect(self.plot_selected_matrix)
        self.plot_amplitude_button.clicked.connect(self.plot_amplitude_noise_diffusion)
        self.show_matrix_action.triggered.connect(self.plot_selected_matrix)
        self.plot_correlation_action.triggered.connect(self.plot_selected_correlation)

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
        self.save_matrix_action.setObjectName("saveMatrixAction")
        self.save_fit_action = QAction("Save the current fit parameters", self)
        self.save_fit_action.setObjectName("saveFitAction")
        self.fit_all_action = QAction("Fit and save all the matrices", self)
        self.fit_all_action.setObjectName("fitAllAction")
        self.save_correlation_action = QAction("Save the correlation functions", self)
        self.save_correlation_action.setObjectName("saveCorrelationAction")
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

        self.error_details = QPlainTextEdit()
        self.error_details.setObjectName("errorDetails")
        self.error_details.setReadOnly(True)
        self.error_details.setPlaceholderText("Worker error details will appear here.")
        self.error_details.setMaximumHeight(180)
        self.error_details.setVisible(False)
        layout.addWidget(self.error_details)

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
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.setObjectName("cancelButton")
        self.cancel_button.setToolTip("Request cooperative cancellation of the active computation")
        self.cancel_button.setEnabled(False)
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
        process_layout.addWidget(self.cancel_button)
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

    def choose_directory(self, *, load_after_select: bool = True) -> str | None:
        """Choose an acquisition directory and optionally load it immediately."""
        selected = QFileDialog.getExistingDirectory(
            self,
            "Open DDM directory",
            self.directory_edit.text() or self._last_directory(),
        )
        if not selected:
            return None
        self.directory_edit.setText(selected)
        if load_after_select:
            self.load_directory_data(selected)
        return selected

    def load_current_directory(self) -> bool:
        """Load the directory currently shown in the directory field."""
        return self.load_directory_data(self.directory_edit.text())

    def start_processing(self) -> None:
        """Start a sequential background computation from current UI values."""
        if self._job_active:
            return
        try:
            videos = self.video_metadata()
        except ValueError as error:
            self.status_label.setText(str(error))
            self.status_label.setToolTip(str(error))
            return
        request = VideoComputationRequest(
            videos=videos,
            max_couples=self.max_couples_spin.value(),
            points_per_decade=self.points_per_decade_spin.value(),
            sectors=self.direction_spin.value(),
            recompute=self.recompute_radio.isChecked(),
        )
        self.error_details.clear()
        self.error_details.setVisible(False)
        self._progress_value = 0
        self.progress_bar.setValue(0)
        self._job_active = True
        self._job_kind = "processing"
        self._set_job_active(True)

        thread = QThread(self)
        worker = ComputationWorker(
            lambda progress, cancel: run_video_computation(request, progress, cancel)
        )
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.progress.connect(self._worker_progress)
        worker.status.connect(self._worker_status)
        worker.result.connect(self._worker_result)
        worker.failure.connect(self._worker_failure)
        worker.cancelled.connect(self._worker_cancelled)
        worker.finished.connect(thread.quit, Qt.ConnectionType.DirectConnection)
        worker.finished.connect(worker.deleteLater)
        thread.finished.connect(self._thread_finished)
        self._thread = thread
        self._worker = worker
        thread.start()

    def edit_initial_guess(self) -> bool:
        """Edit and retain the initial values for the selected model."""
        model_id = self.model_selector.currentData()
        if not isinstance(model_id, str):
            self.status_label.setText("Select a fit model first")
            return False
        model = get_model(model_id)
        values, fixed = self._initial_guess_state(model_id)
        dialog = InitialGuessDialog(model, values, fixed, parent=self)
        self._initial_guess_dialog = dialog
        if dialog.exec() != dialog.DialogCode.Accepted:
            return False
        self._fit_preferences[model_id] = (dialog.initial_values, dialog.fixed_flags)
        self.status_label.setText(f"Initial guess saved for {model.display_name}")
        return True

    def start_fitting(self) -> None:
        """Fit the selected matrix over the current inclusive slider ranges."""
        if self._job_active:
            return
        path = self._selected_matrix_path
        data = self.selected_matrix
        try:
            request = self.current_fit_request()
        except (TypeError, ValueError) as error:
            self.status_label.setText(str(error))
            self.status_label.setToolTip(str(error))
            return
        if path is None or data is None:
            self.status_label.setText("No matrix selected")
            return

        computation_request = FitComputationRequest(path, data, request)
        self._fit_results.pop(path, None)
        self._fit_requests.pop(path, None)
        self.error_details.clear()
        self.error_details.setVisible(False)
        self._progress_value = 0
        self.progress_bar.setValue(0)
        self._job_active = True
        self._job_kind = "fitting"
        self._set_job_active(True)

        thread = QThread(self)
        worker = ComputationWorker(
            lambda progress, cancel: run_fit(computation_request, progress, cancel)
        )
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.progress.connect(self._worker_progress)
        worker.status.connect(self._worker_status)
        worker.result.connect(self._fit_worker_result)
        worker.failure.connect(self._worker_failure)
        worker.cancelled.connect(self._worker_cancelled)
        worker.finished.connect(thread.quit, Qt.ConnectionType.DirectConnection)
        worker.finished.connect(worker.deleteLater)
        thread.finished.connect(self._thread_finished)
        self._thread = thread
        self._worker = worker
        thread.start()

    def current_fit_request(self) -> FitRequest:
        """Build the current model request from the UI state."""
        data = self.selected_matrix
        model_id = self.model_selector.currentData()
        if data is None:
            raise ValueError("no matrix selected")
        if not isinstance(model_id, str):
            raise TypeError("no fit model selected")
        fit_range = FitRange(
            self.q_min_slider.value(),
            self.q_max_slider.value(),
            self.time_min_slider.value(),
            self.time_max_slider.value(),
        )
        values, fixed = self._initial_guess_state(model_id)
        return FitRequest(model_id, values, fixed, fit_range)

    def _initial_guess_state(
        self, model_id: str
    ) -> tuple[tuple[float | None, ...], tuple[bool, ...]]:
        state = self._fit_preferences.get(model_id)
        if state is not None:
            return state
        request = default_fit_request(
            model_id,
            FitRange(0, 0, 0, 0),
        )
        return request.initial_values, request.fixed_flags

    def cancel_processing(self) -> None:
        """Request cooperative cancellation of the active computation."""
        if self._worker is None or not self._job_active:
            return
        self._worker.request_cancel()
        self.cancel_button.setEnabled(False)
        self.status_label.setText("Cancellation requested")

    def _set_job_active(self, active: bool) -> None:
        widgets = (
            self.directory_edit,
            self.browse_button,
            self.load_button,
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
        for widget in widgets:
            widget.setEnabled(not active)
        for action in (
            self.open_directory_action,
            self.concatenate_action,
            self.split_action,
            self.show_matrix_action,
            self.plot_correlation_action,
            self.save_matrix_action,
            self.save_fit_action,
            self.fit_all_action,
            self.save_correlation_action,
            self.contin_action,
            self.about_action,
        ):
            action.setEnabled(not active)
        self.cancel_button.setEnabled(active)

    def _worker_progress(self, stage: str, completed: int, total: int) -> None:
        if total <= 0:
            return
        value = max(0, min(100, int(100 * completed / total)))
        self._progress_value = max(self._progress_value, value)
        self.progress_bar.setValue(self._progress_value)

    def _worker_status(self, status: str) -> None:
        self.status_label.setText(status)

    def _worker_result(self, value: object) -> None:
        if not isinstance(value, ComputationResult):
            self._worker_failure("worker returned an invalid result", repr(value))
            return
        self.progress_bar.setValue(100)
        self._progress_value = 100
        preferred = value.paths[0] if value.paths else None
        self.load_directory_data(
            Path(self.directory_edit.text()), preferred_matrix_path=preferred
        )
        self.status_label.setText(
            f"Completed {len(value.processed_videos)} video(s); "
            f"kept {len(value.kept_videos)} existing set(s)"
        )

    def _fit_worker_result(self, value: object) -> None:
        if not isinstance(value, FitComputationResult):
            self._worker_failure("worker returned an invalid fit result", repr(value))
            return
        if value.matrix_path != self._selected_matrix_path:
            self.status_label.setText("Ignored a fit result for a different matrix")
            return
        self._fit_results[value.matrix_path] = value.fit
        self._fit_requests[value.matrix_path] = value.fit_request
        failed = sum(not status for status in value.fit.convergence_status)
        self.progress_bar.setValue(100)
        self._progress_value = 100
        self.status_label.setText(
            f"Fit complete for {value.matrix_path.name}; {failed} q fit(s) did not converge"
        )
        self._update_matrix_action_state()

    def _worker_failure(self, message: str, details: str) -> None:
        self.error_details.setPlainText(details)
        self.error_details.setVisible(True)
        operation = "Fitting" if self._job_kind == "fitting" else "Processing"
        self.status_label.setText(f"{operation} failed: {message}")
        self.status_label.setToolTip(details)

    def _worker_cancelled(self) -> None:
        operation = "Fitting" if self._job_kind == "fitting" else "Processing"
        self.status_label.setText(f"{operation} cancelled")

    def _thread_finished(self) -> None:
        self._thread = None
        self._worker = None
        self._job_active = False
        self._job_kind = None
        self._set_job_active(False)
        self._update_processing_state()
        self._update_matrix_action_state()

    def load_directory_data(
        self, directory: str | Path, *, preferred_matrix_path: Path | None = None
    ) -> bool:
        """Load metadata and matrix catalog data into the window.

        The method returns whether all video metadata loaded successfully. Matrix
        discovery is performed independently so an invalid metadata file does not
        hide already available matrices.
        """
        root = Path(directory).expanduser()
        self.directory_edit.setText(str(root))
        self._settings.setValue("last_directory", str(root))

        metadata_error: DDMIOError | None = None
        try:
            metadata = load_directory(root)
        except DDMIOError as error:
            metadata = {}
            metadata_error = error
            self._populate_invalid_video_rows(root, error)
        else:
            self._populate_video_table(metadata)

        try:
            matrix_sets = discover_matrix_sets(root, strict=False)
            matrices = load_matrices(root, strict=False)
        except DDMIOError as error:
            self._clear_matrix_catalog()
            matrix_error = error
        else:
            matrix_error = None
            self._populate_matrix_catalog(matrix_sets, matrices, preferred_matrix_path)

        if metadata_error is not None:
            self.status_label.setText(f"Metadata invalid: {metadata_error}")
            self.status_label.setToolTip(str(metadata_error))
        elif matrix_error is not None:
            self.status_label.setText(f"Matrix loading failed: {matrix_error}")
            self.status_label.setToolTip(str(matrix_error))
        else:
            self.status_label.setText(
                f"Loaded {len(metadata)} video(s) and {len(self._matrices)} matrix/matrices"
            )
            self.status_label.setToolTip(str(root))
        self._update_processing_state()
        self._update_matrix_action_state()
        return metadata_error is None

    def _populate_video_table(self, metadata: dict[Path, VideoMetadata]) -> None:
        rows = sorted(metadata.values(), key=lambda item: str(item.path).casefold())
        self._updating_video_table = True
        try:
            self.video_table.setRowCount(0)
            for row, record in enumerate(rows):
                self.video_table.insertRow(row)
                self._set_video_item(row, 0, str(record.path))
                self._set_video_item(row, 1, f"{record.frame_rate:.12g}")
                self._set_video_item(row, 2, f"{record.pixel_size:.12g}")
                self._set_video_item(row, 3, "Valid", editable=False)
        finally:
            self._updating_video_table = False
        for row in range(self.video_table.rowCount()):
            self._validate_video_row(row)

    def _populate_invalid_video_rows(self, directory: Path, error: DDMIOError) -> None:
        try:
            videos = sorted(
                (path for path in directory.iterdir() if path.is_file() and path.suffix.lower() == ".avi"),
                key=lambda path: path.name.casefold(),
            )
        except OSError:
            videos = []
        self._updating_video_table = True
        try:
            self.video_table.setRowCount(0)
            for row, video in enumerate(videos):
                self.video_table.insertRow(row)
                self._set_video_item(row, 0, str(video))
                self._set_video_item(row, 1, "")
                self._set_video_item(row, 2, "")
                self._set_video_item(row, 3, f"Invalid metadata: {error}", editable=False)
                for column in range(1, 4):
                    item = self.video_table.item(row, column)
                    if item is not None:
                        item.setToolTip(str(error))
        finally:
            self._updating_video_table = False

    def _set_video_item(self, row: int, column: int, text: str, *, editable: bool = True) -> None:
        item = QTableWidgetItem(text)
        if not editable:
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
        self.video_table.setItem(row, column, item)

    def _video_cell_changed(self, row: int, column: int) -> None:
        if self._updating_video_table or column == 3:
            return
        valid, field, message = self._validate_video_row(row)
        if valid:
            self._update_processing_state()
            return
        self.status_label.setText(f"Metadata error in row {row + 1}, {field}: {message}")
        self.status_label.setToolTip(message)
        self._update_processing_state()

    def _validate_video_row(
        self, row: int, *, update_state: bool = True
    ) -> tuple[bool, str, str]:
        fields = (
            (0, "path / name"),
            (1, "frame rate"),
            (2, "pixel size"),
        )
        for column, field in fields:
            item = self.video_table.item(row, column)
            text = item.text().strip() if item is not None else ""
            if not text:
                return self._validation_result(
                    row, column, field, "value is required", update_state
                )
            if column == 0:
                if Path(text).suffix.lower() != ".avi":
                    return self._validation_result(
                        row, column, field, "must name an AVI file", update_state
                    )
                continue
            try:
                value = float(text)
            except ValueError:
                return self._validation_result(row, column, field, "must be numeric", update_state)
            if not isfinite(value) or value <= 0:
                return self._validation_result(
                    row, column, field, "must be finite and positive", update_state
                )
        if update_state:
            self._set_row_validation(row, 3, "", "")
        return True, "", ""

    def _validation_result(
        self, row: int, column: int, field: str, message: str, update_state: bool
    ) -> tuple[bool, str, str]:
        if update_state:
            self._set_row_validation(row, column, field, message)
        return False, field, message

    def _set_row_validation(
        self, row: int, column: int, field: str, message: str
    ) -> tuple[bool, str, str]:
        state_item = self.video_table.item(row, 3)
        if state_item is None:
            self._updating_video_table = True
            try:
                self._set_video_item(row, 3, "", editable=False)
            finally:
                self._updating_video_table = False
            state_item = self.video_table.item(row, 3)
        if state_item is not None:
            state_item.setText("Valid" if not message else f"Invalid: {field} ({message})")
            state_item.setToolTip("" if not message else f"{field}: {message}")
        for index in range(3):
            item = self.video_table.item(row, index)
            if item is not None:
                item.setToolTip(f"{field}: {message}" if index == column and message else "")
        return not message, field, message

    def video_metadata(self) -> tuple[VideoMetadata, ...]:
        """Return validated, structured metadata currently shown in the table."""
        records: list[VideoMetadata] = []
        errors: list[str] = []
        for row in range(self.video_table.rowCount()):
            valid, field, message = self._validate_video_row(row)
            if not valid:
                errors.append(f"row {row + 1}, {field}: {message}")
                continue
            path = Path(self.video_table.item(row, 0).text().strip())
            frame_rate = float(self.video_table.item(row, 1).text())
            pixel_size = float(self.video_table.item(row, 2).text())
            records.append(VideoMetadata(path, frame_rate, pixel_size))
        if errors:
            raise ValueError("invalid video metadata: " + "; ".join(errors))
        if not records:
            raise ValueError("no video metadata is loaded")
        return tuple(records)

    def _update_processing_state(self) -> None:
        if self._job_active:
            self.process_button.setEnabled(False)
            self.concatenate_action.setEnabled(False)
            self.split_action.setEnabled(False)
            return
        metadata_valid = self.video_table.rowCount() > 0 and all(
            self._validate_video_row(row, update_state=False)[0]
            for row in range(self.video_table.rowCount())
        )
        self.process_button.setEnabled(metadata_valid)
        self.concatenate_action.setEnabled(metadata_valid and self.video_table.rowCount() > 1)
        self.split_action.setEnabled(metadata_valid and self.video_table.rowCount() == 1)

    def _populate_matrix_catalog(
        self,
        matrix_sets: Sequence[MatrixFileSet],
        matrices: dict[Path, DDMData],
        preferred_matrix_path: Path | None = None,
    ) -> None:
        names = display_names(matrix_sets)
        self._matrices = dict(matrices)
        self._fit_results.clear()
        self._fit_requests.clear()
        self.matrix_selector.blockSignals(True)
        try:
            self.matrix_selector.clear()
            self.matrix_selector.addItem("No matrix loaded", None)
            for name, path in names.items():
                if path in self._matrices:
                    self.matrix_selector.addItem(name, path)
            selected_index = 1
            if preferred_matrix_path is not None:
                for index in range(1, self.matrix_selector.count()):
                    if self.matrix_selector.itemData(index) == preferred_matrix_path:
                        selected_index = index
                        break
            if self.matrix_selector.count() > 1:
                self.matrix_selector.setCurrentIndex(selected_index)
        finally:
            self.matrix_selector.blockSignals(False)
        self._matrix_selection_changed(self.matrix_selector.currentIndex())

    def _clear_matrix_catalog(self) -> None:
        self._matrices = {}
        self._fit_results.clear()
        self._fit_requests.clear()
        self._selected_matrix_path = None
        self.matrix_selector.blockSignals(True)
        try:
            self.matrix_selector.clear()
            self.matrix_selector.addItem("No matrix loaded", None)
        finally:
            self.matrix_selector.blockSignals(False)
        self.set_axis_values((), ())
        self._update_matrix_action_state()

    def _matrix_selection_changed(self, index: int) -> None:
        selected = self.matrix_selector.itemData(index)
        path = Path(selected) if selected is not None else None
        data = self._matrices.get(path) if path is not None else None
        self._selected_matrix_path = path if data is not None else None
        if data is None:
            self.set_axis_values((), ())
            self.status_label.setText("No matrix selected")
        else:
            self.set_axis_values(data.q_values, data.lag_times)
            fit_status = "; fit available" if path in self._fit_results else ""
            self.status_label.setText(f"Selected matrix: {path.name}{fit_status}")
        self._update_matrix_action_state()

    @property
    def selected_matrix(self) -> DDMData | None:
        """Return the selected matrix for later workflow stages."""
        return (
            self._matrices.get(self._selected_matrix_path)
            if self._selected_matrix_path is not None
            else None
        )

    @property
    def selected_matrix_path(self) -> Path | None:
        """Return the full path backing the selected matrix display name."""
        return self._selected_matrix_path

    @property
    def selected_fit(self) -> FitResult | None:
        """Return the fit attached to the selected full matrix path, if any."""
        return (
            self._fit_results.get(self._selected_matrix_path)
            if self._selected_matrix_path is not None
            else None
        )

    @property
    def selected_fit_range(self) -> FitRange | None:
        """Return the inclusive range used by the selected fit, if any."""
        request = (
            self._fit_requests.get(self._selected_matrix_path)
            if self._selected_matrix_path is not None
            else None
        )
        return request.fit_range if request is not None else None

    def plot_selected_matrix(self) -> object | None:
        """Open an independent measured/fitted matrix window."""
        data = self.selected_matrix
        if data is None:
            self.status_label.setText("No matrix selected")
            return None
        fit = self.selected_fit
        fit_range = self.selected_fit_range if fit is not None else self.current_fit_request().fit_range
        controller = MatrixPlotController(
            data,
            fit=fit,
            fit_range=fit_range if fit is not None else None,
            title=f"DDM matrix: {self._selected_matrix_path.name}",
        )
        return self._show_plot(controller)

    def plot_selected_correlation(self) -> object | None:
        """Open an independent measured/fitted correlation window."""
        data = self.selected_matrix
        if data is None:
            self.status_label.setText("No matrix selected")
            return None
        fit = self.selected_fit
        fit_range = self.selected_fit_range if fit is not None else self.current_fit_request().fit_range
        controller = CorrelationPlotController(
            data,
            fit=fit,
            fit_range=fit_range,
            title=f"DDM correlation: {self._selected_matrix_path.name}",
        )
        return self._show_plot(controller)

    def plot_fitted_parameters(self) -> object | None:
        """Open an independent plot of all parameters from the selected fit."""
        fit = self.selected_fit
        if fit is None:
            self.status_label.setText("Fit the selected matrix first")
            return None
        controller = FitParameterPlotController(fit)
        return self._show_plot(controller)

    def plot_amplitude_noise_diffusion(self) -> object | None:
        """Open independent amplitude, noise, and diffusion plots."""
        fit = self.selected_fit
        if fit is None:
            self.status_label.setText("Fit the selected matrix first")
            return None
        controller = AmplitudeNoiseDiffusionPlotController(fit)
        return self._show_plot(controller)

    def _show_plot(self, controller: object) -> object:
        self._plot_controllers.append(controller)
        controller.show()
        return controller

    def _update_matrix_action_state(self) -> None:
        if self._job_active:
            for widget in (
                self.model_selector,
                self.q_min_slider,
                self.q_max_slider,
                self.time_min_slider,
                self.time_max_slider,
                self.initial_guess_button,
                self.fit_button,
                self.fitted_parameters_button,
                self.plot_matrix_button,
                self.plot_amplitude_button,
            ):
                widget.setEnabled(False)
            self.merge_button.setEnabled(False)
            self.average_button.setEnabled(False)
            for action in (
                self.open_directory_action,
                self.concatenate_action,
                self.split_action,
                self.show_matrix_action,
                self.plot_correlation_action,
                self.save_matrix_action,
                self.save_fit_action,
                self.fit_all_action,
                self.save_correlation_action,
                self.contin_action,
                self.about_action,
            ):
                action.setEnabled(False)
            return
        has_matrix = self.selected_matrix is not None
        has_catalog = bool(self._matrices)
        has_multiple = len(self._matrices) > 1
        has_fit = self.selected_fit is not None
        for widget in (
            self.model_selector,
            self.q_min_slider,
            self.q_max_slider,
            self.time_min_slider,
            self.time_max_slider,
            self.initial_guess_button,
            self.fit_button,
            self.plot_matrix_button,
        ):
            widget.setEnabled(has_matrix)
        self.fitted_parameters_button.setEnabled(has_fit)
        self.plot_amplitude_button.setEnabled(has_fit)
        self.merge_button.setEnabled(has_multiple)
        self.average_button.setEnabled(has_multiple)
        self.show_matrix_action.setEnabled(has_matrix)
        self.plot_correlation_action.setEnabled(has_matrix)
        self.save_matrix_action.setEnabled(has_matrix)
        self.save_fit_action.setEnabled(False)
        self.fit_all_action.setEnabled(has_catalog)
        self.save_correlation_action.setEnabled(False)
        self.contin_action.setEnabled(has_matrix)

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
            self.cancel_button,
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
        if self._job_active and self._thread is not None:
            self.cancel_processing()
            if not self._thread.wait(10_000):
                self.status_label.setText("Cancellation is still in progress")
                event.ignore()
                return
            self._job_active = False
            self._thread = None
            self._worker = None
        for controller in tuple(self._plot_controllers):
            controller.close()
        self._plot_controllers.clear()
        self._settings.setValue("geometry", self.saveGeometry())
        self._settings.setValue("last_directory", self.directory_edit.text())
        self._settings.sync()
        super().closeEvent(event)


def create_main_window(*, settings: QSettings | None = None) -> DDMMainWindow:
    """Create a shell window without starting the Qt event loop."""
    return DDMMainWindow(settings=settings)


__all__ = ["DDMMainWindow", "create_main_window"]
