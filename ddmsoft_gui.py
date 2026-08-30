"""Qt desktop interface for DDMSoft.

The numerical modules deliberately remain independent from this module.  This
keeps the processing code usable from notebooks and makes the GUI safe to
import in tests without opening a window.
"""

from __future__ import annotations

import os
import sys
import threading
import traceback
from pathlib import Path

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6.QtCore import QObject, Qt, QThread, Signal, Slot
from PySide6.QtGui import QColor, QFont, QIcon, QPalette
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QSplitter,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from contin import CONTIN
from fitDDM import fitOneDDMmatrix, getRadius, mergeDDM
from generateDDM import FFTStack, concatenateVideos, timeDependantDDM
from namedConcepts import FITMODELS, FITPARAMDEFAULTS, FITPARAMNAMES
from utilities import (
    loadAnalyzedVideos,
    loadDirectory,
    musthaves,
    renameTimestamp,
    saveAutocorrelationCSV,
    saveCONTINfit,
    saveFitTextFile,
    saveMatrixCSV,
    water_viscosity,
    extractCrudef,
)


MODE_NAMES = list(FITMODELS)
APP_DIR = Path(__file__).resolve().parent


class WorkSignals(QObject):
    progress = Signal(int, str)
    finished = Signal(object)
    failed = Signal(str)


class VideoWorker(QObject):
    """Process a directory without blocking Qt's event loop."""

    def __init__(self, jobs, max_couples, points_per_decade, angles, recompute, cancel_event):
        super().__init__()
        self.jobs = jobs
        self.max_couples = max_couples
        self.points_per_decade = points_per_decade
        self.angles = angles
        self.recompute = recompute
        self.cancel_event = cancel_event
        self.signals = WorkSignals()

    @Slot()
    def run(self):
        try:
            total_jobs = len(self.jobs)
            completed = 0
            for video, params in self.jobs.items():
                if self.cancel_event.is_set():
                    break
                try:
                    frequency = float(params["framerate"])
                    pixelsize = float(params["pixelsize"])
                except (KeyError, TypeError, ValueError) as error:
                    raise ValueError(
                        f"Invalid parameters for {Path(video).name}: "
                        "framerate and pixel size must be numbers"
                    ) from error
                if not np.isfinite(frequency) or frequency <= 0 or not np.isfinite(pixelsize) or pixelsize <= 0:
                    raise ValueError(
                        f"Invalid parameters for {Path(video).name}: "
                        "framerate and pixel size must be positive"
                    )

                matrix_dir = Path(video).parent / "ddm_matrices"
                matrix_stem = Path(video).stem
                matrix_exists = all(
                    (matrix_dir / f"{matrix_stem}{suffix}").exists()
                    for suffix in musthaves
                )
                if not self.recompute and matrix_exists:
                    completed += 1
                    self.signals.progress.emit(
                        int(completed / total_jobs * 100),
                        f"Skipped {Path(video).name}: matrix already exists",
                    )
                    continue

                self.signals.progress.emit(
                    int(completed / total_jobs * 100),
                    f"Loading {Path(video).name}",
                )

                def report(progress, maximum):
                    fraction = progress / maximum if maximum else 0
                    overall = (completed + min(fraction, 1)) / total_jobs
                    self.signals.progress.emit(
                        int(overall * 100),
                        f"Processing {Path(video).name}",
                    )

                stack = FFTStack(
                    frequency,
                    pixelsize,
                    self.max_couples,
                    self.points_per_decade,
                    Nangle=self.angles,
                    progress_callback=report,
                )
                stack.loadVideo(video)
                stack.fftVideo()
                if self.cancel_event.is_set():
                    break
                stack.stackToDDM()
                completed += 1
                self.signals.progress.emit(
                    int(completed / total_jobs * 100),
                    f"Finished {Path(video).name}",
                )
            self.signals.finished.emit(not self.cancel_event.is_set())
        except Exception:
            self.signals.failed.emit(traceback.format_exc())


class TimeDependentWorker(QObject):
    def __init__(self, video, params, max_couples, points_per_decade, partitions, cancel_event):
        super().__init__()
        self.video = video
        self.params = params
        self.max_couples = max_couples
        self.points_per_decade = points_per_decade
        self.partitions = partitions
        self.cancel_event = cancel_event
        self.signals = WorkSignals()

    @Slot()
    def run(self):
        try:
            self.signals.progress.emit(5, "Loading video partitions")
            analysis = timeDependantDDM(
                float(self.params["framerate"]),
                float(self.params["pixelsize"]),
                self.max_couples,
                self.points_per_decade,
                self.partitions,
            )
            analysis.loadVideo(self.video)
            if self.cancel_event.is_set():
                self.signals.finished.emit(False)
                return
            self.signals.progress.emit(35, "Computing Fourier transforms")
            analysis.fftAllStacks()
            if self.cancel_event.is_set():
                self.signals.finished.emit(False)
                return
            self.signals.progress.emit(60, "Averaging DDM matrices")
            analysis.ddmAllStacks()
            self.signals.finished.emit(True)
        except Exception:
            self.signals.failed.emit(traceback.format_exc())


class FitWorker(QObject):
    def __init__(self, jobs, model, initial, fixed, qmin, qmax, dtmin, dtmax, cancel_event):
        super().__init__()
        self.jobs = jobs
        self.model = model
        self.initial = initial
        self.fixed = fixed
        self.qmin = qmin
        self.qmax = qmax
        self.dtmin = dtmin
        self.dtmax = dtmax
        self.cancel_event = cancel_event
        self.signals = WorkSignals()

    @Slot()
    def run(self):
        try:
            results = {}
            for index, (path, data) in enumerate(self.jobs.items(), start=1):
                if self.cancel_event.is_set():
                    break
                _, dts, qs = data
                qmin = min(max(self.qmin, 0), max(len(qs) - 2, 0))
                qmax = min(max(self.qmax, qmin + 1), len(qs))
                dtmin = min(max(self.dtmin, 0), max(len(dts) - 2, 0))
                dtmax = min(max(self.dtmax, dtmin + 1), len(dts))
                self.signals.progress.emit(
                    int((index - 1) / len(self.jobs) * 100),
                    f"Fitting {Path(path).name}",
                )
                result = fitOneDDMmatrix(
                    data,
                    model=self.model,
                    ini=list(self.initial),
                    fixed=list(self.fixed),
                    qmin=qmin,
                    qmax=qmax,
                    dtmin=dtmin,
                    dtmax=dtmax,
                )
                results[path] = (result, qmin, qmax, dtmin, dtmax)
            self.signals.finished.emit({"results": results, "model": self.model})
        except Exception:
            self.signals.failed.emit(traceback.format_exc())


class ContinWorker(QObject):
    def __init__(self, tau, ddmdata, gammas, alphas, maxiter, cancel_event):
        super().__init__()
        self.tau = tau
        self.ddmdata = ddmdata
        self.gammas = gammas
        self.alphas = alphas
        self.maxiter = maxiter
        self.cancel_event = cancel_event
        self.signals = WorkSignals()

    @Slot()
    def run(self):
        try:
            solution = None
            run = CONTIN(
                self.tau,
                self.ddmdata,
                self.gammas,
                alpha=self.alphas,
                maxiter=self.maxiter,
            )
            for index, item in enumerate(run):
                if self.cancel_event.is_set():
                    break
                if np.isscalar(item):
                    self.signals.progress.emit(
                        int(index / max(len(self.alphas), 1) * 100),
                        f"CONTIN regularization {index + 1}/{len(self.alphas)}",
                    )
                else:
                    solution = item
            self.signals.finished.emit(solution)
        except Exception:
            self.signals.failed.emit(traceback.format_exc())


class PlotCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None):
        self.figure = Figure(figsize=(8, 5), tight_layout=True)
        super().__init__(self.figure)
        self.setParent(parent)


class PlotDialog(QDialog):
    def __init__(self, title, plotter, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(980, 680)
        self.canvas = PlotCanvas(self)
        plotter(self.canvas.figure)
        toolbar = NavigationToolbar2QT(self.canvas, self)
        layout = QVBoxLayout(self)
        layout.addWidget(toolbar)
        layout.addWidget(self.canvas)


class FitParametersDialog(QDialog):
    def __init__(self, model, values, fixed, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Fit parameter guesses")
        self.setMinimumWidth(540)
        self.model = model
        self.fields = []
        self.checks = []
        layout = QVBoxLayout(self)
        intro = QLabel("Set starting values. A checked parameter remains fixed during fitting.")
        intro.setWordWrap(True)
        layout.addWidget(intro)
        form = QFormLayout()
        for index, (name, value, is_fixed) in enumerate(
            zip(FITPARAMNAMES[model], values, fixed)
        ):
            row = QWidget()
            row_layout = QHBoxLayout(row)
            row_layout.setContentsMargins(0, 0, 0, 0)
            field = QLineEdit("" if value == "" else str(value))
            field.setObjectName(f"parameter_{index}")
            check = QCheckBox("fixed")
            check.setChecked(is_fixed)
            row_layout.addWidget(field, 1)
            row_layout.addWidget(check)
            form.addRow(name, row)
            self.fields.append(field)
            self.checks.append(check)
        layout.addLayout(form)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def values(self):
        result = []
        for index, field in enumerate(self.fields):
            text = field.text().strip()
            if text == "" and index >= len(self.fields) - 2:
                result.append("")
            else:
                result.append(float(text))
        return result, [check.isChecked() for check in self.checks]


class MatrixPickerDialog(QDialog):
    def __init__(self, names, title, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(600, 420)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Select one or more matrices:"))
        self.list = QListWidget()
        self.list.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)
        for path, name in names:
            item = QListWidgetItem(name)
            item.setData(Qt.ItemDataRole.UserRole, path)
            self.list.addItem(item)
        layout.addWidget(self.list)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def selected(self):
        return [item.data(Qt.ItemDataRole.UserRole) for item in self.list.selectedItems()]


class ContinDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("CONTIN relaxation-rate distribution")
        self.resize(1020, 700)
        self.solution = None
        self.canvas = PlotCanvas(self)
        layout = QVBoxLayout(self)
        controls = QGridLayout()
        self.q_index = QSpinBox()
        self.gamma_min = QLineEdit("1e-13")
        self.gamma_max = QLineEdit("1e-11")
        self.gamma_count = QSpinBox()
        self.gamma_count.setValue(10)
        self.alpha_min = QLineEdit("0.01")
        self.alpha_max = QLineEdit("1")
        self.alpha_count = QSpinBox()
        self.alpha_count.setValue(3)
        self.max_iter = QSpinBox()
        self.max_iter.setRange(1, 1000)
        self.max_iter.setValue(10)
        fields = [
            ("q index", self.q_index),
            ("smallest decay rate", self.gamma_min),
            ("biggest decay rate", self.gamma_max),
            ("decay rates", self.gamma_count),
            ("smallest alpha", self.alpha_min),
            ("biggest alpha", self.alpha_max),
            ("alpha values", self.alpha_count),
            ("iterations", self.max_iter),
        ]
        for index, (label, field) in enumerate(fields):
            controls.addWidget(QLabel(label), index // 4, (index % 4) * 2)
            controls.addWidget(field, index // 4, (index % 4) * 2 + 1)
        layout.addLayout(controls)
        buttons = QHBoxLayout()
        self.estimate = QPushButton("Estimate")
        self.plot = QPushButton("Plot selected")
        self.save = QPushButton("Save")
        self.close_button = QPushButton("Close")
        for button in (self.estimate, self.plot, self.save, self.close_button):
            buttons.addWidget(button)
        layout.addLayout(buttons)
        layout.addWidget(self.canvas)
        self.close_button.clicked.connect(self.close)

    def settings(self):
        return {
            "q_index": self.q_index.value(),
            "gamma_min": float(self.gamma_min.text()),
            "gamma_max": float(self.gamma_max.text()),
            "gamma_count": self.gamma_count.value(),
            "alpha_min": float(self.alpha_min.text()),
            "alpha_max": float(self.alpha_max.text()),
            "alpha_count": self.alpha_count.value(),
            "max_iter": self.max_iter.value(),
        }


class DDMWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("DDMSoft | Differential Dynamic Microscopy")
        self.setMinimumSize(1180, 760)
        self.resize(1420, 900)
        icon_path = APP_DIR / "logo" / "logo.png"
        if icon_path.exists():
            self.setWindowIcon(QIcon(str(icon_path)))

        self.params = {}
        self.computed_data = {}
        self.display_data = {}
        self.fitted_data = {}
        self.contin_solutions = {}
        self.current_matrix = None
        self.current_path = None
        self.fit_guesses = {
            model: (list(values), [False] * len(values))
            for model, values in FITPARAMDEFAULTS.items()
        }
        self.worker_thread = None
        self.cancel_event = None
        self.contin_dialog = None

        self._build_ui()
        self._apply_style()

    def _build_ui(self):
        root = QWidget()
        root_layout = QVBoxLayout(root)
        root_layout.setContentsMargins(24, 20, 24, 18)
        root_layout.setSpacing(14)

        header = QHBoxLayout()
        title_block = QVBoxLayout()
        title = QLabel("DDMSoft")
        title.setObjectName("appTitle")
        subtitle = QLabel("Real-time differential dynamic microscopy")
        subtitle.setObjectName("subtitle")
        title_block.addWidget(title)
        title_block.addWidget(subtitle)
        header.addLayout(title_block)
        header.addStretch()
        self.directory = QLineEdit(str(Path.home()))
        self.directory.setPlaceholderText("Choose a directory containing .avi videos")
        self.directory.setMinimumWidth(420)
        browse = QPushButton("Browse")
        load = QPushButton("Load directory")
        browse.clicked.connect(self.browse_directory)
        load.clicked.connect(self.load_directory)
        header.addWidget(self.directory)
        header.addWidget(browse)
        header.addWidget(load)
        root_layout.addLayout(header)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setChildrenCollapsible(False)
        splitter.addWidget(self._build_sidebar())
        splitter.addWidget(self._build_workspace())
        splitter.setSizes([390, 950])
        root_layout.addWidget(splitter, 1)
        self.setCentralWidget(root)

        status = self.statusBar()
        self.status_label = QLabel("Ready. Load a video directory to begin.")
        status.addWidget(self.status_label, 1)
        self.global_progress = QProgressBar()
        self.global_progress.setFixedWidth(220)
        self.global_progress.setRange(0, 100)
        self.global_progress.setValue(0)
        status.addPermanentWidget(self.global_progress)

    def _build_sidebar(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(4, 2, 10, 2)

        process_box = QGroupBox("1  Process videos")
        process = QFormLayout(process_box)
        self.keep_existing = QCheckBox("Keep existing matrices")
        self.keep_existing.setChecked(True)
        self.max_couples = QSpinBox()
        self.max_couples.setRange(0, 1_000_000)
        self.max_couples.setValue(300)
        self.points_per_decade = QSpinBox()
        self.points_per_decade.setRange(1, 200)
        self.points_per_decade.setValue(20)
        self.angle_count = QSpinBox()
        self.angle_count.setRange(1, 360)
        self.angle_count.setValue(1)
        process.addRow(self.keep_existing)
        process.addRow("Max frame couples", self.max_couples)
        process.addRow("Lag times per decade", self.points_per_decade)
        process.addRow("Angular sectors", self.angle_count)
        buttons = QHBoxLayout()
        self.process_button = QPushButton("Process videos")
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.setEnabled(False)
        self.process_button.clicked.connect(self.process_videos)
        self.cancel_button.clicked.connect(self.cancel_work)
        buttons.addWidget(self.process_button)
        buttons.addWidget(self.cancel_button)
        process.addRow(buttons)
        layout.addWidget(process_box)

        matrix_box = QGroupBox("2  Matrix to analyse")
        matrix_form = QFormLayout(matrix_box)
        self.matrix_combo = QComboBox()
        self.matrix_combo.currentIndexChanged.connect(self.select_matrix)
        refresh = QPushButton("Refresh matrices")
        refresh.clicked.connect(self.refresh_matrices)
        matrix_form.addRow(self.matrix_combo)
        matrix_form.addRow(refresh)
        self.q_min = self._range_spin()
        self.q_max = self._range_spin()
        self.t_min = self._range_spin()
        self.t_max = self._range_spin()
        self.q_min.valueChanged.connect(self.update_preview)
        self.q_max.valueChanged.connect(self.update_preview)
        self.t_min.valueChanged.connect(self.update_preview)
        self.t_max.valueChanged.connect(self.update_preview)
        matrix_form.addRow("q range", self._paired_spins(self.q_min, self.q_max))
        matrix_form.addRow("time range", self._paired_spins(self.t_min, self.t_max))
        layout.addWidget(matrix_box)

        fit_box = QGroupBox("3  Fit relaxation models")
        fit_form = QFormLayout(fit_box)
        self.fit_model = QComboBox()
        self.fit_model.addItems(MODE_NAMES)
        self.fit_model.currentTextChanged.connect(self.fit_model_changed)
        fit_form.addRow("Model", self.fit_model)
        guess_button = QPushButton("Edit initial guesses")
        guess_button.clicked.connect(self.edit_fit_guesses)
        fit_form.addRow(guess_button)
        self.fit_button = QPushButton("Fit selected matrix")
        fit_all_button = QPushButton("Fit all matrices")
        self.fit_button.clicked.connect(self.fit_selected)
        fit_all_button.clicked.connect(self.fit_all)
        fit_form.addRow(self.fit_button, fit_all_button)
        self.temperature = QLineEdit()
        self.temperature.setPlaceholderText("optional, deg C")
        self.viscosity = QLineEdit("water")
        self.viscosity.setPlaceholderText("water or Pa s")
        fit_form.addRow("Temperature", self.temperature)
        fit_form.addRow("Viscosity", self.viscosity)
        contin_button = QPushButton("CONTIN distribution...")
        contin_button.clicked.connect(self.open_contin)
        fit_form.addRow(contin_button)
        layout.addWidget(fit_box)

        plot_box = QGroupBox("4  Inspect and export")
        plot_grid = QGridLayout(plot_box)
        plot_actions = [
            ("Inspect curve", self.show_inspector),
            ("Matrix image", self.show_matrix_image),
            ("Amplitude / noise", self.show_amplitude),
            ("Fit parameters", self.show_fit_parameters),
            ("Correlation curves", self.show_correlations),
            ("Save matrix", self.save_matrix),
            ("Save fit", self.save_fit),
            ("Save correlation", self.save_correlation),
        ]
        for index, (label, callback) in enumerate(plot_actions):
            button = QPushButton(label)
            button.clicked.connect(callback)
            plot_grid.addWidget(button, index // 2, index % 2)
        layout.addWidget(plot_box)

        tools_box = QGroupBox("Tools")
        tools = QGridLayout(tools_box)
        merge = QPushButton("Merge matrices...")
        average = QPushButton("Average matrices...")
        split = QPushButton("Time-dependent DDM...")
        concatenate = QPushButton("Concatenate videos...")
        rename = QPushButton("Rename timestamps...")
        merge.clicked.connect(lambda: self.merge_matrices("merge"))
        average.clicked.connect(lambda: self.merge_matrices("average"))
        split.clicked.connect(self.time_dependent)
        concatenate.clicked.connect(self.concatenate_videos)
        rename.clicked.connect(self.rename_timestamps)
        for index, button in enumerate((merge, average, split, concatenate, rename)):
            tools.addWidget(button, index // 2, index % 2)
        layout.addWidget(tools_box)
        layout.addStretch()
        scroll.setWidget(panel)
        return scroll

    @staticmethod
    def _range_spin():
        spin = QSpinBox()
        spin.setRange(0, 0)
        return spin

    @staticmethod
    def _paired_spins(first, second):
        widget = QWidget()
        layout = QHBoxLayout(widget)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(first)
        layout.addWidget(QLabel("to"))
        layout.addWidget(second)
        return widget

    def _build_workspace(self):
        tabs = QTabWidget()

        overview = QWidget()
        overview_layout = QVBoxLayout(overview)
        self.video_table = QTableWidget(0, 3)
        self.video_table.setHorizontalHeaderLabels(["Video", "Frame rate (Hz)", "Pixel size (m)"])
        self.video_table.horizontalHeader().setStretchLastSection(True)
        self.video_table.setAlternatingRowColors(True)
        self.video_table.itemChanged.connect(self.video_parameters_changed)
        overview_layout.addWidget(self.video_table, 1)
        help_text = QLabel(
            "Edit frame rate and pixel size directly in the table when acquisition metadata "
            "is incomplete. Processed matrices appear in the Matrix tab."
        )
        help_text.setWordWrap(True)
        overview_layout.addWidget(help_text)
        tabs.addTab(overview, "Video catalog")

        matrix_tab = QWidget()
        matrix_layout = QVBoxLayout(matrix_tab)
        self.matrix_info = QLabel("No matrix selected")
        self.matrix_info.setObjectName("cardLabel")
        matrix_layout.addWidget(self.matrix_info)
        self.preview_canvas = PlotCanvas()
        matrix_layout.addWidget(self.preview_canvas, 1)
        tabs.addTab(matrix_tab, "Matrix inspector")

        fit_tab = QWidget()
        fit_layout = QVBoxLayout(fit_tab)
        self.fit_info = QLabel("Fit results will appear here")
        self.fit_info.setObjectName("cardLabel")
        fit_layout.addWidget(self.fit_info)
        self.fit_canvas = PlotCanvas()
        fit_layout.addWidget(self.fit_canvas, 1)
        tabs.addTab(fit_tab, "Fit results")
        self.tabs = tabs
        return tabs

    def _apply_style(self):
        app = QApplication.instance()
        app.setStyle("Fusion")
        palette = QPalette()
        palette.setColor(QPalette.ColorRole.Window, QColor("#111827"))
        palette.setColor(QPalette.ColorRole.WindowText, QColor("#e5e7eb"))
        palette.setColor(QPalette.ColorRole.Base, QColor("#0b1220"))
        palette.setColor(QPalette.ColorRole.AlternateBase, QColor("#162235"))
        palette.setColor(QPalette.ColorRole.Text, QColor("#e5e7eb"))
        palette.setColor(QPalette.ColorRole.Button, QColor("#22324a"))
        palette.setColor(QPalette.ColorRole.ButtonText, QColor("#f8fafc"))
        palette.setColor(QPalette.ColorRole.Highlight, QColor("#4f8cff"))
        palette.setColor(QPalette.ColorRole.HighlightedText, QColor("#ffffff"))
        app.setPalette(palette)
        app.setStyleSheet(
            """
            QMainWindow { background: #111827; }
            QLabel#appTitle { font-size: 28px; font-weight: 700; color: #f8fafc; }
            QLabel#subtitle { font-size: 13px; color: #93a4bd; }
            QGroupBox { border: 1px solid #2b3a52; border-radius: 8px; margin-top: 12px; padding: 12px 8px 8px 8px; font-weight: 600; color: #b9c8dc; }
            QGroupBox::title { subcontrol-origin: margin; left: 12px; padding: 0 5px; }
            QPushButton { min-height: 32px; padding: 0 10px; border: 1px solid #3c5270; border-radius: 5px; background: #22324a; }
            QPushButton:hover { background: #2e4667; }
            QPushButton:disabled { color: #718096; background: #172235; }
            QLineEdit, QSpinBox, QComboBox { min-height: 30px; border: 1px solid #3c5270; border-radius: 4px; padding: 0 7px; background: #0b1220; }
            QTableWidget, QTextEdit, QListWidget { border: 1px solid #2b3a52; background: #0b1220; alternate-background-color: #162235; gridline-color: #2b3a52; }
            QTabWidget::pane { border: 1px solid #2b3a52; border-radius: 6px; }
            QTabBar::tab { padding: 10px 18px; color: #93a4bd; }
            QTabBar::tab:selected { color: #f8fafc; border-bottom: 2px solid #4f8cff; }
            QLabel#cardLabel { padding: 10px; border-radius: 6px; background: #162235; color: #b9c8dc; }
            QProgressBar { border: 1px solid #3c5270; border-radius: 4px; text-align: center; }
            QProgressBar::chunk { background: #4f8cff; }
            """
        )
        font = QFont("Sans Serif", 10)
        app.setFont(font)

    def browse_directory(self):
        path = QFileDialog.getExistingDirectory(self, "Choose video directory", self.directory.text())
        if path:
            self.directory.setText(path)

    def load_directory(self):
        path = self.directory.text().strip()
        if not path or not Path(path).is_dir():
            self.show_error("Directory not found", "Choose an existing directory containing AVI videos.")
            return
        try:
            self.params = loadDirectory(path)
            self._populate_video_table()
            self.refresh_matrices()
            self.set_status(f"Loaded {len(self.params)} video(s) from {path}")
        except Exception as error:
            self.show_error("Could not load directory", str(error))

    def _populate_video_table(self):
        self.video_table.blockSignals(True)
        self.video_table.setRowCount(0)
        for row, (path, params) in enumerate(sorted(self.params.items())):
            self.video_table.insertRow(row)
            name = QTableWidgetItem(Path(path).name)
            name.setData(Qt.ItemDataRole.UserRole, path)
            name.setFlags(name.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.video_table.setItem(row, 0, name)
            self.video_table.setItem(row, 1, QTableWidgetItem(str(params.get("framerate", ""))))
            self.video_table.setItem(row, 2, QTableWidgetItem(str(params.get("pixelsize", ""))))
        self.video_table.resizeColumnsToContents()
        self.video_table.blockSignals(False)

    def video_parameters_changed(self, item):
        row = item.row()
        path_item = self.video_table.item(row, 0)
        if path_item is None:
            return
        path = path_item.data(Qt.ItemDataRole.UserRole)
        if path in self.params:
            self.params[path]["framerate"] = self.video_table.item(row, 1).text().strip()
            self.params[path]["pixelsize"] = self.video_table.item(row, 2).text().strip()

    def refresh_matrices(self):
        path = self.directory.text().strip()
        if not Path(path).is_dir():
            return
        try:
            self.computed_data, self.display_data = loadAnalyzedVideos(path)
        except Exception as error:
            self.show_error("Could not load matrices", str(error))
            return
        current = self.current_path
        self.matrix_combo.blockSignals(True)
        self.matrix_combo.clear()
        for display, matrix_path in self.display_data.items():
            self.matrix_combo.addItem(display, matrix_path)
        self.matrix_combo.blockSignals(False)
        if self.matrix_combo.count():
            index = self.matrix_combo.findData(current)
            self.matrix_combo.setCurrentIndex(max(index, 0))
            self.select_matrix(self.matrix_combo.currentIndex())
        else:
            self.current_path = None
            self.current_matrix = None
            self.matrix_info.setText("No computed matrix in this directory")
            self.preview_canvas.figure.clear()
            self.preview_canvas.draw_idle()

    def select_matrix(self, index):
        if index < 0:
            return
        path = self.matrix_combo.itemData(index)
        if path not in self.computed_data:
            return
        self.current_path = path
        self.current_matrix = self.computed_data[path]
        ddm, dts, qs = self.current_matrix
        self.q_min.setRange(0, max(len(qs) - 1, 0))
        self.q_max.setRange(0, max(len(qs) - 1, 0))
        self.t_min.setRange(0, max(len(dts) - 1, 0))
        self.t_max.setRange(0, max(len(dts) - 1, 0))
        self.q_min.setValue(0)
        self.q_max.setValue(max(len(qs) - 1, 0))
        self.t_min.setValue(0)
        self.t_max.setValue(max(len(dts) - 1, 0))
        self.matrix_info.setText(
            f"{Path(path).name}   |   matrix {ddm.shape[0]} x {ddm.shape[1]}   |   "
            f"q {qs.min():.3g} to {qs.max():.3g} m^-1"
        )
        self.update_preview()

    def update_preview(self):
        if self.current_matrix is None:
            return
        self._plot_inspector(self.preview_canvas.figure, self.current_matrix, self.current_path)
        self.preview_canvas.draw_idle()

    def _valid_ranges(self):
        if self.current_matrix is None:
            raise ValueError("Select a computed matrix first")
        qmin, qmax = self.q_min.value(), self.q_max.value()
        tmin, tmax = self.t_min.value(), self.t_max.value()
        if qmax <= qmin or tmax <= tmin:
            raise ValueError("The upper q and time limits must be greater than the lower limits")
        return qmin, qmax + 1, tmin, tmax + 1

    def _jobs_from_table(self):
        for row in range(self.video_table.rowCount()):
            self.video_parameters_changed(self.video_table.item(row, 1))
        return dict(self.params)

    def process_videos(self):
        jobs = self._jobs_from_table()
        if not jobs:
            self.show_error("No videos", "Load a directory containing AVI videos first.")
            return
        if self.angle_count.value() > 1:
            message = "Angular sectors reduce statistics per direction. Continue?"
            if QMessageBox.question(self, "Directional DDM", message) != QMessageBox.StandardButton.Yes:
                return
        self.cancel_event = threading.Event()
        worker = VideoWorker(
            jobs,
            self.max_couples.value(),
            self.points_per_decade.value(),
            self.angle_count.value(),
            not self.keep_existing.isChecked(),
            self.cancel_event,
        )
        self.start_worker(worker, self.processing_finished)

    def processing_finished(self, completed):
        self.refresh_matrices()
        self.set_status("Processing complete" if completed else "Processing cancelled")

    def time_dependent(self):
        jobs = self._jobs_from_table()
        if len(jobs) != 1:
            self.show_error("One video required", "Load a directory containing exactly one video.")
            return
        partitions, accepted = self.get_int("Time-dependent DDM", "Number of portions:", 10, 1, 1000)
        if not accepted:
            return
        video, params = next(iter(jobs.items()))
        self.cancel_event = threading.Event()
        worker = TimeDependentWorker(
            video,
            params,
            self.max_couples.value(),
            self.points_per_decade.value(),
            partitions,
            self.cancel_event,
        )
        self.start_worker(worker, lambda _: self.refresh_matrices())

    def edit_fit_guesses(self):
        model = FITMODELS[self.fit_model.currentText()]
        values, fixed = self.fit_guesses[model]
        dialog = FitParametersDialog(model, values, fixed, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            try:
                self.fit_guesses[model] = dialog.values()
            except ValueError as error:
                self.show_error("Invalid fit parameter", str(error))

    def fit_model_changed(self):
        self.fit_info.setText(f"Selected model: {self.fit_model.currentText()}")

    def fit_selected(self):
        if self.current_path is None:
            self.show_error("No matrix", "Select a computed matrix first.")
            return
        try:
            ranges = self._valid_ranges()
        except ValueError as error:
            self.show_error("Invalid fit range", str(error))
            return
        model = FITMODELS[self.fit_model.currentText()]
        initial, fixed = self.fit_guesses[model]
        self.cancel_event = threading.Event()
        worker = FitWorker(
            {self.current_path: self.current_matrix}, model, initial, fixed, *ranges, self.cancel_event
        )
        self.start_worker(worker, self.fit_finished)

    def fit_all(self):
        if not self.computed_data:
            self.show_error("No matrices", "Process or load matrices first.")
            return
        try:
            ranges = self._valid_ranges()
        except ValueError as error:
            self.show_error("Invalid fit range", str(error))
            return
        model = FITMODELS[self.fit_model.currentText()]
        initial, fixed = self.fit_guesses[model]
        self.cancel_event = threading.Event()
        worker = FitWorker(self.computed_data, model, initial, fixed, *ranges, self.cancel_event)
        self.start_worker(worker, self.fit_finished)

    def fit_finished(self, payload):
        results = payload["results"]
        model = payload["model"]
        for path, (result, qmin, qmax, _, _) in results.items():
            amplitude, noise, model_params, fit, _ = result
            dts = self.computed_data[path][1]
            qs = self.computed_data[path][2]
            ddm_fit = amplitude * (1 - fit) + noise
            self.fitted_data[path] = {
                "dts": dts,
                "qs": qs[qmin:qmax],
                "ddm": ddm_fit,
                "A": amplitude,
                "B": noise,
                "params": model_params,
                "f": fit,
                "model": model,
            }
        self._plot_fit_tab()
        self.set_status(f"Fit complete for {len(results)} matrix/matrices")

    def _plot_inspector(self, figure, data, path=None):
        figure.clear()
        ddm, dts, qs = data
        axes = figure.subplots(1, 2)
        qindex = min(max(self.q_min.value(), 0), len(qs) - 1)
        try:
            raw_f = extractCrudef(ddm)
        except (IndexError, ZeroDivisionError):
            raw_f = np.zeros_like(ddm)
        axes[0].semilogx(dts * qs[qindex] ** 2, raw_f[:, qindex], "o", ms=3, color="#4f8cff")
        axes[1].semilogx(dts, ddm[:, qindex], "o", ms=3, color="#f59e0b")
        fit = self.fitted_data.get(path)
        if fit is not None and len(fit["qs"]):
            fit_index = int(np.argmin(np.abs(fit["qs"] - qs[qindex])))
            if fit_index < fit["f"].shape[1]:
                axes[0].plot(fit["dts"] * fit["qs"][fit_index] ** 2, fit["f"][:, fit_index], color="#f472b6")
                axes[1].plot(fit["dts"], fit["ddm"][:, fit_index], color="#f472b6")
        axes[0].set_title(f"Correlation, q = {qs[qindex] / 1e6:.2f} um^-1")
        axes[0].set_xlabel("tau q^2 [s/m^2]")
        axes[0].set_ylabel("f(q, tau)")
        axes[1].set_title("DDM signal")
        axes[1].set_xlabel("lag time [s]")
        axes[1].set_ylabel("D(q, tau)")
        for axis in axes:
            axis.grid(alpha=0.2)
        figure.tight_layout()

    def _plot_fit_tab(self):
        self._plot_inspector(self.fit_canvas.figure, self.current_matrix, self.current_path)
        self.fit_canvas.draw_idle()
        fit = self.fitted_data.get(self.current_path)
        if fit:
            self.fit_info.setText(
                f"{self.fit_model.currentText()}   |   {len(fit['qs'])} q values fitted"
            )

    def show_inspector(self):
        if self.current_matrix is None:
            self.show_error("No matrix", "Select a computed matrix first.")
            return
        dialog = PlotDialog(
            f"DDM inspector | {Path(self.current_path).name}",
            lambda figure: self._plot_inspector(figure, self.current_matrix, self.current_path),
            self,
        )
        dialog.show()
        self._retain_dialog(dialog)

    def show_matrix_image(self):
        if self.current_matrix is None:
            self.show_error("No matrix", "Select a computed matrix first.")
            return
        def plot(figure):
            figure.clear()
            axes = figure.subplots(1, 2 if self.current_path in self.fitted_data else 1)
            axes = np.atleast_1d(axes)
            ddm, dts, qs = self.current_matrix
            extent = [qs[0] / 1e6, qs[-1] / 1e6, dts[0], dts[-1]]
            image = axes[0].imshow(ddm, aspect="auto", origin="lower", extent=extent, cmap="plasma")
            axes[0].set_title("Measured DDM matrix")
            axes[0].set_xlabel("q [um^-1]")
            axes[0].set_ylabel("lag time [s]")
            figure.colorbar(image, ax=axes[0], shrink=0.85)
            fit = self.fitted_data.get(self.current_path)
            if fit:
                image = axes[1].imshow(
                    fit["ddm"], aspect="auto", origin="lower",
                    extent=[fit["qs"][0] / 1e6, fit["qs"][-1] / 1e6, fit["dts"][0], fit["dts"][-1]],
                    cmap="plasma",
                )
                axes[1].set_title("Fitted region")
                axes[1].set_xlabel("q [um^-1]")
                figure.colorbar(image, ax=axes[1], shrink=0.85)
            figure.tight_layout()
        dialog = PlotDialog("DDM matrix image", plot, self)
        dialog.show()
        self._retain_dialog(dialog)

    def show_amplitude(self):
        fit = self._current_fit()
        if fit is None:
            return
        def plot(figure):
            figure.clear()
            axes = figure.subplots(2, 1, sharex=True)
            qs = fit["qs"] / 1e6
            axes[0].plot(qs, fit["A"], ".", label="amplitude A")
            axes[0].plot(qs, fit["B"], ".", label="noise B")
            axes[0].legend()
            axes[0].set_ylabel("intensity")
            diffusion = fit["params"][0]
            axes[1].plot(qs, diffusion, ".", color="#f472b6")
            axes[1].set_xlabel("q [um^-1]")
            axes[1].set_ylabel("diffusion [m^2/s]")
            figure.tight_layout()
        dialog = PlotDialog("Amplitude, noise and diffusion", plot, self)
        dialog.show()
        self._retain_dialog(dialog)

    def show_fit_parameters(self):
        fit = self._current_fit()
        if fit is None:
            return
        names = FITPARAMNAMES[fit["model"]]
        def plot(figure):
            figure.clear()
            axes = np.atleast_1d(figure.subplots(len(fit["params"]), 1, sharex=True))
            for axis, name, values in zip(axes, names, fit["params"]):
                axis.plot(fit["qs"] / 1e6, values, ".")
                axis.set_ylabel(name)
            axes[-1].set_xlabel("q [um^-1]")
            figure.tight_layout()
        dialog = PlotDialog("Fitted parameters", plot, self)
        dialog.show()
        self._retain_dialog(dialog)

    def show_correlations(self):
        if self.current_matrix is None:
            self.show_error("No matrix", "Select a computed matrix first.")
            return
        def plot(figure):
            figure.clear()
            axis = figure.subplots()
            ddm, dts, qs = self.current_matrix
            raw = extractCrudef(ddm)
            first, last = self.q_min.value(), self.q_max.value() + 1
            indexes = np.linspace(first, max(first, last - 1), min(5, last - first), dtype=int)
            for index in np.unique(indexes):
                axis.semilogx(dts * qs[index] ** 2, raw[:, index], "o", ms=3, label=f"{qs[index] / 1e6:.2f}")
            fit = self.fitted_data.get(self.current_path)
            if fit:
                for index in range(min(5, fit["f"].shape[1])):
                    axis.semilogx(fit["dts"] * fit["qs"][index] ** 2, fit["f"][:, index])
            axis.set_xlabel("tau q^2 [s/m^2]")
            axis.set_ylabel("f(q, tau)")
            axis.legend(title="q [um^-1]", frameon=False)
            figure.tight_layout()
        dialog = PlotDialog("Correlation functions", plot, self)
        dialog.show()
        self._retain_dialog(dialog)

    def _current_fit(self):
        if self.current_path is None or self.current_path not in self.fitted_data:
            self.show_error("No fit", "Fit the selected matrix before opening this view.")
            return None
        return self.fitted_data[self.current_path]

    def save_matrix(self):
        if self.current_matrix is None:
            self.show_error("No matrix", "Select a computed matrix first.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save matrix data", "", "Text files (*.txt *.csv)")
        if path:
            saveMatrixCSV(path, self.current_matrix)
            self.set_status(f"Saved matrix data to {path}")

    def _temperature_viscosity(self):
        temperature = self.temperature.text().strip()
        if not temperature:
            return None, None
        temp_kelvin = float(temperature) + 273.15
        viscosity = self.viscosity.text().strip()
        viscosity_value = water_viscosity(temp_kelvin) if viscosity.lower() == "water" else float(viscosity)
        if not np.isfinite(temp_kelvin) or temp_kelvin <= 0 or not np.isfinite(viscosity_value) or viscosity_value <= 0:
            raise ValueError("Temperature and viscosity must be finite positive values")
        return temp_kelvin, viscosity_value

    def save_fit(self):
        fit = self._current_fit()
        if fit is None:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save fit parameters", "", "Text files (*.txt *.csv)")
        if path:
            temperature, viscosity = self._temperature_viscosity()
            saveFitTextFile(
                path,
                fit["qs"],
                fit["A"],
                fit["B"],
                fit["params"],
                FITPARAMNAMES[fit["model"]][:-2],
                temperature=temperature,
                viscosity=viscosity,
            )
            self.set_status(f"Saved fit parameters to {path}")

    def save_correlation(self):
        fit = self._current_fit()
        if fit is None:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save autocorrelation", "", "CSV files (*.csv)")
        if path:
            saveAutocorrelationCSV(path, self.current_matrix, fit["A"], fit["B"], fit["qs"])
            self.set_status(f"Saved autocorrelation to {path}")

    def merge_matrices(self, mode):
        if len(self.display_data) < 2:
            self.show_error("Not enough matrices", "At least two matrices are required.")
            return
        names = [(path, display) for display, path in self.display_data.items()]
        dialog = MatrixPickerDialog(names, f"{mode.title()} matrices", self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        paths = dialog.selected()
        if not paths:
            return
        try:
            selected = {path: self.computed_data[path] for path in paths}
            output = mergeDDM(selected, mode=mode)
            self.refresh_matrices()
            self.set_status(f"Created {mode} matrix: {output}")
        except Exception as error:
            self.show_error(f"Could not {mode} matrices", str(error))

    def rename_timestamps(self):
        path = QFileDialog.getExistingDirectory(self, "Choose directory", self.directory.text())
        if path:
            try:
                renameTimestamp(path)
                self.set_status("Timestamp names updated")
            except Exception as error:
                self.show_error("Could not rename timestamps", str(error))

    def concatenate_videos(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self,
            "Select videos in order",
            self.directory.text(),
            "AVI videos (*.avi *.AVI)",
        )
        if len(paths) < 2:
            return
        output, _ = QFileDialog.getSaveFileName(
            self, "Save concatenated video", self.directory.text(), "AVI videos (*.avi)"
        )
        if not output:
            return
        try:
            concatenateVideos(paths, output)
            self.set_status(f"Concatenated {len(paths)} videos into {output}")
        except Exception as error:
            self.show_error("Could not concatenate videos", str(error))

    def open_contin(self):
        if self.current_matrix is None:
            self.show_error("No matrix", "Select a computed matrix first.")
            return
        dialog = ContinDialog(self)
        self.contin_dialog = dialog
        _, dts, qs = self.current_matrix
        dialog.q_index.setRange(0, len(qs) - 1)
        dialog.q_index.setValue(self.q_min.value())
        dialog.matrix_path = self.current_path
        dialog.matrix_data = self.current_matrix
        dialog.estimate.clicked.connect(lambda: self.run_contin(dialog))
        dialog.plot.clicked.connect(lambda: self.plot_contin(dialog))
        dialog.save.clicked.connect(lambda: self.save_contin(dialog))
        dialog.show()

    def run_contin(self, dialog):
        try:
            settings = dialog.settings()
            path = dialog.matrix_path
            ddm, dts, qs = dialog.matrix_data
            q_index = settings["q_index"]
            tau = dts * qs[q_index] ** 2
            worker = ContinWorker(
                tau,
                ddm[:, q_index],
                np.linspace(settings["gamma_min"], settings["gamma_max"], settings["gamma_count"]),
                np.linspace(settings["alpha_min"], settings["alpha_max"], settings["alpha_count"]),
                settings["max_iter"],
                threading.Event(),
            )
            dialog.running_path = path
            dialog.running_q_index = q_index
            self.start_worker(worker, lambda solution: self.contin_finished(dialog, path, q_index, solution))
        except Exception as error:
            self.show_error("Invalid CONTIN settings", str(error))

    def contin_finished(self, dialog, path, q_index, solution):
        if solution is None:
            self.show_error("CONTIN failed", "No solution was returned.")
            return
        solution.q_index = q_index
        try:
            temperature, viscosity = self._temperature_viscosity()
            solution.sizes = getRadius(solution.gamma_range, viscosity, temperature) if temperature and viscosity else None
        except Exception:
            solution.sizes = None
        self.contin_solutions.setdefault(path, {})[dialog.matrix_data[2][q_index]] = solution
        dialog.solution = solution
        self.plot_contin(dialog)
        self.set_status("CONTIN estimate complete")

    def plot_contin(self, dialog):
        solution = dialog.solution
        if solution is None:
            self.show_error("No CONTIN result", "Run Estimate first.")
            return
        alpha_index = solution.alphas.index(solution.chosen_alpha)
        def plot(figure):
            figure.clear()
            axes = figure.subplots(1, 2)
            fit = solution.alpha_amplitude[alpha_index] * (1 - solution.alpha_ddmfit[alpha_index]) + solution.alpha_noise[alpha_index]
            axes[0].semilogx(solution.tau, solution.ddmdata, ".", label="data")
            axes[0].semilogx(solution.tau, fit, label="CONTIN")
            axes[0].set_xlabel("tau q^2 [s/m^2]")
            axes[0].legend(frameon=False)
            distribution = solution.alpha_g[alpha_index]
            distribution = distribution / np.nansum(distribution)
            x = solution.sizes * 1e9 if solution.sizes is not None else solution.gamma_range
            axes[1].plot(x, distribution)
            axes[1].set_xlabel("R_h [nm]" if solution.sizes is not None else "gamma [m^2/s]")
            axes[1].set_ylabel("intensity")
            figure.tight_layout()
        dialog.canvas.figure.clear()
        plot(dialog.canvas.figure)
        dialog.canvas.draw_idle()

    def save_contin(self, dialog):
        if dialog.solution is None:
            self.show_error("No CONTIN result", "Run Estimate first.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save CONTIN result", "", "CSV files (*.csv)")
        if path:
            q = dialog.matrix_data[2][dialog.solution.q_index]
            saveCONTINfit(path, self.contin_solutions, dialog.matrix_path, q)
            self.set_status(f"Saved CONTIN result to {path}")

    def start_worker(self, worker, on_finished):
        if self.worker_thread is not None:
            self.show_error("Busy", "Wait for the current operation to finish.")
            return
        thread = QThread(self)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.signals.progress.connect(self.worker_progress)
        worker.signals.finished.connect(on_finished)
        worker.signals.finished.connect(thread.quit)
        worker.signals.failed.connect(self.worker_failed)
        worker.signals.failed.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self.worker_finished)
        self.worker = worker
        self.cancel_event = getattr(worker, "cancel_event", None)
        self.worker_thread = thread
        self.process_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.fit_button.setEnabled(False)
        thread.start()

    def worker_progress(self, value, message):
        self.global_progress.setValue(value)
        self.set_status(message)

    def worker_failed(self, details):
        self.show_error("Operation failed", details)
        self.set_status("Operation failed")

    def worker_finished(self):
        self.worker_thread = None
        self.worker = None
        self.process_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        self.fit_button.setEnabled(True)
        self.global_progress.setValue(0)

    def cancel_work(self):
        if self.cancel_event is not None:
            self.cancel_event.set()
            self.set_status("Cancelling after the current video...")

    def _retain_dialog(self, dialog):
        if not hasattr(self, "_dialogs"):
            self._dialogs = []
        self._dialogs.append(dialog)
        dialog.finished.connect(lambda _: self._dialogs.remove(dialog) if dialog in self._dialogs else None)

    def get_int(self, title, label, value, minimum, maximum):
        from PySide6.QtWidgets import QInputDialog
        return QInputDialog.getInt(self, title, label, value, minimum, maximum)

    def set_status(self, message):
        self.status_label.setText(message)

    def show_error(self, title, message):
        QMessageBox.critical(self, title, message)

    def closeEvent(self, event):
        if self.worker_thread is not None:
            if self.cancel_event is not None:
                self.cancel_event.set()
            if not self.worker_thread.wait(3000):
                self.set_status("Please wait for the current operation to finish")
                event.ignore()
                return
        event.accept()


def configure_application(app):
    app.setApplicationName("DDMSoft")
    app.setOrganizationName("DDMSoft")
    app.setApplicationVersion("2.0")


def main(argv=None):
    app = QApplication.instance() or QApplication(argv or sys.argv)
    configure_application(app)
    window = DDMWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
