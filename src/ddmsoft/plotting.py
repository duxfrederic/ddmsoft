"""Independent Qt/Matplotlib controllers for analytical DDM plots."""

from __future__ import annotations

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QLabel, QMainWindow, QSlider, QVBoxLayout, QWidget

from .contin import CONTINResult
from .fitting import get_model
from .models import DDMData, FitRange, FitResult


def _correlation_from_ddm(data: DDMData) -> np.ndarray:
    """Create the legacy normalized correlation used by the analytical plot."""
    background = data.matrix[0]
    amplitude = (
        data.matrix[-3:].mean(axis=0) - background
        if data.matrix.shape[0] >= 3
        else data.matrix[-1] - background
    )
    return np.divide(
        amplitude[None, :] - (data.matrix - background[None, :]),
        amplitude[None, :],
        out=np.full_like(data.matrix, np.nan, dtype=float),
        where=amplitude[None, :] != 0,
    )


class _PlotController:
    """Common window ownership for modeless plots."""

    def __init__(self, title: str, figure: Figure) -> None:
        self.figure = figure
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, None)
        self.window = QMainWindow()
        self.window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self.window.setWindowTitle(title)
        self.window.setMinimumSize(700, 450)
        central = QWidget(self.window)
        self.layout = QVBoxLayout(central)
        self.layout.addWidget(self.toolbar)
        self.layout.addWidget(self.canvas, 1)
        self.window.setCentralWidget(central)

    def show(self) -> QMainWindow:
        """Show this controller's modeless window and return the window."""
        self.window.show()
        return self.window

    def close(self) -> None:
        self.window.close()


class CorrelationPlotController(_PlotController):
    """Interactive correlation/DDM plot with an independent q slider."""

    def __init__(
        self,
        data: DDMData,
        *,
        fit: FitResult | None = None,
        fit_range: FitRange | None = None,
        title: str = "DDM correlation",
    ) -> None:
        self.data = data
        self.fit = fit
        self.fit_range = fit_range or FitRange(
            0, data.q_values.size - 1, 0, data.lag_times.size - 1
        )
        self.measured_correlation = _correlation_from_ddm(data)
        self._fit_times = data.lag_times[self.fit_range.to_slices()[1]]
        if fit is not None and fit.correlation.shape[0] != self._fit_times.size:
            raise ValueError("fit correlation rows do not match fit_range time values")
        super().__init__(title, Figure(figsize=(11, 5)))
        self.correlation_axis, self.ddm_axis = self.figure.subplots(1, 2)
        self.figure.subplots_adjust(bottom=0.2, wspace=0.28)
        self._q_index = min(data.q_values.size // 3, data.q_values.size - 1)
        self._measured_correlation_line, = self.correlation_axis.semilogx(
            [], [], "s", label="measured"
        )
        self._measured_ddm_line, = self.ddm_axis.semilogx(
            [], [], "s", label="measured"
        )
        self._fit_correlation_line, = self.correlation_axis.semilogx(
            [], [], label="fit"
        )
        self._fit_ddm_line, = self.ddm_axis.semilogx([], [], label="fit")
        self._fit_correlation_line.set_visible(fit is not None)
        self._fit_ddm_line.set_visible(fit is not None)
        initial_q = data.q_values[self._q_index]
        marker_times = data.lag_times[[self.fit_range.time_min, self.fit_range.time_max]]
        positive_times = np.maximum(data.lag_times, np.finfo(float).tiny)
        self.correlation_axis.set_xlim(
            positive_times[0] * initial_q**2,
            positive_times[-1] * initial_q**2,
        )
        self.ddm_axis.set_xlim(positive_times[0], positive_times[-1])
        self._correlation_markers = tuple(
            self.correlation_axis.axvline(
                max(float(marker_time * initial_q**2), np.finfo(float).tiny),
                ls="--",
                color="0.5",
            )
            for marker_time in marker_times
        )
        self._ddm_markers = tuple(
            self.ddm_axis.axvline(
                max(float(marker_time), np.finfo(float).tiny),
                ls="--",
                color="0.5",
            )
            for marker_time in marker_times
        )
        self.correlation_axis.set_xlabel(r"tau q^2 [s/m^2]")
        self.correlation_axis.set_ylabel(r"f(q,tau)")
        self.correlation_axis.set_title("Correlation function")
        self.ddm_axis.set_xlabel("tau [s]")
        self.ddm_axis.set_ylabel("DDM")
        self.ddm_axis.set_title("DDM curve")
        self.correlation_axis.set_ylim(-0.1, 1.1)
        self.correlation_axis.legend(frameon=False)
        self.ddm_axis.legend(frameon=False)
        self.q_slider = QSlider(Qt.Orientation.Horizontal, self.window)
        self.q_slider.setRange(0, data.q_values.size - 1)
        self.q_slider.setValue(self._q_index)
        self.q_slider.setSingleStep(1)
        self.q_slider.setToolTip("Select the q column to inspect")
        self.q_label = QLabel(self.window)
        self.q_label.setToolTip("Selected q index and physical wavenumber")
        self.layout.addWidget(self.q_label)
        self.layout.addWidget(self.q_slider)
        self.q_slider.valueChanged.connect(self.set_q_index)
        self.canvas.mpl_connect("key_press_event", self._on_key)
        self.canvas.mpl_connect("scroll_event", self._on_scroll)
        self.set_q_index(self._q_index)

    @property
    def q_index(self) -> int:
        return self._q_index

    def set_q_index(self, index: int) -> None:
        """Select a q column and update every dependent artist."""
        index = max(0, min(int(index), self.data.q_values.size - 1))
        self._q_index = index
        if self.q_slider.value() != index:
            self.q_slider.blockSignals(True)
            self.q_slider.setValue(index)
            self.q_slider.blockSignals(False)
        q_value = self.data.q_values[index]
        self.q_label.setText(
            f"q index {index} / {self.data.q_values.size - 1}; q = {q_value:.6g} m^-1"
        )
        times = self.data.lag_times
        self._measured_correlation_line.set_data(
            times * q_value**2, self.measured_correlation[:, index]
        )
        self._measured_ddm_line.set_data(times, self.data.matrix[:, index])
        marker_times = self.data.lag_times[
            [self.fit_range.time_min, self.fit_range.time_max]
        ]
        for marker, marker_time in zip(self._correlation_markers, marker_times):
            marker.set_xdata([marker_time * q_value**2, marker_time * q_value**2])
        for marker, marker_time in zip(self._ddm_markers, marker_times):
            marker.set_xdata([marker_time, marker_time])
        if self.fit is not None:
            fit_index = int(np.argmin(np.abs(self.fit.q_values - q_value)))
            if np.isclose(self.fit.q_values[fit_index], q_value, rtol=1e-9, atol=0.0):
                self._fit_correlation_line.set_visible(True)
                self._fit_ddm_line.set_visible(True)
                self._fit_correlation_line.set_data(
                    self._fit_times * self.fit.q_values[fit_index] ** 2,
                    self.fit.correlation[:, fit_index],
                )
                self._fit_ddm_line.set_data(
                    self._fit_times,
                    self.fit.fitted_matrix[:, fit_index],
                )
            else:
                self._fit_correlation_line.set_visible(False)
                self._fit_ddm_line.set_visible(False)
        self.correlation_axis.relim()
        self.correlation_axis.autoscale_view(scalex=True, scaley=False)
        self.ddm_axis.relim()
        self.ddm_axis.autoscale_view()
        self.canvas.draw_idle()

    def _on_key(self, event: object) -> None:
        key = getattr(event, "key", None)
        if key == "right":
            self.set_q_index(self.q_index + 1)
        elif key == "left":
            self.set_q_index(self.q_index - 1)

    def _on_scroll(self, event: object) -> None:
        step = getattr(event, "step", 0)
        button = getattr(event, "button", None)
        if step > 0 or button == "up":
            self.set_q_index(self.q_index + 1)
        elif step < 0 or button == "down":
            self.set_q_index(self.q_index - 1)


class MatrixPlotController(_PlotController):
    """Modeless measured/fitted matrix image window."""

    def __init__(
        self,
        data: DDMData,
        *,
        fit: FitResult | None = None,
        fit_range: FitRange | None = None,
        title: str = "DDM matrix",
    ) -> None:
        self.data = data
        self.fit = fit
        self.fit_range = fit_range
        super().__init__(title, Figure(figsize=(11, 5)))
        self.measured_axis, self.fit_axis = self.figure.subplots(1, 2, sharey=True)
        self.figure.subplots_adjust(wspace=0.12)
        extent = (
            data.q_values[0],
            data.q_values[-1],
            data.lag_times[0],
            data.lag_times[-1],
        )
        self.measured_image = self.measured_axis.imshow(
            data.matrix, origin="lower", aspect="auto", extent=extent
        )
        self.measured_axis.set_title("Measured DDM matrix")
        self.measured_axis.set_xlabel("q [m^-1]")
        self.measured_axis.set_ylabel("lag time [s]")
        self.fit_image = None
        if fit is not None:
            if fit_range is None:
                raise ValueError("fit_range is required when a fitted matrix is supplied")
            q_slice, time_slice = fit_range.to_slices()
            fit_extent = (
                data.q_values[q_slice][0],
                data.q_values[q_slice][-1],
                data.lag_times[time_slice][0],
                data.lag_times[time_slice][-1],
            )
            self.fit_image = self.fit_axis.imshow(
                fit.fitted_matrix,
                origin="lower",
                aspect="auto",
                extent=fit_extent,
            )
        self.fit_axis.set_title("Fitted DDM matrix" if fit is not None else "No fitted matrix")
        self.fit_axis.set_xlabel("q [m^-1]")


class FitParameterPlotController(_PlotController):
    """Plot every fitted parameter against q, including failed-q diagnostics."""

    def __init__(self, result: FitResult, *, title: str = "Fitted parameters") -> None:
        self.result = result
        model = get_model(result.model_id)
        super().__init__(title, Figure(figsize=(9, 5)))
        self.axis = self.figure.subplots()
        self.figure.subplots_adjust(bottom=0.16)
        self.parameter_names = tuple(
            parameter.display_name for parameter in model.parameters[:-2]
        ) + ("Amplitude", "Background")
        values = (*result.model_parameters, result.amplitude, result.noise)
        self.parameter_lines = tuple(
            self.axis.semilogx(result.q_values, values[index], marker="o", label=name)[0]
            for index, name in enumerate(self.parameter_names)
        )
        failed = result.q_values[~np.asarray(result.convergence_status, dtype=bool)]
        self.failed_q_values = failed
        self.failed_markers = tuple(
            self.axis.axvline(q_value, color="#b00020", ls=":", alpha=0.7)
            for q_value in failed
        )
        self.axis.set_xlabel("q [m^-1]")
        self.axis.set_ylabel("parameter value")
        self.axis.set_title("Fitted parameters")
        self.axis.legend(frameon=False)


class AmplitudeNoiseDiffusionPlotController(_PlotController):
    """Plot amplitude, background noise, and the primary diffusion parameter."""

    def __init__(self, result: FitResult, *, title: str = "Amplitude, noise, diffusion") -> None:
        self.result = result
        model = get_model(result.model_id)
        diffusion_index = next(
            (
                index
                for index, parameter in enumerate(model.physical_parameters)
                if parameter.identifier.startswith("diffusion")
            ),
            0,
        )
        diffusion_name = model.physical_parameters[diffusion_index].display_name
        self.diffusion_name = diffusion_name
        super().__init__(title, Figure(figsize=(11, 4)))
        self.axes = tuple(self.figure.subplots(1, 3))
        self.figure.subplots_adjust(wspace=0.3)
        series = (
            ("Amplitude", result.amplitude),
            ("Background noise", result.noise),
            (diffusion_name, result.model_parameters[diffusion_index]),
        )
        self.lines = tuple(
            axis.semilogx(result.q_values, values, marker="o")[0]
            for axis, (_, values) in zip(self.axes, series)
        )
        for axis, (label, _) in zip(self.axes, series):
            axis.set_title(label)
            axis.set_xlabel("q [m^-1]")
            axis.set_ylabel(label)
        failed = result.q_values[~np.asarray(result.convergence_status, dtype=bool)]
        self.failed_q_values = failed
        self.failed_markers = tuple(
            tuple(
                axis.axvline(q_value, color="#b00020", ls=":", alpha=0.7)
                for axis in self.axes
            )
            for q_value in failed
        )


class CONTINPlotController(_PlotController):
    """Interactive CONTIN alpha plot with no module-global state."""

    def __init__(self, result: CONTINResult, *, title: str = "CONTIN") -> None:
        self.result = result
        super().__init__(title, Figure(figsize=(11, 5)))
        self.correlation_axis, self.distribution_axis = self.figure.subplots(1, 2)
        self.figure.subplots_adjust(bottom=0.2, wspace=0.28)
        self._alpha_index = result.selected_index
        self.fit_line, = self.correlation_axis.semilogx([], [], label="CONTIN")
        self.data_line, = self.correlation_axis.semilogx(
            result.tau, result.ddmdata, ".", label="measured"
        )
        self.distribution_line, = self.distribution_axis.plot([], [], label="distribution")
        self.correlation_axis.set_xlabel("tau q^2 [s/m^2]")
        self.correlation_axis.set_ylabel("DDM")
        self.distribution_axis.set_ylabel("intensity")
        self.correlation_axis.legend(frameon=False)
        self.distribution_axis.legend(frameon=False)
        self.alpha_slider = QSlider(Qt.Orientation.Horizontal, self.window)
        self.alpha_slider.setRange(0, result.alphas.size - 1)
        self.alpha_slider.setValue(self._alpha_index)
        self.alpha_slider.setSingleStep(1)
        self.alpha_slider.setToolTip("Inspect each retained CONTIN alpha candidate")
        self.alpha_label = QLabel(self.window)
        self.layout.addWidget(self.alpha_label)
        self.layout.addWidget(self.alpha_slider)
        self.alpha_slider.valueChanged.connect(self.set_alpha_index)
        self.canvas.mpl_connect("key_press_event", self._on_key)
        self.canvas.mpl_connect("scroll_event", self._on_scroll)
        self.set_alpha_index(self._alpha_index)

    @property
    def alpha_index(self) -> int:
        return self._alpha_index

    def set_alpha_index(self, index: int) -> None:
        """Select an alpha candidate and update its fit and distribution."""
        index = max(0, min(int(index), self.result.alphas.size - 1))
        self._alpha_index = index
        if self.alpha_slider.value() != index:
            self.alpha_slider.blockSignals(True)
            self.alpha_slider.setValue(index)
            self.alpha_slider.blockSignals(False)
        amplitude = self.result.alpha_amplitude[index]
        noise = self.result.alpha_noise[index]
        self.fit_line.set_data(
            self.result.tau,
            amplitude * (1.0 - self.result.alpha_ddmfit[index]) + noise,
        )
        distribution = self.result.alpha_g[index]
        distribution = distribution / np.sum(distribution) if np.sum(distribution) else distribution
        x_values = self.result.gamma_range
        x_label = "gamma [m^2/s]"
        if self.result.sizes is not None:
            x_values = self.result.sizes * 1e9
            x_label = "hydrodynamic radius [nm]"
        self.distribution_line.set_data(x_values, distribution)
        self.distribution_axis.set_xlabel(x_label)
        self.alpha_label.setText(
            f"alpha index {index} / {self.result.alphas.size - 1}; "
            f"alpha = {self.result.alphas[index]:.6g}"
        )
        self.correlation_axis.relim()
        self.correlation_axis.autoscale_view()
        self.distribution_axis.relim()
        self.distribution_axis.autoscale_view()
        self.canvas.draw_idle()

    def _on_key(self, event: object) -> None:
        key = getattr(event, "key", None)
        if key == "right":
            self.set_alpha_index(self.alpha_index + 1)
        elif key == "left":
            self.set_alpha_index(self.alpha_index - 1)

    def _on_scroll(self, event: object) -> None:
        step = getattr(event, "step", 0)
        button = getattr(event, "button", None)
        if step > 0 or button == "up":
            self.set_alpha_index(self.alpha_index + 1)
        elif step < 0 or button == "down":
            self.set_alpha_index(self.alpha_index - 1)


# Names that describe the windows in the roadmap and preserve a compact API.
DDMPlotController = CorrelationPlotController
CONTINController = CONTINPlotController


__all__ = [
    "AmplitudeNoiseDiffusionPlotController",
    "CONTINController",
    "CONTINPlotController",
    "CorrelationPlotController",
    "DDMPlotController",
    "FitParameterPlotController",
    "MatrixPlotController",
]
