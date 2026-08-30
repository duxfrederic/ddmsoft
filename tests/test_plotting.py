from __future__ import annotations

import numpy as np
import pytest

from ddmsoft.contin import run_contin
from ddmsoft.fitting import default_fit_request, fit_ddm
from ddmsoft.models import DDMData, FitRange
from ddmsoft.plotting import (
    CONTINPlotController,
    CorrelationPlotController,
    MatrixPlotController,
)


@pytest.fixture
def app():
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    instance = QApplication.instance() or QApplication([])
    yield instance
    for widget in instance.topLevelWidgets():
        widget.close()


def _data() -> DDMData:
    return DDMData(
        np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.8, 0.9, 1.0]]),
        np.array([0.01, 0.1, 1.0]),
        np.array([1e6, 2e6, 3e6]),
    )


def test_correlation_plot_q_callback_updates_owned_artists(app):
    data = _data()
    controller = CorrelationPlotController(data)
    controller.set_q_index(2)
    assert controller.q_index == 2
    assert "q index 2" in controller.q_label.text()
    assert np.array_equal(
        controller._measured_ddm_line.get_ydata(), data.matrix[:, 2]
    )
    controller._on_key(type("Event", (), {"key": "left"})())
    assert controller.q_index == 1
    controller._on_scroll(type("Event", (), {"step": 1, "button": "up"})())
    assert controller.q_index == 2
    controller.close()


def test_main_plot_windows_and_contin_windows_are_independent(app):
    data = _data()
    first = CorrelationPlotController(data)
    second = CorrelationPlotController(data)
    assert first.figure is not second.figure
    first.set_q_index(0)
    second.set_q_index(2)
    first.close()
    assert second.q_index == 2

    tau = np.linspace(0.01, 0.2, 6)
    gamma = np.linspace(1.0, 5.0, 5)
    contin_data = 1.2 * (1 - np.exp(-tau * 3.0)) + 0.1
    result = run_contin(tau, contin_data, gamma, alpha=[0.01, 0.1], maxiter=1)
    contin_first = CONTINPlotController(result)
    contin_second = CONTINPlotController(result)
    contin_first.set_alpha_index(0)
    contin_second.set_alpha_index(1)
    contin_first.close()
    assert contin_second.alpha_index == 1
    contin_second.close()
    second.close()


def test_matrix_plot_contains_measured_and_fitted_artists(app):
    data = _data()
    fit_range = FitRange(0, 2, 0, 2)
    result = fit_ddm(data, default_fit_request("stretch", fit_range))
    controller = MatrixPlotController(data, fit=result, fit_range=fit_range)
    assert controller.measured_image is not None
    assert controller.fit_image is not None
    controller.close()
