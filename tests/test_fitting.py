from __future__ import annotations

import numpy as np
import pytest

from ddmsoft.fitting import (
    MODEL_REGISTRY,
    FitCancelled,
    default_fit_request,
    evaluate_model,
    estimate_amplitude_background,
    fit_ddm,
)
from ddmsoft.models import DDMData, FitRange, FitRequest
from .fixtures import FIT_MODEL_IDS, generate_model_data


def test_registry_contains_every_legacy_model_with_stable_identifiers():
    assert tuple(MODEL_REGISTRY) == FIT_MODEL_IDS
    assert all(model.identifier == model_id for model_id, model in MODEL_REGISTRY.items())
    assert all(len(model.parameters) == len(model.defaults) for model in MODEL_REGISTRY.values())
    assert all(model.parameters[-2].identifier == "amplitude" for model in MODEL_REGISTRY.values())
    assert all(model.parameters[-1].identifier == "background" for model in MODEL_REGISTRY.values())


@pytest.mark.parametrize("model_id", FIT_MODEL_IDS)
def test_every_model_fits_generated_data_and_recreates_returned_curve(model_id):
    data = generate_model_data(model_id)
    fit_range = FitRange(
        q_min=0,
        q_max=data.q_values.size - 1,
        time_min=0,
        time_max=data.lag_times.size - 1,
    )
    result = fit_ddm(data, default_fit_request(model_id, fit_range))
    assert result.correlation.shape == data.matrix.shape
    assert result.fitted_matrix.shape == data.matrix.shape
    assert len(result.convergence_status) == data.q_values.size
    assert all(result.convergence_status)
    parameters = np.column_stack(result.model_parameters)
    for index, q_value in enumerate(result.q_values):
        correlation = evaluate_model(model_id, parameters[index], q_value, data.lag_times)
        expected = result.amplitude[index] * (1.0 - correlation) + result.noise[index]
        assert np.allclose(result.correlation[:, index], correlation)
        assert np.allclose(result.fitted_matrix[:, index], expected)
    assert np.mean((result.fitted_matrix - data.matrix) ** 2) < 1e-5


def test_inclusive_final_q_and_time_positions_are_fitted():
    data = generate_model_data("stretch")
    fit_range = FitRange(1, 3, 2, 6)
    result = fit_ddm(data, default_fit_request("stretch", fit_range))
    assert np.array_equal(result.q_values, data.q_values[1:4])
    assert result.correlation.shape == (5, 3)
    assert np.allclose(
        result.correlation[:, -1],
        evaluate_model(
            "stretch",
            [result.model_parameters[0][-1], result.model_parameters[1][-1]],
            data.q_values[-1],
            data.lag_times[2:7],
        ),
    )


def test_a_b_estimates_are_aligned_to_nonzero_q_min():
    data = generate_model_data("stretch")
    amplitude, background = estimate_amplitude_background(data.matrix)
    assert amplitude[2] != amplitude[3]
    assert background[2] != background[3]

    import ddmsoft.fitting as fitting

    initial_by_q: list[tuple[float, float]] = []

    def capture_initial(objective, initial, **kwargs):
        initial_by_q.append((initial[-2], initial[-1]))
        return type("Result", (), {"x": initial, "success": True, "message": ""})()

    monkeypatch = pytest.MonkeyPatch()
    monkeypatch.setattr(fitting, "minimize", capture_initial)
    try:
        fit_ddm(data, default_fit_request("stretch", FitRange(2, 3, 0, 6)))
    finally:
        monkeypatch.undo()
    assert np.allclose(
        initial_by_q,
        np.column_stack((amplitude[2:4], background[2:4])),
    )


def test_fixed_values_and_caller_owned_inputs_are_preserved():
    data = generate_model_data("stretch")
    initial = [2e-12, 0.82, None, None]
    fixed = [True, True, False, False]
    initial_before = initial.copy()
    fixed_before = fixed.copy()
    request = FitRequest("stretch", initial, fixed, FitRange(0, 3, 0, 6))
    result = fit_ddm(data, request)
    assert initial == initial_before
    assert fixed == fixed_before
    assert np.all(result.model_parameters[0] == 2e-12)
    assert np.all(result.model_parameters[1] == 0.82)


def test_failed_q_fit_is_reported_and_later_qs_continue(monkeypatch):
    data = generate_model_data("stretch")
    import ddmsoft.fitting as fitting

    original = fitting.minimize
    calls = 0

    def fail_once(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 1:
            raise ValueError("synthetic optimizer failure")
        return original(*args, **kwargs)

    monkeypatch.setattr(fitting, "minimize", fail_once)
    result = fit_ddm(data, default_fit_request("stretch", FitRange(0, 3, 0, 6)))
    assert result.convergence_status[0] is False
    assert "synthetic optimizer failure" in result.messages[0]
    assert len(result.convergence_status) == 4
    assert any(result.convergence_status[1:])


def test_fit_progress_and_cancellation_are_cooperative():
    data = generate_model_data("stretch")
    events: list[tuple[int, int]] = []
    state = {"cancel": False}

    def report(completed, total):
        events.append((completed, total))
        if completed == 1:
            state["cancel"] = True

    with pytest.raises(FitCancelled):
        fit_ddm(
            data,
            default_fit_request("stretch", FitRange(0, 3, 0, 6)),
            progress=report,
            cancel=lambda: state["cancel"],
        )
    assert events == [(1, 4)]
