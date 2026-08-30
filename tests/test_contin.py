from __future__ import annotations

import numpy as np
import pytest

from ddmsoft.contin import (
    CONTINCancelled,
    CONTINError,
    contin_ranges,
    export_contin,
    gamma_to_radius,
    run_contin,
)
from ddmsoft.science import stokes_einstein_radius


def _synthetic_data():
    tau = np.linspace(0.01, 0.2, 8)
    gamma = np.linspace(1.0, 5.0, 5)
    distribution = np.array([0.0, 0.2, 0.6, 0.2, 0.0])
    correlation = np.exp(-tau[:, None] * gamma[None, :]) @ distribution
    data = 1.4 * (1.0 - correlation) + 0.08
    return tau, gamma, data


def test_contin_returns_all_candidates_and_minimum_residual_selection():
    tau, gamma, data = _synthetic_data()
    events: list[tuple[int, int]] = []
    result = run_contin(
        tau,
        data,
        gamma,
        alpha=[0.01, 0.1, 1.0],
        maxiter=1,
        progress=lambda completed, total: events.append((completed, total)),
    )
    assert result.alpha_g.shape == (3, 5)
    assert result.alpha_ddmfit.shape == (3, 8)
    assert result.selection_method == "minimum residual"
    assert result.selected_index == int(np.argmin(result.alpha_residuals))
    assert result.chosen_alpha == result.alphas[result.selected_index]
    assert events == [(1, 3), (2, 3), (3, 3)]


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"gamma_range": [1.0, 2.0]}, "gamma_range"),
        ({"gamma_range": [1.0, 2.0, 1.5]}, "increasing"),
        ({"gamma_range": [1.0, 2.0, 3.0], "alpha": []}, "alpha"),
        ({"gamma_range": [1.0, 2.0, 3.0], "maxiter": 0}, "maxiter"),
    ],
)
def test_contin_validates_gamma_alpha_and_iteration_ranges(kwargs, message):
    tau, _, data = _synthetic_data()
    gamma_range = kwargs.pop("gamma_range")
    with pytest.raises(CONTINError, match=message):
        run_contin(tau, data, gamma_range, **kwargs)


def test_contin_cancellation_happens_between_alpha_candidates():
    tau, gamma, data = _synthetic_data()
    state = {"cancel": False}

    def progress(completed, total):
        state["cancel"] = completed == 1

    with pytest.raises(CONTINCancelled):
        run_contin(
            tau,
            data,
            gamma,
            alpha=[0.01, 0.1],
            maxiter=1,
            progress=progress,
            cancel=lambda: state["cancel"],
        )


def test_contin_export_uses_each_candidate_amplitude_and_noise(tmp_path):
    tau, gamma, data = _synthetic_data()
    result = run_contin(tau, data, gamma, alpha=[0.01, 0.1], maxiter=1)
    output = export_contin(tmp_path / "contin.txt", result, video="sample", q=2.0)
    text = output.read_text(encoding="utf-8")
    assert "selection method:\tminimum residual" in text
    assert text.count("\nalpha:") == 2
    assert f"amplitude:\t{result.alpha_amplitude[0]:.03e}" in text
    assert f"amplitude:\t{result.alpha_amplitude[1]:.03e}" in text
    assert f"noise:\t{result.alpha_noise[0]:.03e}" in text
    assert f"noise:\t{result.alpha_noise[1]:.03e}" in text


def test_contin_radius_conversion_is_explicit_si():
    gamma = np.array([1e-12, 2e-12])
    radius = gamma_to_radius(gamma, 298.15, 1e-3)
    assert np.allclose(radius[0], stokes_einstein_radius(1e-12, 298.15, 1e-3))
    assert radius[1] < radius[0]


def test_contin_ranges_validate_counts_and_result_can_carry_sizes():
    gamma, alphas = contin_ranges(1e-12, 5e-12, 5, 0.01, 1.0, 3)
    assert gamma.size == 5
    assert alphas.size == 3
    tau, _, data = _synthetic_data()
    result = run_contin(tau, data, gamma, alpha=alphas, maxiter=1)
    sized = result.with_particle_sizes(298.15, 1e-3)
    assert sized.sizes is not None
    assert sized.sizes.shape == gamma.shape
    with pytest.raises(CONTINError, match="gamma_count"):
        contin_ranges(1, 2, 2, 0.1, 1.0, 2)
