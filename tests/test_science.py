import math

import pytest

from ddmsoft.science import (
    BOLTZMANN_CONSTANT,
    stokes_einstein_diffusion,
    stokes_einstein_radius,
    stokes_einstein_temperature,
    water_viscosity,
)


def test_water_viscosity_uses_kelvin_and_pa_seconds():
    viscosity = water_viscosity(293.15)
    assert 0.0008 < viscosity < 0.0012
    with pytest.raises(ValueError):
        water_viscosity(20.0)


def test_stokes_einstein_round_trips_in_si_units():
    temperature = 293.15
    viscosity = water_viscosity(temperature)
    diffusion = 2.0e-12
    radius = stokes_einstein_radius(diffusion, temperature, viscosity)
    assert math.isclose(stokes_einstein_diffusion(radius, temperature, viscosity), diffusion)
    assert math.isclose(
        stokes_einstein_temperature(diffusion, radius, viscosity), temperature, rel_tol=1e-12
    )
    assert BOLTZMANN_CONSTANT == 1.380649e-23
    with pytest.raises(ValueError):
        stokes_einstein_radius(0.0, temperature, viscosity)
