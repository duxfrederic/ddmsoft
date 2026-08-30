"""Small scientific conversions expressed with explicit SI units."""

from __future__ import annotations

import math

BOLTZMANN_CONSTANT = 1.380649e-23  # J/K


def water_viscosity(temperature_kelvin: float) -> float:
    """Return water's dynamic viscosity in Pa s for an absolute temperature."""
    if not math.isfinite(temperature_kelvin) or temperature_kelvin <= 140.0:
        raise ValueError("temperature_kelvin must be finite and greater than 140 K")
    return 2.414e-5 * math.exp(570.6 / (temperature_kelvin - 140.0))


def stokes_einstein_radius(
    diffusion_m2_s: float, temperature_kelvin: float, viscosity_pa_s: float
) -> float:
    """Convert diffusion (m2/s) to hydrodynamic radius (m)."""
    _positive(diffusion_m2_s, "diffusion_m2_s")
    _positive(temperature_kelvin, "temperature_kelvin")
    _positive(viscosity_pa_s, "viscosity_pa_s")
    return BOLTZMANN_CONSTANT * temperature_kelvin / (
        6.0 * math.pi * viscosity_pa_s * diffusion_m2_s
    )


def stokes_einstein_diffusion(
    radius_m: float, temperature_kelvin: float, viscosity_pa_s: float
) -> float:
    """Convert hydrodynamic radius (m) to diffusion (m2/s)."""
    _positive(radius_m, "radius_m")
    _positive(temperature_kelvin, "temperature_kelvin")
    _positive(viscosity_pa_s, "viscosity_pa_s")
    return BOLTZMANN_CONSTANT * temperature_kelvin / (6.0 * math.pi * viscosity_pa_s * radius_m)


def stokes_einstein_temperature(
    diffusion_m2_s: float, radius_m: float, viscosity_pa_s: float
) -> float:
    """Convert diffusion, radius, and viscosity to absolute temperature (K)."""
    _positive(diffusion_m2_s, "diffusion_m2_s")
    _positive(radius_m, "radius_m")
    _positive(viscosity_pa_s, "viscosity_pa_s")
    return diffusion_m2_s * 6.0 * math.pi * viscosity_pa_s * radius_m / BOLTZMANN_CONSTANT


def _positive(value: float, name: str) -> None:
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
