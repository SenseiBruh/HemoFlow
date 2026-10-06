"""Keep the original circular-pipe references separate from the planar CFD."""
import numpy as np


def calculate_reynolds_number(density, velocity, diameter, viscosity):
    return density * velocity * diameter / viscosity


def calculate_pressure_drop(viscosity, length, mean_velocity, diameter):
    """Original healthy circular pipe reference, Pa; NOT a stenosis CFD result."""
    return 32 * viscosity * length * mean_velocity / diameter**2


def calculate_wall_shear_stress(viscosity, mean_velocity, diameter):
    """Original healthy circular pipe reference, Pa."""
    return 8 * viscosity * mean_velocity / diameter


def calculate_womersley_number(density, heart_rate, diameter, viscosity):
    """Return the circular-vessel Womersley number for a sinusoidal scale.

    This is a dimensionless comparison value, not a replacement for a
    measured waveform or a claim that the planar solver is a circular artery.
    """
    omega = 2 * np.pi * heart_rate / 60.0
    return (diameter / 2.0) * np.sqrt(density * omega / viscosity)


def particle_relaxation_time(density, diameter_m, viscosity):
    """Stokes relaxation time for a dilute spherical one-way particle (s)."""
    return density * diameter_m**2 / (18.0 * viscosity)


def calculate_particle_stokes_number(density, diameter_m, viscosity, heart_rate):
    """Particle response time divided by the cardiac angular-time scale."""
    return particle_relaxation_time(density, diameter_m, viscosity) * (2 * np.pi * heart_rate / 60.0)


def classify_flow(reynolds_number):
    """Retained pipe-flow heuristic, not a local diagnosis of a stenotic flow."""
    if reynolds_number < 2000:
        return "Laminar (pipe reference)"
    if reynolds_number < 4000:
        return "Transitional (pipe reference)"
    return "High-Re / potentially turbulent (pipe reference)"


def pulse_factor(config, time):
    from waveform import pulse_curve
    return pulse_curve(config, time)


def planar_reference(config):
    """Analytic steady, fully developed planar Poiseuille values for verification."""
    return (12 * config.viscosity * config.length * config.mean_velocity / config.diameter**2,
            6 * config.viscosity * config.mean_velocity / config.diameter)
