import numpy as np


def calculate_reynolds_number(density, velocity, diameter, viscosity):
    return density * velocity * diameter / viscosity


def calculate_pressure_drop(viscosity, length, mean_velocity, diameter):
    return 32 * viscosity * length * mean_velocity / diameter**2


def calculate_wall_shear_stress(viscosity, mean_velocity, diameter):
    return 8 * viscosity * mean_velocity / diameter


def classify_flow(reynolds_number):
    if reynolds_number < 2000:
        return "Laminar"
    elif reynolds_number < 4000:
        return "Transitional"
    else:
        return "High-Re / potentially turbulent"


def calculate_local_hemodynamics(
    density,
    viscosity,
    x,
    radius_x,
    local_mean_velocity
):
    diameter_x = 2 * radius_x

    reynolds_x = (
        density
        * local_mean_velocity
        * diameter_x
        / viscosity
    )

    wall_shear_x = (
        8
        * viscosity
        * local_mean_velocity
        / diameter_x
    )

    pressure_gradient_x = (
        32
        * viscosity
        * local_mean_velocity
        / diameter_x**2
    )

    cumulative_pressure_drop = np.zeros_like(x)

    for i in range(1, len(x)):
        dx = x[i] - x[i - 1]

        average_gradient = (
            pressure_gradient_x[i]
            + pressure_gradient_x[i - 1]
        ) / 2

        cumulative_pressure_drop[i] = (
            cumulative_pressure_drop[i - 1]
            + average_gradient * dx
        )

    return (
        diameter_x,
        reynolds_x,
        wall_shear_x,
        pressure_gradient_x,
        cumulative_pressure_drop
    )