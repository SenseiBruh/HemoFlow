import numpy as np


def create_velocity_field(
    length,
    diameter,
    mean_velocity,
    stenosis_percent,
    stenosis_center,
    stenosis_length,
):
    """
    Create an approximate velocity field through a vessel
    containing a smooth stenosis.

    This is currently a quasi-1D analytical approximation,
    not a full Navier-Stokes CFD solution.
    """

    radius = diameter / 2

    # Computational grid
    x = np.linspace(0, length, 500)
    y = np.linspace(-radius, radius, 160)

    X, Y = np.meshgrid(x, y)

    # Convert stenosis percentage to decimal
    severity = stenosis_percent / 100.0

    # Prevent impossible/zero-radius geometry
    severity = np.clip(severity, 0.0, 0.95)

    # Start with a constant-radius vessel
    radius_x = np.full_like(
        x,
        radius,
        dtype=float,
    )

    # ---------------------------------
    # Smooth cosine-shaped stenosis
    # ---------------------------------

    if stenosis_length > 0 and severity > 0:

        half_length = stenosis_length / 2

        normalized_position = (
            (x - stenosis_center)
            / half_length
        )

        inside_stenosis = (
            np.abs(normalized_position) <= 1
        )

        shape = 0.5 * (
            1
            + np.cos(
                np.pi
                * normalized_position[inside_stenosis]
            )
        )

        radius_x[inside_stenosis] = (
            radius
            * (1 - severity * shape)
        )

    # ---------------------------------
    # Conservation of volumetric flow
    #
    # A1 * V1 = A2 * V2
    # ---------------------------------

    local_mean_velocity = (
        mean_velocity
        * (radius / radius_x) ** 2
    )

    # Expand 1D arrays so they can operate
    # across the entire 2D computational grid
    radius_2d = radius_x[np.newaxis, :]
    mean_velocity_2d = (
        local_mean_velocity[np.newaxis, :]
    )

    # ---------------------------------
    # Locally parabolic flow profile
    # ---------------------------------

    velocity = (
        2
        * mean_velocity_2d
        * (
            1
            - (Y / radius_2d) ** 2
        )
    )

    # Mask everything outside vessel walls
    outside_vessel = (
        np.abs(Y) > radius_2d
    )

    velocity = np.ma.masked_where(
        outside_vessel,
        velocity,
    )

    return (
        X,
        Y,
        velocity,
        x,
        radius_x,
        local_mean_velocity,
    )