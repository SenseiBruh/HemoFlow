import numpy as np
import matplotlib.pyplot as plt


def calculate_reynolds_number(density, velocity, diameter, viscosity):
    return density * velocity * diameter / viscosity

def calculate_pressure_drop(viscosity, length, mean_velocity, diameter):
    """
    Hagen-Poiseuille pressure drop for fully developed laminar
    flow through a rigid circular tube.

    Returns pressure drop in Pascals.
    """
    return 32 * viscosity * length * mean_velocity / diameter**2


def calculate_wall_shear_stress(viscosity, mean_velocity, diameter):
    """
    Wall shear stress for fully developed laminar pipe flow.

    Returns shear stress in Pascals.
    """
    return 8 * viscosity * mean_velocity / diameter

def classify_flow(reynolds_number):
    """
    Simple pipe-flow heuristic.
    Real blood flow can be more complicated because vessels are
    pulsatile, elastic, curved, and nonuniform.
    """
    if reynolds_number < 2000:
        return "Laminar"
    elif reynolds_number < 4000:
        return "Transitional"
    else:
        return "High-Re / potentially turbulent"

def create_velocity_field(
    length,
    diameter,
    mean_velocity,
    stenosis_percent,
    stenosis_center,
    stenosis_length
    ):
    radius = diameter / 2

    x = np.linspace(0, length, 500)
    y = np.linspace(-radius, radius, 160)

    X, Y = np.meshgrid(x, y)

    # Convert stenosis percentage to decimal
    severity = stenosis_percent / 100

    # Prevent impossible geometry
    severity = np.clip(severity, 0, 0.95)

    # Start with a normal constant vessel radius
    radius_x = np.full_like(x, radius)

    # Create smooth stenosis
    if stenosis_length > 0 and severity > 0:

        half_length = stenosis_length / 2

        normalized_position = (
            (x - stenosis_center) / half_length
        )

        inside_stenosis = (
            np.abs(normalized_position) <= 1
        )

        # Smooth cosine-shaped narrowing
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

    # Conservation of flow:
    # A1 * V1 = A2 * V2
    local_mean_velocity = (
        mean_velocity
        * (radius / radius_x) ** 2
    )

    # Convert 1D radius arrays into 2D arrays
    radius_2d = radius_x[np.newaxis, :]
    velocity_2d = local_mean_velocity[np.newaxis, :]

    # Locally parabolic velocity profile
    velocity = (
        2
        * velocity_2d
        * (1 - (Y / radius_2d) ** 2)
    )

    # Hide points outside the vessel
    outside_vessel = np.abs(Y) > radius_2d

    velocity = np.ma.masked_where(
        outside_vessel,
        velocity
    )

    return (
        X,
        Y,
        velocity,
        x,
        radius_x,
        local_mean_velocity
    )

def plot_velocity_field(
    X,
    Y,
    velocity,
    x,
    radius_x,
    reynolds_number,
    flow_regime
    ):
    plt.figure(figsize=(12, 4))

    heatmap = plt.pcolormesh(
        X * 1000,
        Y * 1000,
        velocity,
        shading="auto",
        cmap="coolwarm"
    )

    plt.colorbar(
        heatmap,
        label="Blood Velocity (m/s)"
    )

    # Draw vessel walls
    plt.plot(
        x * 1000,
        radius_x * 1000,
        color="black",
        linewidth=2
    )

    plt.plot(
        x * 1000,
        -radius_x * 1000,
        color="black",
        linewidth=2
    )

    plt.xlabel("Distance Along Vessel (mm)")
    plt.ylabel("Radial Position (mm)")

    plt.title(
        f"Blood Flow Velocity | "
        f"Re = {reynolds_number:.0f} | "
        f"{flow_regime}"
    )

    plt.tight_layout()
    plt.show()


def main():
    print("\n--- HemoFlow ---")
    print("Straight Vessel Blood Flow Simulator\n")

    # Blood properties
    density = float(
        input("Blood density (kg/m^3) [1060]: ") or 1060
    )

    viscosity = float(
        input("Blood viscosity (Pa*s) [0.0035]: ") or 0.0035
    )

    # Vessel properties
    diameter_mm = float(
        input("Vessel diameter (mm) [4]: ") or 4
    )

    length_mm = float(
        input("Vessel length (mm) [60]: ") or 60
    )

    mean_velocity = float(
        input("Mean blood velocity (m/s) [0.30]: ") or 0.30
    )
    stenosis_percent = float(
    input("Stenosis severity (%) [50]: ") or 50
    )

    stenosis_center_mm = float(
    input("Stenosis center position (mm) [30]: ") or 30
    )

    stenosis_length_mm = float(
    input("Stenosis length (mm) [15]: ") or 15
    )

    # Convert mm -> meters (Reynolds Number equation assumes consistent SI units, so converting)
    diameter = diameter_mm / 1000
    length = length_mm / 1000
    stenosis_center = stenosis_center_mm / 1000
    stenosis_length = stenosis_length_mm / 1000

    reynolds_number = calculate_reynolds_number(
        density,
        mean_velocity,
        diameter,
        viscosity
    )
    pressure_drop = calculate_pressure_drop(
    viscosity,
    length,
    mean_velocity,
    diameter
    )

    wall_shear_stress = calculate_wall_shear_stress(
    viscosity,
    mean_velocity,
    diameter
    )

    flow_regime = classify_flow(reynolds_number)

    print("\n--- Results ---")
    print(f"Reynolds Number: {reynolds_number:.1f}")
    print(f"Flow Classification: {flow_regime}")
    print(f"Pressure Drop: {pressure_drop:.2f} Pa")
    print(f"Wall Shear Stress: {wall_shear_stress:.3f} Pa")

    (
    X,
    Y,
    velocity,
    x,
    radius_x,
    local_mean_velocity
    ) = create_velocity_field(
    length,
    diameter,
    mean_velocity,
    stenosis_percent,
    stenosis_center,
    stenosis_length
    )

    plot_velocity_field(
    X,
    Y,
    velocity,
    x,
    radius_x,
    reynolds_number,
    flow_regime
    )
    maximum_velocity = np.max(velocity)
    print(
    f"Maximum Local Velocity: "
    f"{maximum_velocity:.3f} m/s"
    )

if __name__ == "__main__":
    main()