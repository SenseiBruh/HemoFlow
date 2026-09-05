import numpy as np
import matplotlib.pyplot as plt


def calculate_reynolds_number(density, velocity, diameter, viscosity):
    return density * velocity * diameter / viscosity


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


def create_velocity_field(length, diameter, mean_velocity):
    radius = diameter / 2

    x = np.linspace(0, length, 400)
    y = np.linspace(-radius, radius, 120)

    X, Y = np.meshgrid(x, y)

    # Poiseuille velocity profile
    velocity = 2 * mean_velocity * (1 - (Y / radius) ** 2)

    return X, Y, velocity


def plot_velocity_field(X, Y, velocity, reynolds_number, flow_regime):
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

    plt.xlabel("Distance Along Vessel (mm)")
    plt.ylabel("Radial Position (mm)")

    plt.title(
        f"Blood Flow Velocity | Re = {reynolds_number:.0f} | {flow_regime}"
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

    # Convert mm -> meters
    diameter = diameter_mm / 1000
    length = length_mm / 1000

    reynolds_number = calculate_reynolds_number(
        density,
        mean_velocity,
        diameter,
        viscosity
    )

    flow_regime = classify_flow(reynolds_number)

    print("\n--- Results ---")
    print(f"Reynolds Number: {reynolds_number:.1f}")
    print(f"Flow Classification: {flow_regime}")

    X, Y, velocity = create_velocity_field(
        length,
        diameter,
        mean_velocity
    )

    plot_velocity_field(
        X,
        Y,
        velocity,
        reynolds_number,
        flow_regime
    )


if __name__ == "__main__":
    main()