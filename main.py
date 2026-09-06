import numpy as np

from physics import (
    calculate_reynolds_number,
    calculate_pressure_drop,
    calculate_wall_shear_stress,
    calculate_local_hemodynamics,
    classify_flow,
)

from solver import create_velocity_field

from visualization import (
    plot_velocity_field,
    plot_hemodynamic_profiles,
    animate_particles,
)


def main():
    print("\n--- HemoFlow ---")
    print("Straight Vessel Blood Flow Simulator\n")

    # -----------------------------
    # Blood properties
    # -----------------------------

    density = float(
        input("Blood density (kg/m^3) [1060]: ") or 1060
    )

    viscosity = float(
        input("Blood viscosity (Pa*s) [0.0035]: ") or 0.0035
    )

    # -----------------------------
    # Vessel properties
    # -----------------------------

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

    # -----------------------------
    # Pulsatile flow properties
    # -----------------------------

    heart_rate = float(
        input("Heart rate (BPM) [72]: ") or 72
    )

    pulsatility_percent = float(
        input("Pulsatility (%) [40]: ") or 40
    )

    # -----------------------------
    # Unit conversions
    # mm -> meters
    # -----------------------------

    diameter = diameter_mm / 1000
    length = length_mm / 1000

    stenosis_center = stenosis_center_mm / 1000
    stenosis_length = stenosis_length_mm / 1000

    # -----------------------------
    # Straight-vessel calculations
    # -----------------------------

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

    flow_regime = classify_flow(
        reynolds_number
    )

    print("\n--- Baseline Results ---")

    print(
        f"Reynolds Number: "
        f"{reynolds_number:.1f}"
    )

    print(
        f"Flow Classification: "
        f"{flow_regime}"
    )

    print(
        f"Straight-Vessel Pressure Drop: "
        f"{pressure_drop:.2f} Pa"
    )

    print(
        f"Straight-Vessel Wall Shear Stress: "
        f"{wall_shear_stress:.3f} Pa"
    )

    # -----------------------------
    # Create stenosed velocity field
    # -----------------------------

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

    # -----------------------------
    # Spatial hemodynamics
    # -----------------------------

    (
        diameter_x,
        reynolds_x,
        wall_shear_x,
        pressure_gradient_x,
        cumulative_pressure_drop
    ) = calculate_local_hemodynamics(
        density,
        viscosity,
        x,
        radius_x,
        local_mean_velocity
    )

    maximum_velocity = np.max(
        velocity
    )

    print("\n--- Local Hemodynamic Results ---")

    print(
        f"Maximum Local Velocity: "
        f"{maximum_velocity:.3f} m/s"
    )

    print(
        f"Maximum Reynolds Number: "
        f"{np.max(reynolds_x):.1f}"
    )

    print(
        f"Maximum Wall Shear Stress: "
        f"{np.max(wall_shear_x):.3f} Pa"
    )

    print(
        f"Total Approx. Pressure Drop: "
        f"{cumulative_pressure_drop[-1]:.2f} Pa"
    )

    print(
        f"Minimum Vessel Diameter: "
        f"{np.min(diameter_x) * 1000:.2f} mm"
    )

    # -----------------------------
    # Spatial profile plots
    # -----------------------------

    plot_hemodynamic_profiles(
        x,
        radius_x,
        local_mean_velocity,
        reynolds_x,
        wall_shear_x,
        cumulative_pressure_drop
    )

    # -----------------------------
    # Animated flow visualization
    # -----------------------------

    animation = animate_particles(
        X,
        Y,
        velocity,
        x,
        radius_x,
        local_mean_velocity,
        reynolds_number,
        flow_regime,
        heart_rate,
        pulsatility_percent
    )


if __name__ == "__main__":
    main()