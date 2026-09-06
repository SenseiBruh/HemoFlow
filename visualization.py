import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

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
        linewidth=2,
        zorder=11
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

def animate_particles(
    X,
    Y,
    velocity,
    x,
    radius_x,
    local_mean_velocity,
    reynolds_number,
    flow_regime,
    heart_rate,
    pulsatility_percent,
    number_of_particles=120
    ):
    fig, ax = plt.subplots(figsize=(12, 4))

    # Heartbeat parameters
    frequency = heart_rate / 60.0
    pulsatility = pulsatility_percent / 100.0

    # Keep pulsatility in a reasonable range for this simplified model
    pulsatility = np.clip(pulsatility, 0.0, 0.95)

    # Draw the velocity heatmap
    max_display_velocity = (
    np.max(velocity)
    * (1 + pulsatility)
    )

    heatmap = ax.pcolormesh(
    X * 1000,
    Y * 1000,
    velocity,
    shading="auto",
    cmap="coolwarm",
    vmin=0,
    vmax=max_display_velocity
    )
    
    status_text = ax.text(
    0.02,
    0.95,
    "",
    transform=ax.transAxes,
    verticalalignment="top",
    bbox=dict(
        boxstyle="round",
        facecolor="white",
        alpha=0.8
    )
    )

    fig.colorbar(
        heatmap,
        ax=ax,
        label="Blood Velocity (m/s)"
    )

    # Draw vessel walls
    ax.plot(
        x * 1000,
        radius_x * 1000,
        color="black",
        linewidth=2
    )

    ax.plot(
        x * 1000,
        -radius_x * 1000,
        color="black",
        linewidth=2
    )

    # ------------------------------------
    # Create particles
    # ------------------------------------

    rng = np.random.default_rng(42)

    particle_x = rng.uniform(
        x[0],
        x[-1],
        number_of_particles
    )

    # Normalized vertical position inside vessel
    # -1 = bottom wall
    #  0 = center
    # +1 = top wall
    particle_eta = rng.uniform(
        -0.85,
        0.85,
        number_of_particles
    )

    local_radius = np.interp(
        particle_x,
        x,
        radius_x
    )

    particle_y = (
        particle_eta * local_radius
    )

    particles = ax.scatter(
    particle_x * 1000,
    particle_y * 1000,
    s=28,
    facecolor="white",
    edgecolor="black",
    linewidth=0.7,
    zorder=10
    )

    ax.set_xlabel(
        "Distance Along Vessel (mm)"
    )

    ax.set_ylabel(
        "Radial Position (mm)"
    )

    ax.set_title(
        f"Animated Blood Flow | "
        f"Re = {reynolds_number:.0f} | "
        f"{flow_regime}"
    )

    # Simulation timestep
    dt = 0.005

    # ------------------------------------
    # Animation update function
    # ------------------------------------

    def update(frame):
        nonlocal particle_x
        nonlocal particle_eta

                # Simulation time
        time = frame * dt

        # Simple pulsatile heartbeat
        pulse_factor = (
            1
            + pulsatility
            * np.sin(
                2
                * np.pi
                * frequency
                * time
            )
        )

        # --------------------------------
        # Update heatmap
        # --------------------------------

        current_velocity_field = (
            velocity * pulse_factor
        )

        heatmap.set_array(
            current_velocity_field.ravel()
        )

        # --------------------------------
        # Update particles
        # --------------------------------

        mean_u = np.interp(
            particle_x,
            x,
            local_mean_velocity
        )

        # Apply heartbeat
        mean_u = mean_u * pulse_factor

        # Local parabolic velocity profile
        particle_velocity = (
            2
            * mean_u
            * (1 - particle_eta**2)
        )

        # Move particles
        particle_x += (
            particle_velocity * dt
        )

        # Recycle particles leaving vessel
        exited = particle_x > x[-1]

        number_exited = np.sum(exited)

        if number_exited > 0:
            particle_x[exited] = x[0]

            particle_eta[exited] = (
                rng.uniform(
                    -0.85,
                    0.85,
                    number_exited
                )
            )

        # Make particles follow vessel geometry
        local_radius = np.interp(
            particle_x,
            x,
            radius_x
        )

        particle_y = (
            particle_eta
            * local_radius
        )

        particles.set_offsets(
            np.column_stack(
                (
                    particle_x * 1000,
                    particle_y * 1000
                )
            )
        )

        # --------------------------------
        # Update status text
        # --------------------------------

        inlet_velocity = (
            local_mean_velocity[0]
            * pulse_factor
        )

        status_text.set_text(
            f"Time: {time:.2f} s\n"
            f"Heart Rate: {heart_rate:.0f} BPM\n"
            f"Inlet Mean Velocity: "
            f"{inlet_velocity:.3f} m/s"
        )

        return (
            particles,
            heatmap,
            status_text
        )

    animation = FuncAnimation(
    fig,
    update,
    interval=20,
    blit=False,
    cache_frame_data=False
    )

    plt.tight_layout()
    plt.show()

    return animation

def plot_hemodynamic_profiles(
    x,
    radius_x,
    local_mean_velocity,
    reynolds_x,
    wall_shear_x,
    cumulative_pressure_drop
):
    x_mm = x * 1000
    diameter_mm = radius_x * 2 * 1000

    # -------------------------
    # Velocity
    # -------------------------

    plt.figure(figsize=(10, 4))

    plt.plot(
        x_mm,
        local_mean_velocity
    )

    plt.xlabel("Position Along Vessel (mm)")
    plt.ylabel("Mean Velocity (m/s)")
    plt.title("Local Mean Blood Velocity")

    plt.grid(alpha=0.3)
    plt.tight_layout()


    # -------------------------
    # Reynolds number
    # -------------------------

    plt.figure(figsize=(10, 4))

    plt.plot(
        x_mm,
        reynolds_x
    )

    plt.xlabel("Position Along Vessel (mm)")
    plt.ylabel("Reynolds Number")
    plt.title("Local Reynolds Number")

    plt.grid(alpha=0.3)
    plt.tight_layout()


    # -------------------------
    # Wall shear stress
    # -------------------------

    plt.figure(figsize=(10, 4))

    plt.plot(
        x_mm,
        wall_shear_x
    )

    plt.xlabel("Position Along Vessel (mm)")
    plt.ylabel("Wall Shear Stress (Pa)")
    plt.title("Approximate Local Wall Shear Stress")

    plt.grid(alpha=0.3)
    plt.tight_layout()


    # -------------------------
    # Pressure loss
    # -------------------------

    plt.figure(figsize=(10, 4))

    plt.plot(
        x_mm,
        cumulative_pressure_drop
    )

    plt.xlabel("Position Along Vessel (mm)")
    plt.ylabel("Cumulative Pressure Drop (Pa)")
    plt.title("Approximate Pressure Loss Along Vessel")

    plt.grid(alpha=0.3)
    plt.tight_layout()


    # -------------------------
    # Vessel diameter
    # -------------------------

    plt.figure(figsize=(10, 4))

    plt.plot(
        x_mm,
        diameter_mm
    )

    plt.xlabel("Position Along Vessel (mm)")
    plt.ylabel("Vessel Diameter (mm)")
    plt.title("Vessel Geometry")

    plt.grid(alpha=0.3)
    plt.tight_layout()
