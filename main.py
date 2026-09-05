import numpy as np
import matplotlib.pyplot as plt

# Blood properties
density = 1060          # kg/m^3
viscosity = 0.0035      # Pa*s

# Vessel properties
diameter = 0.004        # meters
length = 0.060          # meters
mean_velocity = 0.30    # m/s

# Reynolds number
reynolds_number = density * mean_velocity * diameter / viscosity

print(f"Reynolds Number: {reynolds_number:.1f}")

# Vessel radius
radius = diameter / 2

# Create coordinate grid
x = np.linspace(0, length, 400)
y = np.linspace(-radius, radius, 120)

X, Y = np.meshgrid(x, y)

# Poiseuille velocity profile
velocity = 2 * mean_velocity * (1 - (Y / radius) ** 2)

# Plot heatmap
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
    f"Blood Flow Velocity | Re = {reynolds_number:.0f}"
)

plt.tight_layout()
plt.show()