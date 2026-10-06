"""Passive point tracers advected by computed (u,v); fixed memory footprint."""
import math
import numpy as np
from numba import njit


@njit(cache=True, nogil=True)
def _sample(field, x, y, nx, ny):
    x = min(max(x, 0.0), nx - 1.000001)
    y = min(max(y, 0.0), ny - 1.000001)
    ix, iy = int(x), int(y)
    a, b = x - ix, y - iy
    return ((1-a)*(1-b)*field[iy*nx+ix] + a*(1-b)*field[iy*nx+ix+1]
            + (1-a)*b*field[(iy+1)*nx+ix] + a*b*field[(iy+1)*nx+ix+1])


@njit(cache=True, nogil=True)
def _advect(x, y, vx, vy, ux, uy, solid, nx, ny, lattice_steps, max_speed,
            model_code, relaxation_steps):
    exited = np.zeros(len(x), dtype=np.bool_)
    # At most 0.20 cell per substep, with midpoint and end collision checks.
    particle_speed = 0.0
    if model_code == 1:
        for i in range(len(x)):
            particle_speed = max(particle_speed, math.hypot(vx[i], vy[i]))
    substeps = max(1, int(math.ceil(lattice_steps * max(max_speed, particle_speed) / 0.20)))
    h = lattice_steps / substeps
    blend = 1.0
    if model_code == 1:
        # Exponential Stokes relaxation is stable for large particles while
        # retaining a one-way, dilute-particle approximation.
        blend = 1.0 - math.exp(-h / max(relaxation_steps, 1e-12))
    for i in range(len(x)):
        for _ in range(substeps):
            u = _sample(ux, x[i], y[i], nx, ny)
            v = _sample(uy, x[i], y[i], nx, ny)
            if model_code == 0:
                vx[i], vy[i] = u, v
            else:
                # Predict a midpoint using the local fluid velocity, sample
                # again there, then relax the particle toward that velocity.
                mid_vx = vx[i] + 0.5 * blend * (u - vx[i])
                mid_vy = vy[i] + 0.5 * blend * (v - vy[i])
                mx, my = x[i] + 0.5*h*mid_vx, y[i] + 0.5*h*mid_vy
                um = _sample(ux, mx, my, nx, ny)
                vm = _sample(uy, mx, my, nx, ny)
                vx[i] += blend * (um - vx[i])
                vy[i] += blend * (vm - vy[i])
            xx = x[i] + h * vx[i]
            yy = y[i] + h * vy[i]
            if xx >= nx - 1 or xx < 0:
                exited[i] = True
                break
            if (my < 0.5 or my > ny - 1.5 or yy < 0.5 or yy > ny - 1.5):
                break
            if solid[int(round(my)), min(nx-1, max(0, int(round(mx))))] or solid[int(round(yy)), int(round(xx))]:
                # Keep the previous valid point: never jump across a plaque.
                break
            x[i], y[i] = xx, yy
    return exited


class Tracers:
    def __init__(self, solver, trail_length=32):
        self.rng = np.random.default_rng(solver.config.seed + 100_003)
        g = solver.geometry
        count = solver.config.particle_count
        chosen = self.rng.choice(solver.nodes, count)
        # Starting at fluid cell centers also populates pre-existing wake pockets.
        self.x = (chosen % solver.nx).astype(float)
        self.y = (chosen // solver.nx).astype(float)
        self.history = np.full((trail_length, count, 2), np.nan)
        self.history_index = 0
        self.inject_index = 0
        self.recycled = 0
        self.ids = np.zeros(count, dtype=np.uint32)
        self.vx = solver.ux[chosen].astype(float).copy()
        self.vy = solver.uy[chosen].astype(float).copy()
        self.record(solver)

    @property
    def model_code(self):
        return 1 if self._model == "inertial" else 0

    def _update_particle_model(self, solver):
        self._model = solver.config.particle_model
        diameter_m = solver.config.particle_diameter_um * 1e-6
        from physics import particle_relaxation_time
        self.relaxation_time_s = particle_relaxation_time(
            solver.config.particle_density_kg_m3, diameter_m, solver.config.viscosity)
        self.relaxation_steps = max(self.relaxation_time_s / solver.dt, 1e-8)

    def description(self, solver):
        self._update_particle_model(solver)
        return {
            "model": solver.config.particle_model,
            "diameter_um": float(solver.config.particle_diameter_um),
            "density_kg_m3": float(solver.config.particle_density_kg_m3),
            "relaxation_time_s": float(self.relaxation_time_s),
            "relaxation_steps": float(self.relaxation_steps),
        }

    def positions(self, solver):
        return np.column_stack((self.x * solver.dx * 1000,
                                (solver.geometry.y[0] + self.y * solver.dx) * 1000))

    def advance(self, solver, steps):
        self._update_particle_model(solver)
        maximum = float(np.hypot(solver.ux[solver.nodes], solver.uy[solver.nodes]).max())
        exited = _advect(self.x, self.y, self.vx, self.vy, solver.ux, solver.uy,
                         solver.geometry.solid, solver.nx, solver.ny, steps, maximum,
                         self.model_code, self.relaxation_steps)
        n = int(exited.sum())
        if n:
            self.x[exited] = 0.25
            self.y[exited] = self.rng.uniform(0.75, solver.ny - 1.75, n)
            self.vx[exited] = 0.0
            self.vy[exited] = 0.0
            self.history[:, exited] = np.nan
            self.ids[exited] += 1
            self.recycled += n

    def record(self, solver):
        self.history[self.history_index] = self.positions(solver)
        self.history_index = (self.history_index + 1) % len(self.history)

    def trails(self):
        history = np.concatenate((self.history[self.history_index:], self.history[:self.history_index]))
        return history.transpose(1, 0, 2)

    def inject(self, solver, x_mm, y_mm, count=60):
        """User-controlled local dye injection; replaces a bounded group of tracers."""
        x, y = x_mm / 1000 / solver.dx, (y_mm / 1000 - solver.geometry.y[0]) / solver.dx
        candidates_x = self.rng.normal(x, 0.65, count * 5)
        candidates_y = self.rng.normal(y, 0.65, count * 5)
        inserted = 0
        for xx, yy in zip(candidates_x, candidates_y):
            ix, iy = int(round(xx)), int(round(yy))
            if 0 <= ix < solver.nx and 1 <= iy < solver.ny-1 and not solver.geometry.solid[iy, ix]:
                idx = self.inject_index % len(self.x)
                self.x[idx], self.y[idx] = xx, yy
                self.vx[idx], self.vy[idx] = 0.0, 0.0
                self.history[:, idx] = np.nan
                self.ids[idx] += 1
                self.inject_index += 1
                inserted += 1
                if inserted == count:
                    break
        return inserted
