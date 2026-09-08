"""D2Q9 stress-regularized MRT lattice Boltzmann solver in physical units.

Weakly compressible, Newtonian, rigid-wall planar model. Halfway bounce-back
on the rasterized wall; Zou-He velocity inlet and fixed-density outlet.
No synthetic swirl, stochastic body forces, or prescribed particle paths.
"""
import math
import numpy as np
from numba import njit
from geometry import make_geometry, stenosis_report
from kernel import advance as _advance_fast
from kernel import advance_parallel
from physics import calculate_womersley_number, particle_relaxation_time, calculate_particle_stokes_number
from waveform import pulse_value, pulse_curve, phase_label

CX = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1], dtype=np.int32)
CY = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1], dtype=np.int32)
OPP = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6], dtype=np.int32)
W = np.array([4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36])


@njit(cache=True, nogil=True)
def _advance(f, out, nodes, sources, dirs, columns, profile, nx, omega_p,
             rho, ux, uy, steps, iteration, dt, amplitude, frequency):
    g = np.empty(9)
    eq = np.empty(9)
    for step in range(steps):
        time = (iteration + step + 1) * dt
        pulse = 1 + amplitude * math.sin(2 * math.pi * frequency * time)
        for n in range(len(nodes)):
            node = nodes[n]
            for k in range(9):
                g[k] = f[sources[n, k], dirs[n, k]]
            column = columns[n]
            if column == 0:
                u = profile[node // nx] * pulse
                r = (g[0] + g[2] + g[4] + 2 * (g[3] + g[6] + g[7])) / (1 - u)
                g[1] = g[3] + 2 * r * u / 3
                g[5] = g[7] + (g[4] - g[2]) / 2 + r * u / 6
                g[8] = g[6] + (g[2] - g[4]) / 2 + r * u / 6
            elif column == nx - 1:
                r = 1.0
                u = -1 + (g[0] + g[2] + g[4] + 2 * (g[1] + g[5] + g[8])) / r
                g[3] = g[1] - 2 * r * u / 3
                g[6] = g[8] + (g[4] - g[2]) / 2 - r * u / 6
                g[7] = g[5] + (g[2] - g[4]) / 2 - r * u / 6
            r = 0.0
            u = 0.0
            v = 0.0
            for k in range(9):
                r += g[k]
                u += g[k] * CX[k]
                v += g[k] * CY[k]
            if r <= 0 or not math.isfinite(r):
                raise FloatingPointError("Non-finite or non-positive fluid density.")
            u /= r
            v /= r
            rho[node], ux[node], uy[node] = r, u, v
            uu = u * u + v * v
            for k in range(9):
                cu = CX[k] * u + CY[k] * v
                eq[k] = W[k] * r * (1 + 3 * cu + 4.5 * cu * cu - 1.5 * uu)
            pxx, pxy, pyy = 0.0, 0.0, 0.0
            for k in range(9):
                neq = g[k] - eq[k]
                pxx += neq * CX[k] * CX[k]
                pxy += neq * CX[k] * CY[k]
                pyy += neq * CY[k] * CY[k]
            for k in range(9):
                # Relax the bulk/acoustic and higher moments at rate one;
                # retain the physical shear relaxation on the deviatoric stress.
                projected = 4.5 * W[k] * (0.5*(CX[k]*CX[k]-CY[k]*CY[k])*(pxx-pyy) +
                                         2*CX[k]*CY[k]*pxy)
                out[node, k] = eq[k] + (1 - omega_p) * projected
        f, out = out, f
    return f, out


class FlowSolver:
    def __init__(self, config):
        self.config = config.validate()
        self.geometry = make_geometry(config)
        g = self.geometry
        self.ny, self.nx = g.solid.shape
        self.dx = g.dx
        gap = g.min_gap_cells / config.cells_across
        # Nominal throat maximum <= 0.10 lattice units at peak inlet.
        # Keep actual Mach/density variation visible; never silently change viscosity.
        u_lattice = min(0.045, 0.10 * gap / (1.5 * (1 + config.pulsatility_percent / 100)))
        self.dt = u_lattice * self.dx / config.mean_velocity
        self.velocity_scale = self.dx / self.dt
        nu = config.viscosity / config.density * self.dt / self.dx**2
        self.tau = 0.5 + 3 * nu
        if self.tau < 0.501 or self.tau > 2.0:
            raise ValueError(f"Unresolved viscosity/time scale (tau={self.tau:.5f}). Increase cells_across for high Re, or lower it for very viscous flow.")
        self.omega_p = 1 / self.tau
        self.iteration = 0
        self.kernel = _advance_fast
        self.threads_used = 1
        self.nodes = np.flatnonzero(g.fluid.ravel()).astype(np.int32)
        self.columns = self.nodes % self.nx
        rows = self.nodes // self.nx
        self.sources = np.empty((len(self.nodes), 9), dtype=np.int32)
        self.dirs = np.empty_like(self.sources)
        for k in range(9):
            sy, sx = rows - CY[k], self.columns - CX[k]
            valid = (sy >= 0) & (sy < self.ny) & (sx >= 0) & (sx < self.nx)
            source = np.clip(sy, 0, self.ny - 1) * self.nx + np.clip(sx, 0, self.nx - 1)
            bounce = ~valid | g.solid.ravel()[source]
            self.sources[:, k] = np.where(bounce, self.nodes, source)
            self.dirs[:, k] = np.where(bounce, OPP[k], k)
        # Flux-normalized discrete planar inlet: mean input is exact on this grid.
        parabolic = np.maximum(0, 1 - (2 * g.y / config.diameter)**2)
        parabolic[~g.fluid[:, 0]] = 0
        parabolic /= parabolic[g.fluid[:, 0]].mean()
        self.profile = parabolic * u_lattice
        # Initialize from a divergence-free streamfunction through the local gap.
        # This is just a start condition; the solver then evolves momentum itself.
        height = g.upper - g.lower
        eta = np.clip((g.y[:, None] - g.lower) / height, 0, 1)
        initial_pulse = pulse_value(0.0, config.pulsatility_percent/100,
                                    config.heart_rate/60, int(config.pulse_shape == "systolic"))
        q = config.mean_velocity * config.diameter * initial_pulse
        u = 6 * q / height * eta * (1 - eta)
        slope = np.gradient(g.lower, self.dx)[None, :] + eta * np.gradient(height, self.dx)[None, :]
        v = u * slope
        u[g.solid], v[g.solid] = 0, 0
        u[:, 0] = self.profile * self.velocity_scale * initial_pulse
        v[:, 0] = v[:, -1] = 0
        self.ux = (u / self.velocity_scale).ravel().copy()
        self.uy = (v / self.velocity_scale).ravel().copy()
        # A viscous pressure estimate reduces the initial acoustic transient.
        # It is an initial condition, not a pressure field imposed during evolution.
        gradient = 12 * config.viscosity * q / height**3
        pressure_1d = np.zeros(self.nx)
        pressure_1d[:-1] = np.cumsum((0.5 * (gradient[1:] + gradient[:-1]) * self.dx)[::-1])[::-1]
        density_1d = 1 + pressure_1d * 3 / (config.density * self.velocity_scale**2)
        self.rho = np.broadcast_to(density_1d, (self.ny, self.nx)).copy().ravel()
        self.rho[g.solid.ravel()] = 1.0
        self.f = np.empty((self.nx * self.ny, 9))
        speed2 = self.ux**2 + self.uy**2
        for k in range(9):
            cu = CX[k] * self.ux + CY[k] * self.uy
            self.f[:, k] = W[k] * self.rho * (1 + 3 * cu + 4.5 * cu**2 - 1.5 * speed2)
        self.other = self.f.copy()

    @property
    def time(self):
        return self.iteration * self.dt

    def step(self, steps=1):
        if not isinstance(steps, int) or steps < 1:
            raise ValueError("steps must be a positive integer")
        c = self.config
        self.f, self.other = self.kernel(self.f, self.other, self.nodes, self.sources, self.dirs,
                                     self.columns, self.profile, self.nx, self.omega_p,
                                     self.rho, self.ux, self.uy, steps, self.iteration, self.dt,
                                     c.pulsatility_percent / 100, c.heart_rate / 60,
                                     int(c.pulse_shape == "systolic"))
        self.iteration += steps

    def set_threads(self, count):
        import numba
        count = min(max(1,int(count)),numba.config.NUMBA_NUM_THREADS)
        numba.set_num_threads(count)
        self.kernel = _advance_fast if count == 1 else advance_parallel
        self.threads_used = count

    def tune_threads(self):
        """Benchmark short scratch runs; keep the actual fluid state untouched."""
        import numba
        from time import perf_counter
        if self.config.compute_threads:
            self.set_threads(self.config.compute_threads)
            return self.threads_used
        timings = {}
        c = self.config
        for count in (1,2,4,8):
            if count > numba.config.NUMBA_NUM_THREADS:
                continue
            self.set_threads(count)
            f, other = self.f.copy(),self.other.copy()
            r,u,v = self.rho.copy(),self.ux.copy(),self.uy.copy()
            args=(f,other,self.nodes,self.sources,self.dirs,self.columns,self.profile,self.nx,self.omega_p,
                  r,u,v,24,0,self.dt,c.pulsatility_percent/100,c.heart_rate/60,int(c.pulse_shape=="systolic"))
            self.kernel(*args)
            runs=[]
            for _ in range(3):
                t=perf_counter();self.kernel(*args);runs.append(perf_counter()-t)
            timings[count]=float(np.median(runs))
        # More cores must provide a measurable gain to justify their overhead.
        best=min(timings,key=timings.get)
        if timings[best] > .9*timings[1]:
            best=1
        self.set_threads(best)
        return best

    def fields(self):
        shape = (self.ny, self.nx)
        u = self.ux.reshape(shape) * self.velocity_scale
        v = self.uy.reshape(shape) * self.velocity_scale
        p = (self.rho.reshape(shape) - 1) * self.config.density * self.velocity_scale**2 / 3
        vort = np.gradient(v, self.dx, axis=1) - np.gradient(u, self.dx, axis=0)
        return u, v, p, vort

    def diagnostics(self):
        u, v, p, vort = self.fields()
        g = self.geometry
        mask = g.fluid
        count = mask.sum(axis=0)
        mean_u = (u * mask).sum(axis=0) / count
        mean_p = (p * mask).sum(axis=0) / count
        # Approximate wall-tangent shear using first fluid node and half-cell
        # distance to the staircase bounce-back boundary. Curved-wall WSS needs
        # resolution convergence; this estimate is deliberately labeled.
        low = np.argmax(mask, axis=0)
        high = self.ny - 1 - np.argmax(mask[::-1], axis=0)
        ix = np.arange(self.nx)
        lower_slope = np.gradient(g.lower, self.dx)
        upper_slope = np.gradient(g.upper, self.dx)
        tangential_low = (u[low, ix] + lower_slope * v[low, ix]) / np.sqrt(1 + lower_slope**2)
        tangential_high = (u[high, ix] + upper_slope * v[high, ix]) / np.sqrt(1 + upper_slope**2)
        wss_low = 2 * self.config.viscosity * tangential_low / self.dx
        wss_high = 2 * self.config.viscosity * tangential_high / self.dx
        velocity_lat = np.hypot(self.ux[self.nodes], self.uy[self.nodes])
        backflow = mask & (u < -0.01 * self.config.mean_velocity)
        interior = mask.copy()
        interior[:, :2] = interior[:, -2:] = False
        wss_abs = np.concatenate((np.abs(wss_low), np.abs(wss_high)))
        pulse = float(pulse_curve(self.config, np.asarray([self.time]))[0])
        particle_tau = particle_relaxation_time(
            self.config.particle_density_kg_m3, self.config.particle_diameter_um * 1e-6,
            self.config.viscosity)
        return dict(u=u, v=v, p=p, vorticity=vort, mean_u=mean_u, mean_p=mean_p,
                    wss_low=wss_low, wss_high=wss_high,
                    delta_p=float(mean_p[1] - mean_p[-2]),
                    max_velocity=float(np.max(np.hypot(u[mask], v[mask]))),
                    min_u=float(np.min(u[interior])),
                    backflow_percent=float(100 * backflow[interior].mean()),
                    mach=float(velocity_lat.max() * math.sqrt(3)),
                    density_variation=float(np.max(np.abs(self.rho[self.nodes] - 1))),
                    flux_error=float((mean_u[1] * count[1] - mean_u[-2] * count[-2]) /
                                     max(abs(mean_u[1] * count[1]), 1e-12)),
                    pulse_factor=pulse,
                    phase=phase_label(self.config, self.time),
                    inlet_speed=float(self.config.mean_velocity * pulse),
                    max_wss=float(np.nanmax(wss_abs)),
                    mean_abs_wss=float(np.nanmean(wss_abs)),
                    max_vorticity=float(np.nanmax(np.abs(vort[mask]))),
                    womersley=float(calculate_womersley_number(
                        self.config.density, self.config.heart_rate,
                        self.config.diameter, self.config.viscosity)),
                    particle_relaxation_time_s=float(particle_tau),
                    particle_stokes_number=float(calculate_particle_stokes_number(
                        self.config.particle_density_kg_m3, self.config.particle_diameter_um * 1e-6,
                        self.config.viscosity, self.config.heart_rate)),
                    geometry_report=stenosis_report(self.config, g))

    def check_stability(self):
        r = self.rho[self.nodes]
        speed = np.hypot(self.ux[self.nodes], self.uy[self.nodes])
        if not np.isfinite(r).all() or not np.isfinite(speed).all() or r.min() <= 0 or speed.max() > 0.30:
            raise FloatingPointError("The flow left the stable preview range. Increase cells_across, reduce velocity/severity, and restart.")
