"""Unrolled D2Q9 stress-regularized collision, identical shear physics to v1.

The nine directions are fixed. Writing their moments explicitly avoids many
small indexed loops in every cell; there is no fast-math or altered viscosity.
"""
import math
from numba import njit, prange
from waveform import pulse_value


@njit(cache=True, nogil=True)
def advance(f, out, nodes, sources, dirs, columns, profile, nx, omega,
            rho, ux, uy, steps, iteration, dt, amplitude, frequency, shape_code):
    shear_factor = .25 * (1 - omega)
    for step in range(steps):
        pulse = pulse_value((iteration + step + 1)*dt, amplitude, frequency, shape_code)
        for n in range(len(nodes)):
            node = nodes[n]
            g0 = f[sources[n,0],dirs[n,0]]
            g1 = f[sources[n,1],dirs[n,1]]
            g2 = f[sources[n,2],dirs[n,2]]
            g3 = f[sources[n,3],dirs[n,3]]
            g4 = f[sources[n,4],dirs[n,4]]
            g5 = f[sources[n,5],dirs[n,5]]
            g6 = f[sources[n,6],dirs[n,6]]
            g7 = f[sources[n,7],dirs[n,7]]
            g8 = f[sources[n,8],dirs[n,8]]
            column = columns[n]
            if column == 0:
                u = profile[node // nx] * pulse
                r = (g0 + g2 + g4 + 2*(g3 + g6 + g7))/(1-u)
                g1 = g3 + 2*r*u/3
                g5 = g7 + (g4-g2)/2 + r*u/6
                g8 = g6 + (g2-g4)/2 + r*u/6
            elif column == nx-1:
                u = -1 + g0 + g2 + g4 + 2*(g1+g5+g8)
                g3 = g1 - 2*u/3
                g6 = g8 + (g4-g2)/2 - u/6
                g7 = g5 + (g2-g4)/2 - u/6
            r = g0+g1+g2+g3+g4+g5+g6+g7+g8
            if r <= 0 or not math.isfinite(r):
                raise FloatingPointError("Non-finite or non-positive fluid density.")
            u = (g1-g3+g5-g6-g7+g8)/r
            v = (g2-g4+g5+g6-g7-g8)/r
            rho[node], ux[node], uy[node] = r, u, v
            uu, vv = u*u, v*v
            base = 1 - 1.5*(uu+vv)
            # Non-equilibrium deviatoric moments, computed without a 9x9 transform.
            normal = shear_factor * (g1+g3-g2-g4-r*(uu-vv))
            cross = shear_factor * (g5-g6+g7-g8-r*u*v)
            ru, rd = r/9, r/36
            out[node,0] = 4*r/9*base
            out[node,1] = ru*(base+3*u+4.5*uu)+normal
            out[node,2] = ru*(base+3*v+4.5*vv)-normal
            out[node,3] = ru*(base-3*u+4.5*uu)+normal
            out[node,4] = ru*(base-3*v+4.5*vv)-normal
            a, b = u+v, -u+v
            out[node,5] = rd*(base+3*a+4.5*a*a)+cross
            out[node,6] = rd*(base+3*b+4.5*b*b)-cross
            out[node,7] = rd*(base-3*a+4.5*a*a)+cross
            out[node,8] = rd*(base-3*b+4.5*b*b)-cross
        f, out = out, f
    return f, out


@njit(cache=True, nogil=True, parallel=True)
def advance_parallel(f, out, nodes, sources, dirs, columns, profile, nx, omega,
            rho, ux, uy, steps, iteration, dt, amplitude, frequency, shape_code):
    shear_factor = .25 * (1 - omega)
    for step in range(steps):
        pulse = pulse_value((iteration + step + 1)*dt, amplitude, frequency, shape_code)
        for n in prange(len(nodes)):
            node = nodes[n]
            g0 = f[sources[n,0],dirs[n,0]]
            g1 = f[sources[n,1],dirs[n,1]]
            g2 = f[sources[n,2],dirs[n,2]]
            g3 = f[sources[n,3],dirs[n,3]]
            g4 = f[sources[n,4],dirs[n,4]]
            g5 = f[sources[n,5],dirs[n,5]]
            g6 = f[sources[n,6],dirs[n,6]]
            g7 = f[sources[n,7],dirs[n,7]]
            g8 = f[sources[n,8],dirs[n,8]]
            column = columns[n]
            if column == 0:
                u = profile[node // nx] * pulse
                r = (g0 + g2 + g4 + 2*(g3 + g6 + g7))/(1-u)
                g1 = g3 + 2*r*u/3
                g5 = g7 + (g4-g2)/2 + r*u/6
                g8 = g6 + (g2-g4)/2 + r*u/6
            elif column == nx-1:
                u = -1 + g0 + g2 + g4 + 2*(g1+g5+g8)
                g3 = g1 - 2*u/3
                g6 = g8 + (g4-g2)/2 - u/6
                g7 = g5 + (g2-g4)/2 - u/6
            r = g0+g1+g2+g3+g4+g5+g6+g7+g8
            if r <= 0 or not math.isfinite(r):
                rho[node] = math.nan
                ux[node] = math.nan
                uy[node] = math.nan
                for k in range(9):
                    out[node,k] = math.nan
                continue
            u = (g1-g3+g5-g6-g7+g8)/r
            v = (g2-g4+g5+g6-g7-g8)/r
            rho[node], ux[node], uy[node] = r, u, v
            uu, vv = u*u, v*v
            base = 1 - 1.5*(uu+vv)
            # Non-equilibrium deviatoric moments, computed without a 9x9 transform.
            normal = shear_factor * (g1+g3-g2-g4-r*(uu-vv))
            cross = shear_factor * (g5-g6+g7-g8-r*u*v)
            ru, rd = r/9, r/36
            out[node,0] = 4*r/9*base
            out[node,1] = ru*(base+3*u+4.5*uu)+normal
            out[node,2] = ru*(base+3*v+4.5*vv)-normal
            out[node,3] = ru*(base-3*u+4.5*uu)+normal
            out[node,4] = ru*(base-3*v+4.5*vv)-normal
            a, b = u+v, -u+v
            out[node,5] = rd*(base+3*a+4.5*a*a)+cross
            out[node,6] = rd*(base+3*b+4.5*b*b)-cross
            out[node,7] = rd*(base-3*a+4.5*a*a)+cross
            out[node,8] = rd*(base-3*b+4.5*b*b)-cross
        f, out = out, f
    return f, out
