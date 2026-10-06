"""Unrolled D2Q9 stress-regularized collision, identical shear physics to v1.

Serial and parallel dispatchers have distinct function identities so their
Numba disk caches cannot select the other compilation. Regression tests
compare both against the independent loop-based reference.

The nine directions are fixed. Writing their moments explicitly avoids many
small indexed loops in every cell; there is no fast-math or altered viscosity.
"""
import math
from numba import njit, prange
from waveform import pulse_value
from outlet_candidate import characteristic_outlet
from outlet_section import section_density


@njit(cache=True, inline="always")
def extrapolated_boundary(flat, links, adjacent, prescribed_u, is_inlet, shear_factor):
    """Extrapolate the adjacent node's same-time streamed stress to an open end.

    Read only the previous population buffer, never mutable macro fields.
    This stress-regularized non-equilibrium extrapolation is based on
    Guo, Zheng & Shi (2002), DOI 10.1088/1009-1963/11/4/310.
    """
    h0 = flat[links[adjacent,0]]
    h1 = flat[links[adjacent,1]]
    h2 = flat[links[adjacent,2]]
    h3 = flat[links[adjacent,3]]
    h4 = flat[links[adjacent,4]]
    h5 = flat[links[adjacent,5]]
    h6 = flat[links[adjacent,6]]
    h7 = flat[links[adjacent,7]]
    h8 = flat[links[adjacent,8]]
    r = h0+h1+h2+h3+h4+h5+h6+h7+h8
    if r <= 0 or not math.isfinite(r):
        return math.nan, math.nan, math.nan, math.nan, math.nan
    u = (h1-h3+h5-h6-h7+h8)/r
    v = (h2-h4+h5+h6-h7-h8)/r
    normal = shear_factor * (h1+h3-h2-h4-r*(u*u-v*v))
    cross = shear_factor * (h5-h6+h7-h8-r*u*v)
    if is_inlet:
        return r, prescribed_u, 0.0, normal, cross
    return 1.0, u, v, normal, cross


@njit(cache=True, nogil=True)
def advance(f, out, nodes, links, columns, profile, nx, omega,
            rho, ux, uy, steps, iteration, dt, amplitude, frequency, shape_code,
            boundary_code=0):
    shear_factor = .25 * (1 - omega)
    for step in range(steps):
        pulse = pulse_value((iteration + step + 1)*dt, amplitude, frequency, shape_code)
        outlet_r = section_density(f, nx, profile) if boundary_code == 3 else 1.0
        flat = f.ravel()
        for n in range(len(nodes)):
            node = nodes[n]
            g0 = flat[links[n,0]]
            g1 = flat[links[n,1]]
            g2 = flat[links[n,2]]
            g3 = flat[links[n,3]]
            g4 = flat[links[n,4]]
            g5 = flat[links[n,5]]
            g6 = flat[links[n,6]]
            g7 = flat[links[n,7]]
            g8 = flat[links[n,8]]
            column = columns[n]
            if column == 0 and not boundary_code:
                u = profile[node // nx] * pulse
                r = (g0 + g2 + g4 + 2*(g3 + g6 + g7))/(1-u)
                g1 = g3 + 2*r*u/3
                g5 = g7 + (g4-g2)/2 + r*u/6
                g8 = g6 + (g2-g4)/2 + r*u/6
            elif column == nx-1 and not boundary_code:
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
            # Non-equilibrium deviatoric moments, computed without a 9x9 transform.
            normal = shear_factor * (g1+g3-g2-g4-r*(u*u-v*v))
            cross = shear_factor * (g5-g6+g7-g8-r*u*v)
            if boundary_code and (column == 0 or column == nx-1):
                adjacent = n+1 if column == 0 else n-1
                r,u,v,normal,cross = extrapolated_boundary(
                    flat,links,adjacent,profile[node//nx]*pulse,column == 0,shear_factor)
            if boundary_code == 2 and column == nx-1:
                r,u,v,normal,cross = characteristic_outlet(f,node,nx,omega)
            if boundary_code == 3 and column == nx-1:
                r = outlet_r
            rho[node], ux[node], uy[node] = r, u, v
            uu, vv = u*u, v*v
            base = 1 - 1.5*(uu+vv)
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
def advance_parallel(f, out, nodes, links, columns, profile, nx, omega,
            rho, ux, uy, steps, iteration, dt, amplitude, frequency, shape_code,
            boundary_code=0):
    shear_factor = .25 * (1 - omega)
    for step in range(steps):
        pulse = pulse_value((iteration + step + 1)*dt, amplitude, frequency, shape_code)
        outlet_r = section_density(f, nx, profile) if boundary_code == 3 else 1.0
        flat = f.ravel()
        for n in prange(len(nodes)):
            node = nodes[n]
            g0 = flat[links[n,0]]
            g1 = flat[links[n,1]]
            g2 = flat[links[n,2]]
            g3 = flat[links[n,3]]
            g4 = flat[links[n,4]]
            g5 = flat[links[n,5]]
            g6 = flat[links[n,6]]
            g7 = flat[links[n,7]]
            g8 = flat[links[n,8]]
            column = columns[n]
            if column == 0 and not boundary_code:
                u = profile[node // nx] * pulse
                r = (g0 + g2 + g4 + 2*(g3 + g6 + g7))/(1-u)
                g1 = g3 + 2*r*u/3
                g5 = g7 + (g4-g2)/2 + r*u/6
                g8 = g6 + (g2-g4)/2 + r*u/6
            elif column == nx-1 and not boundary_code:
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
            # Non-equilibrium deviatoric moments, computed without a 9x9 transform.
            normal = shear_factor * (g1+g3-g2-g4-r*(u*u-v*v))
            cross = shear_factor * (g5-g6+g7-g8-r*u*v)
            if boundary_code and (column == 0 or column == nx-1):
                adjacent = n+1 if column == 0 else n-1
                r,u,v,normal,cross = extrapolated_boundary(
                    flat,links,adjacent,profile[node//nx]*pulse,column == 0,shear_factor)
            if boundary_code == 2 and column == nx-1:
                r,u,v,normal,cross = characteristic_outlet(f,node,nx,omega)
            if boundary_code == 3 and column == nx-1:
                r = outlet_r
            rho[node], ux[node], uy[node] = r, u, v
            uu, vv = u*u, v*v
            base = 1 - 1.5*(uu+vv)
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


@njit(cache=True, nogil=True)
def mass_inventory(f, previous, nodes):
    """Inventory and last-step change; sum differences before large totals."""
    mass = 0.0
    change = 0.0
    for node in nodes:
        for k in range(9):
            value = f[node, k]
            mass += value
            change += value - previous[node, k]
    return mass, change
