"""Experimental isothermal characteristic outlet, lattice units throughout.

Derivation: characteristic projection of isothermal NS, with transverse and
viscous source terms; cf. Wissocq et al., arXiv:1701.07734, equations 15-30.
This HemoFlow adaptation is not a reproduction of that paper's implementation.
First-order one-sided normal derivatives and explicit Euler boundary update.
Incoming transverse characteristic uses a quiescent reservoir on reverse flow.
Normal velocity is NEVER clipped. No sponge, pressure filtering, or bulk force.
The bulk collision, inlet, raster walls and physical viscosity are unchanged.
"""
import math
from numba import njit

SIGMA = 0.05
RELAX_LENGTH_IN_GAPS = 15.0


@njit(cache=True, inline="always")
def moments(f, node):
    h = f[node]
    r = h[0]+h[1]+h[2]+h[3]+h[4]+h[5]+h[6]+h[7]+h[8]
    if r <= 0 or not math.isfinite(r):
        return math.nan, math.nan, math.nan
    return r, (h[1]-h[3]+h[5]-h[6]-h[7]+h[8])/r, (h[2]-h[4]+h[5]+h[6]-h[7]-h[8])/r


@njit(cache=True, inline="always")
def characteristic_outlet(f, node, nx, omega):
    # All reads are from the OLD population buffer, avoiding parallel races.
    r,u,v = moments(f,node)
    r1,u1,v1 = moments(f,node-1)
    r2,u2,v2 = moments(f,node-2)
    row = node//nx
    ny = len(f)//nx
    rp,up,vp = moments(f,node+nx)
    rm,um,vm = moments(f,node-nx)
    # Half-way no-slip wall ghost values at the two outlet/wall corners.
    if row == 1:
        rm,um,vm = r,-u,-v
    if row == ny-2:
        rp,up,vp = r,-u,-v
    rx,ux,vx = r-r1,u-u1,v-v1
    ry,uy,vy = (rp-rm)*0.5,(up-um)*0.5,(vp-vm)*0.5
    nu = (1/omega-0.5)/3
    # One-sided approximation to normal curvature; centered tangential part.
    fx = nu*(u-2*u1+u2+up-2*u+um)
    fy = nu*(v-2*v1+v2+vp-2*v+vm)
    cs = 1/math.sqrt(3)
    lp = (u+cs)*(rx/3+r*cs*ux)
    tm = -(v*ry/3+r*vy/3-r*cs*v*uy)
    # Weak target-pressure relaxation, not a prescribed constant nodal density.
    # Local Mach is used here; this choice is an explicit adaptation.
    # Hold the reference relaxation length fixed when moving the outlet:
    # 15 nominal gaps = 60 mm for this 4 mm-gap study, on every grid.
    # There are ny-2 fluid cells across the nominal inlet gap.
    relax = SIGMA*(1-(u*u+v*v)*3)*cs/(RELAX_LENGTH_IN_GAPS*(ny-2))
    lm = tm-r*cs*fx+relax*(r-1)/3
    rn = r-1.5*(lp+lm)-v*ry-r*vy
    un = u-(lp-lm)/(2*r*cs)-v*uy+fx
    if u >= 0:
        vn = v-u*vx-v*vy-ry/(3*r)+fy
    else:
        # Re-entering transverse characteristic gets reservoir v=0. Upwind
        # speed -u is positive; no extrapolation against its propagation.
        vn = v+u*v-v*vy-ry/(3*r)+fy
    # Post-collision deviatoric stress from finite-difference strain. The bulk
    # kernel retains only these deviatoric moments, so match its normalization.
    tau = 1/omega
    normal = -(tau-1)*r*(ux-vy)/6
    cross = -(tau-1)*r*(uy+vx)/12
    return rn,un,vn,normal,cross
