"""Experimental section-mean acoustic impedance with regularized extrapolation.

p_out(new) - p_ref = rho_ref * c_s * (U_adjacent(old) - U_ref).
Only the plane/section-mean acoustic mode is matched. This is an explicit
HemoFlow adaptation, not a claim to reproduce a published CBC implementation.
It is not a physiological resistance, a Windkessel, or a general backflow proof.
The node velocities and shear stress come from the ordinary adjacent streamed
state. In particular the first fluid node is not assigned a bulk velocity.
Read only the old population buffer. No clipping, filtering, or bulk forcing.
"""
import math
from numba import njit
from outlet_candidate import moments


@njit(cache=True)
def section_density(f, nx, profile):
    ny = len(f)//nx
    velocity = 0.0
    reference = 0.0
    for row in range(1, ny-1):
        _, u, _ = moments(f, row*nx+nx-2)
        velocity += u
        reference += profile[row]
    # profile is the flux-normalized, time-independent mean inlet profile.
    return 1.0 + math.sqrt(3.0)*(velocity-reference)/(ny-2)
