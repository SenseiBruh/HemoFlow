"""Wall-shear measurement diagnostics; no CFD solver or field smoothing.

SI units throughout. Tangent points toward increasing x; normal points INTO
the fluid on each wall. Signed shear = mu * t @ (grad(U)+grad(U).T) @ n.
The local polynomial is an unvalidated diagnostic on a smooth target wall,
not a replacement for the frozen raster-wall solver's stress measurement.
"""
from dataclasses import dataclass
import numpy as np

VERSION = "wall-audit-1.0"


def frame(slope, wall):
    if wall not in ("lower", "upper"):
        raise ValueError("wall must be lower or upper")
    t = np.array([1.0, slope]) / np.hypot(1, slope)
    n = np.array([-slope, 1.0]) / np.hypot(1, slope)
    return t, n if wall == "lower" else -n


def projected_shear(gradient, slope, wall, mu):
    t, n = frame(slope, wall)
    return float(mu * (t @ (gradient + gradient.T) @ n))


@dataclass
class Stencil:
    rows: np.ndarray
    cols: np.ndarray
    inverse: np.ndarray
    dx: float
    slope: float
    wall: str
    condition: float
    reason: str

    def evaluate(self, u, v, mu):
        result = dict(valid=False, reason=self.reason, nodes=len(self.rows),
                      condition=self.condition, shear_pa=None,
                      wall_tangent_velocity_m_s=None, wall_normal_velocity_m_s=None)
        if self.reason:
            return result
        samples = np.column_stack((u[self.rows, self.cols], v[self.rows, self.cols]))
        if not np.isfinite(samples).all():
            result["reason"] = "nonfinite_fluid_sample"
            return result
        coeff = self.inverse @ samples
        gradient = coeff[1:3].T / self.dx
        t, n = frame(self.slope, self.wall)
        result.update(valid=True, shear_pa=projected_shear(gradient, self.slope, self.wall, mu),
                      wall_tangent_velocity_m_s=float(coeff[0] @ t),
                      wall_normal_velocity_m_s=float(coeff[0] @ n))
        return result


def build_stencil(x, y, solid, wall_x, wall_y, slope, wall, radius_cells=4.0):
    """Unweighted quadratic LS fit using only fluid nodes in a disk.

    Basis [1, X, Y, X^2, XY, Y^2], X/Y scaled by dx. Intercept is FREE:
    no artificial no-slip datum is inserted, allowing a wall residual check.
    No solid-node zeros enter the fit. Singular/undersized stencils fail.
    """
    x, y, solid = np.asarray(x), np.asarray(y), np.asarray(solid, dtype=bool)
    if len(x) < 2 or len(y) < 2 or solid.shape != (len(y), len(x)):
        raise ValueError("Invalid grid/mask shape")
    dx = float(x[1] - x[0])
    if dx <= 0 or not np.allclose(np.diff(x), dx, rtol=1e-8, atol=1e-15) or not np.allclose(np.diff(y), dx, rtol=1e-8, atol=1e-15):
        raise ValueError("A uniform square grid is required")
    if radius_cells <= 0:
        raise ValueError("radius_cells must be positive")
    xx = np.flatnonzero(np.abs(x-wall_x) <= radius_cells*dx*(1+1e-10))
    yy = np.flatnonzero(np.abs(y-wall_y) <= radius_cells*dx*(1+1e-10))
    rr, cc = np.meshgrid(yy, xx, indexing="ij")
    a, b = (x[cc]-wall_x)/dx, (y[rr]-wall_y)/dx
    keep = (~solid[rr, cc]) & (a*a+b*b <= radius_cells**2+1e-9)
    rr, cc, a, b = rr[keep], cc[keep], a[keep], b[keep]
    mat = np.column_stack((np.ones(len(a)), a, b, a*a, a*b, b*b))
    reason, condition = "", float("inf")
    inverse = np.empty((6, len(a)))
    if len(a) < 12:
        reason = "fewer_than_12_fluid_nodes"
    else:
        singular = np.linalg.svd(mat, compute_uv=False)
        condition = float(singular[0]/singular[-1]) if singular[-1] > 0 else float("inf")
        if not np.isfinite(condition) or condition > 1e6:
            reason = "ill_conditioned_stencil"
        else:
            inverse = np.linalg.pinv(mat, rcond=1e-12)
    return Stencil(rr, cc, inverse, dx, float(slope), wall, condition, reason)


def legacy_profile(u, v, solid, lower, upper, dx, mu):
    """Reproduce frozen solver.py diagnostics; independently check saved WSS."""
    fluid = ~solid
    if np.any(fluid.sum(axis=0) < 1):
        raise ValueError("A column has no fluid")
    low = np.argmax(fluid, axis=0)
    high = len(fluid)-1-np.argmax(fluid[::-1], axis=0)
    ix = np.arange(fluid.shape[1])
    out = {}
    for wall, surface, row in (("lower", lower, low), ("upper", upper, high)):
        slope = np.gradient(surface, dx)
        out[wall] = 2*mu*(u[row, ix]+slope*v[row, ix])/np.hypot(1, slope)/dx
    return out


def difference(a, b):
    """Direct B-vs-A diagnostics; never mask invalid values or swap walls."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if a.shape != b.shape or a.size == 0 or not np.isfinite(a).all() or not np.isfinite(b).all():
        return dict(valid=False, reason="missing_or_nonfinite_values", count=int(a.size))
    rms = float(np.sqrt(np.mean((b-a)**2)))
    ref = float(np.sqrt(np.mean(a*a)))
    return dict(valid=True, count=int(a.size), rms_difference_pa=rms,
                reference_rms_pa=ref, relative_l2_percent=100*rms/ref if ref > 1e-12 else None,
                mean_abs_a_pa=float(np.mean(np.abs(a))), mean_abs_b_pa=float(np.mean(np.abs(b))),
                peak_abs_a_pa=float(np.max(np.abs(a))), peak_abs_b_pa=float(np.max(np.abs(b))))


def manufactured_checks():
    """Exact polynomial and curved divergence-free field measurement tests."""
    rows = []
    mu = .0035
    for radius in (4.0, 5.0):
        for wall in ("lower", "upper"):
            sign = 1 if wall == "lower" else -1
            for slope in (0., -.6, .6, 1.2):
                for offset in (.1, .5, .9):
                    for gamma in (-500., 500.):
                        dx = .004/80
                        x = np.arange(-8, 9)*dx
                        y = (np.arange(-10, 11)+offset)*dx
                        X, Y = np.meshgrid(x, y)
                        distance = sign*(Y-slope*X)/np.hypot(1,slope)
                        solid = distance <= 0
                        t, _ = frame(slope, wall)
                        u, v = gamma*distance*t[0], gamma*distance*t[1]
                        got = build_stencil(x,y,solid,0,0,slope,wall,radius).evaluate(u,v,mu)
                        expected = mu*gamma
                        error = abs(got["shear_pa"]-expected) if got["valid"] else float("inf")
                        rows.append(dict(test="inclined_linear_shear", radius_cells=radius, wall=wall,
                            slope=slope, offset_cells=offset, cells=80, expected_pa=expected,
                            measured_pa=got["shear_pa"], absolute_error_pa=error, limit_pa=1e-9,
                            passed=error < 1e-9))
            for cells in (40,80,120,160):
                H, U = .002, .3
                dx = 2*H/cells
                x = np.arange(-8,9)*dx
                y = (np.arange(cells+2)-.5)*dx-H
                X,Y = np.meshgrid(x,y)
                solid = (Y<=-H)|(Y>=H)
                u, v = 1.5*U*(1-(Y/H)**2), np.zeros_like(Y)
                point_y = -H if wall == "lower" else H
                got = build_stencil(x,y,solid,0,point_y,0,wall,radius).evaluate(u,v,mu)
                exact = 3*mu*U/H
                error = abs(got["shear_pa"]-exact) if got["valid"] else float("inf")
                legacy = legacy_profile(u,v,solid,np.full(len(x),-H),np.full(len(x),H),dx,mu)[wall][8]
                rows.append(dict(test="planar_poiseuille", radius_cells=radius, wall=wall,
                    cells=cells, expected_pa=exact, measured_pa=got["shear_pa"],
                    absolute_error_pa=error, limit_pa=1e-9, passed=error<1e-9,
                    legacy_pa=legacy, legacy_expected_pa=exact*(1-1/(2*cells)),
                    legacy_relative_error_percent=100*(legacy-exact)/exact))
            # Streamfunction psi = sign*A/2*(y-f(x))^2 gives an exactly
            # divergence-free no-slip field on y=f(x), with curved f.
            for cells in (40,80,120,160):
                dx, A, slope, curvature = .004/cells, 500., .6, 1000.
                x = np.arange(-8,9)*dx
                y = (np.arange(-10,11)+.35)*dx
                X,Y = np.meshgrid(x,y)
                f = slope*X+.5*curvature*X*X
                d = Y-f
                u, v = sign*A*d, sign*A*d*(slope+curvature*X)
                solid = sign*d<=0
                got = build_stencil(x,y,solid,0,0,slope,wall,radius).evaluate(u,v,mu)
                exact = mu*A*(1+slope*slope)
                error = abs(got["shear_pa"]-exact) if got["valid"] else float("inf")
                # Predeclared diagnostic tolerance: 5% at N40, decreasing as dx^2.
                limit = .05*exact*(40/cells)**2
                rows.append(dict(test="curved_streamfunction", radius_cells=radius, wall=wall,
                    cells=cells, expected_pa=exact, measured_pa=got["shear_pa"],
                    absolute_error_pa=error, limit_pa=limit, passed=error<limit))
    return rows
