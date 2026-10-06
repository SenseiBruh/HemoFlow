"""Physical moments and concurrency checks for the open-end reconstruction."""
from dataclasses import replace
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import Config
from solver import FlowSolver, CX, CY, W


class OpenBoundaryChecks(unittest.TestCase):
    def test_boundary_moments_follow_adjacent_streamed_stress(self):
        s = FlowSolver(Config(geometry_mode='healthy', cells_across=24,
                              boundary_model='regularized'))
        before = s.f.copy()
        s.step()
        for x in (0, s.nx-1):
            node = (s.ny//2)*s.nx + x
            n = int(np.flatnonzero(s.nodes == node)[0])
            adjacent = n+1 if x == 0 else n-1
            streamed = before.ravel()[s.links[adjacent]]
            rn = streamed.sum()
            un, vn = streamed @ CX / rn, streamed @ CY / rn
            rb = rn if x == 0 else 1.0
            ub = s.profile[node//s.nx] * s.diagnostics()['pulse_factor'] if x == 0 else un
            vb = 0.0 if x == 0 else vn
            f = s.f[node]
            self.assertAlmostEqual(float(f.sum()), rb, places=13)
            self.assertAlmostEqual(float(f@CX/rb), ub, places=13)
            self.assertAlmostEqual(float(f@CY/rb), vb, places=13)
            cu = CX*un + CY*vn
            eq = W*rn*(1+3*cu+4.5*cu**2-1.5*(un**2+vn**2))
            n_normal = (streamed-eq) @ (CX*CX-CY*CY)
            n_cross = (streamed-eq) @ (CX*CY)
            b_normal = f @ (CX*CX-CY*CY) - rb*(ub*ub-vb*vb)
            b_cross = f @ (CX*CY) - rb*ub*vb
            self.assertAlmostEqual(float(b_normal), (1-s.omega_p)*n_normal, places=13)
            self.assertAlmostEqual(float(b_cross), (1-s.omega_p)*n_cross, places=13)

    def test_regularized_parallel_matches_serial_with_both_geometries(self):
        for mode in ('healthy', 'single'):
            c = Config(geometry_mode=mode, cells_across=24, boundary_model='regularized')
            a, b = FlowSolver(c), FlowSolver(replace(c, compute_threads=2))
            a.set_threads(1)
            b.set_threads(2)
            a.step(501)
            b.step(501)
            np.testing.assert_array_equal(a.f, b.f)
            self.assertLess(abs(a.mass_balance()['mass_balance_relative']), 1e-9)

    def test_recorded_extrema_match_field_coordinates(self):
        s = FlowSolver(Config(geometry_mode='healthy', cells_across=24))
        s.step(51)
        d = s.diagnostics()
        g = s.geometry
        x = int(np.argmin(np.abs(1000*g.x-d['max_wss_x_mm'])))
        wss = d['wss_low'] if d['max_wss_wall_index'] == 0 else d['wss_high']
        self.assertEqual(abs(float(wss[x])), d['max_wss'])
        x = int(np.argmin(np.abs(1000*g.x-d['min_u_x_mm'])))
        y = int(np.argmin(np.abs(1000*g.y-d['min_u_y_mm'])))
        self.assertEqual(float(d['u'][y,x]), d['min_u'])


if __name__ == '__main__':
    unittest.main()
