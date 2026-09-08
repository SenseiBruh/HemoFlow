"""Physical and behavioral regressions; run with the standard-library unittest."""
from dataclasses import replace
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import Config
from geometry import make_geometry
from particles import Tracers
from physics import calculate_pressure_drop, calculate_wall_shear_stress
from solver import FlowSolver


class PhysicalChecks(unittest.TestCase):
    def test_original_pipe_reference_values(self):
        c = Config()
        self.assertAlmostEqual(c.reynolds, 363.4285714285714)
        self.assertAlmostEqual(calculate_pressure_drop(c.viscosity,c.length,c.mean_velocity,c.diameter), 126.0)
        self.assertAlmostEqual(calculate_wall_shear_stress(c.viscosity,c.mean_velocity,c.diameter), 2.1)

    def test_healthy_channel_matches_planar_solution(self):
        c = Config(geometry_mode="healthy", length_mm=16, stenosis_center_mm=8,
                   stenosis_length_mm=4, mean_velocity=.1, pulsatility_percent=0, cells_across=32)
        s = FlowSolver(c)
        s.step(10000)
        s.check_stability()
        d = s.diagnostics()
        # Match the pressure measurement stations, not the complete domain length.
        measured_length = s.geometry.x[-2] - s.geometry.x[1]
        dp_expected = 12*c.viscosity*measured_length*c.mean_velocity/c.diameter**2
        wss_expected = 6*c.viscosity*c.mean_velocity/c.diameter
        y = s.geometry.y[1:-1]
        expected = 1.5*c.mean_velocity*(1-(2*y/c.diameter)**2)
        profile = d["u"][1:-1, s.nx//2]
        error = np.linalg.norm(profile-expected)/np.linalg.norm(expected)
        self.assertLess(error, .03)
        self.assertLess(abs(d["delta_p"]/dp_expected-1), .05)
        self.assertLess(abs(np.mean(d["wss_low"][3:-3])/wss_expected-1), .05)
        self.assertLess(abs(d["flux_error"]), .01)
        self.assertEqual(d["backflow_percent"], 0)

    def test_seeded_plaque_case_has_computed_backflow(self):
        c = Config(pulsatility_percent=0, cells_across=32)
        s = FlowSolver(c)
        s.step(7000)
        s.check_stability()
        d = s.diagnostics()
        # Inspect downstream of a plaque apex, in fluid cells clear of its surface.
        wake = np.zeros_like(s.geometry.fluid)
        for p in s.geometry.plaques:
            start = p["center_mm"]/1000
            wake |= ((s.geometry.x > start) & (s.geometry.x < start + 3*c.diameter))[None,:]
        wake &= s.geometry.fluid
        wake[:2] = wake[-2:] = False
        self.assertLess(np.min(d["u"][wake]), -.02*c.mean_velocity)
        self.assertGreater(d["backflow_percent"], 1.0)
        self.assertLess(d["mach"], .25)

    def test_seed_reproducibility_and_connected_throats(self):
        c = Config()
        a, b = make_geometry(c), make_geometry(c)
        self.assertTrue(np.array_equal(a.solid,b.solid))
        self.assertEqual(a.plaques,b.plaques)
        self.assertFalse(np.array_equal(a.solid,make_geometry(replace(c,seed=43)).solid))
        for seed in [0,1,7,42,99]:
            g = make_geometry(replace(c,seed=seed,plaque_count=5))
            self.assertGreaterEqual(g.min_gap_cells,6)
            # Every fluid column is a single open interval, with flat open ends.
            transitions = np.count_nonzero(np.diff(g.solid.astype(int),axis=0),axis=0)
            self.assertTrue(np.all(transitions == 2))
            self.assertTrue(np.all(np.any(g.fluid[:,:-1] & g.fluid[:,1:],axis=0)))
            self.assertTrue(np.array_equal(g.solid[:,0],g.solid[:,-1]))
        with self.assertRaises(ValueError):
            make_geometry(replace(c,geometry_mode="single",stenosis_percent=95))

    def test_tracers_stay_in_fluid_and_trails_stay_bounded(self):
        c = Config(cells_across=32,particle_count=400)
        s = FlowSolver(c)
        s.step(4000)
        p = Tracers(s)
        shape = p.history.shape
        for _ in range(100):
            s.step(16)
            s.check_stability()
            p.advance(s,16)
            p.record(s)
        iy,ix = np.rint(p.y).astype(int),np.rint(p.x).astype(int)
        self.assertTrue(np.isfinite(p.x).all() and np.isfinite(p.y).all())
        self.assertFalse(s.geometry.solid[iy,ix].any())
        self.assertEqual(p.history.shape,shape)
        self.assertEqual(p.trails().shape,(400,32,2))
        # Deliberately place tracers at a valid outlet location to test recycling.
        p.x[:10] = s.nx-1.001
        p.y[:10] = s.ny/2
        p.advance(s,16)
        self.assertGreaterEqual(p.recycled,10)
        self.assertTrue(np.all(p.x[:10] < 2))


if __name__ == "__main__":
    unittest.main()
