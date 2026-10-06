"""Independent conservation, indexing, and physical-scaling regressions."""
from dataclasses import replace
from pathlib import Path
import sys
import tempfile
import unittest
import csv
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from config import Config
from geometry import make_geometry,stenosis_report
from solver import FlowSolver, _advance
from recorder import RunRecorder
from analyze_recording import summarize


class NumericalCorrections(unittest.TestCase):
    def test_time_scale_keeps_physical_viscosity_velocity_and_heartbeat(self):
        c=Config(geometry_mode="healthy",time_scale=1)
        a,b=FlowSolver(c),FlowSolver(replace(c,time_scale=.5))
        self.assertAlmostEqual(b.dt,a.dt/2)
        for s in (a,b):
            nu=(s.tau-.5)/3*s.dx**2/s.dt
            self.assertAlmostEqual(nu,c.viscosity/c.density,places=14)
            self.assertAlmostEqual(s.profile[s.geometry.fluid[:,0]].mean()*s.velocity_scale,c.mean_velocity)
            self.assertEqual(s.config.heart_rate,c.heart_rate)

    def test_population_indexing_and_parallel_collision_match_reference(self):
        c=Config(cells_across=24,pulse_shape="sine",time_scale=1,compute_threads=1,boundary_model="zou_he")
        a,b=FlowSolver(c),FlowSolver(c)
        a.set_threads(1);b.set_threads(2)
        f,other=a.f.copy(),a.other.copy()
        r,u,v=a.rho.copy(),a.ux.copy(),a.uy.copy()
        reference,_=_advance(f,other,a.nodes,a.links//9,a.links%9,a.columns,a.profile,
                             a.nx,a.omega_p,r,u,v,51,0,a.dt,
                             c.pulsatility_percent/100,c.heart_rate/60)
        a.step(51);b.step(51)
        np.testing.assert_allclose(a.f[a.nodes],reference[a.nodes],rtol=2e-12,atol=2e-14)
        np.testing.assert_array_equal(a.f,b.f)

    def test_independent_mass_balance_and_injected_defect(self):
        for mode in ("healthy","single"):
            s=FlowSolver(Config(geometry_mode=mode,cells_across=24,time_scale=.5))
            s.step(211)
            d=s.mass_balance()
            self.assertTrue(d["mass_balance_valid"])
            self.assertLess(abs(d["mass_balance_relative"]),1e-9)
            # Changing an interior population creates a real conservation
            # defect. A tautological residual could not detect this.
            node=s.mass_nodes[len(s.mass_nodes)//2]
            s.f[node,0]+=1e-5
            self.assertGreater(abs(s.mass_balance()["mass_balance_relative"]),1e-6)

    def test_overall_raster_report_uses_actual_nodes(self):
        c=Config(geometry_mode="single",stenosis_percent=37,cells_across=40)
        g=make_geometry(c);r=stenosis_report(c,g)
        expected=100*(1-g.min_gap_cells/c.cells_across)
        self.assertAlmostEqual(r["overall_realized_raster_diameter_narrowing_percent"],expected)
        self.assertAlmostEqual(r["overall_realized_continuous_diameter_narrowing_percent"],37)
        self.assertNotAlmostEqual(expected,37)
        self.assertAlmostEqual(r["overall_equivalent_circular_3d_area_reduction_percent"],100*(1-(1-expected/100)**2))

    def test_sample_schedule_does_not_accumulate_rounding_drift(self):
        s=FlowSolver(Config(geometry_mode="healthy",cells_across=24,time_scale=.5))
        interval=s.dt*9.2
        with tempfile.TemporaryDirectory() as directory:
            rec=RunRecorder(directory,every_s=interval)
            rec.record(s,s.diagnostics(),0)
            for _ in range(105):
                s.step(1)
                if rec.due(0,s.time):rec.record(s,s.diagnostics(),0)
            rec.close()
            with rec.csv_path.open() as file:rows=list(csv.DictReader(file))
            for i,row in enumerate(rows):
                self.assertLessEqual(abs(float(row['time_s'])-i*interval),s.dt+1e-12)
            self.assertEqual(float(rows[-1]['time_scale']),.5)
            self.assertIn('mass_storage_rate_kg_m_s',rows[-1])

    def test_cycle_summary_weights_irregular_times_and_marks_partial_cycles(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"samples.csv"
            with path.open("w",newline="") as f:
                writer=csv.DictWriter(f,fieldnames=["epoch","time_s","cycle_index","phase_fraction","delta_p_pa"])
                writer.writeheader()
                for t in [0,.1,.9,1,1.1,1.9,2.1,2.2]:
                    writer.writerow(dict(epoch=0,time_s=t,cycle_index=int(t),phase_fraction=t%1,delta_p_pa=t))
            result=summarize(path)
            self.assertEqual(len(result["cycles"]),3)
            self.assertTrue(result["cycles"][1]["complete_cycle"])
            self.assertFalse(result["cycles"][2]["complete_cycle"])
            self.assertAlmostEqual(result["cycles"][0]["delta_p_pa_mean"],.5)
            self.assertAlmostEqual(result["cycles"][0]["delta_p_pa_std"],1/(12**.5))
            self.assertEqual(len(summarize(path,discard_cycles=1,complete_only=True)["cycles"]),1)


if __name__=="__main__":unittest.main()
