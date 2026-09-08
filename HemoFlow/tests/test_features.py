"""Regression checks for the continuous viewer features."""
from dataclasses import replace
from pathlib import Path
import sys
import csv
import tempfile
import time
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import Config
from engine import FlowWorker, interpolate_frames
from geometry import make_geometry, stenosis_report
from particles import Tracers
from recorder import RunRecorder
from solver import FlowSolver
from analyze_recording import summarize
from waveform import pulse_curve


class GeometryAndViewerFeatures(unittest.TestCase):
    def test_default_random_count_is_small_but_explicit_three_is_available(self):
        counts = []
        layouts = []
        for seed in range(12):
            c = Config(seed=seed, cells_across=32)
            a = make_geometry(c)
            b = make_geometry(c)
            self.assertEqual(a.plaques, b.plaques)
            self.assertLessEqual(len(a.plaques), 2)
            counts.append(len(a.plaques))
            layouts.append(tuple(round(p["center_mm"], 3) for p in a.plaques))
        self.assertGreater(len(set(counts)), 1)
        self.assertGreater(len(set(layouts)), 1)
        three = [len(make_geometry(Config(seed=seed, cells_across=32, max_random_plaques=3)).plaques)
                 for seed in range(24)]
        self.assertIn(3, three)

    def test_custom_plaque_edit_changes_wall_geometry(self):
        base = Config(seed=42, cells_across=32)
        original = make_geometry(base)
        plaque = dict(original.plaques[0])
        plaque["center_mm"] = min(base.length_mm - 8, plaque["center_mm"] + 4)
        plaque["wall"] = "top" if plaque["wall"] == "bottom" else "bottom"
        edited = replace(base, geometry_mode="custom", plaques=[plaque], plaque_count=1,
                         randomize_count=False)
        changed = make_geometry(edited)
        self.assertEqual(changed.plaques[0]["wall"], plaque["wall"])
        self.assertFalse(np.array_equal(original.solid, changed.solid))

    def test_systolic_waveform_reaches_visible_extrema(self):
        c = Config(pulse_shape="systolic")
        period = 60 / c.heart_rate
        t = np.linspace(0, period, 2001)
        speed = c.mean_velocity * pulse_curve(c, t)
        self.assertAlmostEqual(float(speed.min()), .18, places=3)
        self.assertAlmostEqual(float(speed.max()), .42, places=3)
        self.assertLess(float(t[speed.argmax()]), .25 * period)

    def test_worker_publishes_bounded_interpolatable_frames(self):
        c = Config(seed=7, cells_across=24, particle_count=100, compute_threads=1)
        worker = FlowWorker(c, publish_interval=.02)
        worker.start()
        try:
            deadline = time.time() + 6
            while time.time() < deadline:
                frames = worker.status()[2]
                if len(frames) >= 2:
                    break
                time.sleep(.01)
            frames = worker.status()[2]
            self.assertGreaterEqual(len(frames), 2)
            self.assertLessEqual(len(frames), 8)
            latest = frames[-1]
            self.assertGreater(latest.iteration, 0)
            self.assertGreater(latest.dt, 0)
            item = interpolate_frames(frames, time.perf_counter(), delay=.005)
            self.assertIsNotNone(item)
            self.assertEqual(item[3].shape, latest.positions.shape)
        finally:
            worker.stop()
            worker.join(3)
        self.assertEqual(worker.status()[0], "stopped")

    def test_particle_visibility_toggle_keeps_interpolation_safe(self):
        c = Config(seed=9, cells_across=24, particle_count=100, compute_threads=1)
        worker = FlowWorker(c, publish_interval=.01)
        worker.start()
        try:
            deadline = time.time() + 6
            while time.time() < deadline and len(worker.status()[2]) < 2:
                time.sleep(.01)
            worker.set_particles_enabled(False)
            time.sleep(.08)
            frames = worker.status()[2]
            self.assertFalse(frames[-1].particles_enabled)
            self.assertEqual(frames[-1].positions.shape, (0, 2))
            item = interpolate_frames(frames, time.perf_counter(), delay=0)
            self.assertEqual(item[3].shape, (0, 2))
            worker.set_particles_enabled(True)
            time.sleep(.08)
            frames = worker.status()[2]
            self.assertTrue(frames[-1].particles_enabled)
            self.assertEqual(frames[-1].positions.shape, (100, 2))
        finally:
            worker.stop()
            worker.join(3)

    def test_stenosis_report_separates_diameter_and_area_equivalents(self):
        c = Config(geometry_mode="single", stenosis_percent=50, stenosis_center_mm=30,
                   stenosis_length_mm=12, cells_across=64, particle_count=100)
        g = make_geometry(c)
        report = stenosis_report(c, g)
        item = report["plaques"][0]
        self.assertAlmostEqual(item["requested_diameter_narrowing_percent"], 50.0)
        d = item["realized_raster_diameter_narrowing_percent"] / 100.0
        expected_area = 100 * (1 - (1 - d) ** 2)
        self.assertAlmostEqual(item["equivalent_circular_3d_area_reduction_percent"], expected_area)
        self.assertGreater(item["equivalent_circular_3d_area_reduction_percent"],
                           item["realized_raster_diameter_narrowing_percent"])

    def test_old_settings_load_with_v3_defaults(self):
        import json
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "v2_settings.json"
            path.write_text(json.dumps({"density": 1060.0, "viscosity": .0035,
                                        "diameter_mm": 4.0, "length_mm": 60.0,
                                        "mean_velocity": .3, "stenosis_percent": 50.0,
                                        "stenosis_center_mm": 30.0, "stenosis_length_mm": 15.0,
                                        "heart_rate": 72.0, "pulsatility_percent": 40.0,
                                        "plaque_count": 3, "seed": 42,
                                        "geometry_mode": "random", "cells_across": 40,
                                        "particle_count": 1000}), encoding="utf-8")
            loaded = Config.load(path)
            self.assertEqual(loaded.particle_model, "tracer")
            self.assertEqual(loaded.particle_display_limit, 650)

    def test_inertial_particle_properties_have_finite_relaxation(self):
        c = Config(cells_across=24, particle_count=120, particle_model="inertial",
                   particle_diameter_um=100, particle_density_kg_m3=1200,
                   compute_threads=1)
        s = FlowSolver(c)
        s.step(24)
        p = Tracers(s)
        p.advance(s, 24)
        desc = p.description(s)
        self.assertEqual(desc["model"], "inertial")
        self.assertGreater(desc["relaxation_time_s"], 0)
        self.assertTrue(np.isfinite(p.x).all() and np.isfinite(p.y).all())

    def test_recorder_is_append_only_and_rate_limited(self):
        c = Config(geometry_mode="healthy", length_mm=16, stenosis_center_mm=8,
                   stenosis_length_mm=4, cells_across=24, compute_threads=1)
        s = FlowSolver(c)
        with tempfile.TemporaryDirectory() as folder:
            rec = RunRecorder(folder, every_s=s.dt * 10)
            for _ in range(30):
                s.step(2)
                rec.record(s, s.diagnostics(), 0, None)
            rec.close()
            csv_path = next(Path(folder).glob("*_timeseries.csv"))
            with csv_path.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertGreaterEqual(len(rows), 2)
            self.assertLess(len(rows), 30)
            self.assertTrue(next(Path(folder).glob("*_metadata.json")).exists())
            self.assertTrue(next(Path(folder).glob("*_epochs.jsonl")).exists())
            summary = summarize(csv_path)
            self.assertEqual(summary["rows"], len(rows))
            self.assertTrue(summary["cycles"])


if __name__ == "__main__":
    unittest.main()
