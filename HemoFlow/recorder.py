"""Crash-safe, bounded-rate scalar recording for unattended HemoFlow runs.

The recorder deliberately stores compact observables at a simulated-time
interval instead of retaining every field in memory.  Full fields remain
available from the existing snapshot export.  A geometry edit creates a new
epoch entry while the scalar stream stays easy to concatenate and group.
"""
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import json
from pathlib import Path
from threading import Lock
import numpy as np

from geometry import stenosis_report
from physics import calculate_womersley_number
from waveform import phase_label, pulse_curve


SCHEMA_VERSION = 1


class RunRecorder:
    """Append-only time-series recorder safe to attach to ``FlowWorker``."""

    def __init__(self, folder="exports/runs", every_s=0.05, run_name=None):
        every_s = float(every_s)
        if not np.isfinite(every_s) or every_s <= 0:
            raise ValueError("recording interval must be a finite positive number")
        self.folder = Path(folder)
        self.folder.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
        self.run_id = run_name or f"hemoflow_run_{stamp}"
        self.every_s = every_s
        self.csv_path = self.folder / f"{self.run_id}_timeseries.csv"
        self.metadata_path = self.folder / f"{self.run_id}_metadata.json"
        self.epochs_path = self.folder / f"{self.run_id}_epochs.jsonl"
        self._file = None
        self._writer = None
        self._last_time = {}
        self._seen_epochs = set()
        self._lock = Lock()
        self._closed = False

    @property
    def paths(self):
        return {
            "timeseries_csv": str(self.csv_path),
            "metadata_json": str(self.metadata_path),
            "epochs_jsonl": str(self.epochs_path),
        }

    def _ensure_open(self, config, geometry, epoch):
        if self._file is None:
            self._file = self.csv_path.open("w", newline="", encoding="utf-8")
            self._writer = csv.DictWriter(self._file, fieldnames=[
                "epoch", "time_s", "cycle_index", "phase_fraction", "phase",
                "pulse_factor", "inlet_speed_m_s", "mean_inlet_velocity_m_s",
                "mean_outlet_velocity_m_s", "inlet_flow_per_depth_m2_s",
                "outlet_flow_per_depth_m2_s", "delta_p_pa", "max_velocity_m_s",
                "min_axial_velocity_m_s", "backflow_percent", "max_wss_pa",
                "mean_abs_wss_pa", "max_vorticity_s_1", "flux_error",
                "mach", "density_variation", "particle_count", "particles_recycled",
                "particle_relaxation_time_s", "particle_stokes_number",
            ])
            self._writer.writeheader()
            self._file.flush()
        if not self.metadata_path.exists():
            c = asdict(config)
            metadata = {
                "schema_version": SCHEMA_VERSION,
                "run_id": self.run_id,
                "created_utc": datetime.now(timezone.utc).isoformat(),
                "record_every_s": self.every_s,
                "units": {
                    "time": "s", "velocity": "m/s", "pressure": "Pa",
                    "wall_shear": "Pa", "vorticity": "1/s", "flow_per_depth": "m^2/s",
                },
                "scope": "2D planar rigid-wall Newtonian CFD; particle probes are one-way",
                "initial_config": c,
                "womersley_number_initial": float(calculate_womersley_number(
                    config.density, config.heart_rate, config.diameter, config.viscosity)),
            }
            self.metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
        if epoch not in self._seen_epochs:
            entry = {
                "epoch": int(epoch),
                "time_s": 0.0,
                "config": asdict(config),
                "stenosis_report": stenosis_report(config, geometry),
            }
            with self.epochs_path.open("a", encoding="utf-8") as handle:
                handle.write(json.dumps(entry, separators=(",", ":")) + "\n")
            self._seen_epochs.add(epoch)

    def record(self, solver, diagnostics, epoch, tracers=None):
        """Append one row when the simulated-time interval has elapsed."""
        with self._lock:
            if self._closed:
                return False
            time_s = float(solver.time)
            last = self._last_time.get(epoch, -np.inf)
            if time_s < last or time_s - last + 1e-12 < self.every_s:
                return False
            self._ensure_open(solver.config, solver.geometry, epoch)
            c = solver.config
            d = diagnostics
            pulse = float(pulse_curve(c, np.asarray([time_s]))[0])
            phase = (time_s * c.heart_rate / 60.0) % 1.0
            wss_low = np.abs(np.asarray(d["wss_low"], dtype=float))
            wss_high = np.abs(np.asarray(d["wss_high"], dtype=float))
            wss_all = np.concatenate((wss_low, wss_high))
            vort = np.asarray(d["vorticity"], dtype=float)
            mask = solver.geometry.fluid
            gap_in = float(solver.geometry.upper[1] - solver.geometry.lower[1])
            gap_out = float(solver.geometry.upper[-2] - solver.geometry.lower[-2])
            count_in = max(int(mask[:, 1].sum()), 1)
            count_out = max(int(mask[:, -2].sum()), 1)
            inlet_u = float(d["mean_u"][1])
            outlet_u = float(d["mean_u"][-2])
            row = {
                "epoch": int(epoch),
                "time_s": time_s,
                "cycle_index": int(np.floor(time_s * c.heart_rate / 60.0)),
                "phase_fraction": float(phase),
                "phase": phase_label(c, time_s),
                "pulse_factor": pulse,
                "inlet_speed_m_s": float(c.mean_velocity * pulse),
                "mean_inlet_velocity_m_s": inlet_u,
                "mean_outlet_velocity_m_s": outlet_u,
                "inlet_flow_per_depth_m2_s": inlet_u * gap_in,
                "outlet_flow_per_depth_m2_s": outlet_u * gap_out,
                "delta_p_pa": float(d["delta_p"]),
                "max_velocity_m_s": float(d["max_velocity"]),
                "min_axial_velocity_m_s": float(d["min_u"]),
                "backflow_percent": float(d["backflow_percent"]),
                "max_wss_pa": float(np.nanmax(wss_all)),
                "mean_abs_wss_pa": float(np.nanmean(wss_all)),
                "max_vorticity_s_1": float(np.nanmax(np.abs(vort[mask]))),
                "flux_error": float(d["flux_error"]),
                "mach": float(d["mach"]),
                "density_variation": float(d["density_variation"]),
                "particle_count": int(len(tracers.x)) if tracers is not None else 0,
                "particles_recycled": int(tracers.recycled) if tracers is not None else 0,
                "particle_relaxation_time_s": float(d["particle_relaxation_time_s"]),
                "particle_stokes_number": float(d["particle_stokes_number"]),
            }
            self._writer.writerow(row)
            self._file.flush()
            self._last_time[epoch] = time_s
            return True

    def due(self, epoch, time_s):
        """Return whether a simulated-time sample is due for an epoch."""
        return float(time_s) >= self._last_time.get(epoch, -np.inf) + self.every_s - 1e-12

    def close(self):
        with self._lock:
            self._closed = True
            if self._file is not None:
                self._file.flush()
                self._file.close()
                self._file = None
