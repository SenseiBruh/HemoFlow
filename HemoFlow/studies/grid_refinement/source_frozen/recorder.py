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
from provenance import run_manifest


SCHEMA_VERSION = 3
EXTRA_METRICS = (
    "volume_flux_mismatch_relative", "inlet_mass_flux_per_depth_kg_m_s",
    "outlet_mass_flux_per_depth_kg_m_s", "fluid_mass_per_depth_kg_m",
    "control_volume_mass_per_depth_kg_m", "link_inlet_mass_flux_kg_m_s",
    "link_outlet_mass_flux_kg_m_s", "mass_storage_rate_kg_m_s",
    "mass_balance_residual_kg_m_s", "mass_balance_relative", "mass_balance_valid",
)
EXTREMUM_LOCATIONS = (
    "max_wss_x_mm", "max_wss_wall_index", "min_u_x_mm", "min_u_y_mm",
    "max_density_x_mm", "max_density_y_mm",
)


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
        self._next_time = {}
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

    def _ensure_open(self, solver, epoch):
        config, geometry = solver.config, solver.geometry
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
                "dt_s", "tau", "compute_threads", "time_scale", "boundary_model",
                "realized_raster_diameter_narrowing_percent", "equivalent_circular_3d_area_reduction_percent",
            ] + list(EXTRA_METRICS) + list(EXTREMUM_LOCATIONS))
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
                    "mass_per_depth": "kg/m", "mass_flux_per_depth": "kg/(m*s)",
                },
                "scope": "2D planar rigid-wall Newtonian CFD; particle probes are one-way",
                "initial_config": c,
                "manifest": run_manifest(solver),
                "diagnostic_definitions": {
                    "flux_error": "Legacy alias of volume_flux_mismatch_relative: (Qin-Qout)/abs(Qin); NOT mass conservation error.",
                    "mass_balance_relative": "Last lattice step: (link inlet - link outlet - independently measured storage)/(rho_ref*mean_U*nominal_gap). Use only where mass_balance_valid is true.",
                    "link_fluxes": "Population crossings at the faces bounding interior columns 1..nx-2; distinct from nodal rho*u fluxes.",
                    "density_variation": "Maximum |rho/rho_ref - 1|; compressibility diagnostic, not a conservation residual.",
                    "backflow_percent": "Percentage of interior fluid nodes with u < -0.01*mean_U; not reverse volumetric flow.",
                    "sampling": "Samples on fixed simulated-time targets, rounded up to a lattice step. Initial/final endpoint rows may be closer together.",
                    "extremum_locations": "Coordinates of the reported sampled extrema; max_wss_wall_index is 0 for lower and 1 for upper wall. Tied extrema use the first array index. Locations can move between samples.",
                },
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
                "numerics": solver.numerics(),
            }
            with self.epochs_path.open("a", encoding="utf-8") as handle:
                handle.write(json.dumps(entry, separators=(",", ":")) + "\n")
            self._seen_epochs.add(epoch)

    def record(self, solver, diagnostics, epoch, tracers=None, force=False):
        """Append one row when the simulated-time interval has elapsed."""
        with self._lock:
            if self._closed:
                return False
            time_s = float(solver.time)
            last = self._last_time.get(epoch, -np.inf)
            if time_s <= last or (not force and not self.due(epoch, time_s)):
                return False
            self._ensure_open(solver, epoch)
            c = solver.config
            d = diagnostics
            pulse = float(pulse_curve(c, np.asarray([time_s]))[0])
            phase = (time_s * c.heart_rate / 60.0) % 1.0
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
                "inlet_flow_per_depth_m2_s": float(d["inlet_flow_per_depth_m2_s"]),
                "outlet_flow_per_depth_m2_s": float(d["outlet_flow_per_depth_m2_s"]),
                "delta_p_pa": float(d["delta_p"]),
                "max_velocity_m_s": float(d["max_velocity"]),
                "min_axial_velocity_m_s": float(d["min_u"]),
                "backflow_percent": float(d["backflow_percent"]),
                "max_wss_pa": float(d["max_wss"]),
                "mean_abs_wss_pa": float(d["mean_abs_wss"]),
                "max_vorticity_s_1": float(d["max_vorticity"]),
                "flux_error": float(d["flux_error"]),
                "mach": float(d["mach"]),
                "density_variation": float(d["density_variation"]),
                "particle_count": int(len(tracers.x)) if tracers is not None else 0,
                "particles_recycled": int(tracers.recycled) if tracers is not None else 0,
                "particle_relaxation_time_s": float(d["particle_relaxation_time_s"]),
                "particle_stokes_number": float(d["particle_stokes_number"]),
                "dt_s": solver.dt, "tau": solver.tau, "time_scale": c.time_scale,
                "boundary_model": c.boundary_model,
                "compute_threads": solver.threads_used,
                "realized_raster_diameter_narrowing_percent": d["geometry_report"]["overall_realized_raster_diameter_narrowing_percent"],
                "equivalent_circular_3d_area_reduction_percent": d["geometry_report"]["overall_equivalent_circular_3d_area_reduction_percent"],
            }
            row.update({name: float(d[name]) for name in EXTRA_METRICS + EXTREMUM_LOCATIONS})
            self._writer.writerow(row)
            self._file.flush()
            self._last_time[epoch] = time_s
            self._next_time[epoch] = (np.floor((time_s + 1e-12) / self.every_s) + 1) * self.every_s
            return True

    def due(self, epoch, time_s):
        """Return whether a simulated-time sample is due for an epoch."""
        return float(time_s) >= self.next_time(epoch) - 1e-12

    def next_time(self, epoch):
        return float(self._next_time.get(epoch, 0.0))

    def close(self):
        with self._lock:
            self._closed = True
            if self._file is not None:
                self._file.flush()
                self._file.close()
                self._file = None
