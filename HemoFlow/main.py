"""Run ``python main.py`` for the continuous viewer.

Use ``--interactive`` to change the physical inputs before opening it, or
``--headless`` to export a numerical run without creating a window.  Headless
recording can run for a fixed duration, a number of cardiac cycles, or until
Ctrl+C for unattended experiments.
"""
import argparse
from dataclasses import replace
import json
import math
from pathlib import Path
import secrets
import sys
import time
from config import Config
from provenance import VERSION
from physics import calculate_pressure_drop, calculate_wall_shear_stress, calculate_womersley_number, classify_flow


def ask(label, default, integer=False):
    while True:
        answer = input(f"{label} [{default:g}]: ").strip()
        if not answer:
            return default
        try:
            value = int(answer) if integer else float(answer)
            if not math.isfinite(value):
                raise ValueError
            return value
        except ValueError:
            print("Please enter a finite number" + (" without decimals." if integer else "."))


def interactive(config):
    print("\nHemoFlow — press Enter to keep any displayed value.\n")
    fields = [
        ("density", "Blood density (kg/m^3)"), ("viscosity", "Blood viscosity (Pa*s)"),
        ("diameter_mm", "Vessel diameter (mm)"), ("length_mm", "Vessel length (mm)"),
        ("mean_velocity", "Mean blood velocity (m/s)"),
        ("stenosis_percent", "Maximum diameter/gap narrowing (%)"),
        ("stenosis_center_mm", "Single-stenosis center (mm; single mode only)"),
        ("stenosis_length_mm", "Stenosis length (mm; max plaque length in random mode)"),
        ("heart_rate", "Heart rate (BPM)"), ("pulsatility_percent", "Pulsatility (+/- %)"),
        ("plaque_count", "Fixed plaque count (when random count is off)"), ("seed", "Random seed"),
        ("cells_across", "Grid cells across vessel (higher = finer/slower)"),
        ("particle_count", "Number of simulated tracers"),
        ("particle_diameter_um", "Particle diameter (micrometres)"),
        ("particle_density_kg_m3", "Particle density (kg/m^3)"),
        ("particle_display_limit", "Maximum particles drawn (0 = all)"),
        ("max_random_plaques", "Maximum random plaques (1-3)"),
        ("time_scale", "Numerical time scale (0.1-1; lower Mach, more time steps)"),
        ("compute_threads", "Solver threads (0 = automatic)")]
    values = {}
    for key, label in fields:
        values[key] = ask(label, getattr(config, key), key in ("plaque_count", "seed", "cells_across", "particle_count", "particle_display_limit", "max_random_plaques", "compute_threads"))
    shape = input(f"Pulse shape: sine / systolic [{config.pulse_shape}]: ").strip().lower()
    if shape:
        values["pulse_shape"] = shape
    randomize = input(f"Randomize plaque count: yes / no [{'yes' if config.randomize_count else 'no'}]: ").strip().lower()
    if randomize:
        values["randomize_count"] = randomize in ("yes", "y", "true", "1")
    mode = input(f"Geometry: random / single / healthy [{config.geometry_mode}]: ").strip().lower()
    values["geometry_mode"] = mode or config.geometry_mode
    model = input(f"Particle model: tracer / inertial [{config.particle_model}]: ").strip().lower()
    if model:
        values["particle_model"] = model
    visible = input(f"Show particles: yes / no [{'yes' if config.particles_visible else 'no'}]: ").strip().lower()
    if visible:
        values["particles_visible"] = visible in ("yes", "y", "true", "1")
    return replace(config, **values).validate()


def headless(config, seconds, output, record_dir=None, record_every=.05,
             continuous=False, cycles=0):
    import numpy as np
    from solver import FlowSolver
    from recorder import RunRecorder
    from provenance import run_manifest
    if cycles < 0 or (not isinstance(cycles, int)):
        raise ValueError("--cycles must be a non-negative integer")
    if cycles:
        seconds = cycles * 60.0 / config.heart_rate
    if continuous and cycles:
        continuous = False
    if not continuous and (not math.isfinite(seconds) or seconds <= 0):
        raise ValueError("--seconds must be finite and positive")
    solver = FlowSolver(config)
    print("Selecting CPU thread count (or applying --threads)...", flush=True)
    solver.tune_threads()
    print(f"Using {solver.threads_used} CPU thread(s); dt={solver.dt:.6g} s; tau={solver.tau:.6f}; time scale={config.time_scale:g}", flush=True)
    target = math.ceil(seconds / solver.dt) if not continuous else None
    recorder = RunRecorder(record_dir, record_every) if record_dir else None
    started = time.perf_counter()
    next_report = started + 5
    stopped_by_user = False
    try:
        if recorder is not None:
            recorder.record(solver, solver.diagnostics(), 0, None)
        while target is None or solver.iteration < target:
            steps = min(200, target - solver.iteration) if target is not None else 200
            if recorder is not None:
                sample_step = math.ceil((recorder.next_time(0) - 1e-12) / solver.dt)
                steps = min(steps, max(1, sample_step - solver.iteration))
            solver.step(steps)
            solver.check_stability()
            if recorder is not None and recorder.due(0, solver.time):
                recorder.record(solver, solver.diagnostics(), 0, None)
            now = time.perf_counter()
            if now >= next_report:
                suffix = " (continuous; Ctrl+C to stop)" if target is None else ""
                rate = solver.time / max(now - started, 1e-9)
                eta = "" if target is None else f"; approx {(target*solver.dt-solver.time)/max(rate,1e-12)/60:.1f} min remaining"
                print(f"Simulated {solver.time:.3f} s ({rate:.4f} simulated s / wall s){eta}{suffix}", flush=True)
                next_report = now + 5
    except KeyboardInterrupt:
        stopped_by_user = True
        print(f"Stopping after {solver.time:.3f} simulated seconds; writing final exports.", flush=True)
    finally:
        if recorder is not None:
            recorder.record(solver, solver.diagnostics(), 0, None, force=True)
            recorder.close()
    d, g = solver.diagnostics(), solver.geometry
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / "field.npz", x_m=g.x, y_m=g.y, solid=g.solid,
                        u_m_s=d["u"], v_m_s=d["v"], pressure_pa=d["p"], vorticity_s=d["vorticity"],
                        time_s=solver.time, config_json=json.dumps(config.__dict__))
    metrics = {key: value for key, value in d.items() if np.isscalar(value)}
    metrics.update(time_s=solver.time, steps=solver.iteration, dt_s=solver.dt,
                   elapsed_s=time.perf_counter()-started, grid=[solver.nx,solver.ny],
                   stopped_by_user=stopped_by_user,
                   particle_model=config.particle_model,
                   particle_diameter_um=config.particle_diameter_um,
                   particle_density_kg_m3=config.particle_density_kg_m3,
                   geometry_report=d["geometry_report"])
    (output / "metrics.json").write_text(json.dumps(metrics, indent=2)+"\n", encoding="utf-8")
    (output / "stenosis_report.json").write_text(json.dumps(d["geometry_report"], indent=2)+"\n", encoding="utf-8")
    if recorder is not None:
        metrics["recording"] = recorder.paths
        (output / "metrics.json").write_text(json.dumps(metrics, indent=2)+"\n", encoding="utf-8")
    config.save(output / "settings.json")
    (output / "run_manifest.json").write_text(json.dumps(run_manifest(solver), indent=2)+"\n", encoding="utf-8")
    print(json.dumps(metrics, indent=2))


def main():
    parser = argparse.ArgumentParser(description="HemoFlow: continuous 2D blood-flow visualization")
    parser.add_argument("--interactive", action="store_true", help="Prompt for all original and new inputs")
    parser.add_argument("--defaults", action="store_true", help="Ignore saved settings")
    parser.add_argument("--config", type=Path, help="Load a specific settings JSON")
    parser.add_argument("--mode", choices=["random", "single", "healthy"])
    parser.add_argument("--seed", type=int)
    parser.add_argument("--bpm", "--heart-rate", dest="heart_rate", type=float, help="Fixed heart rate in BPM")
    parser.add_argument("--cells", type=int, help="Grid cells across the unobstructed vessel")
    parser.add_argument("--threads", type=int, help="Solver threads (0 = automatic)")
    parser.add_argument("--time-scale", type=float, help="0.1-1: smaller reduces lattice Mach at fixed physical inputs (default 0.5; v3 used 1)")
    parser.add_argument("--boundary-model", choices=["regularized", "zou_he"],
                        help="Open-end boundary reconstruction; zou_he reproduces v3.1")
    parser.add_argument("--fps", type=int, help="Viewer target frame rate, 10-60 (does not change physics)")
    parser.add_argument("--pulse-shape", choices=["sine", "systolic"], help="Inlet waveform shape")
    parser.add_argument("--steady", action="store_true", help="Disable the heartbeat modulation")
    parser.add_argument("--particle-model", choices=["tracer", "inertial"], help="One-way particle model")
    parser.add_argument("--particle-diameter-um", type=float, help="Particle diameter in micrometres")
    parser.add_argument("--particle-density", type=float, help="Particle density in kg/m^3")
    parser.add_argument("--particle-display-limit", type=int, help="Maximum particles drawn (0 = all)")
    parser.add_argument("--particles-off", action="store_true", help="Start without particle advection/rendering")
    parser.add_argument("--headless", action="store_true", help="Compute and export without a window")
    parser.add_argument("--seconds", type=float, default=1.0, help="Simulated duration for --headless")
    parser.add_argument("--cycles", type=int, default=0, help="Record this many cardiac cycles instead of --seconds")
    parser.add_argument("--continuous", action="store_true", help="Run headless until Ctrl+C (or use --cycles)")
    parser.add_argument("--record", action="store_true", help="Record compact time-series observables")
    parser.add_argument("--record-dir", default="exports/runs", help="Folder for time-series recording")
    parser.add_argument("--record-every", type=float, default=.05, help="Recording interval in simulated seconds")
    parser.add_argument("--output", default="exports/headless", help="Output folder for --headless")
    args = parser.parse_args()
    saved = Path(__file__).resolve().parent / "settings.json"
    config = Config()
    if args.config:
        config = Config.load(args.config)
    elif not args.defaults and saved.exists():
        config = Config.load(saved)
    overrides = {}
    for attr, value in (("geometry_mode", args.mode), ("seed", args.seed), ("cells_across", args.cells),
                        ("compute_threads", args.threads), ("pulse_shape", args.pulse_shape),
                        ("time_scale", args.time_scale), ("viewer_fps", args.fps),
                        ("boundary_model", args.boundary_model),
                        ("heart_rate", args.heart_rate), ("particle_model", args.particle_model),
                        ("particle_diameter_um", args.particle_diameter_um),
                        ("particle_density_kg_m3", args.particle_density),
                        ("particle_display_limit", args.particle_display_limit)):
        if value is not None:
            overrides[attr] = value
    if args.steady:
        overrides["pulsatility_percent"] = 0.0
    if args.particles_off:
        overrides["particles_visible"] = False
    config = replace(config, **overrides).validate()
    if args.interactive:
        config = interactive(config)
    # A launch gets a fresh random layout by default. Passing --seed keeps a
    # run reproducible, and custom layouts are never replaced by this step.
    if config.randomize_on_launch and args.seed is None and not args.interactive and config.geometry_mode == "random":
        config = replace(config, seed=secrets.randbelow(2**32)).validate()
    print(f"\nHemoFlow v{VERSION} | continuous 2D flow\n", flush=True)
    print(f"Reynolds number: {config.reynolds:.1f} | {classify_flow(config.reynolds)}")
    print(f"Fixed heart rate: {config.heart_rate:g} BPM | Womersley comparison α: "
          f"{calculate_womersley_number(config.density, config.heart_rate, config.diameter, config.viscosity):.2f}")
    dp = calculate_pressure_drop(config.viscosity,config.length,config.mean_velocity,config.diameter)
    wss = calculate_wall_shear_stress(config.viscosity,config.mean_velocity,config.diameter)
    print(f"Original healthy circular-pipe references: pressure drop {dp:.2f} Pa; WSS {wss:.3f} Pa.")
    print("The live CFD is a planar channel model; its pressure and shear differ from a circular pipe.")
    if args.headless:
        headless(config, args.seconds, args.output,
                 record_dir=(args.record_dir if args.record else None),
                 record_every=args.record_every, continuous=args.continuous,
                 cycles=args.cycles)
    else:
        from visualization import Dashboard
        print("Preparing the solver. The first launch compiles its numerical kernels; later launches are faster.", flush=True)
        dashboard = Dashboard(config)
        print("Space: pause | R: regenerate plaques | V: switch map | click vessel: inject tracers", flush=True)
        dashboard.show()


if __name__ == "__main__":
    try:
        main()
    except (ValueError, FloatingPointError, OSError, TypeError) as exc:
        print(f"\nHemoFlow: {exc}", file=sys.stderr)
        raise SystemExit(1)
    except KeyboardInterrupt:
        print("\nStopped.")
