# HemoFlow

See the [repository overview and completed results](../README.md) for the
verification history and latest quantitative findings. This page covers the
application layout and execution.

HemoFlow is a Python project for exploring pulsatile flow through stenosed
two-dimensional channels. It combines an interactive flow viewer with controlled
experiments on outlet placement, cycle settling, and grid sensitivity.

**Status:** numerical verification and sensitivity analysis in progress.
Patient-specific and independent physical validation have not been established.

![HemoFlow viewer](docs/images/viewer.png)

*Example interface from the supplied project. The image is illustrative; it is
not a verification result.*

## Features

- Editable plaque geometry and reproducible random layouts.
- Fixed heart rate with sinusoidal or illustrative systolic inlet modulation.
- Velocity, pressure, vorticity, and approximate wall-shear displays.
- Passive tracers and a dilute, one-way inertial-particle option.
- Headless recording, geometry reports, and source/parameter provenance.
- Separate experiments with fixed physical probes and phase-matched fields.

## Run the viewer

Python 3.12 or 3.13 and the packages in `requirements.txt` are required.
From this repository folder in PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe main.py --defaults --cells 80 --fps 45
```

On Windows, `START_HemoFlow.bat` also creates a local environment when needed,
installs dependencies on its first launch, and forwards command-line options:

```powershell
.\START_HemoFlow.bat --interactive --fps 45
```

`--fps` sets the target display rate. Simulation speed depends on the grid,
timestep, and CPU. `--particles-off` disables particle advection and rendering.
Without `--defaults`, the viewer loads local `settings.json` when present.

## Record a short run

```powershell
.\.venv\Scripts\python.exe main.py --defaults --mode healthy --steady --cells 40 --headless --seconds 0.1 --record --record-every 0.005 --record-dir exports/demo/raw --output exports/demo/final
```

This is a short execution example. Long-run conclusions require settled flow
and the comparison checks described in [verification](docs/verification.md).
CSV recordings can be opened in Excel; JSON files store settings and provenance.

## Repository layout

| Location | Purpose |
| --- | --- |
| `main.py`, `config.py` | Application entry point and inputs |
| `solver.py`, `kernel.py` | Desktop solver and CPU stepping kernels |
| `geometry.py`, `waveform.py`, `physics.py` | Geometry, forcing, and reference calculations |
| `engine.py`, `visualization.py`, `rendering.py` | Background computation and display |
| `particles.py`, `recorder.py`, `provenance.py` | Particle probes and recording |
| `tests/` | Desktop regression tests |
| `studies/length_confirmation/` | Fixed-grid, 120/150 mm domain comparison |
| `studies/grid_refinement/` | Controlled-wave-speed grid comparison at 150 mm |
| `docs/` | Model assumptions, usage, and verification scope |

The desktop application uses the `regularized` outlet by default, with
`zou_he` available as a baseline. The recent study packages use their own
frozen solver with the experimental `section_impedance` outlet. These are
distinct implementations; study conclusions do not automatically transfer to
the desktop viewer. See [study protocols](studies/README.md).

## Checks

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
```

The tests cover planar flow references, boundary reconstruction, mass accounting,
geometry, particles, recording, and viewer-worker behavior. A passing regression
suite does not establish mesh independence or physical validation.

The model is planar, Newtonian, weakly compressible, and rigid-walled. Particle
markers do not represent resolved red blood cells, hematocrit, platelet adhesion,
or thrombus growth. Further details are in [the model description](docs/model.md).
