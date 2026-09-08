# HemoFlow — continuous flow laboratory

HemoFlow v3 is a 2D educational blood-flow laboratory built from the earlier prototype. It keeps the original blood, vessel, stenosis, heartbeat, Reynolds-number, pressure-drop, wall-shear, heatmap, and tracer inputs while adding a continuously running solver, a smoother threaded viewer, randomized plaque geometry, live plaque editing, explicit stenosis measurement, particle-property sweeps, and unattended time-series recording.

## What is new in v3

- The viewer now reports requested narrowing, realized rasterized 2D diameter narrowing, and the corresponding circular-3D area-equivalent value separately. A 50% diameter reduction therefore appears as approximately 75% area reduction; the area value is a comparison label, not a 3D calculation.
- Particle advection/rendering can be turned off without stopping the CFD heatmaps. Particle diameter and density are adjustable, and `tracer` versus one-way `inertial` (Stokes relaxation) behavior can be compared without changing the blood solution.
- Heart rate is fixed for each run and can be changed with the BPM slider or buttons. The pulse amplitude is still labeled separately as `±%`; it never changes BPM.
- The Record control and headless `--record` mode append compact CSV observables with metadata and geometry epochs, so a PC can collect data for hours without retaining an unbounded history in memory.

## Start on Windows

1. Extract the ZIP to a normal folder such as `C:\HemoFlow Project\HemoFlow` (the folder that contains `main.py`).
2. If the ZIP was downloaded from a browser, open its **Properties**, select **Unblock**, and extract it again. The project contains Python source and batch files; it does not contain a bundled executable.
3. Open the extracted `HemoFlow` folder and double-click **START_HemoFlow.bat**. The launcher tries `py -3.13` first and falls back to the `python` command, creates `.venv`, installs the three pinned packages, and starts the viewer. The first launch compiles the Numba kernels.
4. Flow starts automatically and keeps running until you pause or close the window. The seed is shown in the status line; the default random layout is regenerated on a normal launch.

If Windows says that Smart App Control blocked the launch, install or repair Python from [python.org](https://www.python.org/downloads/windows/) or the Microsoft Store, then open a fresh Command Prompt in the HemoFlow folder and run:

```powershell
cd "C:\HemoFlow Project\HemoFlow"
python --version
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe main.py --defaults
```

On a managed PC, Smart App Control may be controlled by your administrator. Keep it enabled when possible and use the terminal path above; do not replace the Python executable with an unknown download.

Double-click **CONFIGURE_HemoFlow.bat** to enter the physical inputs, heartbeat, random plaque settings, grid resolution, particle count, and solver-thread preference. Press Enter to keep a displayed value. **Save settings** in the dashboard writes the current inputs and editable plaque layout to `settings.json`.

If PowerShell is currently one directory too high (`C:\HemoFlow Project`), use the same explicit folder first:

```powershell
cd "C:\HemoFlow Project\HemoFlow"
.venv\Scripts\python.exe main.py --defaults
```

For an existing Python environment:

```powershell
python -m pip install -r requirements.txt
python main.py
```

Use `python main.py --interactive` for the prompts, or `python main.py --headless --seconds 1 --output exports\headless` to run and export fields without opening a window.

## Viewer controls

| Control | Behavior |
| --- | --- |
| Space / Pause | Freeze and resume the same fluid state and particles. A reset never happens when pausing. |
| R / New layout | Generate a new seeded random layout and restart the flow for that geometry. |
| Click a plaque | Select it. A selected plaque shows square width handles and a circular height handle. |
| Drag plaque body | Move it along the vessel. Drag across the centerline to switch a single-wall plaque from bottom to top or top to bottom. Release to validate the geometry and recompute the fluid. |
| Drag square handles / W− / W+ / `[` `]` | Shorten or lengthen the selected plaque. |
| Drag circular handle / H− / H+ / `{` `}` | Reduce or increase its narrowing height. |
| Add plaque / A | Add a valid plaque in the next available open region. |
| Delete selected / Delete | Remove the selected plaque. |
| Click empty vessel area | Inject a small group of passive tracers at that location. |
| V / Map button | Switch the lower map between signed vorticity and signed axial velocity. Blue axial velocity is backward flow. |
| Pulse button | Toggle the sharp systolic waveform and the older smooth sine waveform. The marker and label show the current simulated phase. |
| BPM slider / BPM − / BPM + | Select a fixed heart rate for the next flow state. Release the slider to rebuild; the value does not oscillate during a run. |
| P / Particles button | Turn tracer advection and drawing off/on while the CFD heatmaps continue. Clicking the vessel while particles are off does not inject hidden tracers. |
| P − / P + | Change displayed particle diameter in micrometres. This updates one-way particle behavior without rebuilding the fluid field. |
| Model button | Toggle `tracer` (follows solved velocity) and `inertial` (dilute spherical Stokes relaxation). |
| Smoothing slider | Increase interpolation latency for smoother particles, or reduce it for lower display latency. It does not change viscosity or physical time. |
| Record button | Append compact scalar observables to `exports/runs`; the CSV is rate-limited by simulated time and survives a normal close. |
| S / Save snapshot | Save a PNG dashboard, full-field NPZ, profile CSV, and reproducible custom plaque JSON under `exports`. |
| J / Save settings | Save current inputs and the current plaque layout to `settings.json`. |

Plaque edits are previews while the mouse is held. Releasing the mouse starts a new flow state for the committed geometry, so the downstream wake and any redirection vortices are recomputed from the velocity solution. A layout stays fixed while that run advances. The selected-plaque readout distinguishes requested and realized 2D diameter narrowing from the circular-3D area-equivalent comparison.

## Random plaques and pulsatility

Random mode uses the seed to place wall-attached plaques with varied center, length, severity, asymmetry, shape exponent, and wall. The default automatic count is one or two; severe cases are weighted toward one plaque. Up to three automatic plaques can be requested with `max_random_plaques`, and a fixed count can be used by setting `randomize_count` to `false`. A given seed and configuration reproduce the same layout. A normal launch chooses a fresh seed when `randomize_on_launch` is true; `--seed 123` disables that replacement for a reproducible run.

The default `systolic` pulse rises quickly to a peak and decays through diastole. The inlet speed is still the requested mean speed multiplied by the waveform, so the solver's velocity field, pressure, shear, particle paths, and phase marker all change with systole and diastole. The pulse plot makes the change visible even when the heatmap is moving slowly. `pulse_shape: "sine"` restores the smooth sinusoid.

If the window is still heavy, start with `--particles-off` or lower `--particle-display-limit 300`; the heatmaps and CFD solver continue independently of that draw setting. A finer `--cells` grid improves geometry resolution but increases cost.

## Stenosis and particle experiments

Use the requested plaque value as a controlled input, then use the exported `stenosis_report.json` (or selected-plaque readout) as the measurement actually realized on the grid. For each plaque it includes the continuous and rasterized 2D diameter equivalents, the minimum lumen, throat location, and a circular-3D area-equivalent value:

```text
area reduction = 100 × (1 − (1 − diameter narrowing/100)^2)
```

That conversion is why 50% diameter narrowing is 75% area reduction in an ideal circular vessel. HemoFlow's straight 2D planar grid does not solve a circular cross-section; report both labels when comparing with carotid publications.

### Carotid scope

The v3 executable remains a straight-channel laboratory. A carotid bifurcation is not just another plaque: it needs a common-carotid bulb, ICA/ECA branch geometry, branch flow split, and separate outlet boundary conditions. Those changes are listed as the next roadmap stage so v3 does not label a single-channel result as a carotid prediction.

The default `tracer` particles are passive point probes. Select `inertial` to add one-way spherical Stokes relaxation; sweep `particle_diameter_um` and `particle_density_kg_m3` while holding the CFD configuration fixed. The viewer reports the particle response time and a cardiac-scale Stokes number so a sweep is documented rather than just cosmetic. This is appropriate for a dilute first study. It does not yet model particle-particle interactions, deposition, platelet adhesion, Brownian diffusion, RBC margination, or two-way feedback on blood. Nanometre-scale materials are below the current grid resolution and should be treated as sub-grid probes, not resolved objects.

## Unattended recording

The GUI Record button writes three files per run under `exports/runs`: an append-only scalar `*_timeseries.csv`, `*_metadata.json` with units and scope, and `*_epochs.jsonl` containing the exact configuration and stenosis report for each reset/edit. Recording is sampled in simulated time, so it remains bounded and comparable across computers.

For a reproducible headless run:

```powershell
.venv\Scripts\python.exe main.py --defaults --seed 42 --mode single --headless --cycles 10 --record --record-dir exports\runs\carotid_proxy_50d --output exports\carotid_proxy_50d\final
```

For unattended collection until Ctrl+C:

```powershell
.venv\Scripts\python.exe main.py --defaults --seed 42 --headless --continuous --record --record-every 0.10 --record-dir exports\runs\overnight_seed42 --output exports\overnight_seed42\final
```

The CSV columns include cardiac phase, inlet speed, inlet/outlet flow-per-depth, CFD pressure drop, maximum/mean wall-shear estimate, maximum vorticity, backflow fraction, flux error, lattice Mach number, and particle recycling count. Group by `epoch` and `cycle_index` before computing cycle means or comparing studies.

To reduce a completed recording to per-cycle means, standard deviations, and extrema:

```powershell
.venv\Scripts\python.exe analyze_recording.py exports\runs\overnight_seed42\hemoflow_run_*_timeseries.csv --output exports\runs\overnight_seed42\cycle_summary.json
```

## What carries forward

| Input / feature | Default or current behavior |
| --- | --- |
| Blood density / dynamic viscosity | 1060 kg/m³ / 0.0035 Pa·s |
| Vessel diameter / length | 4 mm / 60 mm |
| Mean inlet speed | 0.30 m/s |
| Stenosis severity | 50% maximum linear diameter/gap reduction; this is not area reduction |
| Stenosis center / maximum plaque length | 30 mm / 15 mm; center is used in single mode, length bounds random plaques |
| Heart rate / pulsatility | 72 BPM / ±40% |
| Reynolds number | Original inlet calculation, 363.4 with defaults |
| Womersley number | Reported as a circular comparison value; α≈3.0 with defaults |
| Pressure drop / wall shear | Original healthy circular-pipe references retained and labeled; live planar CFD profiles added |
| Velocity map / particles | Both remain; tracers follow solved horizontal and vertical velocity |
| Particle properties | 8 µm, 1060 kg/m³; tracer or one-way inertial model; rendering limit 650 |
| Random geometry | One or two plaques by default, seed 42 in the saved template, variable size/shape/wall |
| Stenosis reporting | Requested %, realized rasterized 2D diameter %, and circular-3D area-equivalent % |
| Solver threads | Automatic benchmark by default; use `--threads N` or `compute_threads` to override |

The recovered upload is preserved byte for byte in `legacy/recovered_original.py.txt`. It is kept as a text reference because that source snapshot ends at an incomplete main guard.

## Reading the maps

The top map shows speed magnitude with tracer dots and trails. A narrow passage produces a jet; downstream separation and recirculation can emerge from the computed velocity. The lower map shows either signed local vorticity or signed axial velocity. Vorticity alone does not prove a vortex because wall shear also contributes; use backward-flow coloring, arrows, and particle paths together.

Pressure is relative CFD pressure with a fixed outlet reference, not arterial blood pressure. The reported CFD drop uses cross-section means near the inlet and outlet. Wall shear is estimated from the first fluid node beside a rasterized wall and is less accurate on steep or curved plaques. The vessel is displayed with an expanded vertical scale so millimeter-scale features are readable.

This is a 2D, rigid-wall, Newtonian educational CFD model. It is not patient-specific or clinically validated and does not claim 3D turbulence, vessel elasticity, non-Newtonian rheology, biological plaque growth, or treatment prediction.

## Performance and numerical method

The physics thread owns the solver and particle probes. The Matplotlib thread only draws bounded snapshots, interpolates particle positions, and handles input. That separation keeps the window responsive while Numba releases the interpreter lock during the D2Q9 updates. The first launch benchmarks one, two, four, and eight threads unless a fixed thread count is supplied. The default 40 cells across is clearer but slower; `--cells 32` is a useful faster preview. The smoothing slider trades display latency for a deeper interpolation buffer, while the particle display limit and particle toggle keep Matplotlib work bounded.

`solver.py` uses a stress-regularized D2Q9 lattice Boltzmann model in physical units with halfway bounce-back walls, a flux-normalized parabolic inlet, and a fixed-density Zou–He outlet. `particles.py` uses midpoint advection with bilinear velocity sampling, collision checks, inlet recycling, and bounded trails. There are no decorative swirl paths or random velocity kicks: any recirculation shown by the viewer comes from the solved velocity field.

Run the regression checks with:

```powershell
python -m unittest discover -s tests -v
```

Useful commands:

```powershell
python main.py --defaults --mode single
python main.py --defaults --mode healthy --steady
python main.py --config settings.json --seed 123 --cells 32 --threads 4
python main.py --headless --seconds 1 --output exports\headless
```

Headless exports include `field.npz`, `metrics.json`, and `settings.json`. Simulated seconds are independent of elapsed wall time; finer grids, severe throats, and long wakes require more compute time.

## Project files

`preview_v3.png` is a static dashboard preview of the new controls; the desktop viewer itself remains live and continuously recomputes the field.

| File | Responsibility |
| --- | --- |
| `main.py` / `config.py` | Launch, prompts, validation, saved settings |
| `physics.py` / `waveform.py` | Reference formulas and inlet waveform |
| `geometry.py` | Seeded plaques, custom edits, connected-lumen checks |
| `solver.py` / `kernel.py` | Time-dependent 2D flow and tuned single/parallel kernels |
| `particles.py` | Tracer/inertial advection, trails, recycling, injection |
| `recorder.py` | Rate-limited CSV observables, metadata, and geometry epochs |
| `analyze_recording.py` | Dependency-light per-cycle summary of a recorded CSV |
| `engine.py` | Worker thread, immutable snapshots, interpolation, reset queue |
| `visualization.py` | Continuous dashboard, maps, controls, exports |
| `tests/` | Physical and behavioral regressions |
| `ROADMAP.md` / `VALIDATION.md` | Delivered work, checks, and remaining studies |
