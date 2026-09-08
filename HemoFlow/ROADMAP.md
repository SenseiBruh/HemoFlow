# HemoFlow roadmap (v3)

Requested outcome: preserve the earlier simulator's inputs and readouts, keep it running continuously, make systole and diastole visible, and show blood-flow redirection and recirculation around editable randomly generated plaques.

## Delivered in this build

| Stage | Work | Result |
| --- | --- | --- |
| 1. Recover and preserve | Kept the density, viscosity, vessel size, inlet speed, stenosis, heartbeat, Reynolds number, pressure/WSS references, heatmap, and particle view. | The original inputs remain in `Config`; the circular-pipe formulas are retained and labeled separately from planar CFD. |
| 2. Seeded plaque geometry | Added variable count, center, width, severity, asymmetry, shape exponent, and wall; retained healthy and single-stenosis modes. | Default random mode makes one or two plaques, keeps the lumen connected, and reproduces a layout from its seed. |
| 3. Pulsatile flow | Added a sharper systolic waveform as well as the original sine option. | Inlet speed, solved velocity, pressure, shear, particles, and the pulse marker vary through simulated systole and diastole. |
| 4. Threaded continuous viewer | Moved the solver and tracer state into a worker thread with a bounded snapshot ring. | The Matplotlib thread renders and interpolates published states; pause holds a state, resume continues it, and the simulation has no preset end time. |
| 5. Plaque interaction | Added selection, body dragging, wall switching, width handles, height handles, add, delete, and keyboard/button controls. | Release commits a validated custom geometry and restarts only the flow state for that geometry. |
| 6. Save and compare | Kept dashboard/field/profile/settings exports and added reproducible custom plaque JSON. | A saved snapshot contains the solver time step and iteration, so it can be compared or replayed from the same geometry. |
| 7. Measurement labels | Added requested, continuous, and rasterized 2D narrowing plus circular-3D area-equivalent reporting. | Every plaque export includes its realized throat, diameter equivalent, area-equivalent value, and grid spacing. |
| 8. Particle experiment controls | Added particle visibility, deterministic display decimation, diameter/density controls, and tracer/inertial one-way Stokes relaxation. | Particle sweeps do not alter the CFD field; their assumptions are explicit in the UI and metadata. |
| 9. Fixed-rate recording | Added GUI and headless recording with rate-limited CSV observables, metadata, and geometry epochs. | A run can stay on for hours without retaining an unbounded frame or plot history. |
| 10. Fixed BPM controls and comparison value | Added BPM slider/buttons and Womersley-number readout. | BPM is constant within a run; pulse amplitude remains a separate ±% modulation. |

## Next stages

| Order | Addition | Completion evidence |
| --- | --- | --- |
| 11 | Convergence and boundary studies | Compare at least three grids and longer entrance/exit domains; report pressure, flow, recirculation length, and wall-shear changes. |
| 12 | Improved wall geometry and stress | Add interpolated or immersed boundaries and a stress-tensor WSS estimate; validate on curved-channel benchmarks. |
| 13 | More faithful blood and vessel behavior | Add a documented Carreau–Yasuda option, measured inlet waveforms, and optional wall compliance with matching benchmarks. |
| 14 | 2D carotid bifurcation | Replace the single vertical interval with a common-carotid bulb and ICA/ECA branches, branch flow split, and outlet boundary conditions. This is a geometry/BC change, not a plaque toggle. |
| 15 | 3D vessel and plaque geometry | Replace the planar approximation with 3D flow, cross-sectional measurements, and 3D particle paths. |
| 16 | Biological growth, if desired | Define a growth law and timescale separately from random layout regeneration. |

## Design rules to keep

- Flow and tracers advance in physical simulation time. Rendering latency and smoothing never alter viscosity or heartbeat frequency.
- Recirculation must emerge from the velocity solution. Do not add decorative circular paths and label them computed vortices.
- Retain the healthy planar benchmark and the original circular-pipe references, and distinguish both from stenotic CFD outputs.
- Keep particles in the fluid and keep frame, particle, and plot-history storage bounded for long sessions.
- Keep a plaque layout fixed during a run and reproducible by seed. A geometry edit starts a new solution for that geometry.
- Report diameter narrowing and circular area-equivalent narrowing as separate quantities; do not call the latter a simulated 3D area.
- Treat particles as one-way probes until a documented concentration-coupled model is added.
- Label approximate wall shear, startup transients, weak compressibility, and numerical limits. A successful display is not clinical validation.
