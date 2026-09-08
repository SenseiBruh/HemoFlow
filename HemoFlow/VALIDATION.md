# Checks performed for this v3 delivery

Run date: September 6, 2026. Runtime: Linux, NumPy 2.3.5, Matplotlib 3.10.8, Numba 0.67.0. The dashboard was rendered with Matplotlib's Agg backend and its worker, pause, reset, map, waveform, plaque, injection, and export paths were exercised. A native Windows window and the batch launcher were not executed in this environment.

## Automated checks

```text
python -m unittest discover -s tests -v

13 tests passed
```

The checks cover the original reference formulas, a healthy planar-channel solution, computed backflow downstream of a plaque, seeded connected geometry, tracer confinement/recycling/bounded trails, small random plaque counts, custom wall edits, systolic extrema, bounded worker snapshot interpolation, separate diameter/area-equivalent stenosis reporting, inertial particle relaxation, visibility toggling, rate-limited recording, and per-cycle CSV summarization.

## v3 additions

For a circular comparison, the identity is explicit:

| Diameter narrowing | Circular area-equivalent reduction |
| ---: | ---: |
| 25% | 43.75% |
| 50% | 75.00% |
| 75% | 93.75% |

The viewer and exports label this as a 3D circular area equivalent. The solved domain remains a 2D planar, rasterized lumen, and each plaque also reports its realized grid-based diameter narrowing.

Particle tests confirm that a 100 µm, 1200 kg/m³ one-way inertial probe receives a finite Stokes relaxation time and remains in the fluid. The solver reports both relaxation time and cardiac-scale Stokes number. Turning particles off publishes zero particle arrays while heatmaps continue; turning them back on re-seeds the bounded tracer set without invalid interpolation.

The recorder writes an append-only scalar CSV, metadata JSON, and epoch JSONL. Rows are sampled by simulated time and include cardiac phase, pressure, shear, vorticity, backflow, flux error, and particle recycling. `analyze_recording.py` groups rows by epoch and cardiac cycle for mean/std/min/max summaries.

## Healthy planar reference

A healthy planar channel was run for 10,000 steps with width 4 mm, length 16 mm, mean velocity 0.1 m/s, density 1060 kg/m³, viscosity 0.0035 Pa·s, steady inlet, and 32 cells across. Pressure comparisons use the same two stations as the numerical diagnostic.

| Quantity | Computed | Analytic / comparison | Relative error |
| --- | ---: | ---: | ---: |
| Pressure drop | 4.17374 Pa | 4.13438 Pa | 0.952% |
| Mean lower-wall shear estimate | 0.51108 Pa | 0.52500 Pa | 2.651% |
| Velocity profile | — | Relative L2 error at the middle section | 0.0771% |
| Inlet/outlet flux mismatch | — | Relative to the inlet | 0.276% |
| Backward-flow cells | 0% | None expected | — |

The original circular-pipe references were also checked: defaults give Re = 363.4286, pressure drop = 126.00 Pa, and wall shear = 2.100 Pa. Those values are retained for comparison and are not stenotic planar CFD results.

## Current seeded plaque probe

The shipped settings use `pulse_shape: "systolic"`, seed 42, 50% maximum narrowing, 72 BPM, ±40% pulsatility, and automatic one-or-two-plaque generation. With 40 cells across and the seed held at 42, the current layout contains two plaques and has a 24-cell minimum throat. A 2,000-step probe reaches 0.01905 simulated seconds without a stability failure.

| Quantity | Value |
| --- | ---: |
| Grid | 601 × 42 including boundary rows |
| Physical time step | 9.52381 μs |
| CFD pressure drop at probe | 179.36 Pa |
| Maximum speed at probe | 0.4212 m/s |
| Minimum axial velocity at probe | −0.0409 m/s |
| Backflow cells at probe | 0.937% |
| Maximum lattice Mach number | 0.0695 |
| Maximum density departure | 0.463% |
| Instantaneous flux mismatch | 4.46% |

This short probe is a startup snapshot, not a converged cycle average. Longer runs develop the downstream wake further and cost more wall time.

## Waveform check

For the default 72 BPM and ±40% modulation, the systolic curve varies inlet speed from 0.18 to 0.42 m/s and reaches its peak during the first quarter of each 0.833 s cycle. The dashboard plots the curve and marks `SYSTOLE` or `DIASTOLE` at the current simulated time, so the change remains visible even when the heatmap is rendered at a slower wall-clock rate.

In v3 the BPM slider/buttons rebuild the next flow state at a fixed selected frequency. The `±40%` value remains an inlet-speed amplitude: it makes speed range from 60% to 140% of the configured mean; it does not make the heart rate vary between 32 and 102 BPM.

## Dashboard behavior

The exercised Agg dashboard checks produced the following results:

```json
{
  "pause_resume": true,
  "map_toggle": true,
  "pulse_toggle": true,
  "plaque_select_drag_move": true,
  "plaque_height_handle": true,
  "plaque_add_delete": true,
  "injection": true,
  "snapshot_exports": true,
  "settings_roundtrip": true,
  "regeneration": true,
  "bounded_snapshot_ring": true,
  "figure_axes_after_resets": 20
}
```

The worker publishes at most eight frames, copies particle positions and trails into each immutable snapshot, and accepts bounded injection/reset requests. Resetting a geometry increments an epoch so an older solver cannot overwrite the new layout. Snapshot NPZ exports include the solver iteration and physical time step.

## Preview

`preview.gif` is a short excerpt rendered from the threaded solver. It loops for convenient inspection, while the desktop viewer continues advancing rather than replaying that recording. The lower map and tracer paths should be read together: signed vorticity highlights local rotation, while signed axial velocity and curved/backward tracer paths show whether a recirculating pocket is present.

## Remaining scientific limits

No grid-convergence or outlet-extension study has been completed for plaque cases. No published stenotic benchmark, 3D artery, patient data, or clinical measurement has been quantitatively matched. Wall-shear accuracy is checked on a flat healthy channel only; staircase plaque estimates are noisier and less reliable. The model is a first 2D numerical visualization, not a clinical predictor.
