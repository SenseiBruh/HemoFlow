# Verification and sensitivity analysis

Completed experiments and their source summaries are recorded in
[results.md](results.md). The latest retained length-confirmation report passes
the N80 120/150 mm comparison for its specified observables and regions.

The project distinguishes software regression tests, analytical verification,
numerical sensitivity, and independent physical validation. None of these
categories is interchangeable with another.

| Check | Evidence and scope |
| --- | --- |
| Software regressions | `tests/`: geometry, boundary moments, parallel/serial kernels, mass accounting, recording, and worker behavior |
| Healthy planar benchmark | Steady channel velocity, pressure difference, and wall shear compared with planar reference values |
| Cycle settling | Adjacent-cycle comparisons of pressure, flow, and signed wall-shear profiles |
| Domain length | 120/150 mm comparison at fixed physical probe locations and fixed N80 grid spacing |
| Grid sensitivity | N80/N120 and conditional N160 comparison at fixed 150 mm length and controlled Δx/Δt |
| Timestep sensitivity | A separate requirement; changing Δt also changes LBM relaxation and numerical wave speed |
| Physical validation | Not established; requires an appropriate independent physical benchmark |

## Controlled severe geometry

Recent studies use a symmetric 75% gap reduction, nominal gap 4 mm, plaque
centre x = 30 mm, and plaque length 15 mm. Density, viscosity, waveform, heart
rate, initialization, and the first-cycle inlet perturbation are controlled.
The reported geometry includes the realized raster gap.

Cross-section pressure, density, and flow are recorded at fixed physical
locations. Signed wall-shear profiles cover x = 20–55 mm. Velocity fields are
saved at selected phases using their actual capture times. Global observables
are retained for continuity with earlier recordings.

## Interpretation

Cycle repeatability applies to the recorded quantities and observation window.
Agreement between domain lengths applies to the compared physical region and
grid. Neither result establishes agreement everywhere in the domain.

Grid comparisons keep Δx/Δt fixed to control numerical wave speed, while
relaxation time and raster wall resolution still vary. Two grids provide a
sensitivity comparison; they do not establish an observed convergence order.
Velocity interpolation excludes stencils that cross solid walls, so signed
wall-shear checks remain necessary.

The numerical-health and comparison thresholds are project engineering
criteria. Display interpolation, reflected field alignment, and pressure
filtering do not substitute for passing direct numerical comparisons.

## Evidence retention

Each study saves the source identity, realized timestep, settings, geometry,
recording times, and analysis results. A report should identify the exact
source version and run for every figure or table. Complete raw recordings and
phase fields belong in a retained experiment archive; selected tables and
figures can be versioned in the repository.

This source distribution excludes large experiment directories and raw
CSV/NPZ recordings. Existing smoke-check files describe package setup tests,
not completed long-run experiments. No new scientific pass is asserted by the
repository reorganization.

## Remaining work

- Review controlled grid and timestep comparisons with complete recording windows.
- Resolve sensitivity in quantities that fail the declared criteria.
- Compare against an independent benchmark with matching geometry and inputs.
- Extend particle physics only with corresponding model assumptions and checks.
