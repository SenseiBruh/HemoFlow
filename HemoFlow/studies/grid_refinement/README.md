# Controlled grid comparison at 150 mm

This study measures spatial sensitivity at fixed geometry and controlled
numerical wave speed. It is not a substitute for a separate timestep study or
an independent physical benchmark.

## Protocol

| Stage | Cells across gap | Severe timestep | Healthy timestep |
| --- | ---: | ---: | ---: |
| Initial comparison | 80 | 0.300 µs | 1.500 µs |
| Initial comparison | 120 | 0.200 µs | 1.000 µs |
| Conditional refinement | 160 | 0.150 µs | 0.750 µs |

The sequence is healthy N80, healthy N120, severe N80, severe N120, then
conditional healthy/severe N160. Healthy cases run for 2 s of steady flow;
severe cases run for 10 cardiac cycles by default. The N160 stage starts only
after the N80/N120 comparison passes the declared criteria. Failure does not
automatically change the solver, tolerances, or inlet.

All domains are 150 mm long. Severe geometry is a symmetric cosine stenosis
with 75% planar-gap reduction, nominal gap 4 mm, centre x = 30 mm, and length
15 mm. This is not a measured circular-area reduction. Density is 1060 kg/m³,
viscosity 0.0035 Pa·s, nominal mean inlet speed 0.3 m/s, heart rate 72 BPM, and
systolic velocity modulation 40%. A positive 0.001 axial skew is flux-normalized
and active for the first cycle only. Initialization is native; particles are off.

The experimental `section_impedance` solver is byte-preserved in
`source_frozen/`, with receipts in `FROZEN_SOURCE.json`. Severe Δx/Δt is
166.6667 m/s at each grid. The realized timestep, relaxation parameter, and
raster geometry are recorded, with 0.501 ≤ τ ≤ 2 enforced. The physical target
shape is fixed, but raster wall resolution and relaxation parameter vary.

## Measurements and criteria

Recordings retain global observables, fixed-section pressure/density/flow,
signed wall shear, jet centroids, and downstream probes. Eight phase fields per
cycle, startup snapshots, and endpoint fields are saved. Wall fits at fixed
physical radii 0.20 and 0.25 mm are diagnostics, not replacements for the
original signed-shear criteria.

Healthy planar analytical and numerical-health checks must pass. Within each
severe case, pressure, flow, mean/peak absolute WSS, and signed WSS must change
by less than 1% across each of the last three adjacent-cycle transitions.
Signed-wall regions are x = 20–55, 28–32, 37.5–55, and 50–55 mm.

Grid comparisons use the coarser grid as reference and require both aggregate
last-three-cycle and individual-cycle checks:

| Quantity | Scalar limit | Waveform limit |
| --- | ---: | ---: |
| Pressure difference | 2% mean | 5% |
| Mean absolute WSS | 5% mean | 5% |
| Peak absolute WSS | 10% maximum | 10% |
| Flow | 1% mean | 2% |

Signed-wall space-time L2 must be within 10% in each region. The direct jet
branch check at x = 40 mm must also pass. Final-cycle velocity-field L2 must
be below 10% at every sampled phase in each region: x = 20–55, 20–22,
22.5–37.5, 37.5–55, and 50–55 mm.

Fine-grid velocity is bilinearly interpolated to coarse physical nodes only
when all four fine-grid corners are fluid. Excluded near-wall points are not
verified by that score; shear checks remain mandatory. Reflected alignment or
pressure filtering does not alter acceptance. These engineering screens do
not establish a convergence order, grid-convergence index, or physical validation.

## Execution and runtime

From the repository root, inspect the plan without starting a batch:

```powershell
.\.venv\Scripts\python.exe studies\grid_refinement\run_grid.py --plan-only
```

To start the study:

```powershell
Set-Location studies\grid_refinement
.\START_TESTS.cmd
```

The launcher searches its own, parent, and grandparent folders for `.venv`.
The default fixed thread count is 16, capped by the environment maximum. A
benchmark estimates solver time including the conditional stage. Fine-grid
runs can take multiple days; there is no wall-clock limit or automatic sleep
prevention. Recording, analysis, and compression add overhead.

`--cycles` accepts 4–30 for a new study. `--smoke` runs short setup checks and
bypasses scientific criteria only for those non-scientific smoke runs. Existing
N80 results are not imported; a new batch supplies its own matched reference.

## Outputs and resume

Results are stored under `experiments/grid_refinement_<timestamp>/`.
Pair-specific analysis contains `grid_comparison.json` and `domain_cycles.csv`;
the latter is an inherited filename with grid-labelled comparison metadata.
Per-case cycle metrics and repeatability CSVs open in Excel.

The exporter creates timestamped `UPLOAD_GRID_SUMMARY_...zip` and
`UPLOAD_GRID_REVIEW_...zip` archives. The summary contains field-comparison
results, not the complete raw fields. Full CSV/NPZ records remain necessary.

From this study folder:

```powershell
.\START_TESTS.cmd --resume "FULL_PATH_TO_STUDY"
.\START_TESTS.cmd --export-only "FULL_PATH_TO_STUDY"
```

Completed cases are hash-checked. Interrupted cases restart in new attempt
folders; there is no timestep checkpoint. Package identity and thread count
must match. A stale `RUNNING.lock` can be removed only after confirming the
old process has stopped. Resume does not override failed comparison criteria.

This reorganized copy has edited documentation and therefore a new package
fingerprint. Existing runs must be resumed from their original unchanged
package. Historical `validation/` receipts describe setup checks, not completed
scientific runs. Current cleanup checks are recorded in
[`../../docs/cleanup_checks.json`](../../docs/cleanup_checks.json).

## Headless capacity

The experiment wrapper raises only the desktop-preview nominal-cell ceiling
from 600,000 to 1,000,000, permitting the N160/L150 grid (960,000 nominal cells).
All other configuration, timestep, and relaxation checks remain active.
The frozen collision, boundary, and wall equations are unchanged.
