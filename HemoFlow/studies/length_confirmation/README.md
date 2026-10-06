# Downstream-length confirmation

This study compares 120 and 150 mm domains after repeated cardiac cycles to
assess outlet-distance sensitivity at fixed physical observation locations.
It is numerical verification and sensitivity analysis, not physical validation.

## Protocol

| Case | Domain length | Cells across gap | Timestep | Duration |
| --- | ---: | ---: | ---: | --- |
| Healthy, steady | 150 mm | 80 | 1.500 µs | 2 s |
| D75, pulsatile | 120 mm | 80 | 0.300 µs | 10 cycles |
| D75, pulsatile | 150 mm | 80 | 0.300 µs | 10 cycles |

The severe geometry is a symmetric cosine stenosis: 75% planar-gap reduction,
4 mm nominal gap, 15 mm lesion length, and centre x = 30 mm. Physical inputs
are density 1060 kg/m³, viscosity 0.0035 Pa·s, nominal mean inlet speed 0.3 m/s,
72 BPM, and 40% systolic velocity modulation. Grid spacing is 50 µm.

The solver uses the experimental `section_impedance` outlet and native
initialization. A positive 0.001 axial inlet skew is flux-normalized and active
only during the first cardiac cycle. Healthy flow has no imposed skew.
The frozen solver and its SHA-256 receipts are retained in `source_frozen/`
and `FROZEN_SOURCE.json`. These runs start from initialization, not checkpoints.

## Measurements and criteria

Recordings include global observables, fixed-section pressure/density/flow,
signed wall-shear profiles, jet centroids, and downstream probes. Eight phase
fields per cycle, startup snapshots, and the endpoint are retained. Scalar
sampling intervals are 50 µs through 20 ms and 250 µs thereafter; wall profiles
are nominally sampled every 1 ms. Actual capture times are recorded.

Adjacent-cycle comparisons use phase-matched, piecewise-linear interpolation
at actual timestamps. Relative L2 norms use the earlier cycle as reference.
Pressure drop, mean/peak absolute WSS, flow, and signed-wall profiles must
change by less than 1% across each of the final three adjacent-cycle transitions.
Signed-shear regions are x = 20–55, 28–32, 37.5–55, and 50–55 mm, using both walls.

The final three cycles are also compared between domains, using 120 mm as
reference. The declared limits are:

| Quantity | Mean difference | Waveform difference |
| --- | ---: | ---: |
| Pressure difference | 2% | 5% |
| Mean absolute WSS | 5% | 5% |
| Peak absolute WSS | — | 10% |
| Flow | 1% | 2% |

Signed-wall space-time L2 must be within 10% in each region. These are project
engineering criteria, not clinical standards. No smoothing, sign removal, or
reflected alignment substitutes for direct comparisons. Scalar repeatability
does not establish full-field periodicity. Failed checks retain their results;
the runner does not change the solver or start a severity sweep automatically.

## Execution

From the repository root, inspect the plan without starting a batch:

```powershell
.\.venv\Scripts\python.exe studies\length_confirmation\run_length_confirmation.py --plan-only
```

To start the default study:

```powershell
Set-Location studies\length_confirmation
.\START_TESTS.cmd
```

The default thread count is 16, capped by the environment maximum. A benchmark
estimates solver runtime; recording and export add overhead. There is no
wall-clock limit or automatic sleep prevention. `--cycles` accepts 4–30 for a
new study. `--smoke` runs setup checks, not the scientific protocol.

## Outputs and reproducibility

Results are written beneath `experiments/length_confirmation_<timestamp>/`.
Per-case `cycle_metrics.csv`, `cycle_repeatability.csv`, and
`cycle_settling.json` describe repeatability. `analysis/domain_cycles.csv` and
`.json` describe matched-cycle domain comparisons. CSV files open in Excel.

The exporter creates timestamped `UPLOAD_LENGTH_SUMMARY_...zip` and
`UPLOAD_LENGTH_REVIEW_...zip` archives. The summary contains source, checks,
and cycle summaries; it does not replace the full raw CSV/NPZ record.

From this study folder:

```powershell
.\START_TESTS.cmd --resume "FULL_PATH_TO_STUDY"
.\START_TESTS.cmd --export-only "FULL_PATH_TO_STUDY"
```

Completed cases are hash-checked; interrupted cases restart in new attempt
folders. Resume is at case boundaries, not solver timesteps, and requires the
same package and thread count. A stale `RUNNING.lock` can be removed only after
confirming that the previous process has stopped.

This reorganized copy has edited documentation and therefore a new package
fingerprint. Existing runs must be resumed from their original unchanged
package. The included `validation/` receipts describe earlier setup checks,
not completed scientific studies. Current cleanup checks are recorded in
[`../../docs/cleanup_checks.json`](../../docs/cleanup_checks.json).
