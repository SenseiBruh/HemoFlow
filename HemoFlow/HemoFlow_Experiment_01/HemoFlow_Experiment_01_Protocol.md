# HemoFlow Experiment 01: Verification and Controlled Stenosis Sensitivity

Prepared 2026-09-07 for the inspected HemoFlow v3 Python source.

**Status: proposed protocol and pilot settings, not completed experimental results.**
The simulator has not been modified. No overnight run has been started. The
settings in this pack are compatible with the inspected configuration schema.
Check that your Windows copy has the same v3 command-line options before use.

## 1. What this study can establish

Research question: **At fixed inlet flow and plaque shape, how does realized
2D lumen-gap narrowing affect pressure drop, jet speed, and reversed-flow area
in a rigid, Newtonian planar channel?**

There are three distinct tasks:

1. **Verification:** compare an unobstructed steady channel with its analytical
   solution, then assess numerical resolution and time dependence.
2. **Sensitivity study:** vary stenosis severity with other inputs fixed.
3. **Validation:** compare against physical measurements with matched geometry,
   fluid properties, boundary conditions, and measurement definitions. This is
   a later step, not something achieved by collecting more simulated cycles.

NASA's verification/validation overview distinguishes mathematical checking
from agreement with physical experiments [1]. A patient CFD paper can provide
context or a computational comparison; numerical similarity alone is not
experimental validation.

The current solver is 2D, planar, rigid-wall, Newtonian, with a prescribed
parabolic inlet profile and fixed-density outlet. It is not a circular 3D
artery, carotid bifurcation, patient model, or clinical diagnostic tool.
The default systolic waveform is illustrative, not patient-measured.

## 2. What v3 actually saves

In the viewer, click **Record: OFF** to switch recording on. The solver must be
running for simulated time and logged samples to advance. Click again to stop.
Files go under `C:\HemoFlow Project\HemoFlow\exports\runs`.

| Output | Contents | Important boundary |
| --- | --- | --- |
| `*_timeseries.csv` | Time, epoch, cycle, phase, inlet/outlet mean velocity and flow per depth, pressure drop, global speed/shear/vorticity summaries, reversed-flow fraction, stability diagnostics, particle counts/response parameters | Compact scalar history, not full flow fields or individual particle trajectories |
| `*_metadata.json` | Run ID, creation time, units, requested recording interval, initial configuration, initial Womersley number, model scope | No Git commit, installed package versions, or automatic code checksum |
| `*_epochs.jsonl` | One JSON object per recorded solver epoch: configuration and stenosis report | A geometry/BPM reset starts a new flow solution; do not pool epochs as one stationary experiment |
| GUI **Save snapshot** | PNG, full-field NPZ, axial-profile CSV, stenosis-report JSON, resolved custom-layout settings JSON | One sampled state, not a continuous spatial history; the displayed PNG can be at a different interpolated display time from the latest numerical snapshot |
| Headless final export | `field.npz`, `metrics.json`, `stenosis_report.json`, `settings.json` | Written at normal completion or handled Ctrl+C, not after every sample or a guaranteed solver failure |

Each plaque's report includes ID, wall, center, length, requested narrowing,
continuous-profile narrowing, raster narrowing, minimum continuous/raster
lumen, throat position, throat-cell count, and circular area-equivalent
**reduction percentage**. The number of report entries gives the resolved
number of plaque objects; `config.plaque_count` alone is not reliable when
random-count mode is enabled. A symmetric `wall: both` object represents one
constriction supported on both walls.

Full shape parameters (`asymmetry`, `shape_exponent`) are explicit in a custom
configuration's `plaques` list. In unedited random mode, the logger can instead
have `plaques: null`: its seed and inputs allow reconstruction with the same
code, but the report does not explicitly store those two shape parameters.
For a random layout worth preserving, use **Save snapshot** and retain the
resulting custom-layout `.json`. The provided study settings avoid this issue
by specifying the shapes explicitly.

The epoch record's `time_s` is currently hardcoded to zero. It is not the
actual time recording was enabled; use the time-series timestamps to determine
coverage. A reset restarts simulated time, so identify samples by run, epoch,
and time together.

Data uses SI units except explicitly named `mm`, percentages, and particle
micrometres. Flow per depth is **m^2/s**, not mL/min or circular-vessel flow.
Pressure is relative CFD pressure, not cuff pressure. The CSV is flushed after
each row, which limits buffered loss, but does not guarantee survival of an
abrupt power failure. Final field exports are not restart checkpoints.

### Numerical meanings of selected columns

| Column | Actual definition in inspected v3 |
| --- | --- |
| `delta_p_pa` | Cross-section mean pressure at grid column 1 minus column -2, near the inlet/outlet; not a local trans-stenotic pressure measurement |
| `max_velocity_m_s` | Maximum speed magnitude over all fluid cells; not specifically the throat velocity |
| `backflow_percent` | Percentage of interior fluid cells with axial velocity below `-0.01 * mean_velocity`; first/last two columns excluded. This is reversed-flow area fraction, not reverse volumetric-flow percentage or recirculation length |
| `max_wss_pa` | Largest absolute approximate wall shear over both walls, including boundary columns; not a lesion-region mean |
| `mean_abs_wss_pa` | Unweighted mean of absolute approximate shear samples along both walls; not an arc-length-weighted lesion mean, local TAWSS, or OSI |
| `flux_error` | Signed fractional inlet/outlet volume-flux mismatch; multiply by 100 for percent. Keep instantaneous and time-integrated checks distinct |
| `mach` | Maximum lattice Mach number, used as a numerical compressibility indicator |
| `density_variation` | Maximum absolute lattice density departure from reference density 1; this is a fraction, not percent |

The recorder does not continuously preserve signed shear at each wall point.
It cannot retrospectively supply local TAWSS, OSI, wall-shear gradients,
recirculation length, or particle residence-time distributions from scalar
logs alone. Those require additional spatial/trajectory recording.
Headless v3 does not advance particles at all; zero particle counts there are
expected. This pack is a fluid experiment, not a biomaterial experiment.

## 3. Correct severity definitions

Use full nominal channel gap `H0` and minimum realized grid gap `Hmin`:

`gap narrowing (%) = 100 * (1 - Hmin / H0)`

For a clearly labeled hypothetical circular lumen:

`circular area reduction (%) = 100 * (1 - (Hmin / H0)^2)`

`circular area equivalent (mm^2) = pi * Hmin_mm^2 / 4`

The absolute hypothetical area is derivable, but not directly saved as an
area field in the current report. Calling the 2D gap a real 3D diameter would
overstate what is modeled. Do not use a circular area to calculate the
simulator's planar flow. Clinical reference-diameter definitions, including
NASCET, also need to be matched before clinical stenosis percentages are
compared; these settings are not a clinical NASCET measurement.

| Case | Nominal gap narrowing | Ideal remaining gap (mm) | Circular area-equivalent reduction | Geometry |
| --- | ---: | ---: | ---: | --- |
| S000 | 0% | 4.00 | 0% | Healthy |
| S025 | 25% | 3.00 | 43.75% | One symmetric constriction |
| S050 | 50% | 2.00 | 75% | Same shape and position |
| S075 | 75% | 1.00 | 93.75% | Same shape and position; severe exploratory case |

Plot measured raster narrowing on the results x-axis, and retain requested
narrowing as a separate column. For overlapping plaques, a local throat can
be caused by more than one plaque: report the combined lumen rather than
adding individual plaque percentages.

### Confirmed reporting issue to fix before the formal study

In `geometry.py`, the top-level field
`overall_realized_raster_diameter_narrowing_percent` is computed from the
continuous wall-profile gap, despite its name. Its associated overall area
equivalent uses that same continuous value. A geometry-only check with 37%
requested narrowing and 40 cells across returned 37% in that top-level field,
but 35% in the per-plaque raster field and actual minimum fluid-node count.

For the supplied centered, isolated constrictions, use the per-plaque raster
measurement and cross-check the minimum gap in the solid/fluid mask. For an
arbitrary geometry, compute the global minimum of fluid-node counts times
grid spacing. Per-plaque sampling is at a continuous-gap minimum and should
also be checked against a raster minimum when geometries are irregular.
Correct the naming/calculation and add regression tests before freezing a
formal-study software version. No correction has been applied in this pack.

## 4. Controlled inputs and run matrix

| Variable | Fixed setting |
| --- | --- |
| Density | 1060 kg/m^3 |
| Dynamic viscosity | 0.0035 Pa*s |
| Nominal channel gap | 4 mm |
| Channel length | 60 mm |
| Mean inlet speed | 0.30 m/s |
| Heart rate | 72 BPM, fixed |
| Pulse amplitude | +/-40% of mean inlet speed; use zero for steady verification |
| Pulse shape | `systolic`, illustrative; identical in every pulsatile case |
| Stenosis center | 30 mm |
| Stenosis length | 15 mm |
| Wall / asymmetry / shape exponent | `both` / 0.5 / 1.0 |
| Geometry randomization | Off; explicit custom plaques |
| Particles | Off; no role in fluid-only results |

The pressure required to sustain the prescribed inlet waveform is the
response. The setup does not predict how a patient's cardiac output changes
in response to stenosis. In a planar channel the throat mean-speed increase
for conserved flow per depth scales with inverse gap, not inverse circular
area. A 50% gap reduction therefore does not imply the fourfold mean-speed
increase of an ideal circular vessel with 75% area loss.

Start with coarse pilot runs. The candidate refinement family is 40, 80, and
160 cells across the unobstructed gap. The 75% case has only 10 throat cells
on the 40-cell grid, then 20 and 40 on the finer grids. These counts are not
proof of adequate WSS resolution. Geometry-only checks confirmed the table's
target minimum gaps on all 12 case/grid combinations.

Recommended sequence:

1. A short healthy steady logging pilot, followed by the existing unit tests.
2. Healthy steady verification at 40, 80, and 160 cells, running until steady
   diagnostics settle, not merely for an arbitrary elapsed PC time.
3. S000 and S050 pulsatile pilots at 40 cells to assess runtime and sampling.
4. After the logging/analysis corrections, all four pulsatile cases on all
   three grids: 12 runs. Run the severe case only if lower-severity checks
   succeed and stability/accuracy remain acceptable.
5. A separate fixed-grid solver-time-step sensitivity test and outlet-length
   sensitivity test. v3 does not expose independent time-step control or
   fixed physical pressure probes, so these need a small implementation step.

Changing `cells_across` also changes solver `dt`. Until time-step sensitivity
is isolated, call the family a **coupled space/time refinement study**, not
pure spatial convergence. Changing the recording interval does not change
solver `dt`. For any later GCI calculation, first assess convergence behavior;
do not assume second-order convergence or force an extrapolation through
non-monotonic results. NASA describes using three refinement levels to
estimate observed convergence and test the asymptotic regime [2].

## 5. Analytical healthy-channel verification

Disable pulsatility and remove all plaques. For fully developed planar
Poiseuille flow between parallel plates, using full gap `H`:

`u(y) = 1.5 * Umean * (1 - (2*y/H)^2)`

`delta_p = 12 * mu * Umean * Lmeasurement / H^2`

`abs(wall shear) = 6 * mu * Umean / H`

These follow from the steady Newtonian momentum equation with no slip at
`y = +/-H/2`. At the specified inputs:

- Centerline speed: 0.450 m/s (the mathematical center; evaluate the profile
  at actual cell centers when computing grid errors).
- Wall shear: 1.575 Pa.
- Pressure drop over the full 60 mm: 47.250 Pa.
- Pressure drop at v3's actual 40-cell-grid stations, separated by 59.8 mm:
  **47.0925 Pa**. Recompute the station separation on each finer grid.
- Steady flow per unit depth: 0.0012 m^2/s.

Do not compare the planar solver with the displayed circular-pipe references
of 126 Pa and 2.1 Pa. Those are different geometries.

Compute velocity-profile relative L2 error at the midpoint and errors in
pressure drop and flat-wall shear away from inlet/outlet columns. Retain
numerical and analytical values, not just a pass/fail flag. The existing
healthy-channel test uses a shorter 16 mm channel, 0.1 m/s inlet, and 32
cells; passing it does not replace verification at the study inputs.

For later pulsatile verification, use an analytical oscillatory **planar**
solution with compatible forcing, inlet profile, and boundaries. Reporting
Womersley alpha near 3.02 is a parameter check, not evidence that that
analytical velocity/pressure solution has been reproduced.

## 6. Sampling, warm-up, and stopping

At 72 BPM, one simulated cycle lasts 60/72 = 0.833333... seconds. Twenty
cycles are 16.6667 simulated seconds, which can take far longer on the PC.
Longer runtime is useful to remove startup effects or characterize persistent
unsteadiness; repeated deterministic cycles are not independent experiments.

The GUI currently requests one sample per 0.05 simulated seconds, about
16-17 samples per beat at most. For phase-resolved pilot collection use
headless `--record-every 0.005`. v3 checks for samples only after solver
batches of up to 200 steps. Actual sample gaps can be larger than requested,
and depend on resolution and severity. The final partial batch can give a
different interval. **Always use the saved `time_s` values.**

For the supplied pulsatile cases at 40 cells, the code predicts ordinary
sample intervals of approximately 0.0060, 0.00714, 0.00635, and 0.00556
seconds for S000, S025, S050, and S075. These are scheduling calculations,
not measured solver results. Requesting 0.005 does not guarantee 0.005.
Start with at least roughly 100 actual samples per cycle for scalar trends,
then test finer output sampling; rapid shear/vorticity extrema may need more.
Phase interpolation cannot recover frequencies that were never sampled.

For the formal pulsatile runs, propose 20 total cycles initially:

1. Record startup too, but initially exclude cycles 0-9 from analysis.
2. Compare complete phase-aligned waveforms and cycle means over consecutive
   late cycles. Extend the run if they are still changing.
3. Analyze at least five complete late cycles after an adequate settling
   window; retain the exact accepted start/end times in the run register.
4. Do not retain a partly recorded last cycle just because it has a
   `cycle_index`. Cycle 20 can contain a single endpoint sample after a
   requested 20-cycle run. Require samples bracketing both cycle boundaries
   for time integration; otherwise exclude that cycle.

**Suggested project screening criteria, not accepted clinical standards:**

- Healthy steady profile/pressure/shear errors below 5%, with errors improving
  under refinement; aim tighter once the numerical behavior is established.
- Last-cycle pressure and speed waveform/mean changes below 1% over three
  successive cycle comparisons, with a fixed nonzero scale to handle zeros.
- Use an absolute tolerance such as 0.1 percentage point for near-zero
  reversed-flow-area changes, rather than division by a zero healthy value.
- Time-integrated inlet/outlet flow-per-depth imbalance below 1% over the
  accepted window; also inspect the instantaneous mismatch and density
  variation. Signed errors must not be hidden by cancellation.
- Monitor lattice Mach and density variation. Prefer low-Mach behavior
  (roughly below 0.1) and small density departures (roughly below 1%) as
  conservative screening targets; if missed, investigate sensitivity to
  numerical scaling. Passing the preview stability stop is not an accuracy
  certificate, and refinement alone need not reduce Mach in this code.
- For preselected pressure and speed outcomes, investigate changes above 5%
  between the two finest grids. Treat raw WSS maxima as exploratory until
  their wall/region definition and convergence are established.

Choose and document criteria before interpreting outcome differences. If a
wake remains nonperiodic, do not keep waiting for forced beat-to-beat identity
or quietly relax criteria. Check numerical causes, then use a documented
longer-window statistical-stationarity approach if physically justified.

The bundled `analyze_recording.py` groups by epoch/cycle and returns JSON.
It currently uses arithmetic sample means, includes startup and incomplete
cycles, and its standard deviation describes variation **within a cycle**.
It is useful for inspection, not a final publication reducer. Formal analysis
needs cycle-boundary interpolation, warm-up/partial-cycle filtering, and
time-weighted integration. Within-beat sample SD is not uncertainty between
independent experiments or patient variability.

## 7. Running pilots on your Windows PC

The ZIP contains a folder named `HemoFlow_Experiment_01`. Extract that folder
inside `C:\HemoFlow Project\HemoFlow`. This adds only study settings and this
protocol; it does not replace the application. These commands use the project
Python directly, so environment activation and VS Code interpreter selection
are unnecessary. Run one command per line in PowerShell.

First confirm the correct directory, environment, and v3 options:

```powershell
cd "C:\HemoFlow Project\HemoFlow"
Test-Path .\main.py
Test-Path .\.venv\Scripts\python.exe
.\.venv\Scripts\python.exe main.py --help
```

Both `Test-Path` commands should print `True`; help should include `--record`,
`--cycles`, `--record-every`, and `--headless`. If a required module is missing,
install the project's requirements using that same Python:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

### Healthy logging pilot, no pack extraction required

```powershell
.\.venv\Scripts\python.exe main.py --defaults --mode healthy --steady --headless --seconds 1 --record --record-every 0.005 --record-dir experiments\pilot_healthy_001\raw --output experiments\pilot_healthy_001\final
```

This is a short functional pilot, not proof that one simulated second is
enough to converge. Check for CSV rows and final exports. Stop with Ctrl+C
in this terminal if needed; wait for final exports before closing it.

### Fixed S050 pulsatile pilot using the provided settings

```powershell
.\.venv\Scripts\python.exe main.py --config HemoFlow_Experiment_01\configs\S050.json --headless --cycles 2 --record --record-every 0.005 --record-dir experiments\pilot_S050_N040_001\raw --output experiments\pilot_S050_N040_001\final
```

Use `S025.json`, `S075.json`, or `S000.json` for other fixed cases and change
the run folder name to match. Append `--cells 80` or `--cells 160` for grid
refinement and update `N040` in the folder name. Append `--steady` for the
steady verification task. Do not use random mode or drag geometry mid-run.
In a custom configuration, the active plaque height is set by the plaque's
own `severity_percent`; changing only top-level `stenosis_percent` will not
resize that custom plaque. Use the supplied separate case files unchanged.

After the readiness corrections and pilot review, use `--cycles 20` for an
initial formal run. Every rerun needs a new folder suffix: time-series names
are timestamped, but reusing `--output` overwrites final field/metrics files.
Launching again starts from the initial condition, not the prior final field.
Twenty cycles is an initial duration choice, not automatic convergence.

For a deliberately open-ended pilot, replace `--cycles 2` with
`--continuous`; still specify `--headless --record` and unique output paths.
Finite runs are easier to audit than unbounded collection. Keep the PC plugged
in, prevent sleep for the session, leave the terminal running, and check disk
space and a short pilot before leaving. Screen-off is fine; PC sleep pauses
computation. Do not disable Defender or other security controls for this.

Important current CLI details: `--record` and `--record-every` are applied by
the headless route; the GUI recording button uses its own hardcoded interval.
The inspected headless route also does not call thread tuning, so `--threads`
is not presently evidence that a headless run used that many solver threads.
Record observed runtime and estimate cost before committing to fine grids.

### Inspect a per-cycle summary

After the S050 pilot has completed, this PowerShell code selects the newest
time-series file in that one pilot folder and invokes the existing reducer:

```powershell
$runCsv = Get-ChildItem "experiments\pilot_S050_N040_001\raw" -Filter "*_timeseries.csv" | Sort-Object LastWriteTime -Descending | Select-Object -First 1
.\.venv\Scripts\python.exe analyze_recording.py "$($runCsv.FullName)" --output experiments\pilot_S050_N040_001\cycle_summary.json
```

This is JSON, not an automatic Excel workbook. Keep the raw CSV as well.

## 8. Data organization and Excel

Use one experiment family folder with one immutable subfolder per run.
Within each run keep `raw/`, `final/`, `analysis/`, and a short notes file.
Keep inputs in `configs/` and final figures/report in their own folders.
Do not put files into `.venv`, overwrite raw CSVs in Excel, or discard failed
runs. Retain failed-run status and reason alongside successful results.
Back up completed runs and their matching code/configuration to another
location. A Git commit alone is local until it is pushed or otherwise backed up.

Recommended run-register columns (one row per run/epoch):

| Group | Columns |
| --- | --- |
| Identity | Experiment ID, run ID, epoch, status, UTC start, operator, Git commit, dirty/clean code status |
| Environment | Python version, package-lock filename, OS, CPU, actual thread count if known, wall runtime |
| Geometry | Config filename, geometry hash if available, seed, resolved plaque count, wall, center, length, asymmetry, exponent, nominal gap |
| Severity | Requested gap reduction, continuous gap reduction, measured raster gap reduction, circular area-equivalent reduction, minimum raster gap, throat cells |
| Numerics | Grid size, dx, dt, actual pressure-station coordinates, requested and measured logging intervals |
| Inclusion | Total simulated duration, excluded warm-up time, accepted cycle IDs, partial cycles excluded, convergence/QC decision |
| Results | Time-mean/peak pressure drop, time-mean/peak global speed, reversed-flow-area statistics, exploratory WSS statistics, flux/Mach/density checks |
| Traceability | Raw CSV, metadata, epochs, final field, analysis file, figures, limitations/failure notes |

For a first Excel workbook, use sheets named `Run Register`, `Geometry`,
`Raw Data`, `Cycle Results`, and `Figures`. Import the finished CSV using
**Data > Get Data > From File > From Text/CSV** [3]. Separate CSV rows from
different runs must retain their run IDs; `epoch = 0` occurs in many files.
The per-run metadata can also be imported with Power Query's JSON connector.
The `.jsonl` file contains one JSON object per line and needs a line-aware
import or preprocessing, not a normal single-object JSON import.

Simple Excel geometry formulas, assuming A2 contains nominal gap in mm and
B2 contains measured remaining gap in mm:

- Gap reduction as a fractional Excel percentage: `=1-B2/A2`
- Circular area-equivalent reduction as an Excel percentage: `=1-(B2/A2)^2`
- Hypothetical remaining circular area in mm^2: `=PI()*B2^2/4`

Format the first two as Percentage. In contrast, raw CSV/JSON fields ending
in `_percent` already contain numbers on the 0-100 scale; divide by 100
before applying Excel percentage formatting, or retain a plain numeric
column with `%` in the heading. Do not accidentally display 25 as 2500%.

Use time-weighted cycle means: integrate by the trapezoidal rule over actual
timestamps with cycle boundaries included, then divide by the exact cycle
duration. Use a PivotTable only after that reduction, or for clearly labeled
exploratory sample summaries. Keep spatially local WSS, whole-domain WSS,
temporal variation, grid differences, and independent geometry variation
separate in column names and plots.

Suggested final figures:

1. Analytical and simulated healthy velocity profiles, plus an error table.
2. Grid/time refinement of pressure drop and speed for each severity.
3. Phase-aligned pressure drop and speed over a representative accepted beat.
4. Cycle-mean pressure drop and reversed-flow area versus realized gap
   narrowing, with circular area equivalents explicitly distinguished.
5. Fixed-scale flow/axial-velocity maps for selected cases at matched phase;
   do not silently autoscale each panel to a different velocity range.

A global peak shear at one moment must not be plotted against a literature
lesion-averaged WSS as though they were the same measurement.

## 9. Git and reproducibility

Missing historical commits does not invalidate your project. It does mean
past changes are harder to reconstruct. Start a truthful baseline now:
`Import existing HemoFlow v3 prototype`. Do not invent/backdate a development
history. Keep recovered v1/v2 ZIPs as labeled archival artifacts, not pretend
they were previously committed. Tags can mark important study snapshots [4].

Since you already have a GitHub repository, first locate its existing clone.
At your project location, inspect:

```powershell
git status
git rev-parse --show-toplevel
git remote -v
```

If Git reports no repository, use your existing clone (or clone the existing
GitHub repository into a separate folder, then bring the current source into
it after reviewing differences). Do not blindly run `git init` in a nested
folder or replace another version. No remote URL has been assumed here.

Review or create a `.gitignore` covering `.venv/`, `__pycache__/`, `*.pyc`,
Numba caches (`*.nbc`, `*.nbi`), and large raw/field-export directories. Keep
study configurations, protocol, tests, analysis code, concise result tables,
and selected figures tracked. Store bulky raw data separately with checksums
and documented locations, or use an appropriate large-file/release workflow.
Do not indiscriminately ignore every CSV/JSON: some are essential study inputs
and small results. Never commit credentials or private patient information.

Within the correct repository, stage and review only the intended files,
then commit and tag the prototype baseline:

```powershell
git add -- *.py requirements.txt README.md ROADMAP.md VALIDATION.md VERSION.txt tests
git diff --cached --stat
git status
git commit -m "Import existing HemoFlow v3 prototype"
git tag -a v3-prototype-baseline -m "Prototype before formal verification study"
git rev-parse HEAD
```

Stage the reviewed ignore file and study configs/protocol as well. If this
repository already tracks HemoFlow differently, adapt the stage paths to its
existing organization; the point is deliberate review, not `git add .` over
a folder containing environments and hours of output. Inspect existing staged
changes before committing. Do not reuse or force-overwrite an existing tag.
Push the intended branch and study tag to the verified remote normally; do
not force-push. After fixing the readiness items, make a separate study-ready
commit/tag and record its commit hash in every run's notes.

Capture environment information in each experiment's provenance notes:

```powershell
.\.venv\Scripts\python.exe --version
.\.venv\Scripts\python.exe -m pip freeze
git rev-parse HEAD
git status --short
```

Keep the exact launch command and settings with these records. The current
logger does not capture them automatically. Do not change the solver midway
through a comparison series; a corrected solver version starts a new series.

## 10. Literature and presentation

There was a citation mismatch earlier: `PMC7012768` is Filonova et al.'s
**Verification of the coupled-momentum method with Womersley's Deformable
Wall analytical solution**, not the carotid stenosis paper [5]. The quoted
107 +/- 73 and 19 +/- 14 dyn/cm^2 values belong to **Computational fluid
dynamic characterization of carotid bifurcation stenosis in patient-based
geometries**, Schirmer and Malek, DOI `10.1002/brb3.25`, `PMC3343298` [6].
Correct this bibliography before presenting the project.

Those WSS values convert to 10.7 +/- 7.3 Pa and 1.9 +/- 1.4 Pa because
1 Pa = 10 dyn/cm^2. They describe different, patient-based 3D geometries;
matching them is context, not a pass criterion for an idealized planar geometry.
Do not adjust free plaque geometry until a number agrees and then call that
independent validation. Record the exact study region, averaging method,
cardiac phase, fluid model, and geometry before making any comparison.

For a future medical-device-focused validation project, the FDA provides
experimental velocity/pressure data and CAD/test conditions for its nozzle
and pump benchmarks [7]. Reproducing an appropriate benchmark would require
matching its geometry, including dimensionality, and boundary conditions;
the current straight planar model is not automatically that benchmark.

An honest report title is: **Verification and stenosis-sensitivity analysis
of a 2D pulsatile flow simulator**. Include objective, equations/assumptions,
boundary conditions, numerical tests, actual geometry measurements, results,
uncertainties, code/data availability, and limitations. A useful professor
discussion can happen before full validation; show what works, what has been
tested, and which question needs laboratory or 3D-CFD support. Use a LinkedIn
post to link the reproducible report/repository after checks are complete,
not to claim clinical prediction from a visual match.

After this first study, a second geometry study could fix realized narrowing
and compare symmetric versus eccentric shapes. A later paired random-shape
study could generate a prespecified set of shapes, save each explicitly, and
scale the same shapes across severities. Grid adjustment may change realized
narrowing, so remeasure it every time. Shape replicates are different from
repeated beats; do not combine their uncertainties as though they were the
same source of variation. Leave particle material/size effects for a separate
protocol with trajectory outputs and a validated particle model.

## 11. Readiness checklist and audit status

Before the formal overnight series:

- [ ] Commit the current prototype baseline and back it up.
- [ ] Correct overall continuous/raster stenosis reporting; add tests for
  grid rounding, shifted/asymmetric plaques, and overlapping walls.
- [ ] Save resolved plaques for every geometry epoch, geometry hashes, source
  version, numerical scales, actual thread count, and environment provenance.
- [ ] Add deterministic sample scheduling and signed per-wall/probe recording
  if local WSS, local pressure loss, or recirculation length will be endpoints.
- [ ] Add complete-cycle filtering, warm-up selection, time-weighted means,
  and tidy cycle/geometry CSV export with tests; do not conflate phase SD with
  between-run uncertainty.
- [ ] Add an independent numerical-time-scale control and fixed physical
  pressure stations for temporal/domain checks.
- [ ] Pass the healthy benchmark and review short stenosis pilots.
- [ ] Freeze/tag the corrected version and then collect the formal matrix.

Checks performed while preparing this pack: inspected the relevant source,
parsed the existing short CSV smoke log, reproduced the 37%/35% reporting
discrepancy, and checked all 12 severity/grid geometries without running the
fluid solver. The current full unit suite could not execute in this workspace
because Numba is unavailable. No claim of a new full solver-test pass or
completed fluid experiment is made. Run the tests in your working Windows
environment with:

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
```

## Sources

Accessed 2026-09-07. Protocol settings and screening tolerances above are
proposals for this project, not prescriptions from these sources.

1. NASA, Overview of CFD Verification & Validation:
   https://www.grc.nasa.gov/www/wind/valid/tutorial/overview.html
2. NASA, Examining Spatial (Grid) Convergence:
   https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html
3. Microsoft, Import data from data sources (Power Query):
   https://support.microsoft.com/en-us/excel/import-data-from-data-sources-power-query
4. Pro Git, Recording Changes and Tagging:
   https://git-scm.com/book/en/v2/Git-Basics-Recording-Changes-to-the-Repository
   https://git-scm.com/book/en/v2/Git-Basics-Tagging
5. Filonova et al., DOI 10.1002/cnm.3266:
   https://pubmed.ncbi.nlm.nih.gov/31617679/
6. Schirmer and Malek, DOI 10.1002/brb3.25:
   https://pmc.ncbi.nlm.nih.gov/articles/PMC3343298/
   https://pubmed.ncbi.nlm.nih.gov/22574273/
7. FDA, Benchmark dataset for validating CFD simulation of blood flow through
   generalized medical device geometries:
   https://cdrh-rst.fda.gov/benchmark-dataset-validating-computational-fluid-dynamic-cfd-simulation-blood-flow-through
