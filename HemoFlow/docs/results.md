# Results record

This record summarizes available experiment reports through September 28, 2026,
and the September 29 software check. Values are taken from retained summaries;
the October 5 documentation update did not rerun the scientific simulations.

## Outlet investigations

The earlier outlet candidate passed steady healthy checks but failed both
pulsatile healthy cases at N80 (0.750 and 0.300 µs). Its severe pilots were
therefore blocked. See [earlier readiness](evidence/earlier_outlet_readiness.json).

The later section-impedance revision passed all six required healthy cases,
the healthy pair comparisons, and the recorded acoustic screen. This allowed
short severe pilots. Those pilots still failed the 1% repeatability criterion:
60.0272% for the baseline and 53.7637% for the revision. See
[revision cases](evidence/outlet_revision_cases.json),
[readiness](evidence/outlet_revision_readiness.json), and
[pilot comparison](evidence/outlet_revision_pilot_comparison.json).

This history does not imply that Zou–He is inherently obsolete or that a
boundary change alone verifies severe flow. In this project, boundary choices
were evaluated within a particular solver, geometry, and forcing setup.

## Wall-shear measurement audit

The [wall audit](evidence/wall_audit_summary.json) reports 128 analytical
measurement checks passed, ten recorded cases processed, no skipped cases,
and no processing errors. Reproduction errors for the legacy shear calculation
were at floating-point roundoff scale. This supports implementation consistency;
agreement with a previous estimator is not independent evidence of physical
accuracy near a rasterized stenosis.

## Cycle settling and domain length

For the ten-cycle N80/0.300 µs severe runs, the worst relative L2 change over
the last three adjacent-cycle transitions was:

| Domain | Worst change | Criterion | Result |
| --- | ---: | ---: | --- |
| 90 mm | 3.176391% | <1% | Fail |
| 120 mm | 0.132516% | <1% | Pass |
| 150 mm | 0.039268% | <1% | Pass |

Sources: [90 mm](evidence/settling_L090.json),
[120 mm](evidence/settling_L120.json), [150 mm](evidence/settling_L150.json).
These compare pressure, flow, mean/peak absolute WSS, and both signed wall
profiles in specified regions; they are not full-field periodicity proofs.

The [90/120 mm comparison](evidence/domain_090_120.json) failed its
last-three-cycle domain screen. Far-region signed shear (x = 50–55 mm) differed
by 12.8828–12.9773% in cycles 8–10, exceeding 10%.

The subsequent [120/150 mm comparison](evidence/domain_120_150.json) passed
all declared final-three-cycle domain limits. Worst differences across cycles
8–10 were 0.021995% in mean pressure difference, 0.122092% in its waveform,
0.111324% in mean absolute WSS, 0.182428% in its waveform, 0.009286% in the
peak-WSS waveform, and 0.006912% in the flow waveform. Signed-shear space-time
L2 differences were at most 0.646179% over x = 20–55 mm and 2.984430% over
x = 50–55 mm. The exact JSON also retains throat and downstream-region scores.

The [study definition](evidence/length_study.json) records ten cycles, N80,
0.300 µs for severe cases, 16 threads, native initialization, and a positive
0.001 inlet skew for the first cycle. Physical geometry and measurement
locations are controlled. D75 means a 75% planar-gap reduction.

## Healthy planar reference

The [healthy L150 report](evidence/healthy_L150.json) passed its recorded
benchmark gate at N80 and Δt = 1.500 µs. Over 1.75–2.00 s, errors were
0.086599% for pressure difference, 1.208990% for mean absolute WSS, and
0.005601% for flow per unit depth. Reference values were 23.625 Pa between
x = 20 and 50 mm, 1.575 Pa WSS, and 0.0012 m²/s flow.

The maximum recorded lattice Mach number was 0.023386; density departure was
0.030687%; the relative link-based mass residual was 1.65 × 10⁻¹². These are
numerical-health diagnostics, not measurement uncertainties or physiological
validation errors. Fixed-radius wall fits are additional diagnostics and do
not replace the original WSS estimator in the acceptance record.

## Software checks and outstanding evidence

The [September 29 cleanup receipt](cleanup_checks.json) records 35 passing
unit/regression tests, a manufactured grid-analysis pipeline, headless recording
and CSV summary checks, offscreen rendering, frozen-source hashes, and Git
line-ending checks. Windows launchers and interactive Windows FPS were not
executed in that Linux check. No numerical code was modified for this
documentation update.

Completed results for the later full N80/N120/N160 grid batch are not included
in the reviewed evidence. A runtime estimate or an implemented test plan is not
a completed convergence study. Spatial refinement controls Δx/Δt; a separate
timestep study remains necessary. Independent physical validation and resolved
hematocrit/platelet physics have not been demonstrated.

The present defensible outcome is a reproducible numerical-sensitivity workflow,
successful healthy checks, and a passed fixed-grid longer-domain comparison
for the monitored observables. A general claim that the entire model is
verified or clinically validated would exceed this evidence.
