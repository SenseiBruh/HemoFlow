# Model definition

## Fluid and geometry

The solver uses a D2Q9 lattice Boltzmann formulation with stress regularization
and separate relaxation of shear and other moments. The fluid is Newtonian and
weakly compressible. Walls are rigid and represented by a rasterized mask with
halfway bounce-back.

The default desktop channel has a 4 mm gap and 60 mm length. Default density
is 1060 kg/m³, dynamic viscosity is 0.0035 Pa·s, and nominal mean inlet speed
is 0.30 m/s. The channel has no resolved out-of-plane dimension. Flow reported
per unit depth has units of m²/s.

## Inlet and outlet

The inlet profile is normalized to the requested discrete mean velocity.
Pulsatility modulates velocity at a fixed heart rate. At 72 BPM the period is
60/72 = 0.8333 s. A 40% modulation scales velocity around its nominal value;
it does not change the heart rate by 40%. The systolic waveform is illustrative,
rather than a patient measurement.

The desktop solver has a velocity inlet and fixed-density outlet, using
regularized stress extrapolation by default. Its selectable Zou–He treatment
is retained for baseline comparisons. The two recent study packages instead
use a frozen solver with a section-impedance outlet. That outlet remains an
experimental numerical treatment; its parameters are not a calibrated
physiological downstream circulation.

## Pressure and wall shear

Pressure is derived from lattice density relative to a numerical reference.
Reported differences must use matching physical stations when comparing
domains of different lengths. The model does not prescribe a patient's
absolute arterial pressure.

Wall shear uses velocity gradients near the rasterized wall. Its accuracy is
sensitive to wall resolution and gradient estimation, particularly near a
stenosis. Signed profiles retain flow direction; the magnitude alone cannot
describe reversal. Numerical studies record the estimator and physical
measurement locations.

For a steady, fully developed planar channel of gap H, mean speed U, viscosity
μ, and measurement separation L:

```math
\Delta p = \frac{12\mu UL}{H^2}, \qquad |\tau_w| = \frac{6\mu U}{H}.
```

The circular-pipe reference functions in `physics.py` use different coefficients.
Those values are reference calculations, not the planar solver's acceptance
targets. Likewise, the reported circular-vessel Womersley number is a comparison
scale and not evidence of agreement with a pulsatile planar analytical solution.

## Stenosis measurement

Requested severity is a reduction of the channel gap. Rasterization can change
the realized minimum gap. Geometry reports retain both the requested and
realized values.

For a circular comparison label only, a linear reduction fraction s corresponds
to an area reduction of 1 − (1 − s)². A 50% linear reduction therefore has a 75%
circular-area equivalent. A 75% linear reduction has a 93.75% equivalent.
The solver does not calculate an actual three-dimensional lumen area or a
clinical carotid stenosis classification.

## Particles and display

Tracer particles follow the computed velocity field. The inertial option adds
a Stokes relaxation time for dilute spherical probes with one-way coupling.
Neither option resolves cell deformation, cell collisions, adhesion, or the
effect of particle concentration on the fluid.

The display interpolates published snapshots to improve animation. Its FPS and
smoothing settings affect rendering. Recorded solver measurements are taken
from numerical states. Changing a plaque in the viewer starts a new flow
solution for the edited geometry.
