# Reassessment using energy and plausible motion

The shared setting **LE, 120 Hz, 16 sweeps, epsilon ratio 10,000** is a useful
practical candidate across all six saved scenes. Every rollout finishes five
seconds, energy error stays below 0.36% of initial energy above the potential
floor, and every scene continues to move. Raising the penalty to 100 million
does not improve that cross-scene result: it produces nearly stationary release
states and persistent pendulum stalls, despite smaller energy errors.

Two concerns remain relevant under an energy-and-plausibility criterion:
small pauses and sharp restarts at spring turning points, and smooth unforced
self-spinup/circularization in the central-force scenes. These are not phase,
frequency, reference-RMS, or Hessian-curvature rejection criteria. A user may
accept the smooth rotating motion while still wanting the artificial conversion
of radial motion into spin called out explicitly.

This is a read-only structural screen of existing saved trajectories, not a new
solve or a direct video inspection. It audits 72 existing high-penalty runs
(both orders and the joint reference, 60/120 Hz, 16 sweeps, ratios 10,000 and
100 million). The reproducible script is `audit_motion.py`; its complete metrics
are in `motion_metrics.json`. No reference trajectory columns are read.

## Shared practical candidate

All energy figures below are maximum absolute physical-energy change from the
run-initial value, divided by initial energy above the common potential floor.
They are approximate conservation over the saved five seconds, not exact
conservation or a long-duration guarantee.

| Scene, all LE/120 Hz/16 sweeps/ratio 10,000 | Maximum energy error | Observed structural behavior |
|---|---:|---|
| Spring released from rest | 0.3545% | Ongoing bounded oscillation, radius 0.75003–1.25; small near-zero-speed endpoint pauses followed by sharp restarts |
| Orbiting spring | 0.00429% | Smooth ongoing rotation and radial breathing; no angular-direction reversal; angular momentum rises 45.9% |
| Spring with gravity | 0.04230% | Continues swinging and stretching; no sustained stall or endpoint jump found |
| Pendulum released from rest | 0.1508% | Continues swinging; length error at most 0.0000794 m on a 1 m rod |
| Pendulum with initial kick | 0.1308% | Continues swinging through repeated reversals; length error at most 0.0000420 m |
| Stiff orbiting spring | 0.01643% | Smooth early radial motion settles into almost circular fast rotation; angular momentum rises 272.6% |

The spring-with-gravity candidate's endpoint speed remains between 0.039 and
2.117 m/s. Both pendulum candidates retain substantial motion in the final
second: maximum speeds approximately 1.999 and 1.728 m/s, respectively. They do
not approach the long quiet tails seen in the frozen alternatives.

The largest accepted endpoint increments for the shared candidates are modest:
approximately 0.0208 m for release, 0.0107 m for the orbit, 0.0176 m for gravity,
0.0167/0.0144 m for the pendulums, and 0.0272 m for the stiff spring. These are
physical endpoint increments; an internal block correction may be larger and
then be undone before the accepted endpoint. An internal excursion alone is not
treated as a visible-motion failure here.

## Release pauses and restarts

For LE/120 Hz/16 sweeps/ratio 10,000, spring release repeatedly reaches an
extreme with almost zero speed and then restarts unusually strongly. For
example, at 0.333333 s its radius is 0.7500307 m and velocity is
+0.0002598 m/s; one 1/120 s interval later its velocity is +0.4135931 m/s.
The finite-difference acceleration is 49.60 m/s².

The same pattern repeats at subsequent turning points. With this spring's
observed energy bound, its physical force magnitude is at most 25.0 N
for unit mass. Therefore the restart increment is stronger than any exact
spring trajectory with that energy could produce over the same interval.
This is a small discrete pause/restart signature, not just a different period.
Its visual severity at 120 Hz should be judged from playback rather than
inferred from the acceleration number alone. The candidate still makes broad,
bounded oscillations and keeps its energy approximately constant.

## Smooth spinup and loss of radial motion

With no gravity, a radial anchored spring exerts no torque about its anchor.
The shared orbiting-spring candidate nevertheless changes angular momentum
from 0.88 to 1.28365 over five seconds. Its angular direction never reverses;
the concern is smooth artificial spinup, not noisy or erratic rotation.

The stiff candidate makes the mechanism visually more distinctive. Angular
momentum rises from 0.88 to 3.27846 while endpoint speed grows from an initial
0.8 m/s to approximately 3.24493 m/s in the final second. It remains smooth:
the overall speed range is 0.8–3.26171 m/s, and the largest early acceleration
is 101.67 m/s², consistent with the strong initial stiff-spring force scale.

Its radius settles near 1.01043 m. From 1.333 s onward the remaining radius
range is below 0.001 m; from 1.925 s onward it is below 0.0001 m. During the
last second the radius varies by only 1.02e-7 m. Stored elastic/radial motion
has therefore become almost entirely steady rotation. The total energy looks
excellent while the central-force angular momentum is not preserved.

The energy multiplier can produce such a change because the stationarity row
contains

\[
\eta\left[\frac{M(x-x_n)}{h^2}+\nabla U(x)\right].
\]

With reconstructed velocity, its kinetic part is proportional to endpoint
velocity. It has a tangential component even when the spring force itself is
purely radial. This explains why conserving the one total-energy number does
not guarantee conservation of the orbit's radial-versus-angular character.

All three sufficiently iterated orders produce effectively the same orbit and
stiff behavior at these high penalties. Choosing EL instead of LE therefore
does not remove this structural feature. It can still be an acceptable-looking
animation if smooth bounded motion is the immediate goal; it should not be
described as faithful central-force angular behavior.

## Persistent stalls hidden by excellent energy numbers

The ratio-100-million, 120 Hz, 16-sweep results show explicit non-equilibrium
stalling:

| Scene and order | Final-second maximum speed | Remaining physical drive |
|---|---:|---:|
| Spring release, JOINT | 1.50e-7 m/s | Approximately 25 m/s² inward acceleration |
| Pendulum release, EL | 3.78e-6 m/s | Approximately 5.94 m/s² tangential gravity |
| Pendulum release, JOINT | 4.04e-7 m/s | Approximately 5.94 m/s² tangential gravity |
| Kicked pendulum, EL | 8.43e-8 m/s | Approximately 5.20 m/s² tangential gravity |
| Kicked pendulum, JOINT | 3.45e-7 m/s | Approximately 5.20 m/s² tangential gravity |

The kicked pendulums climb to angle about 0.55922 radians and then stop near
their first apex. Their speeds remain below 0.0001 m/s from approximately
0.325 s (EL) and 0.292 s (JOINT) onward, even though that angle has a nonzero
gravity torque. This is a physically unsupported rest state, independent of
what period or phase the user prefers.

The same extreme penalty also compromises spring release in the split orders:
at 120 Hz, LE stays almost fixed near its initial stretched radius, with only
0.0000321 m radius variation in the final second; EL shows small jitter around
the stretched state instead of a broad oscillation. At 60 Hz the corresponding
EL and LE release runs hit the production runaway guard. These variants do not
provide a useful all-scene setting.

LE avoids persistent stalling in both saved pendulum scenes, even at the higher
penalty. This is a genuine ordering benefit under the user's criterion. It
does not rescue the high-penalty spring-release scene.

## Pendulum chord velocity and visible length error

These methods derive endpoint velocity from the position chord. Even with
exact endpoint lengths,

\[
x^T\frac{x-x_n}{h}=\frac{\|x-x_n\|^2}{2h}>0
\]

for a nonzero step on the circle. Thus the returned velocity has an outward
radial component. In the shared LE setting, maximum radial components are
0.01667 and 0.01245 m/s for the two pendulum scenes, compared with maximum
speeds approximately 2.00 and 1.73 m/s. This is a small finite-step kinematic
defect; the rendered endpoint positions still trace a circle to the small
length errors above. It is not, by itself, evidence of a visible outward spiral.

The saved trajectories satisfy the general chord identity, including their
small length changes, to numerical precision. No endpoint velocity projection
or independent velocity is present. If the returned velocity later drives
collisions or other velocity constraints, this defect deserves separate
attention; the current saved position motion remains bounded.

## Practical interpretation

The LE/120 Hz/16-sweep/ratio-10,000 setting deserves direct playback as the
current cross-scene candidate. Energy is approximately conserved and none of
the six scenes stalls. The smooth self-spinup and small turning-point
pause/restart signatures are the meaningful remaining plausibility questions.

The extreme-penalty frozen alternatives fail the motion criterion even when
their energy error is tiny. Conversely, smooth orbiting candidates should not
be rejected solely because their phase, period, or trajectory differs from a
continuous reference. This audit makes no claim about longer than the existing
five-second horizon.
