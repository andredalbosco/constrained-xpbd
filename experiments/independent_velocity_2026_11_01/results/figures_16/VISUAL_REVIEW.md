# Static review of the recommended sixteen-sweep setting

Setting: independent velocity, beta_v=1, 60 Hz, 16 sweeps, energy-penalty
ratio 1e8. These figures reuse existing trajectories; no solver was rerun.
The four-sweep movie remains a separate artifact and is not presented as
playback of this setting.

The full path/length/speed and actual/chord-energy/coupling-factor figures,
plus all six contact sheets, were visually inspected. Contact sheets sample
0, 0.25, 0.5, 1, 2, 3, 4, and 5 seconds. This is **sampled static inspection**,
not continuous video viewing or proof that no intermediate artifact occurs.
The energy caption uses the agreed `E(x,v)` notation.

- **LE keeps moving in all six examples.** Its final-second travel is at
  least 1.064 m across these scenes. The sampled paths show no new freeze
  or large spatial discontinuity compared with the four-sweep setting.
  The two pendulums keep swinging, and the gravity spring keeps swinging
  and stretching.
- **Turning corrections remain.** Released spring has recurring sharp
  velocity/chord mismatch and coupling-factor pulses. More sweeps greatly
  reduce actual-energy drift but do not remove that kinematic defect.
  The maximum accepted chord-velocity increment is about 0.672 m/s.
- **Stiff spring still becomes nearly circular.** Radial oscillations fade
  while speed rises toward about 3.24 m/s. This persists at sixteen sweeps;
  better energy closure does not validate the radial-to-spin transfer.
- **JOINT still has delayed motion and apex stalls.** Released spring
  stays near its initial pose for roughly the first 2.3 seconds, makes one
  excursion toward the compressed state around 2.5–3 seconds, then freezes
  there. Total travel is about 0.5 m and final-second travel only
  0.000000074 m. The two JOINT pendulum cases likewise alternate delayed
  excursions with long near-stationary intervals; their longest stalls are
  about 1.58 and 1.62 seconds. Coupling factors cross or approach zero.

LE's maximum actual energy error across the six scenes is approximately
**0.000705%**, attained by released pendulum. Its chord-energy diagnostic
can differ: the stiff spring reaches about **15.887%** and velocity/chord
difference about **0.435 m/s**. The formulation permits this difference;
the traces expose it for judging the accepted motion rather than silently
reconstructing velocity.

The sixteen-sweep setting is a useful next baseline for sustained motion
and stored-velocity energy. It remains a candidate with turning and mode-
transfer limitations, rather than a general stability result.
