# Motion and energy reassessment

This review reads the existing campaign at
`/home/blacknoir/Projects/constrained-xpbd/experiments/block_orders_2026_11_01/results`.
It performs no solver runs. The figures show actual accepted positions,
derived velocity, length from the fixed pivot, and physical mechanical energy.
Continuous reference motion and RMS position comparisons are absent.

The main candidate setting is **120 Hz, 16 sweeps, epsilon / initial energy
above minimum potential = 10000**, shared by all six examples. The stress
comparison is **60 Hz, 16 sweeps, ratio 1e8**, with EL, LE, and JOINT.

## What was inspected

All twelve contact sheets were inspected at times 0, 0.25, 0.5, 1, 2, 3, 4,
and 5 seconds. The position paths, length traces, physical-energy traces,
and the release turning-point detail were also inspected. Movie posters
and decoded sampled movie frames support playback QA. This is **sampled
frame and numerical trace inspection**, not a claim of watching every frame
as a continuous animation.

Both final movies pass a complete `ffmpeg -v error -i <movie> -f null -`
decode with exit status 0. Each contains 151 H.264 frames at 30 fps and
lasts 5.033 seconds. The main movie is 1800 by 1200; the high-penalty movie
is 1800 by 1500. Final one-second posters and independently decoded
three-second frames were visually inspected.

The contact sheets label frozen trajectories and held failed prefixes in
red. The stress movie retains those cases; it does not replace them with
successful runs. Energy traces use accepted states only. Rejected runaway
trial energies remain in the campaign summary, rather than being presented
as accepted motion.

## Main candidate observations

| Example | EL maximum actual energy error | LE maximum actual energy error | Visible motion assessment |
|---|---:|---:|---|
| Released spring | 29.128% | 0.3545% | Both move. EL sticks near turns and restarts; LE has shorter apex pauses and strong restarts. LE is promising but has a remaining turning defect. |
| Orbiting spring | 0.00429% | 0.00429% | Ongoing, smooth-looking accepted paths in sampled frames; the two orders overlap closely. |
| Stiff orbiting spring | 0.01643% | 0.01643% | Ongoing smooth-looking rotation. Radial oscillation fades toward a nearly circular orbit while tangential speed rises. This altered motion deserves further physical checks. |
| Spring with gravity | 0.04230% | 0.04230% | Ongoing swinging and stretching with smooth-looking sampled paths; the orders overlap closely. |
| Released pendulum | 0.5298% | 0.1508% | Both continue swinging with approximately fixed length. LE has the smaller physical-energy change. |
| Kicked pendulum | 0.4711% | 0.1308% | Both continue swinging. LE has the smaller physical-energy change. |

Energy error percentages are `100 * abs(E(t) - E(0)) / (E(0) - minimum U)`.
They measure actual mechanical energy rather than relaxed-constraint closure.
No frequency or phase score enters this assessment.

`release_velocity_detail.png` displays every archived timestep near two
turns. EL has a roughly 0.479 m/s restart following a longer near-stationary
interval. LE also has recurring restart defects: the displayed later turn
has a roughly 0.413 m/s velocity increment between two 1/120-second steps.
Thus the LE release result should not be called completely smooth merely
because its energy drift is small.

## High-penalty observations

- Released spring: EL fails after 0.400 s and LE after 0.383 s. Their final
  accepted prefixes are held and labeled during playback. JOINT remains
  effectively frozen for five seconds, traveling about 0.000000375 m.
- Released pendulum: EL travels only about 0.00000916 m and JOINT about
  0.00000101 m. Both fail the motion criterion despite small energy change.
  LE continues swinging, with maximum actual energy error about 0.1543%.
- Kicked pendulum: EL and JOINT move to an apex, then become nearly
  stationary. Their five-second travel is only about 0.159 m. LE continues
  swinging, with maximum actual energy error about 0.1487%.
- Orbiting, stiff, and gravity springs: all orders have ongoing accepted
  motion and small actual energy changes in these runs. The stiff example
  again loses much of its radial oscillation while spinning up. Smoothness
  alone does not establish that this mode transfer is physically sound.

## Outputs

- `main_motion.mp4`: six examples, EL/LE side by side, shared candidate setting.
- `high_penalty_motion.mp4`: six examples, EL/LE/JOINT, including frozen and failed cases.
- `main_poster.png`, `high_penalty_poster.png`: one-second playback samples.
- `main_*_contacts.png`, `high_penalty_*_contacts.png`: twelve contact sheets.
- `main_paths_length_energy.png`, `high_penalty_paths_length_energy.png`:
  full accepted paths, lengths, and signed actual energy error percentages.
- `release_velocity_detail.png`: timestep-scale view of turning defects.
- `all_selected_runs.csv`: exact source run IDs and measured selection metadata.
- `reassess_visuals.py`: reproduction script, reading the canonical archive by default.

The evidence supports LE at the shared candidate setting as a useful next
baseline under the motion-plus-energy criterion. It does not establish
general robustness, eliminate turning defects, or validate mode transfer
in the stiff orbiting spring.
