# Visual review and movie QA

Reviewed the trajectory traces and sampled static images. This is not a claim of continuous human playback or of testing a native application.

## Final comparison

`comparison.mp4` shows released, orbiting, and stiff orbiting springs at 60 Hz and epsilon ratio 1e8. The columns are full free-body CPD with P/L/E, an angular-row ablation with Kstar unchanged, and the previous independent-velocity LE method with beta_v=1. CPD uses 16 global updates; independent velocity uses 16 Gauss-Seidel sweeps. These are iteration budgets, not matched execution times.

The CPD panels show the relative coordinate r=x0-x1 of the free two-mass system with reduced mass 1. The independent-velocity panels show the internally equivalent anchored mass-1 model. This mapping is explicit in the movie and plot captions. All three primary scenes are unforced, so their plotted actual mechanical energy is kinetic plus elastic energy.

Inspected `unforced_poster.png` at 1 s, the independently decoded `decoded_t3.png` at 3 s, and the complete time traces in `unforced_relative_modes.png`. Text, spring length, speed, angular momentum, actual energy, and the independent-velocity coupling mismatch are readable. The two sampled movie frames match their stated times and the measured trajectories.

The traces and sampled frames distinguish three effects:

- Full CPD leaves the released spring almost motionless.
- Full CPD orbiting cases keep rotating but retain an almost constant stretched radius. This is frozen extension, rather than a frozen whole trajectory.
- Removing the angular row permits angular momentum growth and conversion of radial spring motion into rotation. Independent velocity also shows that conversion in the stiff orbiting case. Small energy error alone does not establish plausible spring dynamics.

See `../MODE_AUDIT.md` for measured radius windows and angular momentum. Phase and frequency agreement with a reference are not used as acceptance criteria.

## Supplemental inspection

Also inspected the free-body, anchored, and hard-rod family trace figures and the compliant released/kicked pendulum contact sheets. The figures identify anchored cases as adaptations and exact rods as extensions. Failed trajectories remain visible as finite prefixes. The gravity cases require the separate physical-energy decomposition and stopping diagnostics in `../../RESULTS.md`; the internal CPD target is not substituted for actual mechanical energy.

Thirteen contact sheets were generated. This receipt does not assert that every contact-sheet frame was individually reviewed.

## Decode receipt

- Final file: `comparison.mp4`.
- H.264, 1680 x 1320 pixels, 30 frames/s, 151 frames, 5.033333 s.
- `ffmpeg -v error -i results/figures/comparison.mp4 -f null -` completed with exit status 0 and no decoder errors.
- The 3 s image was extracted from the encoded movie with FFmpeg and inspected, rather than taken directly from the source animation canvas.

The final source scripts, numerical tables, plots, and movie are stable for archiving. No solver reruns were needed for this visual review.
