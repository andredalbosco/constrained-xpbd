# Symplectic predictor plus an energy multiplier: bounded, but usually worse motion

The present correction is not ready as a general replacement for the
symplectic predictor. A displacement cap controls the large energy updates,
but the corrected trajectories usually lose accuracy and angular momentum.
The near-circular orbit is the main favorable exception in this suite.

## Precisely what was tested

All new runs use compliant springs, including gravity-driven spring pendulums.
The symplectic baseline is

    lambda_force = -h² alpha^-1 C(x_old)
    x_star = x_old + h v_old + h² a_gravity + M^-1 J_old^T lambda_force
    v_star = (x_star-x_old)/h.

For energy-only methods, one energy row is linearized at x_star:

    H=E(x_star)-E_initial
    g=grad E(x_star), D=g^T M^-1 g
    delta_eta=-H/D, delta_x=M^-1 g delta_eta.

The physical energy is the same endpoint energy used previously, including
elastic and gravitational potential, with velocity reconstructed from the
accepted displacement. There is no reservoir, velocity rescaling, nonlinear
root solve, or repeated projection. The same initial energy target is used
throughout each run.

The next step's force cache is recomputed from the accepted position. This is
essential for the energy-only family: carrying an unchanged lambda after the
energy correction would reuse a stale elastic force. Cache refresh is not an
additional position correction or a claim that the full KKT system is solved.

Separate paired EL, LE, and joint controls carry their solved lambda across
steps, exactly as in the preceding warm-start study. Their post-predictor
material correction is not symplectic. They are labeled `paired_*`; they
must not be conflated with the energy-only method.

## Numerical variants and scope

Sixteen methods: ordinary XPBD; symplectic predictor; energy-only raw,
half-strength, tenth-strength, cap 0.1, half plus cap 0.1, half plus cap 0.01,
reject oversized (0.1), objective-direction rejection plus cap; and paired
EL/LE/joint with raw or fractional-repair-plus-cap acceptance.

The cap bounds the accepted extra displacement by

    cap * (h |v_old| + h² |a_net(x_old)|).

The energy-only half-plus-cap method scales the linear multiplier increment
by at most 1/2 and scales it further if necessary to satisfy the cap. Paired
repair+cap uses the previous whole-proposal filter relative to material-only
XPBD. Their common numeric cap does not make those different filters identical.

There are eight cases: released spring, vertical gravity spring, small swing,
large swing from rest, swing with tangential velocity, coupled radial and
angular motion, eccentric orbit, and a true circular orbit. Each has one moving
bob. No chain, cloth, grid, or hard rod was added: the symplectic elastic-force
proof does not cover an arbitrary stored hard-constraint reaction.

- 384 runs at 10 seconds: 16 methods x 8 cases x 240/480/960 Hz.
- 192 runs at 100 seconds: 8 selected baseline/bounded methods x the same cases
  and rates. Unprotected raw methods were not repeated for 100 seconds.
- One frozen correction pass per small substep, no convergence sweeps.
- RK4 references checked by subdivision doubling for every case and step size.
- Seven explicit stops among the ten-second runs; zero among the selected
  hundred-second runs. Several completed raw runs still have enormous spikes.

## Results at 480 Hz, ten seconds

The energy-only candidate below is fraction=0.5, cap=0.1. Energy error is the
maximum absolute error divided by initial energy above minimum potential.
Position error is RMS distance to the independently refined reference.

| Case | Predictor energy error | With energy correction | Predictor position RMS | With energy correction |
|---|---:|---:|---:|---:|
| Released spring | 1.0526% | 0.1450% | 0.002009 m | 0.004343 m |
| Vertical gravity spring | 1.0526% | 0.1450% | 0.001607 m | 0.003475 m |
| Small swing | 0.3121% | 0.0058% | 0.000365 m | 0.001799 m |
| Large swing from rest | 0.6741% | 0.0291% | 0.002075 m | 0.065508 m |
| Tangential swing | 0.5605% | 0.0052% | 0.002950 m | 0.005187 m |
| Coupled motion | 1.0481% | 0.0323% | 0.004965 m | 0.235189 m |
| Eccentric orbit | 0.6192% | 0.0245% | 0.000814 m | 0.253860 m |
| Circular orbit | 0.00362% | about 0.00004% | 0.000951 m | 0.000854 m |

The eccentric orbit's angular momentum begins at 0.88. The baseline ends at
0.88000000000035; the energy-only correction ends at 0.990261. Paired LE with
repair+cap ends at 1.233131 and has 0.923417 m RMS error. The energy-only update
therefore improves over paired warm LE, but is much worse than its own
symplectic baseline.

## One hundred seconds

Clipping preserves boundedness in these tests; phase and angular errors still
accumulate. At 480 Hz, position RMS changes as follows:

- Released spring: baseline 0.00356 m; energy-only half+cap 0.07802 m.
- Large spring swing: baseline 0.00306 m; corrected 0.72917 m.
- Eccentric orbit: baseline 0.00139 m; corrected 1.37759 m.
- Circular orbit: baseline 0.00175 m; corrected 0.000705 m.

Eccentric-orbit angular momentum ends at 0.88000000000188 for the baseline and
1.210123 for the corrected trajectory. The small maximum energy error persists
despite this drift. A smaller cap or rejecting large corrections does not
resolve the orbital failure: those filters are mostly inactive on its ordinary
steps. Objective-direction rejection often preserves more baseline-like
behavior by rejecting most updates, but does not give a uniformly useful method.

Refining h helps the ten-second eccentric orbit: corrected RMS at 240/480/960 Hz
is 0.43411/0.25386/0.13838 m. Baseline is 0.001749/0.000814/0.000392 m. Thus the
data support a large error coefficient, not a claim of proven nonconvergence.

## Why raw correction is dangerous

For a released scalar spring with extension q_old and a=h²k/m, the symplectic
predictor is q_star=(1-a)q_old. Its endpoint-energy gradient is -k a q_old.
One linear energy correction asks for

    delta_q=-(1-a)q_old/2.

This stays order q_old as h shrinks. With q_old=0.25, k=100, m=1 and h=1/480,
the first raw endpoint has roughly 576.75 times the initial energy. Such a run
can remain finite afterward, so counting only NaNs or safety stops is inadequate.
The independent first-step oracle confirms this in the new implementation.

## Verification and files

The tests verify exact symplectic-predictor parity, its central-force angular
invariant over 4,800 steps, the analytical raw-energy first step, the correction
cap, constitutive force-cache refresh, and equilibrium singular-gradient handling.
All five standalone CTests pass. Plot auditing independently reconstructs the
saved mechanical energies and position RMS errors for 48 full trajectory groups,
checks baseline completion, and checks all recorded energy-only correction caps.

Source: `symplectic_energy.hpp`, `symplectic_energy_study.cpp`.
CSV: `symplectic_energy_results/` and `symplectic_energy_long/`.
Plots and audit: `symplectic_energy_plots/`.
Reproducer: `plot_symplectic_energy.py` (NumPy and Matplotlib).

    cmake --build /tmp/xpbd-correction-build -j4
    /tmp/xpbd-correction-build/xpbd-symplectic-energy-study OUTPUT_DIR 10
    /tmp/xpbd-correction-build/xpbd-symplectic-energy-study LONG_OUTPUT_DIR 100

The next issue is the direction and branch chosen by the energy correction,
particularly at turning points and when radial and angular modes coexist.
These results do not justify moving this correction into a general solver or
extending it to a chain yet.
