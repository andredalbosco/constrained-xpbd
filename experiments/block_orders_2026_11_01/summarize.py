"""Build measured comparison tables without treating failed prefixes as runs."""
import csv
from collections import Counter
from pathlib import Path
import numpy as np
from solver import scenes
from reference import reference_trajectory

ROOT=Path(__file__).resolve().parent


def number(value):
    return f'{float(value):.5g}'


def main():
    result=ROOT/'results'
    rows=list(csv.DictReader((result/'summary.csv').open()))
    def pick(case,order,hz=60,sweeps=4,ratio=100.):
        return next(r for r in rows if r['case']==case and r['order']==order
                    and int(r['hz'])==hz and int(r['sweeps'])==sweeps
                    and float(r['epsilon_ratio'])==ratio)
    text=['# Measured block-order results','',
          'These are refreshed nonlinear block-Gauss-Seidel runs of the agreed '
          'position-only objective. The material block solves position/material; '
          'the energy block solves position/relaxation/energy multiplier together. '
          'Both retain position stationarity and immediately update their primal variables.','',
          'The campaign contains 468 five-second configurations: six scenes, 60/120 Hz, '
          '1/4/16 sweeps, and four energy-penalty ratios (1, 100, 10000, 100000000). '
          'Each EL/LE configuration has an approximate joint-system control; ordinary '
          'XPBD is an additional control. JOINT uses the same Gauss-Newton/omitted-energy-curvature '
          'system; it is not the original published CPD solver.','',
          '## Completion and physical quality','',
          '| Method | Finished five seconds | Guard failures | Runs with max normalized energy error > 10 |',
          '|---|---:|---:|---:|']
    for order in ('EL','LE','JOINT','XPBD'):
        subset=[r for r in rows if r['order']==order]
        complete=sum(int(r['completed']) for r in subset)
        large=sum(float(r['max_rel_energy'])>10 for r in subset)
        text.append(f'| {order} | {complete}/{len(subset)} | {len(subset)-complete} | {large} |')
    text+=['',
           'Finished means the requested time budget was reached within the explicit guards. '
           'It does not imply plausible motion or convergence. Of the 463 completed configurations, '
           "297 contain at least one step above that method's residual tolerance. EL/LE/JOINT "
           'test all four stationary equations. XPBD convergence uses its material/position '
           'equations; energy remains a diagnostic for that control. These criteria establish '
           'first-order stationarity, not a certified constrained minimum.','',
           'Maxima and final energy include an offending rejected candidate. RMS and path length '
           'use accepted states only; a failed-prefix RMS is not a full-run accuracy result. '
           'The normalized energy error is absolute error from run-initial mechanical energy '
           'divided by initial energy above the declared potential minimum.','',
           '## Representative matched settings','',
           'Five seconds, 60 Hz, four sweeps; epsilon is 100 times initial energy above '
           'the potential minimum, held fixed for the run. All entries below finish.','',
           '| Scene | Order | Max energy error (%) | Position RMS (m) | Converged steps (%) |',
           '|---|---|---:|---:|---:|']
    for scene in scenes():
        for order in ('EL','LE'):
            r=pick(scene.name,order)
            text.append(f"| {scene.name} | {order} | {100*float(r['max_rel_energy']):.4g} | "
                        f"{number(r['rms_position'])} | {100*float(r['converged_fraction']):.4g} |")
    text+=['',
           'At four sweeps LE ends with the energy correction, while EL ends with the material '
           'correction. The resting-pendulum length errors illustrate the tradeoff: '
           f"EL {number(pick('pendulum_release','EL')['max_geometry_error'])} m versus "
           f"LE {number(pick('pendulum_release','LE')['max_geometry_error'])} m. "
           'The true nonlinear rows can reopen after the other block.','',
           '## Small energy error can hide incorrect motion','',
           'For the unforced orbit at 60 Hz, 16 sweeps and penalty ratio 100000000, '
           'EL, LE and JOINT all report convergence at every step and agree to numerical precision. '
           'Their measured quantities are:','',
           '| Order | Max normalized energy error | Position RMS (m) | Max angular-momentum error |',
           '|---|---:|---:|---:|']
    for order in ('EL','LE','JOINT'):
        r=pick('spring_orbit',order,sweeps=16,ratio=1e8)
        text.append(f"| {order} | {number(r['max_rel_energy'])} | {number(r['rms_position'])} | {number(r['max_angular_error'])} |")
    release=pick('spring_release','JOINT',sweeps=16,ratio=1e8)
    scene=next(s for s in scenes() if s.name=='spring_release')
    reference,_=reference_trajectory(scene,np.arange(301)/60.)
    reference_path=float(np.sum(np.linalg.norm(np.diff(reference,axis=0),axis=1)))
    text+=['',
           'Angular momentum is measured but is not constrained in this energy-only adaptation. '
           'Reaching the stationary equations does not establish correct motion. The same observed '
           'root is reached by both orders and the joint control; changing EL to LE cannot resolve '
           'this case. The separate curvature audit examines whether selected stationary roots '
           'also satisfy the second-order condition.','',
           f"The high-penalty joint spring release travels only {number(release['path_length'])} m "
           f"over five seconds, versus {number(reference_path)} m in the continuous reference. "
           f"Its max normalized energy error is {number(release['max_rel_energy'])}, and every step "
           'meets the residual tolerance. This is the nearly frozen energy-level branch.','',
           '## A concrete EL failure mechanism','',
           'The second step of the 60 Hz, one-sweep spring-release run with penalty ratio 10000 '
           'isolates the order effect from different trajectory histories. Starting both orders '
           'from the same incoming state, the E block gives an energy-multiplier increment near '
           '-383.46 and a position correction of -0.07197 m. Its frozen position row closes, '
           'but the refreshed position residual is about 102112.92. The following L block moves '
           '27.60 m, producing physical energy about 1.41e6. From exactly that incoming state, '
           'LE with one sweep has max block movement 0.00639 m and endpoint energy about 3.06143.','',
           'The diagnostic verifies the exact nonlinear remainder in both position stationarity '
           'and energy. With 64 sweeps EL, LE and JOINT reach the same root for this one step. '
           'JOINT with one iteration also overshoots. Binding variables does not eliminate '
           'the omitted-curvature/nonlinear-step problem. See [DIAGNOSTIC.md](DIAGNOSTIC.md).','',
           '## Stationarity versus a minimum','',
           'A selected-root audit evaluates the full primal Lagrangian Hessian in position and '
           'rho, restricted to the physical-energy constraint tangent. It uses the analytical '
           'spring Hessian rather than the approximation used by the solver.','',
           '| Selected roots | Count | Tangent-space curvature |',
           '|---|---:|---|',
           '| Orbit, EL/LE/JOINT, 60 Hz, 16 sweeps, ratio 1e8 | 900 | Positive |',
           '| Released spring, JOINT, same settings | 300 | Near-zero / inconclusive |',
           '| Released spring, LE, first step, 60 Hz, one sweep, ratio 1e4 | 1 | Negative definite |','',
           'The inaccurate orbit therefore survives this local second-order check; global '
           'optimality and trajectory accuracy are separate questions. The first LE release step '
           'is especially instructive: eta is about -270.27 and the tangent Hessian eigenvalues '
           'are about -996300 and -974641, despite a scaled stationarity residual of 3.27e-12. '
           'It reaches a stationary maximum, rather than a minimum. Full residual convergence '
           'alone is insufficient for accepting an optimization solution. See '
           '`curvature_audit.py` and `diagnostics/curvature_audit.json`.','',
           '## Scope and interpretation','',
           '- This tests the two requested block orders with the current positive material compliance '
           'and hard-rod extension, not an independently optimized velocity.',
           '- EL and LE have matched numbers of block visits. A JOINT iteration does different work '
           'from a scalar-block visit; iteration counts are not equal compute budgets.',
           '- No line search, caps, hidden substeps, force clipping, energy rescale, or fallback integrator '
           'is used in the primary campaign.',
           '- Energy targets use incoming physical energy and are fixed within each step. Targets can '
           'therefore inherit prior finite-budget energy errors. Fixed run-initial targets are not tested here.',
           '- epsilon is fixed per run; the induced energy compliance changes as incoming energy changes.',
           '- Endpoint velocity is reconstructed. No endpoint-tangency constraint is added to the hard rod; '
           'the radial-velocity defect is recorded.',
           '- There are no networks or multi-pendulums in this focused comparison. The global-energy '
           'support/local-material support problem remains to be tested there.','',
           'The useful finding is twofold: block order materially affects truncated nonlinear solves, '
           'and strong energy enforcement still admits nearly frozen or inaccurate converged motion. '
           'These trials do not establish a generally stable replacement for EL.','',
           'Validation: 11,600 equation/update checks, an independent audit of all 468 runs and '
           '208,735 trajectory rows, and 5,292 replayed serialized visits. Continuous references '
           'pass tighter-tolerance refinement and comparison with the previously verified reference implementation.','']
    (ROOT/'RESULTS.md').write_text('\n'.join(text))
    print('Wrote RESULTS.md')


if __name__=='__main__':
    main()
