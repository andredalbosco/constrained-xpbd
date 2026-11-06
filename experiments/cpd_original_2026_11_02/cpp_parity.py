"""Cross-check the new planar solver against unchanged local 3D C++ CPD.

The frozen reference CSV is shipped. --rebuild-reference recompiles the old local
baseline with the shipped harness; old sources are read but never modified.
This is not comparison with the paper authors' unavailable implementation.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np

import cpd_solver

ROOT = Path(__file__).resolve().parent
OLD = Path('/home/blacknoir/Projects/dynamics-research')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rebuild-reference', action='store_true')
    args = parser.parse_args()
    commands = []
    if args.rebuild_reference:
        eigen = (Path('/usr/include/eigen3') if Path('/usr/include/eigen3/Eigen/Core').is_file()
                 else OLD/'external/eigen')
        with tempfile.TemporaryDirectory(prefix='xpbd-cpd-cpp-parity-') as directory:
            executable = Path(directory)/'parity'
            command = ['g++', '-std=c++17', '-O2', '-I'+str(eigen), '-I'+str(OLD/'src/core'),
                       '-I'+str(OLD/'experiments/work_balance'), str(ROOT/'cpp_parity.cpp'),
                       str(OLD/'experiments/work_balance/baselines.cpp'), '-o', str(executable)]
            subprocess.run(command, check=True)
            commands.append(command)
            reference = subprocess.run([str(executable)], check=True, capture_output=True, text=True)
            (ROOT/'CPP_REFERENCE.csv').write_text(reference.stdout)
    scenes = {scene.name: scene for scene in cpd_solver.scenes()}
    with (ROOT/'CPP_REFERENCE.csv').open() as source:
        reference_rows = list(csv.DictReader(source))
    rows = []
    for reference in reference_rows:
        scene = scenes[reference['case']]
        cap, epsilon = int(reference['cap']), float(reference['epsilon'])
        actual = cpd_solver.step(scene, scene.x0, scene.v0, 1/60, epsilon,
                                 policy='fixed', iterations=cap)
        x = np.array([[float(reference['x0']), float(reference['y0'])],
                      [float(reference['x1']), float(reference['y1'])]])
        v = np.array([[float(reference['vx0']), float(reference['vy0'])],
                      [float(reference['vx1']), float(reference['vy1'])]])
        row = dict(case=scene.name, cap=cap, epsilon=epsilon,
                   old_iterations=int(reference['iterations']), new_iterations=actual.iterations,
                   new_accepted=actual.accepted,
                   position_error=float(np.max(abs(actual.x-x))),
                   velocity_error=float(np.max(abs(actual.v-v))),
                   rho_error=abs(actual.rho-float(reference['rho'])))
        assert actual.accepted and row['old_iterations'] == actual.iterations == cap
        assert max(row['position_error'], row['velocity_error'], row['rho_error']) < 1e-10
        rows.append(row)
    with (ROOT/'CPP_PARITY.csv').open('w', newline='') as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    inputs = [OLD/'experiments/work_balance/baselines.cpp',
              OLD/'experiments/work_balance/baselines.hpp', OLD/'src/core/system.hpp',
              ROOT/'cpp_parity.cpp', ROOT/'CPP_REFERENCE.csv']
    receipt = dict(status='passed', cases=len(rows),
                   max_position_error=max(r['position_error'] for r in rows),
                   max_velocity_error=max(r['velocity_error'] for r in rows),
                   max_rho_error=max(r['rho_error'] for r in rows),
                   comparison='unchanged earlier local baseline, not authors source code',
                   source_sha256={str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},
                   rebuild_commands=commands)
    (ROOT/'CPP_PARITY.json').write_text(json.dumps(receipt, indent=2, sort_keys=True))
    (ROOT/'CPP_PARITY.md').write_text('''# Earlier full-CPD C++ source cross-check

The unchanged local three-dimensional CPD baseline is compared with the new planar
implementation for free two-mass springs. Each mass is two, COM is zero, and
relative initial conditions match the anchored unit-mass examples. The cases are
release, orbit, and stiff spring at h=1/60, epsilon=1e4/1e8, and caps 1/4/16.
Old tolerance 1e-30 makes iteration counts agree with each fixed cap in all cases.

The old 3D implementation uses its disclosed rank-revealing COD solve for the
dependent angular row. The new planar implementation omits the unused dimensions
and uses a direct full-row solve. Agreement here therefore also checks that the
planar reduction preserves the relevant full free-body equations.

This comparison establishes parity with the earlier local implementation for
these 18 one-step cases. It does not establish parity with authors' source code
or general equivalence of their different rank-loss policies.

```bash
OPENBLAS_NUM_THREADS=1 python3 cpp_parity.py
# Recompile unchanged old sources when they are available:
OPENBLAS_NUM_THREADS=1 python3 cpp_parity.py --rebuild-reference
```

```json
'''+json.dumps(receipt, indent=2, sort_keys=True)+'\n```\n')
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
