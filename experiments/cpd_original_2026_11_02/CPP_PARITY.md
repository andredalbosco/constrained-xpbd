# Earlier full-CPD C++ source cross-check

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
{
  "cases": 18,
  "comparison": "unchanged earlier local baseline, not authors source code",
  "max_position_error": 4.336808689942018e-18,
  "max_rho_error": 1.6805133673525319e-18,
  "max_velocity_error": 2.7755575615628914e-16,
  "rebuild_commands": [],
  "source_sha256": {
    "/home/blacknoir/Global/constrained-xpbd-cpd-original-2026-11-02/CPP_REFERENCE.csv": "50fe527a77bd9f04bf0b22be3c48bc1e629376555ce728f5c662a308ea2c5fcc",
    "/home/blacknoir/Global/constrained-xpbd-cpd-original-2026-11-02/cpp_parity.cpp": "6d0407790841054121f8f005e74a14703c5c2b2c40aa1528bf9c274341f5cdc6",
    "/home/blacknoir/Projects/dynamics-research/experiments/work_balance/baselines.cpp": "14df774bc0e687bc761eb52715f29054b793338aaf1deb8a9b9f15a536a7fa03",
    "/home/blacknoir/Projects/dynamics-research/experiments/work_balance/baselines.hpp": "cfea006162e84227497b74994907d295c51406a910cb65736be41084fa5e8ab9",
    "/home/blacknoir/Projects/dynamics-research/src/core/system.hpp": "49cae1558316474439051e446d5745d835b222d1490c2910482533c64988337a"
  },
  "status": "passed"
}
```
