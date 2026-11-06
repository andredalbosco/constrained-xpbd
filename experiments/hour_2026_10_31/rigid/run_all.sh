#!/usr/bin/env bash
set -euo pipefail
study_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$study_dir"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release "$@"
cmake --build build -j2
ctest --test-dir build --output-on-failure
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-6}"
for pair in \
  short:results order:order_results accuracy:accuracy_results \
  stress:stress_results long:long_results momentum:momentum_results \
  metric:metric_results objectives:objectives_results \
  objectives_order:objectives_order_results reaction_order:reaction_order_results \
  equilibria:equilibria_results equilibria_stable:equilibria_stable_results \
  momentum_stable:momentum_stable_results; do
  mode="${pair%%:*}"
  directory="${pair#*:}"
  ./build/rigid_study "$mode" "$directory"
done
./build/rigid_probes
./build/rigid_modes
./build/rigid_modal_transfer
./build/rigid_controller_probe
"${XPBD_PYTHON:-python3}" build_viewer.py
MPLCONFIGDIR="${TMPDIR:-/tmp}/xpbd-rigid-plots" "${XPBD_PYTHON:-python3}" analyze.py
