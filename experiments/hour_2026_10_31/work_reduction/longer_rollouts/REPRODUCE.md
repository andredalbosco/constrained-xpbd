# Five-second follow-up

Archive location: `hour_2026_10_31/work_reduction/longer_rollouts`. This folder depends on the frozen parent hour package and its first/second-round headers. It does not alter them. Sources were copied from the canonical hour package, with changes retained in `source_changes.diff`.

The numerical experiment changes requested duration from one to five seconds. It preserves the same old states, mass ratio10, N8/N32, stiffness factors1/10, h=1/30 and1/60, tolerances, 128-evaluation total budget, and stop-at-first-failure rule. Methods run separately so GN-polish has priority. Added failure diagnostics record the inherited solver reason and rejected residuals; they do not change steps. The inherited `kkt_line_search_stalled` label can also occur when the evaluation budget ends inside line search. `failed_evaluations=128` identifies an exhausted total budget; it is not proof that extra budget would converge. The copied numerical solver header is unchanged. GN-polish permits64 GN evaluations then exact Newton with the remaining budget. The physical reference used by trajectory errors remains RK4 with256 substeps, checked against128. A separate512/1024 pair probes the two uncertified hard references without changing that baseline.

Exact compilation command used (replace `SOURCE` with `rollouts`, `reference_audit`, or `reference_fine`):

```bash
clang++ -std=c++17 -O2 -DNDEBUG \
 -I/home/blacknoir/Projects/constrained-xpbd/external/eigen \
 -I/home/blacknoir/Projects/constrained-xpbd/experiments/composition_2026_10_31 \
 -I/home/blacknoir/Projects/constrained-xpbd/experiments/round2_2026_10_31/local_work \
 -I/home/blacknoir/Projects/constrained-xpbd/experiments/hour_2026_10_31/work_reduction \
 -I/home/blacknoir/Projects/constrained-xpbd/experiments/hour_2026_10_31/local_solver \
 SOURCE.cpp -o SOURCE
./rollouts gn_polish reduced_GN_polish > gn_polish.log
./rollouts full full > full.log
./rollouts exact_parallel reduced_exact > exact_parallel.log
./reference_audit references > reference.log
./reference_fine references_fine > reference_fine.log
python3 analyze.py /home/blacknoir/Projects/constrained-xpbd/experiments/hour_2026_10_31/work_reduction/rollout_results
```

`finalize.py` records the observed wall-clock interruption and appends its status report after analysis; it is specific to this bounded run.

The rollout compiler needed only the `work_reduction` include directory, and the reference compilers needed only `local_solver`; listing both is equivalent. CMake provides portable path overrides. After archiving, its default `HOUR_ROOT=../..` resolves to the parent hour package. During Global staging pass `-DHOUR_ROOT=/home/blacknoir/Projects/constrained-xpbd/experiments/hour_2026_10_31`.

Recorded times overlap other work and are not performance evidence. A trace ends at its last converged accepted step; summary failure residuals describe the rejected next step. A zero-step failure has no physical RMS, not zero error. CSV physical RMS on uncertified cases must not be presented as validated accuracy. Prefix checks compare each accepted first-second prefix to the original saved experiment; no original experiment is overwritten or rerun.
