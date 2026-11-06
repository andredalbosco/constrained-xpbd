# Validation from the saved repository package

On 2026-10-31, the saved source at `/home/blacknoir/Projects/constrained-xpbd/experiments/round2_2026_10_31` configured and built successfully into the fresh directory `/tmp/xpbd-round2-canonical-build` with Clang 20.1.8 and the Release configuration. All 10 executable targets linked. The build used the package's default dependency paths in the repository, without Global source-path overrides.

All 4 CTest checks passed. Their output is retained in [canonical_ctest.log](canonical_ctest.log). Check output and generated validation CSVs went into the build directory, leaving saved research data intact.

Reproduction:

```bash
cmake -S /home/blacknoir/Projects/constrained-xpbd/experiments/round2_2026_10_31 -B /tmp/xpbd-round2-canonical-build -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/xpbd-round2-canonical-build -j 4
ctest --test-dir /tmp/xpbd-round2-canonical-build --output-on-failure
```

Archive SHA-256 verification and local report-link checks also passed. The original main solver and scene files retained their pre-archive hashes. This validates saved-source portability on this workstation; it does not establish portability to other machines or compilers.
