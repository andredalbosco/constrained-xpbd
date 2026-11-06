# Validation from the saved repository package

On 2026-10-31, the saved source at `/home/blacknoir/Projects/constrained-xpbd/experiments/hour_2026_10_31` configured and built successfully into the fresh directory `/tmp/xpbd-hour-canonical-build` with Clang 20.1.8 and the Release configuration. All 28 executable targets linked. The build used the package's default dependency paths in the repository, without Global source-path overrides.

All 11 CTest checks passed. Their output is retained in [canonical_ctest.log](canonical_ctest.log). Check output and generated validation CSVs went into the build directory, leaving saved research data intact.

Reproduction:

```bash
cmake -S /home/blacknoir/Projects/constrained-xpbd/experiments/hour_2026_10_31 -B /tmp/xpbd-hour-canonical-build -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/xpbd-hour-canonical-build -j 4
ctest --test-dir /tmp/xpbd-hour-canonical-build --output-on-failure
```

Archive SHA-256 verification and local report-link checks also passed. The original main solver and scene files retained their pre-archive hashes. This validates saved-source portability on this workstation; it does not establish portability to other machines or compilers.

A final read-only browser connection retry still returned `No browser is available`, and supported discovery returned an empty list. Interactive HTML control verification therefore remains unavailable; the native MP4 checks are independent.
