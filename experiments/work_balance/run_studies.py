#!/usr/bin/env python3
"""Run the independent C++ work studies and collect their CSVs and logs."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess


def main():
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=repo / "build")
    parser.add_argument("--output", type=Path, default=repo / "results/work_balance")
    args = parser.parse_args()
    build, output = args.build_dir.resolve(), args.output.resolve()
    studies = ("kinetic-work", "midpoint-gs", "spring", "pendulum", "coupled-momentum")
    for study in studies:
        binary = build / f"xpbd-{study}-study"
        if not binary.is_file():
            parser.error(f"missing {binary}; configure and build with XPBD_BUILD_WORK_STUDIES=ON")

    def run(study):
        directory = output / study
        directory.mkdir(parents=True, exist_ok=True)
        # These three programs write CSVs themselves; the other studies print CSV.
        filename = "run.log" if study in ("spring", "pendulum", "coupled-momentum") else "summary.csv"
        with (directory / filename).open("w") as stdout, (directory / "stderr.log").open("w") as stderr:
            process = subprocess.run([str(build / f"xpbd-{study}-study")], cwd=directory,
                                     stdout=stdout, stderr=stderr, check=False)
        return study, process.returncode

    with ThreadPoolExecutor(max_workers=len(studies)) as pool:
        results = list(pool.map(run, studies))
    for study, code in results:
        print(f"{'PASS' if code == 0 else 'FAIL'} {study}: {output / study}")
    return int(any(code != 0 for _, code in results))


if __name__ == "__main__":
    raise SystemExit(main())
