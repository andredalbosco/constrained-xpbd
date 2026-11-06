#!/usr/bin/env python3
"""Run the direct-energy coupling studies, retaining recorded failure cases."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess


def main():
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=repo / "build")
    parser.add_argument("--output", type=Path, default=repo / "results/energy_coupling")
    args = parser.parse_args()
    build, output = args.build_dir.resolve(), args.output.resolve()
    studies = ("block", "spring", "pendulum")
    for study in studies:
        binary = build / f"xpbd-energy-{study}-study"
        if not binary.is_file():
            parser.error(f"missing {binary}; build with XPBD_BUILD_ENERGY_STUDIES=ON")

    def run(study):
        directory = output / study
        directory.mkdir(parents=True, exist_ok=True)
        filename = "summary.csv" if study == "block" else "run.log"
        with (directory / filename).open("w") as stdout, (directory / "stderr.log").open("w") as stderr:
            process = subprocess.run([str(build / f"xpbd-energy-{study}-study")], cwd=directory,
                                     stdout=stdout, stderr=stderr, check=False)
        return study, process.returncode

    with ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(run, studies))
    for study, code in results:
        print(f"{'PASS' if code == 0 else 'FAIL'} study checks: {study}: {output / study}")
    print("PASS means the study's validation checks passed; unstable variants remain recorded in the CSVs.")
    return int(any(code != 0 for _, code in results))


if __name__ == "__main__":
    raise SystemExit(main())
