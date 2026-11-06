#!/usr/bin/env python3
"""Rebuild and execute the recorded experiment design, preserving existing results by default."""
import argparse, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
TIERS=['core','parameters','solver','long','stress','traces','convergence','stress_reference','serial_cost']
def main():
    p=argparse.ArgumentParser();p.add_argument('--threads',type=int,default=8);p.add_argument('--tiers',nargs='+',choices=TIERS+['smoke'],default=TIERS);p.add_argument('--output',type=Path,default=ROOT/'rerun_results');p.add_argument('--xpbd-source',type=Path,default=Path('/home/blacknoir/Projects/constrained-xpbd'));a=p.parse_args()
    if a.output.exists() and any(a.output.iterdir()):p.error('output must be empty; choose a new --output directory')
    def run(cmd):subprocess.run([str(x) for x in cmd],cwd=ROOT,check=True)
    run(['cmake','-S',ROOT,'-B',ROOT/'build','-DCMAKE_BUILD_TYPE=Release',f'-DXPBD_SOURCE={a.xpbd_source}']);run(['cmake','--build',ROOT/'build','-j4']);run(['ctest','--test-dir',ROOT/'build','--output-on-failure'])
    a.output.mkdir(parents=True,exist_ok=True)
    for tier in a.tiers:
        run([sys.executable,'design.py','--tier',tier]);run([ROOT/'build/study','--manifest',ROOT/'manifests'/f'{tier}.csv','--out',a.output/tier,'--threads',1 if tier=='serial_cost' else a.threads])
    run([ROOT/'build/landscape',a.output/'landscape.csv'])
    print('Simulation outputs:',a.output)
if __name__=='__main__':main()
