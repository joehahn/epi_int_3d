#!/usr/bin/env python3
"""Run the test simulations that the notebooks analyze, in parallel.

Each folder under the given paths that holds an inputs.py is a simulation; it is run from that
folder with source/epi_int_lite.py (or, with --fast, fast/epi_int_fast.py, which gives the same results
bit for bit), writing output/ (and output.log). Simulations that already have an output/ folder are
skipped unless --force is given.

Run:  uv run python tools/run_simulations.py                     #all 22 test simulations
      uv run python tools/run_simulations.py tests/spreading       #one test
      uv run python tools/run_simulations.py tests/libration -j 8  #the 20 libration runs
      uv run python tools/run_simulations.py tests/nominal --fast  #16 min instead of ~2 h
"""
import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DRIVER = REPO/'source'/'epi_int_lite.py'
FAST_DRIVER = REPO/'fast'/'epi_int_fast.py'


def run(folder, driver=DRIVER):
    t0 = time.time()
    with open(folder/'output.log', 'w') as log:
        rc = subprocess.run([sys.executable, str(driver)], cwd=folder, stdout=log,
                            stderr=subprocess.STDOUT).returncode
    return folder, rc, time.time() - t0


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('paths', nargs='*', default=['tests'])
    parser.add_argument('-j', '--jobs', type=int, default=os.cpu_count())
    parser.add_argument('--force', action='store_true', help='rerun simulations that have output/')
    parser.add_argument('--fast', action='store_true', help='use the faster driver, fast/epi_int_fast.py')
    args = parser.parse_args()
    folders = sorted({f.parent for p in args.paths for f in (REPO/p).resolve().rglob('inputs.py')})
    todo = [f for f in folders if args.force or not (f/'output').is_dir()]
    for f in sorted(set(folders) - set(todo)):
        print(f'skip (has output/): {f.relative_to(REPO)}')
    with ThreadPoolExecutor(args.jobs) as pool:
        for folder, rc, dt in pool.map(lambda f: run(f, FAST_DRIVER if args.fast else DRIVER), todo):
            status = 'done' if rc == 0 else f'FAILED (exit {rc}, see output.log)'
            print(f'{status}: {folder.relative_to(REPO)}  ({dt/60:.1f} min)', flush=True)


if __name__ == '__main__':
    main()
