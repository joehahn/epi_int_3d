#!/usr/bin/env python3
"""Run the test simulations that the notebooks analyze, in parallel.

Each folder under the given paths that holds an inputs.py is a simulation; it is run from that
folder, writing output/ (and output.log). Simulations that already have an output/ folder are
skipped unless --force is given.
- The tests ported from epi_int_lite (unperturbed, spreading, libration) run with
  source/epi_int_lite.py, or with --fast with fast/epi_int_3d.py and epi_int_lite_gravity = True,
  which gives the same results bit for bit.
- epi_int_3d's own tests (EPI_INT_3D_TESTS: tests/nominal, tests/inclined, and any test with an
  `inclination` in inputs.py) always run with fast/epi_int_3d.py and its default ring gravity, as
  described in the paper.

Run:  uv run python tools/run_simulations.py                     #all 24 test simulations
      uv run python tools/run_simulations.py tests/spreading       #one test
      uv run python tools/run_simulations.py tests/libration -j 8  #the 20 libration runs
      uv run python tools/run_simulations.py tests/spreading --fast  #the faster driver, same results
"""
import argparse
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DRIVER = REPO/'source'/'epi_int_lite.py'
FAST_DRIVER = REPO/'fast'/'epi_int_3d.py'


EPI_INT_3D_TESTS = ('nominal', 'inclined')    #epi_int_3d's own tests, under tests/


def own_test(folder):
    """Is it one of epi_int_3d's own tests, rather than one ported from epi_int_lite? Those run with
    epi_int_3d's own ring gravity; so does any test that needs vertical motion."""
    rel = folder.relative_to(REPO).parts
    if rel[0] == 'tests' and len(rel) > 1 and rel[1] in EPI_INT_3D_TESTS:
        return True
    return re.search(r'^inclination\s*=', (folder/'inputs.py').read_text(), re.M) is not None


def run(folder, driver=DRIVER):
    t0 = time.time()
    cmd = [sys.executable, str(driver)]
    if driver == FAST_DRIVER and not own_test(folder):
        cmd += ['-m', json.dumps({'epi_int_lite_gravity': True})]     #so --fast gives the reference's results
    with open(folder/'output.log', 'w') as log:
        rc = subprocess.run(cmd, cwd=folder, stdout=log,
                            stderr=subprocess.STDOUT).returncode
    return folder, rc, time.time() - t0


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('paths', nargs='*', default=['tests'])
    parser.add_argument('-j', '--jobs', type=int, default=os.cpu_count())
    parser.add_argument('--force', action='store_true', help='rerun simulations that have output/')
    parser.add_argument('--fast', action='store_true', help='use the faster driver, fast/epi_int_3d.py')
    args = parser.parse_args()
    folders = sorted({f.parent for p in args.paths for f in (REPO/p).resolve().rglob('inputs.py')})
    todo = [f for f in folders if args.force or not (f/'output').is_dir()]
    for f in sorted(set(folders) - set(todo)):
        print(f'skip (has output/): {f.relative_to(REPO)}')
    with ThreadPoolExecutor(args.jobs) as pool:
        for folder, rc, dt in pool.map(lambda f: run(f, FAST_DRIVER if args.fast or own_test(f) else DRIVER), todo):
            status = 'done' if rc == 0 else f'FAILED (exit {rc}, see output.log)'
            print(f'{status}: {folder.relative_to(REPO)}  ({dt/60:.1f} min)', flush=True)


if __name__ == '__main__':
    main()
