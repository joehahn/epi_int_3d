#!/usr/bin/env python3
"""Run the test simulations that the notebooks analyze, in parallel.

Each folder under the given paths that holds an inputs.py is a simulation; it is run from that
folder, writing output/ (and output.log). Simulations that already have an output/ folder are
skipped unless --force is given.
- The tests ported from epi_int_lite (unperturbed, spreading, libration) run with
  source/epi_int_lite.py, or with --fast with fast/epi_int_3d.py and epi_int_lite_gravity = True,
  which gives the same results bit for bit.
- epi_int_3d's own tests (EPI_INT_3D_TESTS: tests/nominal, tests/inclined, tests/self_confining (epi_int_lite's
  nominal ringlet, compared with that run's figures), tests/vertical_resonance (test particles at a
  satellite's vertical resonance), tests/bending_waves (a self-gravitating ring's bending wave), tests/lindblad_resonance
  and tests/density_waves (their in-plane counterparts at a Lindblad resonance), and any test with an
  `inclination` in inputs.py) always run with fast/epi_int_3d.py and its default ring gravity, as
  described in the paper.
- VARIANTS: extra runs of a test with some inputs changed, each in its own output_<name>/ (and output_<name>.log):
  tests/bending_waves at 400 and 1600 streamlines, beside its 800, for the paper's wavenumber convergence (Fig. 4);
  tests/density_waves at 400 and 1600 streamlines and with viscosity (runs A, C, D), for the viscous damping (Fig. 9).

Run:  uv run python tools/run_simulations.py                     #all 29 test simulations
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


#extra runs: test folder -> {output folder: modified inputs}
sys.path.insert(0, str(REPO/'tests'))
from vertical_analysis import RESOLUTION_RUNS
from density_analysis import VISCOUS_RUNS, RESOLUTION_VARIANTS
VARIANTS = {'bending_waves': {k: v for k, v in RESOLUTION_RUNS.items() if k != 'output'},
            'density_waves': dict(RESOLUTION_VARIANTS, **VISCOUS_RUNS)}


EPI_INT_3D_TESTS = ('nominal', 'inclined', 'self_confining', 'vertical_resonance', 'bending_waves', 'lindblad_resonance',
                    'density_waves')    #epi_int_3d's own tests, under tests/


def own_test(folder):
    """Is it one of epi_int_3d's own tests, rather than one ported from epi_int_lite? Those run with
    epi_int_3d's own ring gravity; so does any test that needs vertical motion."""
    rel = folder.relative_to(REPO).parts
    if rel[0] == 'tests' and len(rel) > 1 and rel[1] in EPI_INT_3D_TESTS:
        return True
    return re.search(r'^inclination\s*=', (folder/'inputs.py').read_text(), re.M) is not None


def run(folder, driver=DRIVER, output='output', modified=None):
    t0 = time.time()
    cmd = [sys.executable, str(driver)]
    modified = dict(modified or {})
    if driver == FAST_DRIVER and not own_test(folder):
        modified['epi_int_lite_gravity'] = True                        #so --fast gives the reference's results
    if output != 'output':
        modified['output_folder'] = output
    if modified:
        cmd += ['-m', json.dumps(modified)]
    with open(folder/f'{output}.log', 'w') as log:
        rc = subprocess.run(cmd, cwd=folder, stdout=log,
                            stderr=subprocess.STDOUT).returncode
    return folder/output if output != 'output' else folder, rc, time.time() - t0


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('paths', nargs='*', default=['tests'])
    parser.add_argument('-j', '--jobs', type=int, default=os.cpu_count())
    parser.add_argument('--force', action='store_true', help='rerun simulations that have output/')
    parser.add_argument('--fast', action='store_true', help='use the faster driver, fast/epi_int_3d.py')
    args = parser.parse_args()
    folders = sorted({f.parent for p in args.paths for f in (REPO/p).resolve().rglob('inputs.py')})
    runs = [(f, 'output', None) for f in folders]
    runs += [(f, out, mod) for f in folders for out, mod in VARIANTS.get(f.name, {}).items()]
    todo = [r for r in runs if args.force or not (r[0]/r[1]).is_dir()]
    for f, out, mod in runs:
        if (f, out, mod) not in todo:
            print(f'skip (has {out}/): {f.relative_to(REPO)}')
    driver = lambda f: FAST_DRIVER if args.fast or own_test(f) else DRIVER
    with ThreadPoolExecutor(args.jobs) as pool:
        for folder, rc, dt in pool.map(lambda r: run(r[0], driver(r[0]), r[1], r[2]), todo):
            status = 'done' if rc == 0 else f'FAILED (exit {rc}, see output.log)'
            print(f'{status}: {folder.relative_to(REPO)}  ({dt/60:.1f} min)', flush=True)


if __name__ == '__main__':
    main()
