#!/usr/bin/env python3
"""Line-by-line profile of the reference implementation's hot functions on one benchmark case.

The driver source/epi_int_lite.py runs unchanged in a scratch folder; helper_fns' functions are
wrapped by line_profiler before the driver's `from helper_fns import *` picks them up.

Run:  uv run python benchmark/line_profile.py [Nr] [Nt] [steps]     #default 1000 20 20
"""
import json
import os
import runpy
import shutil
import sys
from pathlib import Path

from line_profiler import LineProfiler

BENCH = Path(__file__).resolve().parent
REPO = BENCH.parent
sys.path.insert(0, str(REPO/'source'))
sys.path.insert(0, str(BENCH))
import helper_fns
from benchmark import params, case_name

FUNCTIONS = ['velocity_kick', 'accelerations', 'ring_gravity', 'interpolate_fn',
             'lagrange_poly_fit', 'ring_shear_viscosity', 'ring_bulk_viscosity', 'get_lambda',
             'coords2elem', 'elem2coords', 'sort_particles']


def main():
    Nr, Nt, steps = (int(x) for x in (sys.argv[1:] + ['1000', '20', '20'][len(sys.argv) - 1:])[:3])
    folder = BENCH/'runs'/f'{case_name(Nr, Nt)}_lineprof'
    shutil.rmtree(folder, ignore_errors=True)
    folder.mkdir(parents=True)
    shutil.copy(BENCH/'inputs.py', folder)
    lp = LineProfiler()
    for name in FUNCTIONS:
        setattr(helper_fns, name, lp(getattr(helper_fns, name)))
    os.chdir(folder)
    sys.argv = ['epi_int_lite.py', '-m', json.dumps(params(Nr, Nt, steps))]
    with open('output.log', 'w') as log:
        saved, sys.stdout = sys.stdout, log
        try:
            runpy.run_path(str(REPO/'source'/'epi_int_lite.py'), run_name='__main__')
        finally:
            sys.stdout = saved
    lp.print_stats(output_unit=1e-3, stripzeros=True)


if __name__ == '__main__':
    main()
