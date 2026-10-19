#!/usr/bin/env python3
"""epi_int_3d's faster driver (Phase 0b): the reference driver, source/epi_int_lite.py, run unchanged
with ring self-gravity computed by the Numba kernel in fast/kernels.py. Its results are bit-for-bit
the reference implementation's.

Run it like the reference driver, from a folder holding inputs.py:
    cd tests/nominal && uv run python ../../fast/epi_int_fast.py [-m '{"total_number_of_outputs": 20}']
The kernel uses all cores; set NUMBA_NUM_THREADS to use fewer.
"""
import runpy
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO/'source'))
sys.path.insert(0, str(REPO/'fast'))
import helper_fns
import kernels


def ring_gravity(lambda_, G_ring, r, t, vr, vt, fast_gravity):
    """Drop-in for helper_fns.ring_gravity; epi_int_3d always has fast_gravity = False."""
    if fast_gravity:
        raise ValueError('fast_gravity is not supported: set fast_gravity = False')
    return kernels.ring_gravity_parallel(lambda_, G_ring, r, t, vr, vt)


helper_fns.ring_gravity = ring_gravity
sys.argv[0] = str(REPO/'source'/'epi_int_lite.py')
runpy.run_path(sys.argv[0], run_name='__main__')
