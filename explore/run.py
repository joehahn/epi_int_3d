#!/usr/bin/env python3
"""Run the reference driver, source/epi_int_lite.py, unchanged, with the run-time patches of
tools/patches.py: always the longitude-wrapped streamline-crossing check (issue 1), and optionally

    --taper f        taper the ringlet's edges over a fraction f of the width at each edge
    --fix-edge-own   outer-edge flux = the outer streamline's own, at its own longitudes (issue 4)
    --fix-edge       outer-edge flux = hahn13's p(r - Delta)

Run from a folder holding inputs.py, e.g.
    cd explore/ladder/nr0010 && uv run python ../../run.py [--taper 0.2] [-m '{"dt": 0.25}']
"""
import runpy
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO/'tools'))
import patches


def flag(name, value=False):
    """Remove --name (and its value, if value) from sys.argv; return the value, True, or None."""
    if name not in sys.argv:
        return None
    k = sys.argv.index(name)
    v = sys.argv[k + 1] if value else True
    del sys.argv[k:k + 1 + value]
    return v


fix_edge = 'own' if flag('--fix-edge-own') else ('hahn13' if flag('--fix-edge') else None)
taper = flag('--taper', value=True)
patches.apply(taper=None if taper is None else float(taper), fix_edge=fix_edge)
sys.argv[0] = str(REPO/'source'/'epi_int_lite.py')
runpy.run_path(sys.argv[0], run_name='__main__')
