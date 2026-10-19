#test_regression.py
#
#The 2D regression suite: shortened versions of epi_int_lite's unperturbed, spreading and
#libration tests, run through the driver source/epi_int_lite.py (or the driver named by the
#environment variable EPI_INT_DRIVER, e.g. fast/epi_int_fast.py). Each test checks
#  1. physics: the diagnostics of epi_int_lite's test notebooks (see tests/analysis.py), and
#  2. regression: the final state matches tests/golden/<test>.npz, written by the Python 3 port
#     of epi_int_lite (the reference implementation) to relative precision GOLDEN_RTOL.
#
#Run:  uv run pytest
#Regenerate the golden files (only when the reference implementation itself changes):
#      UPDATE_GOLDEN=1 uv run pytest
#Check the faster code against the same golden files:
#      EPI_INT_DRIVER=fast/epi_int_fast.py uv run pytest

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from analysis import (REPO, load, unperturbed_metrics, spreading_metrics, libration_metrics)

TESTS = REPO/'tests'
GOLDEN = TESTS/'golden'
GOLDEN_RTOL = 1.0e-11
DRIVER = REPO/os.environ.get('EPI_INT_DRIVER', 'source/epi_int_lite.py')


def run_sim(test_folder, output_folder, **params):
    """Run the driver (source/epi_int_lite.py by default) in test_folder, with inputs.py modified by params.

    Warnings (e.g. NumPy deprecations) are raised as errors."""
    params['output_folder'] = str(output_folder)
    cmd = [sys.executable, '-W', 'error', str(DRIVER),
           '-m', json.dumps(params)]
    result = subprocess.run(cmd, cwd=test_folder, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return load(test_folder, output_folder, params)


def check_golden(name, o):
    """Compare the final state with the golden file, or write it if UPDATE_GOLDEN is set."""
    path = GOLDEN/f'{name}.npz'
    state = {k: o[k][-1] for k in ('r', 't', 'vr', 'vt')}
    if os.environ.get('UPDATE_GOLDEN'):
        assert DRIVER == REPO/'source'/'epi_int_lite.py', 'golden files come from the reference implementation'
        GOLDEN.mkdir(exist_ok=True)
        np.savez(path, times=o['times'], **state)
        pytest.skip(f'wrote {path.name}')
    golden = np.load(path)
    np.testing.assert_array_equal(o['times'], golden['times'])
    for k, x in state.items():
        #scale by each array's typical size; angles t are compared modulo 2 pi
        diff = np.angle(np.exp(1j*(x - golden[k]))) if k == 't' else x - golden[k]
        scale = 1.0 if k == 't' else np.abs(golden[k]).max()
        assert np.abs(diff).max() <= GOLDEN_RTOL*scale, f'{k} differs from {path.name}'


def test_unperturbed(tmp_path):
    """Unperturbed particles keep their a, e and angular momentum h (check_orbits.ipynb)."""
    #2000 orbits of the full test's 10^5, 51 streamlines with 1e-8 < e < 1e-2
    p, o = run_sim(TESTS/'unperturbed', tmp_path, total_number_of_outputs=20)
    m = unperturbed_metrics(p, o)
    assert m['e_min'] < 1.0e-6 and m['e_max'] > 1.0e-3
    #the port gives about 7e-13, 2e-13 and 3e-13 here; over the full 1e5 orbits the notebook's
    #figures show changes growing to about 1e-11
    assert m['max_delta_a'] < 5.0e-12
    assert m['max_delta_e'] < 5.0e-12
    assert m['max_delta_h'] < 5.0e-12
    check_golden('unperturbed', o)


def test_spreading(tmp_path):
    """A narrow viscous ringlet spreads like the analytic solution, conserving L (spreading.ipynb)."""
    #run to 2.5 viscous timescales of the full test's 50
    p, o = run_sim(TESTS/'spreading', tmp_path, total_number_of_outputs=50)
    m = spreading_metrics(p, o, t_visc=(1.0, 2.5))
    #the notebook's full run conserved L to 1.1e-7. The 15 streamlines resolve the analytic
    #profile only coarsely, and it assumes an initially infinitely narrow ringlet, so the port's
    #interior streamlines agree with it to about 12% here, as in the notebook's figure
    assert m['max_delta_L'] < 1.0e-7
    assert m['sigma_rel_err_t1'] < 0.15
    assert m['sigma_rel_err_t2.5'] < 0.15
    check_golden('spreading', o)


def test_libration(tmp_path):
    """A narrow eccentric ringlet librates at the BGT frequency (libration.ipynb)."""
    #the full test: 2 streamlines of 241 particles, ring mass 1e-9, 1.6e4 orbits
    p, o = run_sim(TESTS/'libration'/'narrow_ring'/'ring_mass_1.0e-9', tmp_path)
    m = libration_metrics(p, o)
    #the notebook's Python 2 run gave omega/omega_bgt = 1.00554 and H = 0.50032
    assert m['omega_over_omega_bgt'] == pytest.approx(1.0055419800, rel=1.0e-6)
    assert m['H_mean'] == pytest.approx(0.5003165495, rel=1.0e-6)
    assert m['e_prime'] == pytest.approx(0.0289001766, rel=1.0e-6)
    check_golden('libration', o)
