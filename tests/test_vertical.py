#test_vertical.py
#
#Regression tests of epi_int_3d's vertical motion, run through fast/epi_int_3d.py (the reference
#implementation is planar). The z = 0 path is tested by test_regression.py with
#EPI_INT_DRIVER=fast/epi_int_3d.py, which must reproduce the reference bit for bit.
#
#Run:  uv run pytest tests/test_vertical.py

import json
import subprocess
import sys

import numpy as np
import pytest

from analysis import REPO, load
from vertical_analysis import load3d, inclined_metrics, exact_orbits, forced_vertical

TESTS = REPO/'tests'
DRIVER = REPO/'fast'/'epi_int_3d.py'


def run_3d(test_folder, output_folder, **params):
    params['output_folder'] = str(output_folder)
    cmd = [sys.executable, '-W', 'error', str(DRIVER), '-m', json.dumps(params)]
    result = subprocess.run(cmd, cwd=test_folder, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return params


def test_inclined(tmp_path):
    """Inclined test particles keep their inclination, and their nodes regress at Omega_a - nu_a
    (tests/inclined/inclined.ipynb); the in-plane motion is unchanged."""
    #200 orbits of the full test's 1000: 51 streamlines, 1e-8 < e < 1e-2, i = 0.01
    params = run_3d(TESTS/'inclined', tmp_path/'inc', total_number_of_outputs=20)
    p, o = load3d(TESTS/'inclined', tmp_path/'inc', params)
    m = inclined_metrics(p, o)
    #the drift is exact for its own equations (V8): about 1e-13 and 1e-11 here
    assert m['max_delta_inc_rel'] < 1.0e-11
    assert m['max_delta_node'] < 1.0e-9
    #the same particles integrated exactly in the 3D J2 potential, those with e < 1e-4: their
    #nodal rates are slower by a factor cos^8(i) ~ 1 - 4 i^2, because the code's a comes from
    #h_z = h cos(i) (tests/inclined/inclined.ipynb)
    small_e = np.flatnonzero(o['e'][0][:, 0] < 1.0e-4)[:8]
    idx = (small_e, np.zeros_like(small_e))
    node_exact, node_code = exact_orbits(p, o, idx)
    T = o['times'] - o['times'][0]
    d = np.polyfit(T, node_exact, 1)[0]/np.polyfit(T, node_code, 1)[0] - 1
    #(each particle's estimate scatters by ~10% over these 200 orbits; their mean by ~1%)
    assert np.mean(d) == pytest.approx(-4*p['inclination']**2, rel=0.03)
    assert np.allclose(d, -4*p['inclination']**2, rtol=0.2)
    #inclination does not change the in-plane motion: compare with the same run uninclined
    run_3d(TESTS/'inclined', tmp_path/'flat', total_number_of_outputs=20, inclination=0.0)
    flat = load(TESTS/'inclined', tmp_path/'flat', params)[1]
    inc = load(TESTS/'inclined', tmp_path/'inc', params)[1]
    for k in ('r', 't', 'vr', 'vt'):
        np.testing.assert_array_equal(inc[k], flat[k])
    assert not np.load(tmp_path/'flat'/'z.npy').any() and not np.load(tmp_path/'flat'/'vz.npy').any()


def test_vertical_resonance(tmp_path):
    """A zero-mass ring under a growing m = 2 vertical Fourier satellite (fast/satellite.py): each particle's z(t)
    is the exact solution B11, settling into the forced oscillation B6 (tests/vertical_resonance)."""
    #the first 60 of the full test's 200 outputs, t = 1.8e4, 1.8 growth times
    params = run_3d(TESTS/'vertical_resonance', tmp_path/'vr', total_number_of_outputs=60)
    p, o = load3d(TESTS/'vertical_resonance', tmp_path/'vr', params)
    sat, f = forced_vertical(p, o)
    assert sat.r_V == pytest.approx(1.02, rel=1e-12)
    err = np.abs(o['z'] - f['z']).max(axis=(0, 2))/np.abs(f['z']).max(axis=(0, 2))
    near = np.abs(f["x"]) < 1.0e-6                     #the streamline at r_V
    #away from resonance the leapfrog's O(dt^2) error, ~3e-11 absolutely: up to 1.3e-3 of max|z| at |x| = 0.02
    assert err[~near].max() < 2.0e-3
    #the streamline at x = 0 (nu - omega_D ~ 1e-14), growing secularly: 8e-7 in the full run
    assert err[near].max() < 1.0e-5
    #no satellite, no vertical motion
    run_3d(TESTS/'vertical_resonance', tmp_path/'off', total_number_of_outputs=2,
           satellite=dict(p['satellite'], mass_final=-1.0))
    assert not np.load(tmp_path/'off'/'z.npy').any()
