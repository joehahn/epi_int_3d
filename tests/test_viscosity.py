#test_viscosity.py
#
#Unit test of tools/patches.py's shear-torque fix (epi_int_lite issue 12, notes/epi_int_lite_issues.md), which
#fast/epi_int_3d.py applies by default. In a circular Keplerian ring with uniform surface density sigma, the shear
#viscosity's angular momentum flux is F = -nu_s sigma r^2 d(omega)/dr = (3/2) nu_s sigma r^(-1/2), so hahn13's
#Eq. (27) gives the tangential acceleration A = -(1/sigma r) dF/dr = (3/4) nu_s r^(-5/2) inside the ring, and
#-F/lambda r at the inner edge (notes/equations.md W1). The reference's -(1/sigma) d(F/r)/dr gives three times
#that inside the ring.
#
#Run:  uv run pytest tests/test_viscosity.py

import sys

import numpy as np

from analysis import REPO

sys.path.insert(0, str(REPO/'source'))
sys.path.insert(0, str(REPO/'tools'))
import helper_fns
import patches

REFERENCE = helper_fns.ring_shear_viscosity              #as imported, before any patch


def uniform_ring(Nr=101, Nt=64, width=0.01, nu=1.0e-12, sigma=1.0e-6):
    """Circular Keplerian streamlines (G = M = 1) spanning 1 <= r <= 1 + width, with uniform surface density."""
    r1 = np.linspace(1.0, 1.0 + width, Nr)
    t1 = np.linspace(-np.pi, np.pi, Nt, endpoint=False)
    r = np.tile(r1[:, None], (1, Nt))
    t = np.tile(t1, (Nr, 1))
    vr = 0.0*r
    vt = r**-0.5
    lam = sigma*(r1[1] - r1[0]) + 0.0*r
    return nu, sigma, lam, r, t, vr, vt


def tangential(shear, nu, lam, r, t, vr, vt):
    """At of helper_fns.accelerations, viscosity only, with ring_shear_viscosity set to shear."""
    saved = helper_fns.ring_shear_viscosity
    helper_fns.ring_shear_viscosity = shear
    try:
        Ar, At = helper_fns.accelerations(lam, -1.0, nu, -1.0, -1.0, r, t, vr, vt, False, False, False)
    finally:
        helper_fns.ring_shear_viscosity = saved
    return At


def test_shear_torque():
    nu, sigma, lam, r, t, vr, vt = uniform_ring()
    r1 = r[:, 0]
    At_fixed = tangential(patches.ring_shear_viscosity, nu, lam, r, t, vr, vt)
    At_reference = tangential(REFERENCE, nu, lam, r, t, vr, vt)
    interior = slice(2, -2)
    #the fix gives hahn13's Eq. (27) inside the ring, to the finite differences' precision (atol = 0: A ~ 1e-12)
    hahn13 = 0.75*nu*r1**-2.5
    assert np.allclose(At_fixed[interior], hahn13[interior, None], rtol=1.0e-5, atol=0.0)
    #the reference's interior value is three times that
    assert np.allclose(At_reference[interior], 3.0*hahn13[interior, None], rtol=1.0e-5, atol=0.0)
    #at the inner edge both give -F/lambda r
    F = 1.5*nu*sigma*r1**-0.5
    edge = -F[0]/(lam[0, 0]*r1[0])
    assert np.allclose(At_fixed[0], edge, rtol=1.0e-3, atol=0.0)
    assert np.allclose(At_reference[0], edge, rtol=1.0e-3, atol=0.0)
    #the patch is applied, and undone, as fast/epi_int_3d.py and this test expect
    assert helper_fns.ring_shear_viscosity is REFERENCE
