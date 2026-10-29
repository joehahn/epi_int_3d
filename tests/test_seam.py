#test_seam.py
#
#Unit tests of tools/patches.py's seam fix (epi_int_lite issue 10, notes/epi_int_lite_issues.md): pressure and
#viscosity computed on the ring wrapped with two ghost columns per side. The accelerations of a ring cannot depend
#on where its longitudes start, so rotating the whole ring by alpha must leave each particle's acceleration
#unchanged. Without the fix it does not for the particles beside theta = +-pi; with it, it does, and the particles
#away from the seam get the reference's accelerations bit for bit.
#
#Run:  uv run pytest tests/test_seam.py

import sys

import numpy as np

from analysis import REPO

sys.path.insert(0, str(REPO/'source'))
sys.path.insert(0, str(REPO/'tools'))
import helper_fns
import patches

REFERENCE = helper_fns.accelerations                   #as imported, before any patch


def eccentric_ringlet(Nr=4, Nt=60, a0=1.0, width=1.0e-3, e=0.01, e_prime=0.4, J2=0.01, Rp=0.5):
    """Nested eccentric streamlines with an eccentricity gradient, and each streamline's particles offset in
    longitude from the others', so that interpolating between streamlines is not trivial."""
    a = a0 + width*np.arange(Nr)/(Nr - 1)
    ecc = e + e_prime*(a - a.mean())
    M = -np.pi + 2*np.pi*(np.arange(Nt)[None, :] + 0.5 + 0.3*np.arange(Nr)[:, None]/Nr)/Nt
    r = a[:, None]*(1 - ecc[:, None]*np.cos(M))
    t = helper_fns.adjust_angle(M + 2*ecc[:, None]*np.sin(M))
    Omg = np.sqrt((1.0/a**3)*(1 + 1.5*J2*(Rp/a)**2))
    vr = (a*Omg*ecc)[:, None]*np.sin(M)
    vt = (a*Omg)[:, None]*(1 + ecc[:, None]*np.cos(M))
    order = np.argsort(t, axis=1)
    r, t, vr, vt = (np.take_along_axis(x, order, axis=1) for x in (r, t, vr, vt))
    lam = helper_fns.get_lambda(1.0e-10, Nr, J2, Rp, r, t, vr, vt)
    return lam, r, t, vr, vt


def rotated(accelerations, alpha, lam, r, t, vr, vt, nu=1.0e-13):
    """accelerations of the ring rotated by alpha (viscosity only), mapped back to the original particles."""
    t2 = helper_fns.adjust_angle(t + alpha)
    o = np.argsort(t2, axis=1)
    inv = np.argsort(o, axis=1)
    tk = lambda x: np.take_along_axis(x, o, axis=1)
    Ar, At = accelerations(tk(lam), -1.0, nu, nu, -1.0, tk(r), tk(t2), tk(vr), tk(vt), False, False, False)
    return np.take_along_axis(Ar, inv, axis=1), np.take_along_axis(At, inv, axis=1)


def test_seam():
    ring = eccentric_ringlet()
    Nt = ring[1].shape[1]
    seam = np.zeros(Nt, dtype=bool)
    seam[[0, 1, -2, -1]] = True                        #the two particles beside theta = +-pi, either side
    for fix_viscosity in (False, True):
        old = patches.accelerations if fix_viscosity else REFERENCE
        new = patches.seam_accelerations(fix_viscosity)
        A_old, A_new = rotated(old, 0.0, *ring), rotated(new, 0.0, *ring)
        for alpha in (1.234, -2.5):
            A_old_rot, A_new_rot = rotated(old, alpha, *ring), rotated(new, alpha, *ring)
            for k in range(2):
                scale = np.abs(A_new[k]).max()
                #without the fix, rotating the ring changes the forces beside the old seam
                assert np.abs(A_old_rot[k] - A_old[k])[:, seam].max() > 1e-6*scale
                #with it, the forces are independent of where the longitudes start
                assert np.abs(A_new_rot[k] - A_new[k]).max() < 1e-9*scale
        #and away from the seam the fix changes nothing, bit for bit
        for k in range(2):
            assert np.array_equal(A_new[k][:, ~seam], A_old[k][:, ~seam])
            assert not np.array_equal(A_new[k][:, seam], A_old[k][:, seam])


def test_wrap():
    x = np.arange(12.0).reshape(2, 6)
    assert np.array_equal(patches.wrap(x, 1), helper_fns.wrap_ring(x))
    assert np.array_equal(patches.wrap(x, 1, longitude=True), helper_fns.wrap_ring(x, longitude=True))
    assert np.array_equal(patches.wrap(x, 2, longitude=True)[:, 1:-1], helper_fns.wrap_ring(x, longitude=True))
