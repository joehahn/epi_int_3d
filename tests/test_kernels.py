#test_kernels.py
#
#Unit tests of fast/kernels.py's ring gravity on a circular ringlet, where the wire sums are known in
#closed form: the reference kernel pulls with the perturbed particle's lambda (epi_int_lite issue 7),
#the fixed kernel with each perturbing streamline's own.
#
#Run:  uv run pytest tests/test_kernels.py

import sys

import numpy as np

from analysis import REPO

sys.path.insert(0, str(REPO/'source'))
sys.path.insert(0, str(REPO/'fast'))
from helper_fns import get_lambda, wrap_ring
import kernels


def circular_ringlet(Nr=7, Nt=12, a0=1.0, width=1.0e-3, m_ring=1.0e-9, J2=0.01, Rp=0.5):
    """Circular streamlines, equal masses, particles at the same longitudes in every streamline."""
    a = a0 + width*np.arange(Nr)/(Nr - 1)
    t = np.tile(-np.pi + 2*np.pi*(np.arange(Nt) + 0.5)/Nt, (Nr, 1))
    r = np.tile(a[:, None], (1, Nt))
    vr = np.zeros_like(r)
    Omg = np.sqrt((1.0/a**3)*(1 + 1.5*J2*(Rp/a)**2))
    vt = np.tile((a*Omg)[:, None], (1, Nt))
    lam = get_lambda(m_ring, Nr, J2, Rp, r, t, vr, vt)
    return a, lam, [wrap_ring(x, longitude=(k == 1)) for k, x in enumerate((lam, t, r, vr, vt))]


def test_perturber_lambda():
    a, lam, (lw, tw, rw, vrw, vtw) = circular_ringlet()
    G = 1.0
    da = a[None, :] - a[:, None]                       #a_s - a_i
    off = ~np.eye(len(a), dtype=bool)
    lam_s = lam[:, 0]                                  #each streamline's lambda (all equal along it)
    expect_ref = np.where(off, 2*G*lam_s[:, None]/np.where(off, da, 1), 0).sum(axis=1)
    expect_fix = np.where(off, 2*G*lam_s[None, :]/np.where(off, da, 1), 0).sum(axis=1)
    tol = dict(rtol=1e-12, atol=1e-12*np.abs(2*G*lam_s/(a[1] - a[0])).max())   #sums cancel mid-ringlet
    for ref, fix in ((kernels.ring_gravity_serial, kernels.ring_gravity_lam_serial),
                     (kernels.ring_gravity_parallel, kernels.ring_gravity_lam_parallel)):
        Ar, At = ref(lw, G, rw, tw, vrw, vtw)
        Ar_fix, At_fix = fix(lw, G, rw, tw, vrw, vtw)
        np.testing.assert_allclose(Ar[:, 1:-1], np.tile(expect_ref[:, None], (1, Ar.shape[1] - 2)), **tol)
        np.testing.assert_allclose(Ar_fix[:, 1:-1], np.tile(expect_fix[:, None], (1, Ar.shape[1] - 2)), **tol)
        assert not At.any() and not At_fix.any()
    #the fix's correction, sum_s 2 G (lambda_s - lambda_i)/(a_s - a_i): for circular equal-mass
    #streamlines lambda = m_1/(2 pi a), so it is -sum_s G m_1/(pi a_s a_i), a nearly uniform inward pull
    m1 = lam_s*2*np.pi*a
    corr = -np.where(off, G*m1[None, :]/(np.pi*a[None, :]*a[:, None]), 0).sum(axis=1)
    np.testing.assert_allclose(Ar_fix[:, 1] - Ar[:, 1], corr, rtol=1e-6)
