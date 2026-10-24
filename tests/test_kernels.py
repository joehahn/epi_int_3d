#test_kernels.py
#
#Unit tests of fast/kernels.py's ring gravity on a circular ringlet, where the wire sums are known in
#closed form: the reference kernel pulls with the perturbed particle's lambda (epi_int_lite issue 7),
#the fixed kernel with each perturbing streamline's own. And the nearpoint and shape_slope options
#(notes/equations.md G14) on m = 2 streamlines whose velocities are not along their shapes, as for a
#pattern rotating at Omega_p = Omega/2.
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


def m2_ringlet(Nr=6, Nt=2000, m=2, a0=1.0, width=0.5, b0=0.02, slowdown=0.5):
    """Streamlines r_s(theta) = a_s + b_s cos(m theta), with velocities along slowdown x the shape's slope:
    v_r/v_theta = slowdown x s_r, as for a pattern rotating at (1 - slowdown) Omega (G6). The slopes,
    s_r ~ 0.1, are far steeper than a ring's, so that cos(phi) = 1 - O(s_r^2) is resolved."""
    a = a0 + width*np.arange(Nr)/(Nr - 1)
    b = b0*(1 + np.arange(Nr)/Nr)
    t = np.tile(-np.pi + 2*np.pi*(np.arange(Nt) + 0.5)/Nt, (Nr, 1))
    r = a[:, None] + b[:, None]*np.cos(m*t)
    s_r = -m*b[:, None]*np.sin(m*t)/r                  #(dr/dtheta)/r
    vt = np.sqrt(1.0/r)
    vr = slowdown*s_r*vt
    lam = np.full_like(r, 1.0e-10)
    return r, s_r, vr, vt, [wrap_ring(x, longitude=(k == 1)) for k, x in enumerate((lam, t, r, vr, vt))]


def test_wire_geometry_options():
    r, s_r, vr, vt, (lw, tw, rw, vrw, vtw) = m2_ringlet()
    G, lam = 1.0, 1.0e-10
    dr = r[None, :, :] - r[:, None, :]                 #dr[i, s] = r_s - r_i at the particles' longitudes
    off = ~np.eye(r.shape[0], dtype=bool)[:, :, None]
    inv = np.where(off, 1/np.where(off, dr, 1), 0)
    pull = 2*G*lam*inv                                 #2 G lambda/d_r from each streamline s
    slope = {True: s_r, False: vr/vt}                  #tan(phi) from the shape, or from the velocity
    for nearpoint in (False, True):
        for shape in (False, True):
            tph = slope[shape][None, :, :]
            cph = 1.0 if nearpoint else 1/np.sqrt(1 + tph**2)
            Ar_expect = (pull*cph).sum(axis=1)
            At_expect = -(pull*tph*cph).sum(axis=1)
            kern = kernels.get_kernel(False, nearpoint, shape, parallel=False)
            Ar, At = kern(lw, G, rw, tw, vrw, vtw)
            #3-point Lagrange errors: r' to O(b (m dtheta)^3), the slope to O((m dtheta)^2)
            np.testing.assert_allclose(Ar[:, 1:-1], Ar_expect, rtol=1e-6)
            np.testing.assert_allclose(At[:, 1:-1], At_expect, rtol=1e-3, atol=1e-6*np.abs(At_expect).max())
            Ar_p, At_p = kernels.get_kernel(False, nearpoint, shape, parallel=True)(lw, G, rw, tw, vrw, vtw)
            assert np.array_equal(Ar, Ar_p) and np.array_equal(At, At_p)
    #nearpoint drops the cos(phi) factor: here a few parts in 10^3
    Ar_c = kernels.get_kernel(False, False, True, parallel=False)(lw, G, rw, tw, vrw, vtw)[0]
    Ar_n = kernels.get_kernel(False, True, True, parallel=False)(lw, G, rw, tw, vrw, vtw)[0]
    assert np.abs(Ar_n[:, 1:-1]/Ar_c[:, 1:-1] - 1).max() > 1e-3
    #the velocity's slope is half the shape's here, and so is the torque it gives
    At_v = kernels.get_kernel(False, True, False, parallel=False)(lw, G, rw, tw, vrw, vtw)[1]
    At_s = kernels.get_kernel(False, True, True, parallel=False)(lw, G, rw, tw, vrw, vtw)[1]
    np.testing.assert_allclose(At_v[:, 1:-1], 0.5*At_s[:, 1:-1], rtol=1e-3, atol=1e-6*np.abs(At_s).max())
    #with no options, get_kernel is the reference kernel
    assert kernels.get_kernel() is kernels.ring_gravity_parallel
    assert kernels.get_kernel(True, parallel=False) is kernels.ring_gravity_lam_serial
