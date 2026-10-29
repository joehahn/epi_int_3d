#test_kernels.py
#
#Unit tests of fast/kernels.py's ring gravity on a circular ringlet, where the wire sums are known in
#closed form: the reference kernel pulls with the perturbed particle's lambda (epi_int_lite issue 7),
#the fixed kernel with each perturbing streamline's own. And the nearpoint and shape_slope options
#(notes/equations.md G14) on m = 2 streamlines whose velocities are not along their shapes, as for a
#pattern rotating at Omega_p = Omega/2. And the 3D kernel (G16, G17; the paper's Eqn. 23): at z = 0 it is
#the planar nearpoint and shape_slope kernel bit for bit; out of the plane it matches the near point
#computed as a vector, Delta = d0 - (d0.t)t; and a rigidly tilted flat ringlet is pulled within its plane.
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


def test_3d_in_plane_is_planar():
    r, s_r, vr, vt, (lw, tw, rw, vrw, vtw) = m2_ringlet()
    zw = np.zeros_like(rw)
    for fix in (False, True):
        for parallel in (False, True):
            Ar, At = kernels.get_kernel(fix, True, True, parallel=parallel)(lw, 1.0, rw, tw, vrw, vtw)
            Ar3, At3, Az3 = kernels.get_kernel_3d(fix, parallel=parallel)(lw, 1.0, rw, tw, zw)
            assert np.array_equal(Ar, Ar3) and np.array_equal(At, At3) and not Az3.any()


def warped_ringlet(Nr=6, Nt=2000, m=2, a0=1.0, width=0.5, b0=0.02, h0=0.03):
    """m2_ringlet's streamlines r_s = a_s + b_s cos(m theta), also warped: z_s = h_s sin(m theta + s)."""
    a = a0 + width*np.arange(Nr)/(Nr - 1)
    b = b0*(1 + np.arange(Nr)/Nr)
    h = h0*(1 + 0.5*np.arange(Nr)/Nr)
    ph = np.arange(Nr)[:, None]
    t = np.tile(-np.pi + 2*np.pi*(np.arange(Nt) + 0.5)/Nt, (Nr, 1))
    r = a[:, None] + b[:, None]*np.cos(m*t)
    z = h[:, None]*np.sin(m*t + ph)
    s_r = -m*b[:, None]*np.sin(m*t)/r                  #(dr/dtheta)/r
    s_z = m*h[:, None]*np.cos(m*t + ph)/r              #(dz/dtheta)/r
    lam = np.full_like(r, 1.0e-10)
    return r, z, s_r, s_z, [wrap_ring(x, longitude=(k == 1)) for k, x in enumerate((lam, t, r, z))]


def wire_pull(G, lam, dr, dz, s_r, s_z):
    """2 G lam Delta/Delta^2 with Delta = d0 - (d0.t)t, d0 = (dr, 0, dz), t = (s_r, 1, s_z)/N: the near point
    as a vector (G16), not the kernel's components; arrays broadcast, the last axis is (r, theta, z)."""
    d0 = np.stack([dr, np.zeros_like(dr), dz], axis=-1)
    tt = np.stack([s_r, np.ones_like(s_r), s_z], axis=-1)
    tt /= np.linalg.norm(tt, axis=-1, keepdims=True)
    D = d0 - (d0*tt).sum(axis=-1, keepdims=True)*tt
    with np.errstate(invalid='ignore', divide='ignore'):     #a streamline's own term, D = 0, is masked later
        return 2*G*lam*D/(D*D).sum(axis=-1, keepdims=True)


def test_3d_wire():
    r, z, s_r, s_z, (lw, tw, rw, zw) = warped_ringlet()
    G, lam = 1.0, 1.0e-10
    Nr = r.shape[0]
    off = ~np.eye(Nr, dtype=bool)[:, :, None, None]
    #pull on particle (i, j) from streamline s, all at the particle's longitude: [i, s, j, component]
    A = wire_pull(G, lam, r[None, :, :] - r[:, None, :], z[None, :, :] - z[:, None, :], s_r[None, :, :] + 0*r[:, None, :],
                  s_z[None, :, :] + 0*r[:, None, :])
    expect = np.where(off, A, 0).sum(axis=1)
    Ar, At, Az = kernels.get_kernel_3d(False, parallel=False)(lw, G, rw, tw, zw)
    scale = np.abs(expect).max()
    for got, k in ((Ar, 0), (At, 1), (Az, 2)):
        #3-point Lagrange errors: measured 4e-9 (r), 3e-7 (theta), 3e-8 (z) of the largest pull
        np.testing.assert_allclose(got[:, 1:-1], expect[:, :, k], rtol=0, atol=1e-5*scale)
    #the out-of-plane terms matter here: the planar kernel misses them by far more than that
    Ar2, At2 = kernels.get_kernel(False, True, True, parallel=False)(lw, G, rw, tw, rw*0, rw*0 + 1)
    assert np.abs(Ar2[:, 1:-1] - expect[:, :, 0]).max() > 1e-2*np.abs(expect[:, :, 0]).max()
    Ar_p, At_p, Az_p = kernels.get_kernel_3d(False, parallel=True)(lw, G, rw, tw, zw)
    assert np.array_equal(Ar, Ar_p) and np.array_equal(At, At_p) and np.array_equal(Az, Az_p)


def test_3d_tilted_ringlet():
    """A flat circular ringlet tilted by inc about the x axis, z = r tan(inc) sin(theta): every wire and every
    particle lies in that plane, so each pull does too, A.n = 0 with the plane's normal
    n = z-hat - tan(inc)(sin(theta) r-hat + cos(theta) theta-hat), up to the Lagrange fit's errors."""
    a, lam, (lw, tw, rw, vrw, vtw) = circular_ringlet(Nt=400)
    inc = 0.1
    zw = rw*np.tan(inc)*np.sin(tw)
    Ar, At, Az = kernels.get_kernel_3d(True, parallel=False)(lw, 1.0, rw, tw, zw)
    An = Az - np.tan(inc)*(np.sin(tw)*Ar + np.cos(tw)*At)
    A = np.sqrt(Ar**2 + At**2 + Az**2)
    assert np.abs(An[:, 1:-1]).max() < 5e-7*A[:, 1:-1].max()         #measured 1.6e-8
    assert np.abs(Az[:, 1:-1]).max() > 0.05*A[:, 1:-1].max()        #and the pull does have a z part
