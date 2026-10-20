"""Numba kernels for the bookkeeping of epi_int_3d's driver (Phase 0b, step 4).

Each timestep of epi_int_lite converts the ring's coordinates many times: polar to Cartesian and
back for the central body's motion (coordinate_kick), to orbit elements and back for the drift,
between planetocentric and barycentric frames around every velocity kick (mixed2planeto,
planeto2mixed), and to elements again for each streamline's linear density (get_lambda). The
reference does each as whole-array NumPy expressions with many temporary arrays. Here the
elementwise work runs as compiled loops that write into preallocated arrays, and the driver does
the few whole-ring sums between them with NumPy, so that every sum is NumPy's own (pairwise)
summation, as in the reference.

Every expression is written with the reference's operations in the reference's order (e.g.
a*a*Omg/r is ((a*a)*Omg)/r), and there is no fastmath, so the results are bit-for-bit the
reference's. One pitfall: when sin and cos of the same angle are computed in one loop, LLVM fuses
them into a single sincos call, whose sin differs from NumPy's (and libm's) sin in the last bit for
~0.3% of angles. So sin and cos of an angle are always taken in separate loops.

The vertical drift (vertical_drift, V8 of notes/equations.md) has no reference to match; it
follows the same rules.

The kernels work on flattened (1D) views of the (Nr, Nt) arrays, except sort_streamlines, wrap and
the crossing check, which need the streamlines. Each is compiled twice, serial (SERIAL) and with its
loops spread over all cores (PARALLEL); the loops' iterations are independent, so both give the
same results. Threads pay off only for large rings (see fast/epi_int_3d.py).
"""
import numba
import numpy as np
from numba import prange

TWOPI = 2.0*np.pi


@numba.njit
def _adjust(angle):
    """helper_fns.adjust_angle for one angle: into [-pi, pi]."""
    if angle > np.pi:
        angle -= TWOPI
    if angle < -np.pi:
        angle += TWOPI
    return angle


@numba.njit
def _omega_kappa(a, J2_15, Rp2):
    """helper_fns.Omega and Kappa at a, with Ar = 0; J2_15 = 1.5*J2, Rp2 = Rp*Rp."""
    a2 = a*a
    Ra2 = Rp2/a2
    f = 1.0/a2/a
    return np.sqrt(f*(1.0 + J2_15*Ra2)), np.sqrt(f*(1.0 - J2_15*Ra2))


def rt2xy(r, t, vr, vt, x, y, vx, vy, s, c):
    """helper_fns.rt2xy into x, y, vx, vy (s, c are work arrays)."""
    for k in prange(t.size):              #sin and cos in separate loops, see the module docstring
        s[k] = np.sin(t[k])
    for k in prange(t.size):
        c[k] = np.cos(t[k])
    for k in prange(r.size):
        x[k] = r[k]*c[k]
        y[k] = r[k]*s[k]
        vx[k] = vr[k]*c[k] - vt[k]*s[k]
        vy[k] = vr[k]*s[k] + vt[k]*c[k]


def xy2rt(x, y, vx, vy, r, t, vr, vt, s, c):
    """helper_fns.xy2rt into r, t, vr, vt (s, c are work arrays)."""
    for k in prange(x.size):
        r[k] = np.sqrt(x[k]*x[k] + y[k]*y[k])
        t[k] = np.arctan2(y[k], x[k])
    for k in prange(t.size):              #sin and cos in separate loops, see the module docstring
        s[k] = np.sin(t[k])
    for k in prange(t.size):
        c[k] = np.cos(t[k])
    for k in prange(x.size):
        vr[k] = vx[k]*c[k] + vy[k]*s[k]
        vt[k] = -vx[k]*s[k] + vy[k]*c[k]


def polar(x, y, r, t):
    """The r, t of helper_fns.xy2rt."""
    for k in prange(x.size):
        r[k] = np.sqrt(x[k]*x[k] + y[k]*y[k])
        t[k] = np.arctan2(y[k], x[k])


def shift(x, y, dx, dy):
    """x += dx, y += dy."""
    for k in prange(x.size):
        x[k] = x[k] + dx
        y[k] = y[k] + dy


def shift4(x, y, vx, vy, dx, dy, dvx, dvy):
    """x += dx, ... vy += dvy."""
    for k in prange(x.size):
        x[k] = x[k] + dx
        y[k] = y[k] + dy
        vx[k] = vx[k] + dvx
        vy[k] = vy[k] + dvy


def drift(r, t, vr, vt, dt, J2_15, Rp2, Rp, two_GM_Rp, a, e, wt, M, s, c):
    """coords2elem, drift (M += kappa dt) and elem2coords, in place on r, t, vr, vt."""
    for k in prange(r.size):
        #coords2elem
        h = r[k]*vt[k]
        cc = (h*h)/two_GM_Rp
        ak = Rp*(cc + np.sqrt(cc*cc - J2_15))
        Omg, Kap = _omega_kappa(ak, J2_15, Rp2)
        ar = ak/r[k]
        ar2 = ar*ar
        ar3 = ar2*ar
        e_cos_M = (ar2 - 1.0)/2.0
        aK = ak*Kap
        e_sin_M = vr[k]*(ar3/aK)
        e[k] = np.sqrt(e_sin_M*e_sin_M + e_cos_M*e_cos_M)
        Mk = np.arctan2(e_sin_M, e_cos_M)
        wt[k] = _adjust(t[k] - (Omg/Kap)*(Mk + 2.0*e_sin_M))
        #drift
        M[k] = Mk + Kap*dt
        a[k] = ak
    for k in prange(M.size):
        s[k] = np.sin(M[k])
    for k in prange(M.size):
        c[k] = np.cos(M[k])
    for k in prange(r.size):
        #elem2coords
        ak = a[k]
        e_sin_M = e[k]*s[k]
        e_cos_M = e[k]*c[k]
        rk = ak/np.sqrt(1.0 + 2.0*e_cos_M)
        Omg, Kap = _omega_kappa(ak, J2_15, Rp2)
        t[k] = _adjust((Omg/Kap)*(M[k] + 2.0*e_sin_M) + wt[k])
        ra = rk/ak
        ra3 = ra*ra*ra
        vr[k] = (ak*Kap*ra3)*e_sin_M
        vt[k] = ak*ak*Omg/rk
        r[k] = rk


def shift1(x, dx):
    """x += dx."""
    for k in prange(x.size):
        x[k] = x[k] + dx


def vertical_drift(a, z, vz, dt, J2_45, Rp2, nu, s, c):
    """The vertical part of the drift (notes/equations.md V8): rotate (z, vz/nu) by nu*dt, with
    nu = nu_a (V1) at the semimajor axes a that drift() just computed; J2_45 = 4.5*J2. This is
    V2-V4 exactly: the node regresses at Omega_a - nu_a while the mean longitude advances at
    Omega_a. z = vz = 0 stays exactly 0. (nu, s, c are work arrays.)"""
    for k in prange(a.size):
        ak = a[k]
        a2 = ak*ak
        nu[k] = np.sqrt((1.0/a2/ak)*(1.0 + J2_45*(Rp2/a2)))
    for k in prange(a.size):              #sin and cos in separate loops, see the module docstring
        s[k] = np.sin(nu[k]*dt)
    for k in prange(a.size):
        c[k] = np.cos(nu[k]*dt)
    for k in prange(a.size):
        zk = z[k]
        vzk = vz[k]
        z[k] = zk*c[k] + (vzk/nu[k])*s[k]
        vz[k] = vzk*c[k] - nu[k]*zk*s[k]


def sort_streamlines(r, t, vr, vt, z, vz):
    """helper_fns.sort_particles, carrying z and vz along: order each streamline's particles by
    longitude, by insertion sort, which is nearly free because particles barely move in longitude
    per step. Longitudes are distinct, so the order is the reference's argsort order."""
    Nr, Nt = t.shape
    for i in prange(Nr):
        for j in range(1, Nt):
            tj = t[i, j]
            if tj >= t[i, j - 1]:
                continue
            rj, vrj, vtj, zj, vzj = r[i, j], vr[i, j], vt[i, j], z[i, j], vz[i, j]
            k = j - 1
            while k >= 0 and t[i, k] > tj:
                t[i, k + 1] = t[i, k]
                r[i, k + 1] = r[i, k]
                vr[i, k + 1] = vr[i, k]
                vt[i, k + 1] = vt[i, k]
                z[i, k + 1] = z[i, k]
                vz[i, k + 1] = vz[i, k]
                k -= 1
            t[i, k + 1] = tj
            r[i, k + 1] = rj
            vr[i, k + 1] = vrj
            vt[i, k + 1] = vtj
            z[i, k + 1] = zj
            vz[i, k + 1] = vzj


def linear_density(r, vr, vt, lam, factor, J2_15, Rp2, Rp, two_GM_Rp):
    """helper_fns.get_lambda: lambda = factor*(Omega(a)/v), factor = m/(Nr*2*pi)."""
    for k in prange(r.size):
        h = r[k]*vt[k]
        cc = (h*h)/two_GM_Rp
        ak = Rp*(cc + np.sqrt(cc*cc - J2_15))
        Omg, Kap = _omega_kappa(ak, J2_15, Rp2)
        v = np.sqrt(vr[k]*vr[k] + vt[k]*vt[k])
        lam[k] = factor*(Omg/v)


def kick(vr, vt, Ar, At, dt):
    """vr += Ar*dt, vt += At*dt."""
    for k in prange(vr.size):
        vr[k] += Ar[k]*dt
        vt[k] += At[k]*dt


def wrap(c, cw, offset):
    """helper_fns.wrap_ring into cw, shape (Nr, Nt + 2); offset = 2 pi for longitudes, else 0."""
    Nr, Nt = c.shape
    for i in prange(Nr):
        cw[i, 0] = c[i, Nt - 1] - offset
        for j in range(Nt):
            cw[i, j + 1] = c[i, j]
        cw[i, Nt + 1] = c[i, 0] + offset


def crossed(rw, tw):
    """The wrapped crossing check of tools/patches.py: is any particle of streamline i inside
    streamline i - 1, interpolated to its longitude? rw, tw are wrapped by wrap()."""
    Nr, W = rw.shape
    n = 0
    for i in prange(1, Nr):
        s = i - 1
        for j in range(1, W - 1):
            x = tw[i, j]
            x0 = tw[s, j - 1]
            x1 = tw[s, j]
            x2 = tw[s, j + 1]
            dx0 = x - x0
            dx1 = x - x1
            dx2 = x - x2
            dx10 = x1 - x0
            dx20 = x2 - x0
            dx21 = x2 - x1
            l0 = (dx1*dx2)/(dx10*dx20)
            l1 = -(dx0*dx2)/(dx10*dx21)
            l2 = (dx0*dx1)/(dx20*dx21)
            if rw[i, j] - (rw[s, j - 1]*l0 + rw[s, j]*l1 + rw[s, j + 1]*l2) < 0:
                n += 1
    return n > 0


KERNELS = ['rt2xy', 'xy2rt', 'polar', 'shift', 'shift4', 'shift1', 'drift', 'vertical_drift', 'sort_streamlines', 'linear_density', 'kick', 'wrap', 'crossed']
SERIAL = {k: numba.njit(globals()[k]) for k in KERNELS}
PARALLEL = {k: numba.njit(globals()[k], parallel=True) for k in KERNELS}
