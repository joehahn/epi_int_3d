"""Numba kernels for epi_int_3d's faster code (Phase 0b, step 2).

ring_gravity(lambda_, G_ring, r, t, vr, vt) does what source/helper_fns.py's ring_gravity does with
fast_gravity = False, as one compiled loop instead of a Python loop over streamlines that builds
about ten temporary arrays per pass. The arrays are the ring wrapped in longitude by wrap_ring,
shape (Nr, Nt + 2). The acceleration on particle (i, j) due to streamline s = i - shft (mod Nr),
for shft = 1, ..., Nr - 1, is 2 G lambda_ij/dr along the unit normal to streamline s, where dr and
streamline s's velocity are interpolated to longitude t_ij by a 2nd-order Lagrange polynomial
through particles j - 1, j, j + 1 of streamline s.

Each particle's sum runs over shft in the reference's order, and each term is evaluated with the
reference's operations in the reference's order, without fastmath, so the result is bit-for-bit
the reference's. Only the interior columns 1, ..., Nt are computed: the wrapped edge columns 0 and
Nt + 1 are dropped by accelerations() and are left zero here.

The loop over shft sits outside the loop over j, so the innermost loop runs along contiguous
memory over independent particles, which LLVM can vectorize without reordering any sum. The
outer loop over target streamlines i is a prange: with parallel=True each thread handles whole
streamlines, which also leaves every sum's order unchanged.
"""
import numba
import numpy as np


def _ring_gravity(lambda_, G_ring, r, t, vr, vt):
    Nr, W = r.shape
    Ar = np.zeros((Nr, W))
    At = np.zeros((Nr, W))
    two_G = 2.0*G_ring
    for i in numba.prange(Nr):
        for shft in range(1, Nr):
            s = i - shft
            if s < 0:
                s += Nr
            for j in range(1, W - 1):
                x = t[i, j]
                x0 = t[s, j - 1]
                x1 = t[s, j]
                x2 = t[s, j + 1]
                #lagrange_poly_fit's weights
                dx0 = x - x0
                dx1 = x - x1
                dx2 = x - x2
                dx10 = x1 - x0
                dx20 = x2 - x0
                dx21 = x2 - x1
                l0 = (dx1*dx2)/(dx10*dx20)
                l1 = -(dx0*dx2)/(dx10*dx21)
                l2 = (dx0*dx1)/(dx20*dx21)
                ri = r[s, j - 1]*l0 + r[s, j]*l1 + r[s, j + 1]*l2
                vri = vr[s, j - 1]*l0 + vr[s, j]*l1 + vr[s, j + 1]*l2
                vti = vt[s, j - 1]*l0 + vt[s, j]*l1 + vt[s, j + 1]*l2
                dr = ri - r[i, j]
                A = (two_G*lambda_[i, j])/dr
                vi = np.sqrt(vri*vri + vti*vti)
                cos_phi = vti/vi
                sin_phi = vri/vi
                Ar[i, j] += A*cos_phi
                At[i, j] -= A*sin_phi
    return Ar, At


ring_gravity_serial = numba.njit(_ring_gravity)
ring_gravity_parallel = numba.njit(_ring_gravity, parallel=True)
