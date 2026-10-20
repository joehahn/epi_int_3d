#viscous_stability.py
#
#Linear stability of a circular ringlet as the code discretizes it: the axisymmetric modes of the
#reference implementation's own ring forces. The Jacobian of helper_fns.accelerations() (gravity,
#viscosity: the code's finite differences delta_f, df_dr, surface_density, A_P) is taken numerically on
#a circular ringlet whose streamlines all share one set of particle longitudes, so the code's Lagrange
#interpolation is exact and only the radial discretization is tested. The planet's force and the
#inertial terms are added analytically:
#    dr/dt = v_r,  dv_r/dt = v_t^2/r - g_p(r) + A_r,  dv_t/dt = -v_r v_t/r + A_t,
#with g_p = (GM/r^2)[1 + (3/2)J2(Rp/r)^2]. The eigenvalues' real parts are the growth rates. Edge
#streamlines use the code's one-sided stencils.
#
#Run:  uv run python notes/numerics/viscous_stability.py   (about 6 min)
#
#Result (2026-10-19), fastest mode, 333 streamlines unless noted:
#  case                              growth/orbit  dphi/pi   simulation (eccentric, Nt = 40)
#  both viscosities                  3.3e-3        0.82      crossed, 145 orbits
#  both, gravity off                 0 (neutral)             crossed, 168 orbits
#  shear only, full stress (fixed)   0.21          0.52      crossed, 57 orbits
#  shear torque only (unfixed)       0.52          0.52      crossed, 23 orbits
#  bulk only (fixed)                 0 (neutral)             no crossing, 1e3 orbits
#  no viscosity                      0 (neutral)             crossed, 849 orbits
#  bulk doubled                      0 (neutral)             crossed, 186 orbits
#  both viscosities, Nr = 100        0 (neutral)             stable, 5e4 orbits
#So the code's axisymmetric modes explain the shear-only crossings (about 24 e-folds are needed, from
#raggedness ~1e-11 to ~1 spacing), but not the others: those need the streamlines' eccentricity or the
#interpolation in longitude, which this calculation leaves out.

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO/'source'))
sys.path.insert(0, str(REPO/'tools'))
import helper_fns
import patches

J2, Rp, GM = 0.01, 0.5, 1.0
MASS, WIDTH, TAPER = 3.0e-8, 0.0125, 0.2          #the nominal ringlet
NT = 8                                              #particles per streamline (exact interpolation)


def ringlet(Nr):
    """Circular streamlines at the nominal ringlet's tapered spacing, all at the same longitudes."""
    a = 1.0 + WIDTH*patches.tapered_spacing(Nr, TAPER)
    t = np.tile(np.linspace(-np.pi, np.pi, NT, endpoint=False) + np.pi/NT, (Nr, 1))
    return a, t


def ring_forces(r, vr, vt, t, nu_s, nu_b, G_ring, fixed):
    """The code's ring accelerations, per streamline (all particles of a streamline are identical)."""
    R, VR, VT = (np.repeat(x[:, None], NT, axis=1) for x in (r, vr, vt))
    lam = helper_fns.get_lambda(MASS, len(r), J2, Rp, R, t, VR, VT)
    acc = patches.accelerations if fixed else helper_fns.accelerations
    Ar, At = acc(lam, G_ring, nu_s, nu_b, 0.0, R, t, VR, VT, False, False, False)
    Ar = Ar if np.ndim(Ar) else np.zeros_like(R)
    At = At if np.ndim(At) else np.zeros_like(R)
    return Ar.mean(axis=1), At.mean(axis=1)


def growth_rates(Nr, nu_s, nu_b, G_ring, fixed=True, eps=1.0e-10):
    """Eigenvalues and eigenvectors of the linearized code, state (r, v_r, v_t) per streamline."""
    a, t = ringlet(Nr)
    r0, vr0 = a.copy(), np.zeros(Nr)
    vt0 = a*helper_fns.Omega(J2, Rp, a)
    J = np.zeros((3*Nr, 3*Nr))
    #orbital part, exact
    gp_prime = -2.0*GM/r0**3 - 6.0*J2*Rp**2*GM/r0**5
    i = np.arange(Nr)
    J[i, Nr + i] = 1.0                                          #dr/dt = v_r
    J[Nr + i, i] = -vt0**2/r0**2 - gp_prime                     #d(v_t^2/r - g_p)/dr
    J[Nr + i, 2*Nr + i] = 2.0*vt0/r0                            #d(v_t^2/r)/dv_t
    J[2*Nr + i, Nr + i] = -vt0/r0                               #d(-v_r v_t/r)/dv_r
    #ring forces, central differences
    for k in range(3*Nr):
        x = [r0.copy(), vr0.copy(), vt0.copy()]
        x[k//Nr][k % Nr] += eps
        Ap, Tp = ring_forces(*x, t, nu_s, nu_b, G_ring, fixed)
        x[k//Nr][k % Nr] -= 2*eps
        Am, Tm = ring_forces(*x, t, nu_s, nu_b, G_ring, fixed)
        J[Nr:2*Nr, k] += (Ap - Am)/(2*eps)
        J[2*Nr:, k] += (Tp - Tm)/(2*eps)
    w, V = np.linalg.eig(J)
    return a, w, V


def describe(a, w, V):
    """Fastest-growing mode: growth per orbit, where it sits, and how zigzag it is."""
    Nr = len(a)
    k = np.argmax(w.real)
    xi = V[:Nr, k]
    amp = np.abs(xi)
    x = (a - a.mean())/(a[-1] - a[0])
    centre = (x*amp**2).sum()/(amp**2).sum()
    #phase advance between adjacent streamlines, weighted by amplitude: pi for a zigzag
    dphi = np.abs(np.angle(xi[1:]*np.conj(xi[:-1])))
    theta = (dphi*amp[1:]*amp[:-1]).sum()/(amp[1:]*amp[:-1]).sum()
    period = 2*np.pi/helper_fns.Omega(J2, Rp, a.mean())
    return w[k].real*period, centre, theta/np.pi


CASES = [  #(label, Nr, nu_s, nu_b, G_ring, fixed, simulated result at Nt = 40)
    ('both viscosities', 333, 7e-11, 7e-11, 1.0, False, 'crossed, 145 orbits'),
    ('both, gravity off', 333, 7e-11, 7e-11, -1.0, False, 'crossed, 168 orbits'),
    ('shear only, full stress (fixed)', 333, 7e-11, 0.0, 1.0, True, 'crossed, 57 orbits'),
    ('shear torque only (unfixed code)', 333, 7e-11, 0.0, 1.0, False, 'crossed, 23 orbits'),
    ('bulk only (fixed)', 333, 0.0, 7e-11, 1.0, True, 'no crossing, 1e3 orbits'),
    ('no viscosity', 333, 0.0, 0.0, 1.0, False, 'crossed, 849 orbits'),
    ('bulk doubled', 333, 7e-11, 1.4e-10, 1.0, False, 'crossed, 186 orbits'),
    ('both viscosities, Nr = 100', 100, 7e-11, 7e-11, 1.0, False, 'stable, 5e4 orbits'),
]

if __name__ == '__main__':
    print(f'{"case":36s} {"growth/orbit":>12s} {"e-folds in 145 orbits":>22s} {"centre":>7s} {"dphi/pi":>8s}   simulation')
    for label, Nr, nu_s, nu_b, G_ring, fixed, sim in CASES:
        a, w, V = growth_rates(Nr, nu_s, nu_b, G_ring, fixed)
        g, centre, th = describe(a, w, V)
        print(f'{label:36s} {g:12.2e} {g*145:22.1f} {centre:7.2f} {th:8.2f}   {sim}')
