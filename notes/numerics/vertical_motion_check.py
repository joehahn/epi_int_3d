#!/usr/bin/env python3
# Numerical check of the first-order, exactly reversible vertical epicyclic equations V1-V4
# (notes/equations.md), extending hahn25 Eqs. A9-A10 to 3D. Checks the coords <-> elements
# round trip, and the error after one orbit versus an exact integration in the 3D J2 potential.
# Units: GM = R = 1.
#
# Run:  uv run python notes/numerics/vertical_motion_check.py
import numpy as np
from scipy.integrate import solve_ivp

GM, Rp, J2 = 1.0, 1.0, 0.016

def freqs(a):
    n2 = GM/a**3
    Omega = np.sqrt(n2*(1 + 1.5*J2*(Rp/a)**2))
    kappa = np.sqrt(n2*(1 - 1.5*J2*(Rp/a)**2))
    nu = np.sqrt(n2*(1 + 4.5*J2*(Rp/a)**2))     # vertical frequency d^2Phi/dz^2 at z=0
    return Omega, kappa, nu

def elem2coords(a, e, M, wt, inc, node):
    Omega, kappa, nu = freqs(a)
    D = 1 + 2*e*np.cos(M)
    r = a/np.sqrt(D)                                  # A9
    th = (Omega/kappa)*(M + 2*e*np.sin(M)) + wt
    vr = e*a*kappa*np.sin(M)*(r/a)**3
    vt = a*a*Omega/r
    lam = wt + (Omega/kappa)*M                         # mean longitude
    u = lam - node                                     # vertical phase, advances at rate nu
    z = a*inc*np.sin(u)
    vz = a*inc*nu*np.cos(u)
    return r, th, z, vr, vt, vz

def coords2elem(r, th, z, vr, vt, vz):
    h = r*vt
    C = h*h/(2*GM*Rp)
    a = Rp*(C + np.sqrt(C*C - 1.5*J2))                 # A10
    Omega, kappa, nu = freqs(a)
    esM = (vr/(a*kappa))*(a/r)**3
    ecM = 0.5*((a/r)**2 - 1)
    e = np.hypot(esM, ecM)
    M = np.arctan2(esM, ecM)
    wt = th - (Omega/kappa)*(M + 2*esM)
    lam = wt + (Omega/kappa)*M
    inc = np.hypot(z/a, vz/(a*nu))
    u = np.arctan2(z, vz/nu)
    node = lam - u
    return a, e, M, wt, inc, node

def drift(coords, dt):
    a, e, M, wt, inc, node = coords2elem(*coords)
    Omega, kappa, nu = freqs(a)
    M = M + kappa*dt
    node = node + (Omega - nu)*dt                      # nodal regression
    return elem2coords(a, e, M, wt, inc, node)

# --- exact 3D equations of motion in the J2 potential, Cartesian ---
def rhs(t, s):
    x, y, z, vx, vy, vz = s
    rho2 = x*x + y*y + z*z
    rho = np.sqrt(rho2)
    f = GM/rho**3
    q = 1.5*J2*Rp*Rp/rho2
    zz = 5*z*z/rho2
    return [vx, vy, vz, -f*x*(1 + q*(1 - zz)), -f*y*(1 + q*(1 - zz)), -f*z*(1 + q*(3 - zz))]

def cyl2cart(r, th, z, vr, vt, vz):
    c, s = np.cos(th), np.sin(th)
    return [r*c, r*s, z, vr*c - vt*s, vr*s + vt*c, vz]

def cart2cyl(x, y, z, vx, vy, vz):
    r = np.hypot(x, y)
    th = np.arctan2(y, x)
    return r, th, z, (x*vx + y*vy)/r, (x*vy - y*vx)/r, vz

if __name__ == '__main__':
    rng = np.random.default_rng(1)

    # 1. round-trip coords -> elements -> coords, many random particles
    n = 100000
    a0 = rng.uniform(1.5, 3, n)
    el = (a0, rng.uniform(0, 0.05, n), rng.uniform(-np.pi, np.pi, n), rng.uniform(-np.pi, np.pi, n),
          rng.uniform(0, 0.05, n), rng.uniform(-np.pi, np.pi, n))
    c0 = elem2coords(*el)
    c1 = elem2coords(*coords2elem(*c0))
    print('round trip, max |relative error| per coordinate (r, theta, z, vr, vt, vz):')
    scale = [a0, 1, a0, a0**-0.5, a0**-0.5, a0**-0.5]
    print('   ', ['%.1e' % np.max(np.abs((x1 - x0)/s)) for x0, x1, s in zip(c0, c1, scale)])
    # repeated drift steps of dt=0 should not wander
    c = c0
    for _ in range(10000):
        c = elem2coords(*coords2elem(*c))
    print('after 10^4 round trips, max |dr/a| = %.1e, max |dz/a| = %.1e' %
          (np.max(np.abs(c[0] - c0[0])/a0), np.max(np.abs(c[2] - c0[2])/a0)))

    # 2. accuracy vs exact integration: error after one orbit should scale as eps^2 when e = i = eps
    print('\naccuracy after 1 orbit, a=2 (position error / a):')
    print('      e        i      |dx|/a    |dz|/a')
    a = 2.0
    T = 2*np.pi/freqs(a)[0]
    for e, inc in [(0.02, 0), (0.01, 0), (0, 0.02), (0, 0.01), (0.02, 0.02), (0.01, 0.01), (0.005, 0.005)]:
        c0 = elem2coords(a, e, 0.7, 0.3, inc, 1.1)
        sol = solve_ivp(rhs, (0, T), cyl2cart(*c0), rtol=1e-12, atol=1e-14, method='DOP853')
        exact = np.array(sol.y[:, -1])
        approx = np.array(cyl2cart(*drift(c0, T)))
        dx = np.linalg.norm(exact[:3] - approx[:3])/a
        dz = abs(exact[2] - approx[2])/a
        print('   %6.3f  %6.3f   %.2e  %.2e' % (e, inc, dx, dz))
