#!/usr/bin/env python3
# Sympy verification of the first-order vertical epicyclic equations V1-V6, and of the
# in-plane solution of hahn25 Appendix A (H1-H3) as transcribed in the paper's Appendix A
# (see notes/equations.md).
#
# Run:  uv run python notes/sympy/vertical_motion.py
#
# Each check prints PASS or raises AssertionError.
#
# Notation: Omega, kappa, nu here are the frequencies Omega_a, kappa_a, nu_a of the paper;
# the longitude of ascending node (Omega in the paper) is called node. M_p and R are the
# planet's mass and radius (M_p and R_p in the paper).
import sympy as sp

G, Mp, R, J2 = sp.symbols('G M_p R J_2', positive=True)
r, a, rho = sp.symbols('r a rho', positive=True)
z = sp.symbols('z', real=True)


def check(name, expr):
    """Assert that expr simplifies to zero."""
    assert sp.simplify(expr) == 0, f'{name}: residual {sp.simplify(expr)}'
    print(f'PASS  {name}')


# 3D potential of a J2-oblate planet, with rho^2 = r^2 + z^2 and P2 the Legendre polynomial
rho_ = sp.sqrt(r**2 + z**2)
P2 = lambda mu: (3*mu**2 - 1)/2
Phi = -(G*Mp/rho_)*(1 - J2*(R/rho_)**2*P2(z/rho_))

# --- the 3D potential reduces to hahn25 Eq. A2 in the midplane ---
Phi_A2 = -G*Mp/r - J2*G*Mp*R**2/(2*r**3)
check('Phi(r, z=0) equals hahn25 Eq. A2', Phi.subs(z, 0) - Phi_A2)

# --- V5: Phi is even in z, and Phi(r, z) = Phi(r, 0) + nu^2(r) z^2/2 + O(z^4) ---
nu2_r = (G*Mp/r**3)*(1 + sp.Rational(9, 2)*J2*(R/r)**2)
check('V5: Phi is even in z', Phi - Phi.subs(z, -z))
check('V5: Phi = Phi(r,0) + nu^2(r) z^2/2 + O(z^4)',
      sp.series(Phi, z, 0, 4).removeO() - (Phi_A2 + nu2_r*z**2/2))

# --- V1: frequencies, exact in J2 ---
n2 = G*Mp/a**3
Omega2 = n2*(1 + sp.Rational(3, 2)*J2*(R/a)**2)    # hahn25 Eq. A7, squared
kappa2 = n2*(1 - sp.Rational(3, 2)*J2*(R/a)**2)    # hahn25 Eq. A8, squared
nu2 = n2*(1 + sp.Rational(9, 2)*J2*(R/a)**2)       # V1

dPhi_dr = sp.diff(Phi, r).subs(z, 0)
check('Omega^2 = (1/r) dPhi/dr at z=0', (dPhi_dr/r).subs(r, a) - Omega2)
check('kappa^2 = r^-3 d(r^4 Omega^2)/dr at z=0',
      (sp.diff(r**3*dPhi_dr, r)/r**3).subs(r, a) - kappa2)
check('V1: nu^2 = d^2Phi/dz^2 at z=0', sp.diff(Phi, z, 2).subs(z, 0).subs(r, a) - nu2)
check('V1: kappa^2 + nu^2 = 2 Omega^2', kappa2 + nu2 - 2*Omega2)

# --- V2: elements -> coordinates; V3: coordinates -> elements ---
inc = sp.symbols('i', positive=True)
u = sp.symbols('u', real=True)                     # vertical phase u = lambda - node
nu = sp.sqrt(nu2)
z_V2 = a*inc*sp.sin(u)
vz_V2 = a*nu*inc*sp.cos(u)

# V3 recovers i, and (sin u, cos u), hence u = atan2(z, v_z/nu) exactly (mod 2 pi)
i_V3 = sp.sqrt((z_V2/a)**2 + (vz_V2/(a*nu))**2)
check('V3: i(V2) = i', i_V3 - inc)
check('V3: z/(a i) = sin u', z_V2/(a*i_V3) - sp.sin(u))
check('V3: (v_z/nu)/(a i) = cos u', (vz_V2/nu)/(a*i_V3) - sp.cos(u))

# --- V2 + V4: the motion is an exact harmonic oscillator at frequency nu ---
# lambda = wt + (Omega/kappa) M advances at Omega since dM/dt = kappa, and
# the node advances at Omega - nu (V4), so du/dt = Omega - (Omega - nu) = nu.
t, wt, M0, node0 = sp.symbols('t omegatilde M_0 node_0', real=True)
Om, ka = sp.sqrt(Omega2), sp.sqrt(kappa2)
M_t = M0 + ka*t
lam_t = wt + (Om/ka)*M_t
node_t = node0 + (Om - nu)*t                           # V4
u_t = lam_t - node_t
check('V4: du/dt = nu', sp.diff(u_t, t) - nu)
z_t, vz_t = z_V2.subs(u, u_t), vz_V2.subs(u, u_t)
check('V2: dz/dt = v_z', sp.diff(z_t, t) - vz_t)
check('V2: dv_z/dt = -nu^2 z', sp.diff(vz_t, t) + nu2*z_t)

# --- V2 solves the exact vertical equation of motion to first order ---
# z'' = -dPhi/dz evaluated at r = a(1 + d), where d = O(e) and z = O(i).
# The residual -dPhi/dz + nu^2 z has no terms of first order in (d, z).
eps, d, zh = sp.symbols('epsilon d zhat', real=True)
Az = -sp.diff(Phi, z)
resid = (Az + nu2*z).subs({r: a*(1 + eps*d), z: eps*a*zh})
ser = sp.series(resid, eps, 0, 3).removeO()
check('V2: residual of z equation has no O(eps^1) term', ser.coeff(eps, 1))
print('      leading residual, O(e*i):', sp.factor(sp.simplify(ser.coeff(eps, 2))))

# --- V6: the O(e i) residual equals the radial variation of nu^2(r), a d(nu^2)/da * delta * z ---
check('V6: residual = -a d(nu^2)/da * d * z',
      ser.coeff(eps, 2) - (-a*sp.diff(nu2_r, r).subs(r, a))*d*(a*zh))
check('V6: a d(nu^2)/da = -(3GM/a^3)[1 + (15/2) J2 (R/a)^2]',
      a*sp.diff(nu2_r, r).subs(r, a) + 3*G*Mp/a**3*(1 + sp.Rational(15, 2)*J2*(R/a)**2))

# --- V4: nodal regression rate, to first order in J2 ---
x = sp.symbols('x', positive=True)                 # x = J2 (R/a)^2
n = sp.sqrt(n2)
rate = n*(sp.sqrt(1 + sp.Rational(3, 2)*x) - sp.sqrt(1 + sp.Rational(9, 2)*x))
check('V4: Omega - nu = -(3/2) J2 (R/a)^2 n + O(J2^2)',
      sp.series(rate, x, 0, 2).removeO() + sp.Rational(3, 2)*x*n)

# --- H1-H3: the in-plane solution of hahn25 Appendix A, as transcribed in the paper's Appendix A ---
e, M = sp.symbols('e M', positive=True)
Om_a, ka_a = sp.sqrt(Omega2), sp.sqrt(kappa2)
r_A = a/sp.sqrt(1 + 2*e*sp.cos(M))
vr_A = e*a*ka_a*sp.sin(M)*(r_A/a)**3
vt_A = a**2*Om_a/r_A
check('H1: h = r v_theta = a^2 Omega_a exactly', r_A*vt_A - a**2*Om_a)
C = sp.symbols('C', positive=True)
a_A = R*(C + sp.sqrt(C**2 - sp.Rational(3, 2)*J2))
lhs = sp.expand((a**4*Omega2).subs(a, a_A))
check('H2: a = R(C + sqrt(C^2 - 3J2/2)) solves h^2 = a^4 Omega_a^2, with h^2 = 2 G M R C',
      lhs - 2*G*Mp*R*C)
check('H3: e sin M = (v_r/a kappa_a)(a/r)^3', (vr_A/(a*ka_a))*(a/r_A)**3 - e*sp.sin(M))
check('H3: e cos M = [(a/r)^2 - 1]/2', ((a/r_A)**2 - 1)/2 - e*sp.cos(M))

print('\nAll checks of V1-V6 and H1-H3 passed.')
