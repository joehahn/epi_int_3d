#!/usr/bin/env python3
# Sympy verification of the ring-pressure equations of Hahn & Spitale (2013, hahn13) Section 2.3.3
# as encoded in source/helper_fns.py (ring_pressure, A_P, delta_f, surface_density), with z = 0,
# plus the axisymmetric stability of the discrete streamline scheme (P1-P8, see
# notes/equations.md).
#
# Run:  uv run python notes/sympy/ring_pressure.py
#
# Each check prints PASS or raises AssertionError.
import mpmath
import sympy as sp

G, lam, sigma, c, Delta, kappa, k = sp.symbols('G lambda sigma c Delta kappa k', positive=True)
r, eps, phi = sp.symbols('r epsilon phi', real=True)
vr, vt = sp.symbols('v_r v_theta', real=True)
p = sp.Function('p')


def check(name, expr):
    """Assert that expr simplifies to zero."""
    assert sp.simplify(expr) == 0, f'{name}: residual {sp.simplify(expr)}'
    print(f'PASS  {name}')


def taylor(expr, n):
    """Series of expr in Delta about 0 through O(Delta^(n-1))."""
    return sp.series(expr, Delta, 0, n).removeO()


# --- P1: hahn13 Eq. 23, A_p = [p(r - Delta) - p(r)]/lambda = -(1/sigma) dp/dr + O(Delta), sigma = lambda/Delta ---
A23 = (p(r - Delta) - p(r))/(sigma*Delta)
check('P1: Eq. 23, [p(r-D) - p(r)]/lambda -> -(1/sigma) dp/dr',
      taylor(A23, 1) + sp.diff(p(r), r)/sigma)

# --- P2: hahn13 Eqs. 24-25 as coded: sigma_j = 2 lambda/(r_{j+1} - r_{j-1}) (surface_density, delta_f),
#     A_j = -[(p_{j+1} - p_{j-1})/(r_{j+1} - r_{j-1})]/sigma_j (A_P), = -(p_{j+1} - p_{j-1})/(2 lambda),
#     which is -(1/sigma) dp/dr + O(Delta^2) on a uniform grid ---
sig_j = 2*lam/(2*Delta)
A_code = -((p(r + Delta) - p(r - Delta))/(2*Delta))/sig_j
check('P2: coded interior A_j = -(p_{j+1} - p_{j-1})/(2 lambda)', A_code + (p(r + Delta) - p(r - Delta))/(2*lam))
check('P2: coded interior A_j = -(1/sigma) dp/dr + O(Delta^2)',
      taylor(A_code.subs(lam, sigma*Delta), 2) + sp.diff(p(r), r)/sigma)

# --- P3: momentum. With the coded edges A_0 = -p_0/lambda, A_{N-1} = p_{N-2}/lambda (hahn13's p(r - Delta)),
#     the net radial force sum_j lambda A_j of N streamlines is not zero: it is
#     (-p_0 + p_1 + p_{N-2} - p_{N-1})/2, which vanishes only if p_0 = p_1 and p_{N-2} = p_{N-1}.
#     So the scheme violates Newton's third law at O(Delta dp/dr) at each edge. ---
N = 7
ps = sp.symbols(f'p0:{N}')
F = -ps[0] + ps[N - 2] + sum(-(ps[j + 1] - ps[j - 1])/2 for j in range(1, N - 1))
check('P3: net radial pressure force of the coded scheme = (-p_0 + p_1 + p_{N-2} - p_{N-1})/2',
      F - (-ps[0] + ps[1] + ps[N - 2] - ps[N - 1])/2)

# --- P4: the coded components A_r = A v_theta/v, A_theta = -A v_r/v are perpendicular to the velocity,
#     i.e. pressure acts normal to the streamline (hahn13: A_p,r ~ A_p, A_p,theta ~ -A_p v_r/v_theta) ---
v = sp.sqrt(vr**2 + vt**2)
check('P4: (A_r, A_theta) . (v_r, v_theta) = 0', (vt/v)*vr + (-vr/v)*vt)

# --- P5: Toomre. A pressure-supported (p = c^2 sigma) gravitating sheet has
#     omega^2 = kappa^2 - 2 pi G sigma k + c^2 k^2, whose minimum over k is kappa^2 - (pi G sigma/c)^2,
#     so it is stable iff Q = c kappa/(pi G sigma) >= 1. The code sets c = Q pi G sigma/Omega
#     (initialize_streamline), i.e. kappa -> Omega, which is exact to O(J2) ---
w2 = kappa**2 - 2*sp.pi*G*sigma*k + c**2*k**2
kmin = sp.solve(sp.diff(w2, k), k)[0]
check('P5: k at minimum omega^2 = pi G sigma/c^2', kmin - sp.pi*G*sigma/c**2)
check('P5: minimum omega^2 = kappa^2 (1 - 1/Q^2), Q = c kappa/(pi G sigma)',
      w2.subs(k, kmin) - kappa**2*(1 - (sp.pi*G*sigma/(c*kappa))**2))

# --- P6: the coded pressure on an axisymmetric radial mode xi_j = eps cos(j phi), phi = k Delta, of a
#     uniform lattice of streamlines spaced Delta: linearizing A_0 gives A_0 = -omega_p^2 xi_0 with
#     omega_p^2 = (c/Delta)^2 sin^2(phi). It -> c^2 k^2 as phi -> 0 but VANISHES at phi = pi: central
#     differences do not see the grid-scale zigzag (odd-even decoupling), so pressure cannot resist it ---
xi = lambda j: eps*sp.cos(j*phi)
rr = lambda j: j*Delta + xi(j)
sig = lambda j: 2*lam/(rr(j + 1) - rr(j - 1))
pp = lambda j: c**2*sig(j)
A0 = -(pp(1) - pp(-1))/(2*lam)
A0_lin = sp.diff(A0, eps).subs(eps, 0)*eps
wp2 = (c/Delta)**2*sp.sin(phi)**2
check('P6: coded pressure on mode xi_j: A_0 = -(c/Delta)^2 sin^2(phi) xi_0', A0_lin + wp2*eps)
check('P6: omega_p^2 -> c^2 k^2 as Delta -> 0', sp.limit(wp2.subs(phi, k*Delta), Delta, 0) - c**2*k**2)
check('P6: omega_p^2 = 0 for the zigzag phi = pi', wp2.subs(phi, sp.pi))

# --- P7: the coded gravity (ring_gravity, 2 G lambda/Delta per streamline) on the same mode:
#     each pair of streamlines at +-j adds 4 G lambda (1 - cos j phi)/(j Delta)^2 xi_0 (destabilizing), so
#     X(phi) = (4 G sigma/Delta) sum_j (1 - cos j phi)/j^2 = (4 G sigma/Delta)(pi phi/2 - phi^2/4), 0 <= phi <= 2 pi.
#     Its maximum, at the zigzag phi = pi, is pi^2 G sigma/Delta: half the continuum 2 pi G sigma k there ---
j = sp.symbols('j', positive=True, integer=True)
Fj = 2*G*lam/(rr(j) - rr(0)) - 2*G*lam/(rr(0) - rr(-j))
Fj_lin = sp.diff(Fj, eps).subs(eps, 0)
check('P7: pair +-j adds 4 G lambda (1 - cos j phi)/(j Delta)^2 per unit xi_0',
      Fj_lin - 4*G*lam*(1 - sp.cos(j*phi))/(j*Delta)**2)
mpmath.mp.dps = 30
for ph in (mpmath.mpf('0.3'), mpmath.mpf(1), mpmath.mpf(2), mpmath.pi):
    s = mpmath.zeta(2) - mpmath.re(mpmath.polylog(2, mpmath.exp(1j*ph)))   # sum (1 - cos j phi)/j^2
    assert abs(s - (mpmath.pi*ph/2 - ph**2/4)) < 1e-25, ph
print('PASS  P7: sum_j (1 - cos j phi)/j^2 = pi phi/2 - phi^2/4 (numerically via Li_2, phi = 0.3, 1, 2, pi)')
X = (4*G*sigma/Delta)*(sp.pi*phi/2 - phi**2/4)
check('P7: X(pi) = pi^2 G sigma/Delta', X.subs(phi, sp.pi) - sp.pi**2*G*sigma/Delta)
check('P7: X -> 2 pi G sigma k as Delta -> 0', sp.limit(X.subs(phi, k*Delta), Delta, 0) - 2*sp.pi*G*sigma*k)

# --- P8: resolution parameter R = lambda_crit/(2 Delta), lambda_crit = 4 pi^2 G sigma/kappa^2. The cold lattice
#     (c = 0) is axisymmetrically stable iff max X = X(pi) < kappa^2, i.e. iff R < 2; for a ringlet of mass m_r,
#     width Da and N_r streamlines, sigma = m_r/(2 pi a Da), Delta = Da/N_r, R = pi G m_r N_r/(kappa^2 a Da^2) ---
m_r, Da, N_r, a = sp.symbols('m_r Delta_a N_r a', positive=True)
lam_crit = 4*sp.pi**2*G*sigma/kappa**2
R = lam_crit/(2*Delta)
check('P8: R = 2 pi^2 G sigma/(kappa^2 Delta)', R - 2*sp.pi**2*G*sigma/(kappa**2*Delta))
check('P8: X(pi) < kappa^2  <=>  R < 2', (X.subs(phi, sp.pi)/kappa**2) - R/2)
check('P8: R = pi G m_r N_r/(kappa^2 a Da^2)',
      R.subs({sigma: m_r/(2*sp.pi*a*Da), Delta: Da/N_r}) - sp.pi*G*m_r*N_r/(kappa**2*a*Da**2))

print('\nAll checks of P1-P8 passed.')
