#!/usr/bin/env python3
# Sympy (and numerical) verification of the linear theory of a spiral density wave launched at a satellite's inner
# Lindblad resonance (L1-L9; notes/equations.md), the in-plane counterpart of notes/sympy/bending_wave.py (B1-B21).
#
# Run:  uv run python notes/sympy/density_wave.py
#
# Each check prints PASS or raises AssertionError.
#
# Notation (GT): the satellite, of mass M_s, is on a circular orbit of radius a_s in the ring plane, with mean longitude
# lambda_s = Omega_{a_s} t. Omega_a and kappa_a are the angular and epicyclic frequencies at radius a. The m-armed
# pattern rotates at Omega_p, omega_m = m(Omega_a - Omega_p) is the frequency a ring particle feels, and
# x = (r - r_L)/r_L is the distance from the inner Lindblad resonance r_L. A perturbation is Re[f e^{i chi}],
# chi = m(theta - Omega_p t) (the paper's convention, the complex conjugate of Shu's).

import numpy as np
import sympy as sp
from scipy.integrate import quad, dblquad
from scipy.special import fresnel


def check(name, expr):
    expr = sp.simplify(expr)
    assert expr == 0, f'{name}: residual {expr}'
    print('PASS', name)


def close(name, x, y, tol):
    err = abs(x - y)/abs(y)
    assert err < tol, f'{name}: {x} vs {y}, relative error {err:.3g}'
    print('PASS', name, f'(relative error {err:.2g})')


def laplace(beta, m, s):
    f = lambda phi: np.cos(m*phi)/(1.0 + beta**2 - 2.0*beta*np.cos(phi))**s
    return (2.0/np.pi)*quad(f, 0.0, np.pi, limit=200, epsabs=0.0, epsrel=1e-13)[0]


def dlaplace(beta, m, s):
    """d b_s^(m)/d beta, analytically: d/dbeta (1 + beta^2 - 2 beta cos phi)^-s = -2s(beta - cos phi)(...)^-(s+1)."""
    f = lambda phi: -2.0*s*(beta - np.cos(phi))*np.cos(m*phi)/(1.0 + beta**2 - 2.0*beta*np.cos(phi))**(s + 1)
    return (2.0/np.pi)*quad(f, 0.0, np.pi, limit=200, epsabs=0.0, epsrel=1e-13)[0]


# --- L1: the satellite's m-th term. 1/|r - r_s| = (1/a_s) sum_k (1/2) b_{1/2}^(k)(beta) cos k(theta - lambda_s), b^(-k) = b^(k),
# so for m >= 1 the direct potential's m-th term (k = +-m) is Phi_m = -(GM_s/a_s) b_{1/2}^(m)(beta) cos m(theta - lambda_s);
# the indirect term is m = 1 only. Numerically: the m-th cosine coefficient of -GM_s/|r - r_s| in theta - lambda_s ---
G, Ms, a_s = 1.0, 1.0, 1.0
for m in (2, 3, 4):
    for beta in (0.5, 0.63, 0.8):
        r = beta*a_s
        phi_exact = lambda psi: -G*Ms/np.sqrt(r*r + a_s*a_s - 2*r*a_s*np.cos(psi))
        coef = (1/np.pi)*quad(lambda psi: phi_exact(psi)*np.cos(m*psi), 0, 2*np.pi, limit=400, epsabs=0, epsrel=1e-13)[0]
        close(f'L1: m = {m}, beta = {beta}: cos coefficient = -(GM_s/a_s) b_1/2^(m)', coef, -G*Ms/a_s*laplace(beta, m, 0.5), 1e-10)
# the analytic beta-derivative used in fast/satellite.py, against a finite difference
for m in (2, 3):
    d_fd = (laplace(0.63 + 1e-5, m, 0.5) - laplace(0.63 - 1e-5, m, 0.5))/2e-5
    close(f'L1: m = {m}: analytic db/dbeta', dlaplace(0.63, m, 0.5), d_fd, 1e-8)

# --- L2: pattern speed and resonance. Phi_m rotates with the satellite, Omega_p = Omega_{a_s}. The inner Lindblad
# resonance is where omega_m = kappa_a; Keplerian (kappa = Omega): (m - 1) Omega_L = m Omega_{a_s}, r_L = a_s[(m-1)/m]^{2/3} ---
m_, r_, a_ = sp.symbols('m r a_s', positive=True)
Om = lambda x: x**sp.Rational(-3, 2)
rL = a_*((m_ - 1)/m_)**sp.Rational(2, 3)
check('L2: Keplerian m(Omega - Omega_p) = kappa at r_L', sp.simplify(m_*(Om(rL) - Om(a_)) - Om(rL)))

# --- L3: a test particle's forced epicycle. Linearizing r'' - r theta'^2 = -Phi_0' - dPhi_1/dr and
# d(r^2 theta')/dt = -dPhi_1/dtheta about a circular orbit (xi_theta the azimuthal displacement times r):
#   xi_r'' - 2 Omega xi_theta' + (kappa^2 - 4 Omega^2) xi_r = -dPhi_1/dr,   xi_theta'' + 2 Omega xi_r' = -(1/r) dPhi_1/dtheta.
# Along the orbit, chi = m(theta - Omega_p t) advances at omega_m, so d/dt = i omega_m on Re[f e^{i chi}]. The forced
# solution is xi_r = -Psi/(r D), D = kappa^2 - omega_m^2, Psi = r dPhi/dr + (2 Omega/(Omega - Omega_p)) Phi
# (Binney & Tremaine 2008, Eqn. 3.146a; the potential of Goldreich & Tremaine 1979) ---
W, K, w, rr, mm = sp.symbols('Omega kappa omega_m r m', positive=True)
Phi, dPhi = sp.symbols('Phi dPhi')                                 #Phi_1's amplitude and its r-derivative at r
xr, xt = sp.symbols('xi_r xi_theta')
eqs = [-w**2*xr - 2*W*sp.I*w*xt + (K**2 - 4*W**2)*xr + dPhi, -w**2*xt + 2*W*sp.I*w*xr + sp.I*mm*Phi/rr]
sol = sp.solve(eqs, [xr, xt], dict=True)[0]
Wp = W - w/mm                                                      #Omega_p, from omega_m = m(Omega - Omega_p)
Psi = rr*dPhi + 2*W/(W - Wp)*Phi
check('L3: xi_r = -Psi/(r D)', sol[xr] + Psi/(rr*(K**2 - w**2)))
# at r_L (omega_m = kappa) and Keplerian (kappa = Omega), 2 Omega/(Omega - Omega_p) = 2m
check('L3: 2 Omega/(Omega - Omega_p) = 2m Omega/kappa', 2*W/(W - Wp) - 2*mm*W/w)

# --- L4: near r_L, D = kappa^2 - omega_m^2 ~= Dcal x, Dcal = r dD/dr at r_L = 2 kappa r (dkappa/dr - m dOmega/dr);
# Keplerian Dcal = 3(m-1) Omega_L^2, the same as the vertical resonance's (B5) ---
D_kep = Om(r_)**2 - (m_*(Om(r_) - Om(a_)))**2
Dcal = sp.simplify((r_*sp.diff(D_kep, r_)).subs(a_, r_/((m_ - 1)/m_)**sp.Rational(2, 3)))
check('L4: Keplerian Dcal = 3(m-1) Omega_L^2', Dcal - 3*(m_ - 1)*Om(r_)**2)
# with J2 (GM = 1): numerically, r dD/dr against 2 kappa r (kappa' - m Omega')
J2, Rp = 0.01, 0.5
Omf = lambda r: np.sqrt((1/r**3)*(1 + 1.5*J2*(Rp/r)**2))
kaf = lambda r: np.sqrt((1/r**3)*(1 - 1.5*J2*(Rp/r)**2))
Op, mv = Omf(1.6), 2
h = 1e-6
Df = lambda r: kaf(r)**2 - (mv*(Omf(r) - Op))**2
d = lambda f: (f(rl + h) - f(rl - h))/(2*h)
from scipy.optimize import brentq
rl = brentq(lambda r: mv*(Omf(r) - Op) - kaf(r), 0.5, 1.5, xtol=1e-15)
close('L4: at r_L, r dD/dr = 2 kappa r (dkappa/dr - m dOmega/dr) (J2 = 0.01)', rl*d(Df), 2*kaf(rl)*rl*(d(kaf) - mv*d(Omf)), 1e-7)

# --- L5: a self-gravitating sheet. A radial displacement xi = Re[xi e^{ikr}] makes sigma_1 = -sigma dxi/dr = -i k sigma xi,
# whose potential (WKB) is -2 pi G sigma_1/|k|, so its radial pull -dPhi_s/dr = -2 pi G sigma |k| xi pushes the
# displacement on (self-gravity is destabilizing in the plane). Adding it to L3's radial equation (the azimuthal pull is
# smaller by 1/(kr)): D xi = 2 pi G sigma |k| xi - Psi/r. Free waves: D = kappa^2 - omega_m^2 = 2 pi G sigma |k|, so they
# exist only where D > 0, outside r_L (x > 0); |k| = Dcal x/(2 pi G sigma), lambda = 4 pi^2 G sigma/(Dcal x) ---
k, sg, Gs, xi = sp.symbols('k sigma G xi', positive=True)
s1 = -sp.I*k*sg*xi
phis = -2*sp.pi*Gs*s1/k
check('L5: radial self-pull -dPhi_s/dr = +2 pi G sigma |k| xi', -sp.I*k*phis - 2*sp.pi*Gs*sg*k*xi)

# --- L6: group velocity. omega = m Omega - omega_m, omega_m = (kappa^2 - 2 pi G sigma|k|)^{1/2}, so
# c_g = d omega/dk = sgn(k) pi G sigma/omega_m. The wave carries its angular momentum away from r_L, outward (x > 0), so c_g > 0
# and k > 0: crests m theta + int k dr = const have d theta/dr = -k/m < 0, trailing ---
kk = sp.symbols('k', real=True)
mO = sp.symbols('mOmega', positive=True)
omega = mO - sp.sqrt(K**2 - 2*sp.pi*Gs*sg*sp.Abs(kk))
cg = sp.diff(omega, kk)
wm = sp.sqrt(K**2 - 2*sp.pi*Gs*sg*sp.Abs(kk))
check('L6: c_g = sgn(k) pi G sigma/omega_m', sp.simplify((cg - sp.sign(kk)*sp.pi*Gs*sg/wm).subs(kk, sp.Symbol('kp', positive=True))))

# --- L7: the forced wave. For the trailing wave k > 0, |k| xi = -i dxi/dr, so with D = Dcal x, L5 is
#   i eps dxi/dx + x xi = L_L,  eps = 2 pi G sigma/(r_L Dcal),  L_L = -Psi/(r_L Dcal)
# (the bending wave's Eqn. with eps -> -eps: B18). Its solution with no wave at x < 0 is
#   xi = -(i L_L/eps) e^{i x^2/2 eps} int_{-inf}^x e^{-i s^2/2 eps} ds,
# which tends to L_L/x (the test particle's forced response, L3) for x << -eps^{1/2} and to
# L_L/x + Z_w e^{i x^2/2 eps} for x >> eps^{1/2}, |Z_w| = (2 pi/eps)^{1/2}|L_L| = (|Psi|/r_L)(r_L/(G sigma Dcal))^{1/2} ---
x, s, eps, LL = sp.symbols('x s epsilon L_L')
I0 = sp.Function('I0')(x)                                           #int_{-inf}^x e^{-i s^2/2 eps} ds, dI0/dx = e^{-i x^2/2 eps}
xi7 = -(sp.I*LL/eps)*sp.exp(sp.I*x**2/(2*eps))*I0
dxi7 = sp.diff(xi7, x).subs(sp.Derivative(I0, x), sp.exp(-sp.I*x**2/(2*eps)))
check('L7: i eps xi\' + x xi = L_L', sp.I*eps*dxi7 + x*xi7 - LL)
check('L7: e^{i x^2/2 eps} is the free wave, phase increasing outward', sp.I*eps*sp.diff(sp.exp(sp.I*x**2/(2*eps)), x)
      + x*sp.exp(sp.I*x**2/(2*eps)))


def xi_num(xv, epsv, Lv=1.0):
    """L7's xi with the Fresnel integrals: int_{-inf}^x e^{-i s^2/2 eps} ds = (pi eps)^{1/2}[(1 - i)/2 + C(u) - i S(u)],
    u = x/(pi eps)^{1/2}."""
    u = xv/np.sqrt(np.pi*epsv)
    S, C = fresnel(u)
    integral = np.sqrt(np.pi*epsv)*((1 - 1j)/2 + C - 1j*S)
    return -(1j*Lv/epsv)*np.exp(1j*xv**2/(2*epsv))*integral


epsv = 1e-6
for xv in (-0.02, -0.005):
    #to next order (integration by parts), xi = (L_L/x)(1 + i eps/x^2 + O(eps^2/x^4))
    close(f'L7: x = {xv}: xi -> L_L/x inside r_L', xi_num(xv, epsv), (1.0/xv)*(1 + 1j*epsv/xv**2), 5e-3)
for xv in (0.005, 0.02):
    close(f'L7: x = {xv}: |xi - L_L/x| -> (2 pi/eps)^{{1/2}}|L_L| outside r_L', abs(xi_num(xv, epsv) - (1.0/xv)*(1 + 1j*epsv/xv**2)),
          np.sqrt(2*np.pi/epsv), 1e-3)
# the numerical xi satisfies the ODE (finite difference)
xv, hx = 0.003, 1e-9
dxi = (xi_num(xv + hx, epsv) - xi_num(xv - hx, epsv))/(2*hx)
close('L7: Fresnel form satisfies i eps xi\' + x xi = L_L', 1j*epsv*dxi + xv*xi_num(xv, epsv), 1.0, 1e-5)
# its phase winds with dphase/dx = x/eps, i.e. k = Dcal x/(2 pi G sigma) > 0 (L5), trailing
xv = 0.01
dph = np.angle((xi_num(xv + 1e-7, epsv) - 1/(xv + 1e-7))/(xi_num(xv - 1e-7, epsv) - 1/(xv - 1e-7)))/2e-7
close('L7: the wave part\'s phase gradient is x/eps', dph, xv/epsv, 1e-3)

# --- L8: the torque on the ring. The satellite's torque on a ring element displaced by xi is -dPhi/dtheta at its displaced
# place, whose first-order part is -xi_r d2Phi/drdtheta - (xi_theta/r) d2Phi/dtheta2. With xi_theta from L3's azimuthal
# equation (self-gravity's azimuthal pull is smaller by 1/(kr)), its longitude average is -(m/2r) Psi Im(xi_r), for real
# Phi (the Fourier satellite's, with its phase in chi). So T = -pi m sigma r_L Psi int Im(xi_r) dx, and L7's xi gives
# int Im(xi_r) dx = -pi L_L = pi Psi/(r_L Dcal): T = -pi^2 m sigma Psi^2/Dcal (Goldreich & Tremaine 1979), < 0 at an
# inner Lindblad resonance and independent of the ring's self-gravity, like B9 ---
P1 = sp.symbols('Phi', real=True)
dP1 = sp.symbols('dPhi', real=True)
xr_c = sp.symbols('xi_r')
xt_c = (2*W*sp.I*w*xr_c + sp.I*mm*P1/rr)/w**2                       #L3's azimuthal equation solved for xi_theta
avg = lambda a, b: sp.re(a*sp.conjugate(b))/2                       #theta average of Re[a e^{i chi}] Re[b e^{i chi}]
xre, xim = sp.symbols('xi_re xi_im', real=True)
dens = -avg(xr_c, sp.I*mm*dP1) - avg(xt_c/rr, -mm**2*P1)
dens = dens.subs(xr_c, xre + sp.I*xim)
Psi1 = rr*dP1 + 2*W/(W - Wp)*P1
check('L8: torque per unit mass = -(m/2r) Psi Im(xi_r)', sp.expand_complex(dens) + mm/(2*rr)*Psi1*xim)
# int Im(xi) dx = -pi L_L for L7's xi (xi = L_L g, eps = 1): with a Gaussian window e^{-delta x^2}, delta -> 0
xs = np.linspace(-300.0, 300.0, 3000001)
g = xi_num(xs, 1.0)
for delta in (1e-3, 3e-4):
    close(f'L8: int Im(xi) e^(-delta x^2) dx = -pi L_L (delta = {delta})', np.trapezoid(g.imag*np.exp(-delta*xs**2), xs),
          -np.pi, 3e-2)
Gs_, sgs, Dcs, ms, rls, Pss = sp.symbols('G sigma Dcal m r_L Psi', positive=True)
T8 = -sp.pi*ms*sgs*rls*Pss*(sp.pi*Pss/(rls*Dcs))
check('L8: T = -pi^2 m sigma Psi^2/Dcal', T8 + sp.pi**2*ms*sgs*Pss**2/Dcs)

# --- L9: test particles. L3's forced displacement is an eccentricity e_f = |xi_r|/r_L = |Psi|/(r_L^2 |D|) ~= |Psi|/(r_L^2 Dcal |x|)
# = |L_L/x| (L7's x << -eps^{1/2} limit); Keplerian, with Phi = -(GM_s/a_s) b_{1/2}^(m) and Psi = -(GM_s/a_s)(beta db/dbeta + 2m b):
#   e_f = (M_s/M_p) beta_L |beta_L db/dbeta + 2m b_{1/2}^(m)(beta_L)|/[3(m-1)|x|] ---
Ms_, Mp_, bb, db, beta_, xx = sp.symbols('M_s M_p b db beta x', positive=True)
rL_, aS = sp.symbols('r_L a_s', positive=True)
Psi9 = -(Gs_*Ms_/aS)*(beta_*db + 2*ms*bb)
ef = sp.Abs(Psi9)/(rL_**2*3*(ms - 1)*(Gs_*Mp_/rL_**3)*xx)
check('L9: Keplerian e_f', (ef - (Ms_/Mp_)*beta_*(beta_*db + 2*ms*bb)/(3*(ms - 1)*xx)).subs(rL_, beta_*aS))
