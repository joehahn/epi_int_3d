#!/usr/bin/env python3
# Sympy (and numerical) verification of the linear theory of a bending wave launched at a satellite's inner
# vertical resonance (B1-B10; notes/equations.md, notes/phase2c/bending_wave_theory.md).
#
# Run:  uv run python notes/sympy/bending_wave.py
#
# Each check prints PASS or raises AssertionError.
#
# Notation (GT): the satellite, of mass M_s, is on a circular orbit of radius a_s, inclination i_s, node Omega_s and
# mean longitude lambda_s; Omega_a and nu_a are the angular and vertical frequencies at radius a. The ring is
# in the planet's equator, z = 0. The m-armed pattern rotates at Omega_p, and omega_m = m(Omega_a - Omega_p) is the
# forcing frequency a ring particle feels; x = (r - r_V)/r_V is the distance from the vertical resonance r_V.
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import fresnel


def check(name, expr):
    """Assert that expr simplifies to zero."""
    res = sp.simplify(expr)
    assert res == 0, f'{name}: residual {res}'
    print(f'PASS  {name}')


def close(name, x, y, tol):
    """Assert that x equals y to relative precision tol."""
    err = abs(x - y)/abs(y)
    assert err < tol, f'{name}: {x} vs {y}, relative error {err:.2e}'
    print(f'PASS  {name}  (relative error {err:.1e})')


def laplace(beta, m, s):
    """Laplace coefficient b_s^(m)(beta) = (2/pi) int_0^pi cos(m phi)/(1 + beta^2 - 2 beta cos phi)^s dphi."""
    f = lambda phi: np.cos(m*phi)/(1 + beta**2 - 2*beta*np.cos(phi))**s
    return (2/np.pi)*quad(f, 0, np.pi, limit=200, epsabs=0, epsrel=1e-13)[0]


r, th, z, a_s, i_s, u, node, GM_s = sp.symbols('r theta z a_s i_s u Omega_s GM_s', real=True)
m = sp.symbols('m', positive=True, integer=True)

# --- B1: an inclined circular satellite's vertical pull, to first order in i_s ---
# satellite at argument of latitude u = lambda_s - Omega_s, rotated by i_s about the line of nodes
Xs = a_s*(sp.cos(node)*sp.cos(u) - sp.sin(node)*sp.sin(u)*sp.cos(i_s))
Ys = a_s*(sp.sin(node)*sp.cos(u) + sp.cos(node)*sp.sin(u)*sp.cos(i_s))
Zs = a_s*sp.sin(u)*sp.sin(i_s)
d2 = (r*sp.cos(th) - Xs)**2 + (r*sp.sin(th) - Ys)**2 + (z - Zs)**2
lam_s = node + u
d2_claim = r**2 + a_s**2 - 2*r*a_s*sp.cos(th - lam_s) + (z - a_s*i_s*sp.sin(u))**2
check('B1: |r - r_s|^2 = r^2 + a_s^2 - 2 r a_s cos(theta - lambda_s) + (z - a_s i_s sin(lambda_s - Omega_s))^2 + O(i_s^2)',
      sp.series(d2 - d2_claim, i_s, 0, 2).removeO())
check('B1: z_s = a_s i_s sin(lambda_s - Omega_s) + O(i_s^3)', sp.series(Zs - a_s*i_s*sp.sin(u), i_s, 0, 3).removeO())
# F_z = -d/dz(-GM_s/|r - r_s|) = GM_s (z_s - z)/|r - r_s|^3; the z_s term is the forcing, the z term only a
# secular shift of nu (k = 0) and parametric (second-order) terms, which the Fourier satellite leaves out
F = -sp.diff(-GM_s/sp.sqrt(d2_claim), z)
check('B1: F_z = GM_s (z_s - z)/|r - r_s|^3', F - GM_s*(a_s*i_s*sp.sin(u) - z)/d2_claim**sp.Rational(3, 2))

# --- B2: the Laplace expansion 1/D0^3 = (1/2 a_s^3) sum_k b_{3/2}^(k) cos(k psi), psi = theta - lambda_s, b^(-k) = b^(k),
# turns sin(A) cos(k psi) into a sum of sin(A + k psi); the k = -m term has the resonant argument ---
A, psi, k = sp.symbols('A psi k', real=True)
check('B2: sin A [cos(k psi) + cos(-k psi)] = sin(A + k psi) + sin(A - k psi)',
      sp.expand_trig(sp.sin(A)*(sp.cos(k*psi) + sp.cos(-k*psi)) - sp.sin(A + k*psi) - sp.sin(A - k*psi)))
lam = sp.symbols('lambda_s', real=True)
arg = (lam - node) + (-m)*(th - lam)
check('B2: k = -m: (lambda_s - Omega_s) - m(theta - lambda_s) = -[m theta - (m+1) lambda_s + Omega_s]',
      arg + (m*th - (m + 1)*lam + node))
# so F_z contains -A_V sin(m theta - (m+1) lambda_s + Omega_s), A_V = GM_s i_s b_{3/2}^(m)(beta)/(2 a_s^2), beta = r/a_s

# B3 (numerical): the exact force of an inclined satellite (all orders in i_s), Fourier-analyzed in theta and
# lambda_s on a grid, has the B2 coefficient for the (m, -(m+1)) term, and nothing in its cosine
N, inc = 256, 1.0e-5
g = 2*np.pi*np.arange(N)/N
TH, LS = np.meshgrid(g, g, indexing='ij')
for mm in (2, 3, 5):
    beta = ((mm - 1)/(mm + 1))**(2/3)                  # at the Keplerian resonance (B4); a_s = 1, GM_s = 1
    xs, ys, zs = np.cos(LS), np.sin(LS)*np.cos(inc), np.sin(LS)*np.sin(inc)          # node 0, u = lambda_s
    Fz = zs/((beta*np.cos(TH) - xs)**2 + (beta*np.sin(TH) - ys)**2 + zs**2)**1.5      # particle at z = 0
    chi = mm*TH - (mm + 1)*LS
    AV = inc*laplace(beta, mm, 1.5)/2
    close(f'B3: m = {mm}: exact satellite F_z has -A_V sin(m theta - (m+1) lambda_s)', -2*np.mean(Fz*np.sin(chi)), AV, 1e-8)
    assert abs(2*np.mean(Fz*np.cos(chi))) < 1e-12*AV

# --- B4: pattern speed and resonance. lambda_s advances at Omega_{a_s}, the satellite's argument of latitude
# lambda_s - Omega_s at nu_{a_s} (V8), so m Omega_p = (m+1) Omega_{a_s} - (Omega_{a_s} - nu_{a_s}) = m Omega_{a_s} + nu_{a_s};
# the inner vertical resonance m(Omega_a - Omega_p) = nu_a; Keplerian: (m-1) Omega_V = (m+1) Omega_{a_s} ---
W, nus, Wa, nua = sp.symbols('Omega_as nu_as Omega_a nu_a', positive=True)
Wp = W + nus/m
check('B4: m Omega_p = (m+1) Omega_as - (Omega_as - nu_as)', m*Wp - ((m + 1)*W - (W - nus)))
GM, rr, aa = sp.symbols('GM r a', positive=True)
kep = lambda x: sp.sqrt(GM/x**3)
bV = ((m - 1)/(m + 1))**sp.Rational(2, 3)
rV = bV*aa
check('B4: Keplerian r_V = a_s [(m-1)/(m+1)]^(2/3)',
      sp.simplify((m*(kep(rV) - (kep(aa) + kep(aa)/m)) - kep(rV)).subs(m, 2)) + sp.simplify(
          (m*(kep(rV) - (kep(aa) + kep(aa)/m)) - kep(rV)).subs(m, 5)))

# --- B5: near r_V, nu^2 - omega_m^2 = D_V x + O(x^2), D_V = -r d(omega_m^2 - nu^2)/dr = 2 nu_V r_V (dnu/dr - m dOmega/dr);
# with J2 below; Keplerian D_V = 3(m-1) Omega_V^2 ---
Om = sp.Function('Omega')(rr)
nu = sp.Function('nu')(rr)
Wp_ = sp.symbols('Omega_p', positive=True)
DV = -rr*sp.diff((m*(Om - Wp_))**2 - nu**2, rr)
DV_res = DV.subs(Wp_, Om - nu/m)                       # at r_V, m(Omega - Omega_p) = nu
check('B5: D_V = 2 nu r (dnu/dr - m dOmega/dr) at r_V', DV_res - 2*nu*rr*(sp.diff(nu, rr) - m*sp.diff(Om, rr)))
J2, Rp = sp.symbols('J2 R_p', positive=True)
Om_J2 = sp.sqrt(GM/rr**3*(1 + sp.Rational(3, 2)*J2*(Rp/rr)**2))
nu_J2 = sp.sqrt(GM/rr**3*(1 + sp.Rational(9, 2)*J2*(Rp/rr)**2))
DV_J2 = 2*nu_J2*rr*(sp.diff(nu_J2, rr) - m*sp.diff(Om_J2, rr))
check('B5: Keplerian D_V = 3(m-1) Omega_V^2', DV_J2.subs(J2, 0) - 3*(m - 1)*GM/rr**3)
D_J2_series = sp.series(DV_J2, J2, 0, 2).removeO()
check('B5: with J2, D_V = (GM/r^3)[3(m-1) + (3/2)(8m - 15) J2 (R_p/r)^2] + O(J2^2)',
      D_J2_series - GM/rr**3*(3*(m - 1) + sp.Rational(3, 2)*(8*m - 15)*J2*(Rp/rr)**2))

# --- B6: test particles: z'' + nu^2 z = F cos(omega_m t) is forced at z = F cos(omega_m t)/(nu^2 - omega_m^2) ---
t, Fh, w, n_ = sp.symbols('t F omega nu', positive=True)
zf = Fh*sp.cos(w*t)/(n_**2 - w**2)
check('B6: forced vertical oscillation z = F/(nu^2 - omega_m^2) ~ F/(D_V x)', sp.diff(zf, t, 2) + n_**2*zf - Fh*sp.cos(w*t))

# --- B7: a sheet of surface density sigma, corrugated by Z = h cos(kx), pulls itself back with F_z = -2 pi G sigma |k| Z.
# The displaced sheet is the flat one plus the layer -sigma Z delta'(z), whose potential is -d/dz of psi, the potential of
# sigma h cos(kx) delta(z) ---
Gs, sg, kk, hh, xx, zz = sp.symbols('G sigma k h x z', positive=True)
psi_up = -(2*sp.pi*Gs*sg*hh/kk)*sp.cos(kk*xx)*sp.exp(-kk*zz)          # z > 0 (psi is even in z)
psi_dn = psi_up.subs(zz, -zz)                                        # z < 0
check('B7: psi solves Laplace off the sheet', sp.diff(psi_up, xx, 2) + sp.diff(psi_up, zz, 2))
check('B7: psi jumps in dpsi/dz by 4 pi G sigma h cos(kx) (Poisson)',
      (sp.diff(psi_up, zz) - sp.diff(psi_dn, zz)).subs(zz, 0) - 4*sp.pi*Gs*sg*hh*sp.cos(kk*xx))
Phi1_up, Phi1_dn = -sp.diff(psi_up, zz), -sp.diff(psi_dn, zz)
Fz_sheet = -(sp.diff(Phi1_up, zz) + sp.diff(Phi1_dn, zz)).subs(zz, 0)/2
check('B7: F_z = -2 pi G sigma k Z on the sheet', Fz_sheet + 2*sp.pi*Gs*sg*kk*hh*sp.cos(kk*xx))
# so (omega_m^2 - nu^2) = 2 pi G sigma |k| (the bending-wave dispersion relation); near r_V, -D_V x = 2 pi G sigma |k|:
# |k| = D_V |x|/(2 pi G sigma) for x < 0, wavelength 4 pi^2 G sigma/(D_V |x|); no waves for x > 0

# --- B8: the forced wave. A trailing wave (k > 0) has |k| h = -i dh/dr, so (nu^2 - omega_m^2) h + 2 pi G sigma |k| h = F
# becomes -i (2 pi G sigma/r_V) dh/dx + D_V x h = F, solved by
#   h = -i K F exp(-i L x^2/2) int_x^oo exp(i L s^2/2) ds,  K = r_V/(2 pi G sigma), L = D_V K ---
x, s = sp.symbols('x s', real=True)
K, L, DVs = sp.symbols('K L D_V', positive=True)
Fc = sp.symbols('F')                                   # complex forcing amplitude
h = -sp.I*K*Fc*sp.exp(-sp.I*L*x**2/2)*sp.Integral(sp.exp(sp.I*L*s**2/2), (s, x, sp.oo))
ode = -sp.I/K*sp.diff(h, x) + DVs*x*h - Fc
check('B8: h solves -i(2 pi G sigma/r_V) h\' + D_V x h = F', ode.subs(L, DVs*K).doit())


def h_num(xv, Lv, Kv=1.0, Fv=1.0):
    """B8's h(x) with int_x^oo exp(i L s^2/2) ds = sqrt(pi/L) [(1/2 - C(w)) + i (1/2 - S(w))], w = x sqrt(L/pi)."""
    S, C = fresnel(xv*np.sqrt(Lv/np.pi))
    I = np.sqrt(np.pi/Lv)*((0.5 - C) + 1j*(0.5 - S))
    return -1j*Kv*Fv*np.exp(-1j*Lv*xv**2/2)*I


Lv, Kv = 2.0e6, 4.0e5                                  # e.g. L = D_V r_V/(2 pi G sigma), D_V = L/K = 5
DVv = Lv/Kv
for xv in (3e-2, 1e-1):
    close(f'B8: x = +{xv}: h -> F/(D_V x), the test-particle response (B6)', (h_num(xv, Lv, Kv)*DVv*xv).real, 1.0, 1e-3)
for xv in (-3e-2, -1e-1):
    # inside r_V, h is the free wave plus the local forced response F/(D_V x), of relative size 1/(|x| sqrt(2 pi L))
    wave = lambda xv: h_num(xv, Lv, Kv) - 1/(DVv*xv)
    close(f'B8: x = {xv}: |h - F/(D_V x)| -> |F| K sqrt(2 pi/L) = |F| sqrt(r_V/(G sigma D_V))',
          abs(wave(xv)), Kv*np.sqrt(2*np.pi/Lv), 1e-5)
    dx = 1e-9
    kloc = np.angle(wave(xv + dx)/wave(xv - dx))/(2*dx)                             # d(phase)/dx = k r_V
    close(f'B8: x = {xv}: local wavenumber k r_V = -L x (B7)', kloc, -Lv*xv, 1e-3)
check('B8: K sqrt(2 pi/L) = sqrt(r_V/(G sigma D_V))',
      (K*sp.sqrt(2*sp.pi/L)).subs({K: rr/(2*sp.pi*Gs*sg), L: DVs*rr/(2*sp.pi*Gs*sg)}) - sp.sqrt(rr/(Gs*sg*DVs)))

# --- B9: torque on the ring. The forcing potential is -F_z z, so the torque per area is sigma Z dF_z/dtheta, whose
# longitude average is (m sigma/2) Im(h F*). Then T = pi m sigma r_V^2 int Im(h F*) dx = -pi m sigma r_V^2 K |F|^2 int Re G dx,
# G(x) = int_x^oo exp(i L (s^2 - x^2)/2) ds. Swapping s and x conjugates the integrand, so 2 Re int G dx is the
# integral over the whole (x, s) plane, |sqrt(2 pi/L)|^2: Re int G dx = pi/L, and T = -pi^2 m sigma r_V^2 |F|^2/D_V,
# independent of the ring's self-gravity ---
eps = 0.05                                             # a convergence factor exp(-eps (s^2 + x^2)); L = 1
n = 4001
gx = np.linspace(-30, 30, n)
X, S2 = np.meshgrid(gx, gx, indexing='ij')
f = np.exp((0.5j - eps)*S2**2 - (0.5j + eps)*X**2)*((S2 > X) + 0.5*(S2 == X))           # half weight on the diagonal
half = f.sum().real*(gx[1] - gx[0])**2
close('B9: Re of the integral over s > x is half that over the whole plane, pi/|i L/2 - eps|/2',
      half, 0.5*np.pi/abs(0.5j - eps), 2e-3)
Ttot = sp.pi*m*sg*rr**2*(-K*sp.Abs(Fc)**2)*sp.pi/L
check('B9: T = -pi^2 m sigma r_V^2 |F|^2/D_V', Ttot.subs(L, DVs*K) + sp.pi**2*m*sg*rr**2*sp.Abs(Fc)**2/DVs)

# --- B10: group velocity. Inside r_V the wave has omega - m Omega = -omega_m, so omega = m Omega - sqrt(nu^2 + 2 pi G sigma k),
# c_g = d omega/dk = -pi G sigma/omega_m ~ -pi G sigma/nu_V: inward, at a constant speed ---
om = sp.symbols('Omega_a')                             # the ring's angular velocity
wk = m*om - sp.sqrt(n_**2 + 2*sp.pi*Gs*sg*kk)
check('B10: c_g = -pi G sigma/omega_m', sp.diff(wk, kk) + sp.pi*Gs*sg/sp.sqrt(n_**2 + 2*sp.pi*Gs*sg*kk))

# --- B11: a test particle under a satellite whose mass grows as 1 - exp(-t/tau). A circular particle at
# theta = theta_0 + Omega_a t feels F_z = Re[F e^{i(omega_m t + phi)}] (1 - e^{-gamma t}), gamma = 1/tau, omega_m = m(Omega_a - Omega_p).
# Starting from z = v_z = 0 (complex form; z is its real part):
#   z = F e^{i phi}[e^{i omega t}/(nu^2 - omega^2) - e^{(i omega - gamma) t}/((i omega - gamma)^2 + nu^2)] + C e^{i nu t} + D e^{-i nu t}
# with C, D fixed by z(0) = z'(0) = 0. The free part C, D never decays; for gamma << |nu - omega| its size,
# relative to the forced part, is ~ gamma/|nu - omega| ---
gam, phi = sp.symbols('gamma phi', positive=True)
Fz0 = sp.symbols('F', positive=True)
zp = Fz0*sp.exp(sp.I*phi)*(sp.exp(sp.I*w*t)/(n_**2 - w**2)
                           - sp.exp((sp.I*w - gam)*t)/((sp.I*w - gam)**2 + n_**2))
C_, D_ = sp.symbols('C D')
z11 = zp + C_*sp.exp(sp.I*n_*t) + D_*sp.exp(-sp.I*n_*t)
sol = sp.solve([z11.subs(t, 0), sp.diff(z11, t).subs(t, 0)], [C_, D_], dict=True)[0]
z11 = z11.subs(sol)
check('B11: z solves z\'\' + nu^2 z = F e^{i(omega t + phi)}(1 - e^{-gamma t})',
      sp.diff(z11, t, 2) + n_**2*z11 - Fz0*sp.exp(sp.I*(w*t + phi))*(1 - sp.exp(-gam*t)))
check('B11: z(0) = 0', z11.subs(t, 0))
check('B11: z\'(0) = 0', sp.diff(z11, t).subs(t, 0))
free = sp.lambdify((n_, w, gam), (abs(sol[C_]) + abs(sol[D_]))*(n_**2 - w**2)/Fz0)   # free/forced, |e^{i phi}| = 1
for dw in (1e-2, 1e-3):
    ratio = free(1.0, 1.0 - dw, 1e-2*dw)
    close(f'B11: |nu - omega| = {dw}, gamma = 0.01|nu - omega|: free/forced ~ gamma/|nu - omega|', ratio, 1e-2, 2e-2)

# --- B12: the forced inclination of a test particle near r_V. Its forced z (B6) has amplitude A_V/|nu^2 - omega_m^2| ~ A_V/(D_V |x|),
# so i_f = |z|/a; in the Keplerian limit (D_V = 3(m-1) GM_p/r_V^3, B5; a = r_V, beta_V = r_V/a_s)
#   i_f = (M_s/M_p) i_s beta_V^2 b_{3/2}^(m)(beta_V)/[6(m-1)|x|]   (this is Shu 1984's |L_V/x|/r_V) ---
Ms, Mp, isat, bb, xx_, rV_, as_ = sp.symbols('M_s M_p i_s b x r_V a_s', positive=True)
A_V = Gs*Ms*isat*bb/(2*as_**2)
D_V = 3*(m - 1)*Gs*Mp/rV_**3
i_f = A_V/(D_V*xx_)/rV_
check('B12: i_f = A_V/(D_V |x| r_V) = (M_s/M_p) i_s beta_V^2 b/(6(m-1)|x|)',
      i_f - (Ms/Mp)*isat*(rV_/as_)**2*bb/(6*(m - 1)*xx_))
# --- B13: the restoring force of a bending wave in a finite lattice of streamlines. In an infinite lattice of wires
# spaced Delta (lambda = sigma Delta), the corrugation z_j = Z cos(j phi), phi = k Delta, pulls wire 0 back with
# X(phi) Z, X = (2 G sigma/Delta) sum_{n != 0} (1 - cos n phi)/n^2 = (4 G sigma/Delta)(pi phi/2 - phi^2/4)  (P7).
# A ring that ends at its edges lacks the wires beyond them. Their "1" terms, the pull of each wire on the displaced
# particle, sum to (2 G sigma/Delta) sum_{n >= n_e} 1/n^2 ~ 2 G sigma/((n_e - 1/2) Delta) = 2 G sigma/d, with d = the
# distance to the edge streamline plus Delta/2; their cos terms oscillate and sum to O(1/(n_e^2 phi)). So
#   nu^2 - omega_m^2 = -X(k Delta) + 2 G sigma (1/d_in + 1/d_out),  i.e.  D_V |x| = X(k Delta) - 2 G sigma (1/d_in + 1/d_out) ---
Nn = sp.symbols('N', positive=True)
tail = sp.polygamma(1, Nn)                             # sum_{n >= N} 1/n^2
lead = sp.series(tail - 1/(Nn - sp.Rational(1, 2)), Nn, sp.oo, 4).removeO()
check('B13: sum_{n >= N} 1/n^2 - 1/(N - 1/2) = O(N^-3)', sp.limit(lead*Nn**2, Nn, sp.oo))
for Nv in (5, 50):
    close(f'B13: sum_{{n >= {Nv}}} 1/n^2 vs 1/(N - 1/2)', float(sp.polygamma(1, Nv)), 1/(Nv - 0.5), 0.12/Nv**2)    # ~1/(12 N^2)
Nr = 270                                               # the finite lattice: wire i of Nr, G = sigma = Delta = 1
for i, per in ((135, 50), (135, 11), (60, 21), (230, 21)):
    ph = 2*np.pi/per
    nn = np.arange(-i, Nr - i)
    nn = nn[nn != 0]
    finite = 2*np.sum((1 - np.cos(nn*ph))/nn**2)
    X = 4*(np.pi*ph/2 - ph**2/4)
    edges = 2*(1/(i + 0.5) + 1/(Nr - 1 - i + 0.5))
    cos_tail = 2*(1/(i**2*ph) + 1/((Nr - 1 - i)**2*ph))  # the dropped cos terms' size
    close(f'B13: finite lattice, wire {i} of {Nr}, {per} wires per wavelength: X - 2 G sigma (1/d_in + 1/d_out)',
          finite, X - edges, 2*cos_tail/(X - edges))
# --- B14: the whole series. Summed over all integers m (b^(-m) = b^(m)), the z_s part of F_z is exactly
# -(GM_s i_s/2 a_s^2) sum_m b_{3/2}^(m)(beta) sin(m theta - (m+1) lambda_s + Omega_s), at first order in i_s, for any
# beta = r/a_s != 1 (the Laplace coefficient's integral holds for beta > 1 too). Then which m resonate: the argument's
# rate seen by a particle at r is m(Omega_a - Omega_as) - nu_as, resonant where its size is nu_a; with J2 = 0.01, R_p = 0.5
# (the tests), m >= 2 resonate only inside a_s, m <= -2 only outside, m = -1 (the secular term, independent of lambda_s)
# where the nodes precess together and m = 0, 1 nowhere, all three only at r = a_s itself ---
for beta in (0.7, 1.4):
    Ms = np.arange(-120, 121)
    th0, ls0, node0 = 0.7, 2.1, 0.4
    exact = inc*np.sin(ls0 - node0)/(1 + beta**2 - 2*beta*np.cos(th0 - ls0))**1.5     # GM_s = a_s = 1
    series = -inc/2*sum(laplace(beta, abs(mm), 1.5)*np.sin(mm*th0 - (mm + 1)*ls0 + node0) for mm in Ms)
    close(f'B14: beta = {beta}: z_s/Delta_s^3 = -(i_s/2) sum over all m of b^(m) sin(m theta - (m+1) lambda_s + Omega_s)',
          series, exact, 1e-9)
J2_, Rp_ = 0.01, 0.5
freq = lambda r, c: np.sqrt((1 + c*J2_*(Rp_/r)**2)/r**3)
Om, nu = (lambda r: freq(r, 1.5)), (lambda r: freq(r, 4.5))
a_s = 2.0
rs = np.concatenate([np.linspace(0.6, a_s - 1e-3, 20000), np.linspace(a_s + 1e-3, 6.0, 20000)])
for mm in range(-6, 7):
    g_ = np.abs(mm*(Om(rs) - Om(a_s)) - nu(a_s)) - nu(rs)
    flip = (np.sign(g_[:-1]) != np.sign(g_[1:])) & ((rs[:-1] - a_s)*(rs[1:] - a_s) > 0)     # not across r = a_s
    roots = rs[:-1][flip]
    expect = 'inside' if mm >= 2 else 'outside' if mm <= -2 else 'none'
    found = 'none' if not len(roots) else 'inside' if roots.max() < a_s else 'outside' if roots.min() > a_s else 'both'
    assert found == expect and len(roots) <= 1, (mm, roots)
# the paper's form: with the satellite's longitude theta_s = Omega_as t and node Omega_s = (Omega_as - nu_as) t, both zero at
# t = 0 (fast/satellite.py's defaults, which the tests use), the argument is m theta - (m Omega_as + nu_as) t
tt, ths, Wa_, nua_ = sp.symbols('t theta Omega_as nu_as', real=True)
check('B14: m theta - (m+1) theta_s + Omega_s = m theta - (m Omega_as + nu_as) t',
      (m*ths - (m + 1)*Wa_*tt + (Wa_ - nua_)*tt) - (m*ths - (m*Wa_ + nua_)*tt))
print('PASS  B14: m >= 2 resonate inside a_s, m <= -2 outside, m = -1, 0, 1 not away from a_s (J2 = 0.01)')
# --- B15: B8 in the paper's notation, Eqn. (vres_h): with K = r_V/(2 pi G sigma) = 1/(eps D_V), L = D_V K = 1/eps and
# F = i A_V = D_V L_V, h = -i K F exp(-i L x^2/2) int = -(i L_V/eps) exp(-i x^2/2 eps) int_x^oo exp(i s^2/2 eps) ds ---
ep, LVs = sp.symbols('epsilon L_V')
h_B8 = -sp.I*K*Fc*sp.exp(-sp.I*L*x**2/2)*sp.Integral(sp.exp(sp.I*L*s**2/2), (s, x, sp.oo))
h_paper = -sp.I*LVs/ep*sp.exp(-sp.I*x**2/(2*ep))*sp.Integral(sp.exp(sp.I*s**2/(2*ep)), (s, x, sp.oo))
check("B15: B8's h = -(i L_V/eps) e^{-i x^2/2 eps} int_x^oo e^{i s^2/2 eps} ds",
      sp.simplify(h_B8.subs({K: 1/(ep*DVs), L: 1/ep, Fc: DVs*LVs}) - h_paper))
check('B15: K = 1/(eps D_V) for eps = 2 pi G sigma/(r_V D_V)',
      (rr/(2*sp.pi*Gs*sg)) - 1/((2*sp.pi*Gs*sg/(rr*DVs))*DVs))
# --- B16: the massless limit (paper Appendix D). Integrating by parts with e^{i s^2/2 eps} = (eps/i s) d/ds e^{i s^2/2 eps},
# int_x^oo e^{i s^2/2 eps} ds = -(eps/i x) e^{i x^2/2 eps} [1 + O(eps/x^2)] for x > 0; for x < 0 add the whole Fresnel
# integral (2 pi eps)^{1/2} e^{i pi/4}. So h = L_V/x [1 + O(eps/x^2)] for x > 0, h -> h_forced as eps -> 0, and for x < 0
# h = L_V/x - i e^{i pi/4} (2 pi/eps)^{1/2} L_V e^{-i x^2/2 eps}, whose wave term has |Z_w| of B8 ---
sv = sp.symbols('s', real=True)
ep_ = sp.symbols('epsilon', positive=True)
E = sp.exp(sp.I*sv**2/(2*ep_))
check('B16: d/ds[(eps/i s) e^{i s^2/2 eps}] = e^{i s^2/2 eps} - (eps/i s^2) e^{i s^2/2 eps}',
      sp.diff(ep_/(sp.I*sv)*E, sv) - (E - ep_/(sp.I*sv**2)*E))
fres = sp.integrate(E, (sv, -sp.oo, sp.oo))
check('B16: int_-oo^oo e^{i s^2/2 eps} ds = (2 pi eps)^{1/2} e^{i pi/4}',
      sp.simplify(sp.expand_complex(fres - sp.sqrt(2*sp.pi*ep_)*sp.exp(sp.I*sp.pi/4))))
xs_ = sp.symbols('x', real=True, nonzero=True)
check('B16: -(i L_V/eps) e^{-i x^2/2 eps} * (-(eps/i x) e^{i x^2/2 eps}) = L_V/x',
      -sp.I*LVs/ep_*sp.exp(-sp.I*xs_**2/(2*ep_))*(-(ep_/(sp.I*xs_))*sp.exp(sp.I*xs_**2/(2*ep_))) - LVs/xs_)
check('B16: |(i L_V/eps)(2 pi eps)^{1/2}| = (2 pi/eps)^{1/2} |L_V|, the |Z_w| of B8',
      sp.Abs(sp.I*LVs/ep_*sp.sqrt(2*sp.pi*ep_)) - sp.sqrt(2*sp.pi/ep_)*sp.Abs(LVs))
sg_, G_, rV_, D_, AV_, nu_, m_ = sp.symbols('sigma G r_V D A_V nu m', positive=True)
Zw2 = AV_**2*rV_/(G_*sg_*D_)                                   # |Z_w|^2, paper Eqn. (vres_amp)
check('B16: T_V/c_g and sigma |Z_w|^2 are both independent of sigma',
      sp.diff((sp.pi**2*m_*sg_*rV_**2*AV_**2/D_)/(sp.pi*G_*sg_/nu_), sg_) + sp.diff(sg_*Zw2, sg_))
def h_eps(xv, ev):
    """Eqn. (vres_h) with L_V = 1, from B8's Fresnel form (K = 1/eps, L = 1/eps, F = 1)."""
    return h_num(xv, 1/ev, 1/ev, 1.0)
for ev in (1e-6, 1e-8):
    for xv in (0.01, -0.01):
        wave = -1j*np.exp(1j*np.pi/4)*np.sqrt(2*np.pi/ev)*np.exp(-1j*xv**2/(2*ev)) if xv < 0 else 0
        err = abs(h_eps(xv, ev) - 1/xv - wave)*abs(xv)
        assert err < 2*ev/xv**2, (ev, xv, err)
        print(f'PASS  B16: eps = {ev:.0e}, x = {xv:+.2f}: |h - L_V/x - h_wave|/|L_V/x| = {err:.1e} <= ~eps/x^2 = {ev/xv**2:.0e}')
# --- B17: h_forced = L_V/x (paper Eqn. 32b) is the massless response iA_V/D (Section 4.1.2) to first order in x,
# since L_V = iA_V/calD and D = calD x + O(x^2) near r_V ---
x17, A17, cD17, c2 = sp.symbols('x A_V calD c_2', real=True, nonzero=True)
D17 = cD17*x17 + c2*x17**2                                      # D expanded about r_V
ratio = sp.series((sp.I*A17/D17)/((sp.I*A17/cD17)/x17), x17, 0, 1).removeO()
check('B17: (iA_V/D)/(L_V/x) = 1 + O(x), with L_V = iA_V/calD', ratio - 1)
# --- B18: Shu's Eqn. 88 in the paper's convention, -i eps dh/dx + x h = L_V (paper Eqn. vres_ode), which
# Eqn. (vres_h) solves; its homogeneous solution e^{-i x^2/2 eps} is the free wave, absent at x > 0 since
# h -> L_V/x there (B16). The wave's phase x^2/2 eps reaches 2 pi at |x| = (4 pi eps)^{1/2}: the first wavelength ---
x18 = sp.symbols('x', real=True)
e18 = sp.symbols('epsilon', positive=True)
L18 = sp.symbols('L_V')
s18 = sp.symbols('s', real=True)
h18 = -sp.I*L18/e18*sp.exp(-sp.I*x18**2/(2*e18))*sp.Integral(sp.exp(sp.I*s18**2/(2*e18)), (s18, x18, sp.oo))
check('B18: Eqn. (vres_h) solves -i eps dh/dx + x h = L_V', sp.expand(-sp.I*e18*h18.diff(x18).doit() + x18*h18 - L18))
check('B18: e^{-i x^2/2 eps} solves the homogeneous equation',
      -sp.I*e18*sp.diff(sp.exp(-sp.I*x18**2/(2*e18)), x18) + x18*sp.exp(-sp.I*x18**2/(2*e18)))
check('B18: x^2/2 eps = 2 pi at |x| = (4 pi eps)^{1/2}', (sp.sqrt(4*sp.pi*e18))**2/(2*e18) - 2*sp.pi)
# --- B19: free bending waves (paper Section 4.1.1). Z = Re[h e^{i(m theta - omega t)}], h ~ e^{i int k dr},
# omega_m = m(Omega_a - Omega_p) = m Omega_a - omega. Dispersion relation omega_m^2 = nu^2 + 2 pi G sigma |k| gives
# c_g = d omega/dk = -sgn(k) pi G sigma/omega_m; crests m theta + k r = const have d theta/dr = -k/m (trailing iff
# k/m > 0); and with D = calD x, |k| = k > 0 and k -> -i d/dr, D h = -2 pi G sigma k h is -i eps h' + x h = 0 ---
k19 = sp.symbols('k', real=True, nonzero=True)
G19, sg19, nu19, Oa19, m19, r19 = sp.symbols('G sigma nu Omega_a m r', positive=True)
wm19 = sp.sqrt(nu19**2 + 2*sp.pi*G19*sg19*sp.Abs(k19))          # omega_m > 0 near a vertical resonance
om19 = m19*Oa19 - wm19
check('B19: c_g = d omega/dk = -sgn(k) pi G sigma/omega_m', sp.diff(om19, k19) + sp.sign(k19)*sp.pi*G19*sg19/wm19)
th19 = sp.symbols('theta0')
mm = sp.symbols('m_', nonzero=True, real=True)
check('B19: along a crest m theta + k r = const, d theta/dr = -k/m', sp.diff((th19 - k19*r19)/mm, r19) + k19/mm)
for sk, label in ((1, 'inner (m > 0): k > 0 gives c_g < 0'), (-1, 'outer (m < 0): k < 0 gives c_g > 0')):
    kk = sk*sp.Rational(3, 7)
    cg = (-sp.sign(kk)*sp.pi*G19*sg19/wm19.subs(k19, kk))
    assert sp.sign(cg) == -sk, label
    print(f'PASS  B19: {label}')
x19 = sp.symbols('x', real=True)
e19, cD19, rV19 = sp.symbols('epsilon calD r_V', positive=True)
hf = sp.Function('h')(x19)
# D h = -2 pi G sigma k h with k h = -i dh/dr = -(i/r_V) dh/dx, D = calD x, eps = 2 pi G sigma/(r_V calD)
lhs = cD19*x19*hf + 2*sp.pi*G19*sg19*(-sp.I/rV19)*hf.diff(x19)
check('B19: (D h + 2 pi G sigma k h)/calD = -i eps dh/dx + x h, the left side of Eqn. (vres_ode)',
      (lhs/cD19 - (-sp.I*e19*hf.diff(x19) + x19*hf)).subs(e19, 2*sp.pi*G19*sg19/(rV19*cD19)))
print('all bending-wave checks passed')
