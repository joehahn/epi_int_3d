#!/usr/bin/env python3
# Sympy verification of the 3D gravity of a perturbing streamline treated as a straight, infinitely
# long wire (G1-G9; notes/equations.md, notes/phase2/streamline_gravity.md). This is the 3D
# generalization of hahn13 Fig. 3 and Eqs. 19-20.
#
# Run:  uv run python notes/sympy/streamline_gravity.py
#
# Each check prints PASS or raises AssertionError.
#
# Geometry (local frame at the perturbed particle, unit vectors r, theta, z): the perturbing
# streamline crosses the particle's meridional plane (its longitude) at separation
# d = (dr, 0, dz) from the particle, where dr = r' - r and dz = z' - z. Near that point the
# streamline is a straight wire of linear density lam, along the unit tangent
# T = (eps, 1, eta)/sqrt(1 + eps^2 + eta^2), with slopes eps = v_r'/v_theta', eta = v_z'/v_theta'
# (hahn13 takes the tangent from the streamline's velocity; see G6).
import sympy as sp

G, lam, sigma, Delta = sp.symbols('G lambda sigma Delta', positive=True)
dr, dz, eps, eta, u, s = sp.symbols('Delta_r Delta_z epsilon eta u s', real=True)
b = sp.symbols('b', positive=True)


def check(name, expr):
    """Assert that expr simplifies to zero (each component, for a Matrix)."""
    expr = sp.Matrix(expr) if isinstance(expr, (list, tuple, sp.MatrixBase)) else expr
    res = sp.simplify(expr)
    zero = res.is_zero_matrix if isinstance(res, sp.MatrixBase) else res == 0
    assert zero, f'{name}: residual {res}'
    print(f'PASS  {name}')


d = sp.Matrix([dr, 0, dz])
norm = sp.sqrt(1 + eps**2 + eta**2)
T = sp.Matrix([eps, 1, eta])/norm
c = (d.T*T)[0]                       # distance along the wire from its crossing point to the foot
D = d - c*T                          # perpendicular from the particle to the wire
D2 = sp.simplify((D.T*D)[0])

# --- G1: an infinite straight wire pulls with 2 G lam D/|D|^2, D the perpendicular to it ---
check('G1: D is perpendicular to the wire', (D.T*T)[0])
check('G1: |D|^2 = |d|^2 - (d.T)^2', D2 - (dr**2 + dz**2 - c**2))
# Points on the wire are d + s T = D + u T with u = s + c, so the wire's pull on the particle is
# G lam * integral (D + u T)/(|D|^2 + u^2)^(3/2) du; b stands for |D|.
I_perp = sp.integrate(1/(b**2 + u**2)**sp.Rational(3, 2), (u, -sp.oo, sp.oo))
I_par = sp.integrate(u/(b**2 + u**2)**sp.Rational(3, 2), (u, -sp.oo, sp.oo))
check('G1: integral of the point-mass pull along the wire, perpendicular part = 2/|D|^2', I_perp - 2/b**2)
check('G1: integral of the point-mass pull along the wire, parallel part = 0', I_par)
A_exact = 2*G*lam*D/D2               # components (A_r, A_theta, A_z)

# --- G2: to first order in the slopes eps, eta, at any dz/dr ---
rho2 = dr**2 + dz**2
A_G2 = sp.Matrix([2*G*lam*dr/rho2, -2*G*lam*(eps*dr + eta*dz)/rho2, 2*G*lam*dz/rho2])
small = sp.symbols('t', positive=True)   # bookkeeping: eps -> t eps, eta -> t eta
A_t = A_exact.subs({eps: small*eps, eta: small*eta}, simultaneous=True)
A_G2t = A_G2.subs({eps: small*eps, eta: small*eta}, simultaneous=True)
check('G2: exact = first-order form + O(eps^2, eta^2, eps eta)',
      [sp.series(A_t[k] - A_G2t[k], small, 0, 2).removeO() for k in range(3)])
# the second-order terms, for the record: A_r and A_z change at O(slope^2), A_theta at O(slope^3)
err2 = [sp.factor(sp.series(A_t[k], small, 0, 3).removeO().coeff(small, 2)) for k in range(3)]
check('G2: A_theta has no second-order term', err2[1])
print('      second-order parts of A_r, A_z:', err2[0], ',', err2[2])

# --- G3: the planar limit dz = eta = 0 is hahn13 Eq. 20 ---
A_2D = A_exact.subs({dz: 0, eta: 0})
check('G3: planar wire, exact: A = (2 G lam/dr)(1, -eps, 0)',
      A_2D - sp.Matrix([2*G*lam/dr, -2*G*lam*eps/dr, 0]))
# hahn13 (and epi_int_lite, fast/kernels.py): A = (2 G lam/dr)(cos phi, -sin phi), with
# sin phi = v_r'/v' = eps/sqrt(1 + eps^2). That is the exact wire times 1/sqrt(1 + eps^2), i.e. it
# agrees to first order, as hahn13 states.
A_hs13 = sp.Matrix([2*G*lam/dr/sp.sqrt(1 + eps**2), -2*G*lam*eps/dr/sp.sqrt(1 + eps**2), 0])
check('G3: hahn13 Eq. 20 = exact planar wire/sqrt(1 + eps^2)', A_hs13 - A_2D/sp.sqrt(1 + eps**2))
check('G3: hahn13 Eq. 20 = exact planar wire + O(eps^2)',
      [sp.series((A_hs13 - A_2D)[k].subs(eps, small*eps), small, 0, 2).removeO() for k in range(3)])

# --- G4: continuum limit, linear: the bending-wave restoring force -2 pi G sigma |k| Z ---
# A sheet of surface density sigma (streamlines lam = sigma dr') with height Z(r) = Z0 cos(k r),
# |dz| << |dr|, so A_z = sum 2 G lam dz/dr^2 -> 2 G sigma Z0 integral (cos(k x) - 1)/x^2 dx at r = 0
# (the cos(kr) sin(kx)/x^2 part is odd in x and drops out). Same for both signs of k.
x = sp.symbols('x', real=True)
k = sp.symbols('k', positive=True)
I_bend = sp.integrate((sp.cos(k*x) - 1)/x**2, (x, -sp.oo, sp.oo))
check('G4: integral (cos kx - 1)/x^2 dx = -pi |k|', I_bend + sp.pi*k)
print('      so A_z = -2 pi G sigma |k| Z, the restoring term of the bending-wave dispersion relation')

# --- G5: the same on the code's lattice: A_z = -X(phi) Z, X of P7 ---
# Streamlines spaced Delta, Z_j = Z0 cos(j phi): A_z on streamline 0 = sum_{n != 0} 2 G lam Z0
# (cos(n phi) - 1)/(n Delta)^2 = (4 G sigma Z0/Delta) sum_{n >= 1} (cos(n phi) - 1)/n^2.
# The Fourier series sum_{n >= 1} cos(n phi)/n^2 = pi^2/6 - pi phi/2 + phi^2/4 on [0, 2 pi]:
phi = sp.symbols('phi', real=True)
n = sp.symbols('n', positive=True, integer=True)
F = sp.pi**2/6 - sp.pi*phi/2 + phi**2/4
check('G5: mean of F over [0, 2 pi] is 0', sp.integrate(F, (phi, 0, 2*sp.pi)))
check('G5: cosine coefficients of F are 1/n^2', sp.integrate(F*sp.cos(n*phi), (phi, 0, 2*sp.pi))/sp.pi - 1/n**2)
check('G5: sine coefficients of F are 0 (F is even about pi)', sp.integrate(F*sp.sin(n*phi), (phi, 0, 2*sp.pi)))
A_z_lattice = 4*G*sigma/Delta*(F - sp.pi**2/6)
X = (4*G*sigma/Delta)*(sp.pi*phi/2 - phi**2/4)       # P7
check('G5: A_z/Z0 = -X(phi), the X of P7', A_z_lattice + X)
check('G5: X -> 2 pi G sigma k as phi = k Delta -> 0',
      sp.limit(X.subs(phi, k*Delta)/(2*sp.pi*G*sigma*k), Delta, 0) - 1)
# the same 4-term check numerically, partial sum to n = 20000
import math
for ph in (0.3, 1.7, math.pi):
    S = sum((math.cos(m*ph) - 1)/m**2 for m in range(1, 20001))
    assert abs(S - float((F - sp.pi**2/6).subs(phi, ph))) < 1e-4, ph
print('PASS  G5: partial sums of sum (cos n phi - 1)/n^2 match')

# --- G6: velocity slope versus shape slope, for a pattern rotating at Omega_p ---
# A streamline whose shape r'(theta, t) = a + f(theta - Omega_p t), z'(theta, t) = g(theta - Omega_p t)
# rotates rigidly at pattern speed Omega_p, while its particles move at theta = Omega t (to lowest
# order). Then v_r'/v_theta' = (1 - Omega_p/Omega) (1/r) dr'/dtheta, and likewise for z.
t_, Om, Omp, a = sp.symbols('t Omega Omega_p a', positive=True)
f, g = sp.Function('f'), sp.Function('g')
th = Om*t_
rr = a + f(th - Omp*t_)
zz = g(th - Omp*t_)
slope_r = sp.diff(f(s), s).subs(s, th - Omp*t_)/rr            # (1/r) dr'/dtheta at fixed t
slope_z = sp.diff(g(s), s).subs(s, th - Omp*t_)/rr
check('G6: v_r/v_theta = (1 - Omega_p/Omega) (1/r) dr/dtheta', sp.diff(rr, t_)/(rr*Om) - (1 - Omp/Om)*slope_r)
check('G6: v_z/v_theta = (1 - Omega_p/Omega) (1/r) dz/dtheta', sp.diff(zz, t_)/(rr*Om) - (1 - Omp/Om)*slope_z)

# --- G7: hahn13's form plus two small modifications, to the same accuracy as G2 ---
# A wire along theta-hat through (r', z') pulls with A_g = 2 G lam/delta along the meridional
# direction (cos psi, 0, sin psi), delta^2 = dr^2 + dz^2, tan psi = dz/dr. Tilting the wire by phi in
# the r-theta plane (hahn13: sin phi = v_r'/v') rotates the radial part toward -theta-hat by phi;
# tilting it by chi in the theta-z plane (sin chi = v_z'/sqrt(v_theta'^2 + v_z'^2)) rotates the
# vertical part toward -theta-hat by chi:
#   A_r = A_g cos psi cos phi,  A_theta = -A_g (cos psi sin phi + sin psi sin chi),  A_z = A_g sin psi.
delta = sp.sqrt(dr**2 + dz**2)
cphi, sphi = 1/sp.sqrt(1 + eps**2), eps/sp.sqrt(1 + eps**2)
schi = eta/sp.sqrt(1 + eta**2)
Ag = 2*G*lam/delta
A_G7 = sp.Matrix([Ag*(dr/delta)*cphi, -Ag*((dr/delta)*sphi + (dz/delta)*schi), Ag*dz/delta])
A_G7t = A_G7.subs({eps: small*eps, eta: small*eta}, simultaneous=True)
check('G7: hahn13 form + tilts = exact wire + O(eps^2, eta^2, eps eta)',
      [sp.series(A_t[k] - A_G7t[k], small, 0, 2).removeO() for k in range(3)])
check('G7: at dz = 0 it is hahn13 Eq. 20, cos phi and sin phi as coded', A_G7.subs(dz, 0) - A_hs13)
# the coded form: q = dz/dr, A = 2 G lam/(dr (1 + q^2)), A_r = A cos phi, A_theta = -A (sin phi + q sin chi),
# A_z = A q (valid for either sign of dr)
q = dz/dr
Acode = 2*G*lam/(dr*(1 + q**2))
check('G7: = A (cos phi, -(sin phi + q sin chi), q), A = 2 G lam/(dr (1 + q^2))',
      A_G7 - sp.Matrix([Acode*cphi, -Acode*(sphi + q*schi), Acode*q]))

# --- G8: the exact wire in vector form: D = Delta - (Delta.v)v = v x (Delta x v), with v the wire's
# unit direction (the paper's Eqns. for D and A_g); here v = T and Delta = d ---
cross = T.cross(d.cross(T))
check('G8: v x (Delta x v) = Delta - (Delta.v) v = D', cross - D)
check('G8: A_g = 2 G lam (Delta - (Delta.v)v)/(Delta^2 - (Delta.v)^2) = G1',
      2*G*lam*(d - c*T)/(dr**2 + dz**2 - c**2) - A_exact)
check('G8: in the ring plane, G7 (as coded) = exact x cos(phi)', (A_G7 - A_exact*cphi).subs({dz: 0, eta: 0}))

# --- G9: hahn13's notation in 3D. Delta = D is the perpendicular from the particle to the wire
# (hahn13 Fig. 3), A_g = 2 G lam/|Delta| along it. With the wire's unit direction v = (v_r, v_theta,
# v_z)/v', to first order its components are Delta_r = r' - r, Delta_z = z' - z and
# Delta_theta = -(Delta_r v_r' + Delta_z v_z')/v', where (r', z') is where the wire crosses the
# particle's longitude, and |Delta| = sqrt(Delta_r^2 + Delta_z^2). ---
vr_, vt_, vz_ = sp.symbols("v_r' v_theta' v_z'", real=True)
vhat = sp.Matrix([vr_, vt_, vz_])/sp.sqrt(vr_**2 + vt_**2 + vz_**2)
Dv = d - (d.T*vhat)[0]*vhat
Dv_t = Dv.subs({vr_: small*vr_, vz_: small*vz_}, simultaneous=True)
first = sp.Matrix([dr, -(dr*vr_ + dz*vz_)/vt_, dz])          # v' = v_theta' to first order
check('G9: Delta = (Delta_r, -(Delta_r v_r + Delta_z v_z)/v, Delta_z) + O(slope^2)',
      [sp.series(Dv_t[k] - first.subs({vr_: small*vr_, vz_: small*vz_}, simultaneous=True)[k], small, 0, 2).removeO()
       for k in range(3)])
check('G9: |Delta|^2 = Delta_r^2 + Delta_z^2 + O(slope^2)',
      sp.series((Dv_t.T*Dv_t)[0] - dr**2 - dz**2, small, 0, 2).removeO())
# in the plane, Delta = Delta_r cos(phi) (cos(phi), -sin(phi), 0) exactly, cos(phi) = v_theta'/v',
# sin(phi) = v_r'/v': hahn13's perpendicular distance is |Delta| = Delta_r cos(phi), and A_g along
# Delta has hahn13 Eq. 20's components A_g cos(phi), -A_g sin(phi)
Dp = Dv.subs({dz: 0, vz_: 0})
cph, sph = vt_/sp.sqrt(vr_**2 + vt_**2), vr_/sp.sqrt(vr_**2 + vt_**2)
check('G9: planar Delta = Delta_r cos(phi) (cos(phi), -sin(phi), 0) (hahn13 Fig. 3, Eq. 20)',
      Dp - dr*cph*sp.Matrix([cph, -sph, 0]))
