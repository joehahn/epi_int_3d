#!/usr/bin/env python3
# Sympy verification of the 3D gravity of a perturbing streamline treated as a straight, infinitely
# long wire (G1-G14; notes/equations.md, notes/phase2/streamline_gravity.md). This is the 3D
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

# --- G10: the wire along the streamline's shape, not its velocity (for m >= 2 patterns, G6). The
# Lagrange polynomial gives r'(theta), z'(theta); a step dtheta along the streamline moves
# (dr', r dtheta, dz') in the local frame, so the tangent is t = (s_r, 1, s_z)/sqrt(1 + s_r^2 + s_z^2)
# with slopes s_r = (1/r) dr'/dtheta, s_z = (1/r) dz'/dtheta. Then, to first order in the slopes,
# Delta = (Delta_r, -(Delta_r s_r + Delta_z s_z), Delta_z), and in the plane Delta = Delta_r cos(phi)
# (cos(phi), -sin(phi), 0) exactly with tan(phi) = s_r. ---
s_r, s_z = sp.symbols('s_r s_z', real=True)
tsh = sp.Matrix([s_r, 1, s_z])/sp.sqrt(1 + s_r**2 + s_z**2)
Dsh = d - (d.T*tsh)[0]*tsh
Dsh_t = Dsh.subs({s_r: small*s_r, s_z: small*s_z}, simultaneous=True)
first_sh = sp.Matrix([dr, -small*(dr*s_r + dz*s_z), dz])
check('G10: Delta = (Delta_r, -(Delta_r s_r + Delta_z s_z), Delta_z) + O(slope^2)',
      [sp.series(Dsh_t[k] - first_sh[k], small, 0, 2).removeO() for k in range(3)])
check('G10: Delta.t = 0 for that first-order Delta, to first order',
      sp.series((first_sh.T*tsh.subs({s_r: small*s_r, s_z: small*s_z}, simultaneous=True))[0], small, 0, 2).removeO())
Dsh_p = Dsh.subs({dz: 0, s_z: 0})
cps, sps = 1/sp.sqrt(1 + s_r**2), s_r/sp.sqrt(1 + s_r**2)
check('G10: planar Delta = Delta_r cos(phi) (cos(phi), -sin(phi), 0), tan(phi) = s_r',
      Dsh_p - dr*cps*sp.Matrix([cps, -sps, 0]))
# near an m-armed resonance, m (Omega - Omega_p) = +-kappa or +-nu, and kappa, nu = Omega to O(J2), so the
# velocity's slope is (1 - Omega_p/Omega) = +-1/m times the shape's (G6)
mm = sp.symbols('m', positive=True)
for sgn in (1, -1):
    Omp_res = sp.solve(sp.Eq(mm*(Om - Omp), sgn*Om), Omp)[0]
    check(f'G10: at an m-armed resonance (sign {sgn:+d}), 1 - Omega_p/Omega = {sgn:+d}/m', 1 - Omp_res/Om - sgn/mm)

# --- G11: the closed form used in the paper (exact for a straight wire, any slopes). The wire passes
# through the crossing point at separation d = (d_r, 0, d_z) from the particle, along t = (s_r, 1, s_z)
# (not normalized). With c = (d.t)/(t.t), the perpendicular from the particle to the wire is
# Delta = d - c t = (d_r - c s_r, -c, d_z - c s_z), |Delta|^2 = d_r^2 + d_z^2 - c (d_r s_r + d_z s_z),
# and A_g = 2 G lam Delta/|Delta|^2 (G1). In the plane, A_g = (2 G lam/d_r)(1, -s_r, 0) exactly. ---
d_r, d_z = sp.symbols('d_r d_z', real=True)
dd = sp.Matrix([d_r, 0, d_z])
tt = sp.Matrix([s_r, 1, s_z])
cc = (d_r*s_r + d_z*s_z)/(1 + s_r**2 + s_z**2)
Dl = sp.Matrix([d_r - cc*s_r, -cc, d_z - cc*s_z])
check('G11: c = (d.t)/(t.t)', cc - (dd.T*tt)[0]/(tt.T*tt)[0])
check('G11: Delta = d - c t', Dl - (dd - cc*tt))
check('G11: Delta is perpendicular to the wire', (Dl.T*tt)[0])
check('G11: |Delta|^2 = d_r^2 + d_z^2 - c (d_r s_r + d_z s_z)', (Dl.T*Dl)[0] - (d_r**2 + d_z**2 - cc*(d_r*s_r + d_z*s_z)))
check('G11: = G1 (the exact wire) with T = t/|t|', Dl - Dsh.subs({dr: d_r, dz: d_z}))
A_pl = (2*G*lam*Dl/(Dl.T*Dl)[0]).subs({d_z: 0, s_z: 0})
check('G11: planar A_g = (2 G lam/d_r)(1, -s_r, 0)', A_pl - sp.Matrix([2*G*lam/d_r, -2*G*lam*s_r/d_r, 0]))

# --- G12: the paper's derivation, hahn13 Fig. 3 moved and tilted. (1) Raising particle and wire by z
# changes nothing. (2) Raising the wire by d_z, still horizontal and tilted by phi (tan phi = s_r) in its
# plane: the perpendicular is hahn13's in-plane one, length Delta_h = d_r cos(phi) along
# (cos phi, -sin phi), plus the vertical leg d_z -- exact. (3) Tilting the wire by its vertical slope s_z
# adds -d_z s_z to Delta_theta; the remaining errors are O(s_r s_z, s_z^2). ---
cph2, sph2 = 1/sp.sqrt(1 + s_r**2), s_r/sp.sqrt(1 + s_r**2)
Dh = d_r*cph2
D12 = sp.Matrix([Dh*cph2, -Dh*sph2, d_z])
check('G12: step 2 (horizontal raised wire) is exact: = G11 at s_z = 0', D12 - Dl.subs(s_z, 0))
check('G12: step 2, |Delta|^2 = d_r^2 cos^2(phi) + d_z^2', (D12.T*D12)[0] - (d_r**2*cph2**2 + d_z**2))
D123 = D12 + sp.Matrix([0, -d_z*s_z, 0])
err12 = (Dl - D123).subs(s_z, small*s_z)
check('G12: step 3 matches G11 to first order in s_z, up to O(s_r s_z)',
      [sp.series(sp.series(err12[k], small, 0, 2).removeO().coeff(small, 1), s_r, 0, 1).removeO() for k in range(3)])
check('G12: and exactly at O(s_z^0)', [sp.series(err12[k], small, 0, 1).removeO() for k in range(3)])
# the accelerations of the paper, 2 G lam Delta/Delta^2 with |Delta|^2 to the same order
Dsq = d_r**2*cph2**2 + d_z**2
A12 = 2*G*lam*sp.Matrix([d_r*cph2**2, -(d_r*sph2*cph2 + d_z*s_z), d_z])/Dsq
check('G12: A = 2 G lam (d_r cos^2 phi, -(d_r sin phi cos phi + d_z s_z), d_z)/Delta^2', A12 - 2*G*lam*D123/Dsq)
check('G12: in the plane, A = (2 G lam/d_r)(1, -s_r, 0), Eqn. 14 with Delta = d_r cos(phi)',
      A12.subs(d_z, 0) - sp.Matrix([2*G*lam/d_r, -2*G*lam*s_r/d_r, 0]))

# --- G13: the planned kernel, with the distance taken at the crossing point:
# Delta^2 = d_r^2 + d_z^2, A = 2 G lam (d_r cos phi, -(d_r sin phi + d_z s_z), d_z)/Delta^2. It equals
# G12 (the paper's Eqn. 17) to first order in the slopes, and in the plane it is epi_int_lite's
# (2 G lam/d_r)(cos phi, -sin phi). ---
Akern = 2*G*lam*sp.Matrix([d_r*cph2, -(d_r*sph2 + d_z*s_z), d_z])/(d_r**2 + d_z**2)
dk = (Akern - A12).subs({s_r: small*s_r, s_z: small*s_z}, simultaneous=True)
check('G13: kernel = G12 + O(slope^2)', [sp.series(dk[k], small, 0, 2).removeO() for k in range(3)])
check("G13: in the plane, kernel = epi_int_lite's (2 G lam/d_r)(cos phi, -sin phi, 0)",
      Akern.subs(d_z, 0) - 2*G*lam/d_r*sp.Matrix([cph2, -sph2, 0]))

# --- G14: the 2D kernel options (fast/kernels.py, nearpoint and shape_slope). The shape slope
# s_r = (dr'/dtheta)/r' at the crossing point is the derivative of lagrange_poly_fit's 2nd-order
# polynomial through (x0, y0), (x1, y1), (x2, y2), whose weights' derivatives are
# dl0 = (dx1 + dx2)/(dx10 dx20), dl1 = -(dx0 + dx2)/(dx10 dx21), dl2 = (dx0 + dx1)/(dx20 dx21),
# with dxk = x - xk and dxkl = xk - xl as coded. (r' rather than r, since tan(phi) = dr'/(r' dtheta)
# is the polar curve's own tangent; they differ at O(d_r/r), which the wire model drops.) ---
X0, X1, X2, Xx, Y0, Y1, Y2 = sp.symbols('x0 x1 x2 x y0 y1 y2', real=True)
dx0, dx1, dx2 = Xx - X0, Xx - X1, Xx - X2
dx10, dx20, dx21 = X1 - X0, X2 - X0, X2 - X1
l0, l1, l2 = (dx1*dx2)/(dx10*dx20), -(dx0*dx2)/(dx10*dx21), (dx0*dx1)/(dx20*dx21)
poly = Y0*l0 + Y1*l1 + Y2*l2
check('G14: lagrange_poly_fit passes through its three points',
      [poly.subs(Xx, X0) - Y0, poly.subs(Xx, X1) - Y1, poly.subs(Xx, X2) - Y2])
dl0, dl1, dl2 = (dx1 + dx2)/(dx10*dx20), -(dx0 + dx2)/(dx10*dx21), (dx0 + dx1)/(dx20*dx21)
check('G14: d/dx of lagrange_poly_fit = y0 dl0 + y1 dl1 + y2 dl2', sp.diff(poly, Xx) - (Y0*dl0 + Y1*dl1 + Y2*dl2))
# The options' planar pulls, A = (2 G lam/d_r) x (as coded):
#   reference (epi_int_lite):          (cos phi, -sin phi),  tan phi = v_r'/v_theta'
#   shape_slope:                       (cos phi, -sin phi),  tan phi = s_r
#   nearpoint:                         (1, -tan phi),        tan phi = v_r'/v_theta' (hahn13 Eq. 20)
#   nearpoint + shape_slope:           (1, -s_r)             (= G11, G12 in the plane)
# nearpoint is G3's exact planar wire; with tan phi = s_r it is G12's planar limit:
tph = sp.symbols('tau', real=True)
check('G14: nearpoint = exact planar wire (G3) with eps = tan phi',
      sp.Matrix([2*G*lam/dr, -2*G*lam*tph/dr, 0]) - A_2D.subs(eps, tph))
check('G14: the cos phi, sin phi options = nearpoint x cos phi',
      sp.Matrix([1/sp.sqrt(1 + tph**2), -tph/sp.sqrt(1 + tph**2)]) - sp.Matrix([1, -tph])/sp.sqrt(1 + tph**2))
