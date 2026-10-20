#viscous_overstability.py
#
#Sympy checks O1-O5 (notes/equations.md): axisymmetric viscous overstability of a cold (pressureless),
#self-gravitating ringlet, in the continuum and as epi_int_lite discretizes it on a lattice of circular
#streamlines spaced Delta apart. Keplerian (kappa = Omega), constant kinematic viscosities nu_s (shear)
#and nu_b (bulk), so beta = dln(nu)/dln(sigma) = 0. Perturbations ~ exp(i k x + s t).
#
#Run:  uv run python notes/sympy/viscous_overstability.py

import sympy as sp

s, k, W, D = sp.symbols('s k Omega Delta', positive=True)
nu, nub, G, sig, x = sp.symbols('nu_s nu_b G sigma x', positive=True)
phi = sp.symbols('phi', real=True)
n = sp.symbols('n', integer=True)
xi, y = sp.symbols('xi y')                              #radial and azimuthal Lagrangian displacements
I = sp.I
nu1 = nub + sp.Rational(4, 3)*nu                         #radial (compressional) viscosity, nu_b + 4 nu_s/3

#--- O1: continuum dispersion relation, Eulerian, cold, kappa = Omega -----------------------------
#continuity  s sigma1 + i k sigma u = 0
#radial      s u - 2 Omega v = -k^2 nu1 u - i k phi1,   phi1 = -2 pi G sigma1/k
#azimuthal   s v + (Omega/2) u = -k^2 nu v - (3/2) Omega nu (1 + beta) i k sigma1/sigma,  beta = 0
u, v, s1 = sp.symbols('u v sigma1')
eqs = [s*s1 + I*k*sig*u,
       s*u - 2*W*v + k**2*nu1*u + I*k*(-2*sp.pi*G*s1/k),
       s*v + W/2*u + k**2*nu*v + sp.Rational(3, 2)*W*nu*I*k*s1/sig]
M = sp.Matrix([[sp.diff(e, q) for q in (s1, u, v)] for e in eqs])
disp_E = sp.expand(M.det()/sig)
print('O1 continuum dispersion relation (x sigma):', sp.collect(sp.expand(disp_E*sig), s))

#--- O2: the same system in Lagrangian (Hill) form, with general K for k^2 and g for 2 pi G sigma k --
#xi'' - 2 Omega y' - 3 Omega^2 xi = g xi - nu1 K xi'
#y''  + 2 Omega xi'               = -nu K (y' + 3 Omega xi)
K, g = sp.symbols('K g', positive=True)
H = sp.Matrix([[s**2 - 3*W**2 - g + nu1*K*s, -2*W*s],
               [2*W*s + 3*W*nu*K, s**2 + nu*K*s]])
disp_L = sp.expand(H.det())
#Lagrangian and Eulerian forms agree when K = k^2 and g = 2 pi G sigma k: the Lagrangian determinant
#carries an extra factor s (the neutral azimuthal drift), the Eulerian one a factor 1/sigma
lag = sp.expand(sp.cancel(disp_L.subs({K: k**2, g: 2*sp.pi*G*sig*k})/s))
print('O2 Lagrangian = Eulerian dispersion relation:', sp.expand(lag - disp_E*sig) == 0)

#--- O3: growth rate to first order in viscosity, with g = x Omega^2 -------------------------------
eps, d = sp.symbols('epsilon d')
w0 = W*sp.sqrt(1 - x)
Dl = disp_L.subs(g, x*W**2).subs({nu: eps*nu, nub: eps*nub}).subs(s, I*w0 + eps*d)
d1 = sp.solve(sp.expand(sp.series(Dl, eps, 0, 2).removeO()).coeff(eps, 1), d)[0]
gamma = sp.simplify(sp.re(sp.expand_complex(d1)))
claim = K*(nu*(4*x + 2) - 3*nub*(1 - x))/(6*(1 - x))
print('O3 growth rate  gamma = K [nu_s(4x + 2) - 3 nu_b(1 - x)]/[6(1 - x)]:', sp.simplify(gamma - claim) == 0)
#stability criterion and the thresholds quoted in the text
crit = sp.solve(sp.Eq(nu*(4*x + 2), 3*nub*(1 - x)), x)[0]
print('    marginal x for nu_b = nu_s:', sp.nsimplify(crit.subs(nub, nu)), ' for nu_b = 2 nu_s:', sp.nsimplify(crit.subs(nub, 2*nu)))
print('    no gravity (x = 0): stable iff nu_b/nu_s >', sp.Rational(2, 3))

#--- O3b: literature check, beta = dln(nu)/dln(sigma) general, no self-gravity -----------------------
#Lehmann, Schmidt & Salo (2017, ApJ 851, 125, Sec. 2) give the onset of overstability as beta > beta_c =
#(1/3)(gamma - 2/3 - F3/nu0), gamma = nu_b/nu_s, which with F3 = 0 (no thermal effects) is Schmit &
#Tscharnuter's (1995) value. Our cold Eulerian system with (1 + beta) kept, to first order in nu:
beta = sp.symbols('beta', real=True)
eqsb = [s*s1 + I*k*sig*u,
        s*u - 2*W*v + k**2*nu1*u,
        s*v + W/2*u + k**2*nu*v + sp.Rational(3, 2)*W*nu*(1 + beta)*I*k*s1/sig]
Mb = sp.Matrix([[sp.diff(e, q) for q in (s1, u, v)] for e in eqsb])
Db = sp.expand(Mb.det()/sig).subs({nu: eps*nu, nub: eps*nub}).subs(s, I*W + eps*d)
db = sp.solve(sp.expand(sp.series(Db, eps, 0, 2).removeO()).coeff(eps, 1), d)[0]
gb = sp.simplify(sp.re(sp.expand_complex(db)))
beta_c = sp.solve(sp.Eq(gb, 0), beta)[0]
print('O3b growth rate without gravity:', sp.factor(gb), ';  beta_c =', sp.simplify(beta_c))
print('    equals Lehmann et al. (2017) (gamma - 2/3)/3 with gamma = nu_b/nu_s:',
      sp.simplify(beta_c - (nub/nu - sp.Rational(2, 3))/3) == 0)

#--- O4: the code's viscous stencils on a lattice, Bloch wave f_j = f exp(i j phi) ------------------
#delta_f(f)_j = (f_{j+1} - f_{j-1})/2 -> i sin(phi) f;  sigma_j = lambda/delta_r_j, delta_r_j = Delta + delta_f(xi)_j
#shear:  w_j = Omega_j + y'_j/a;  P_j = -nu sigma_j a delta_f(w)_j/delta_r_j;  A_t = -delta_f(P)_j/lambda
#bulk:   P_j = -nu1 sigma_j delta_f(v_r)_j/delta_r_j (the v_r/r term is O(Delta/a) smaller);  A_r = -delta_f(P)_j/lambda
lam, a, eta = sp.symbols('lambda a eta', positive=True)
yd, xid, ud = sp.symbols('ydot xi_amp u_amp')
Dl_ = lambda f: I*sp.sin(phi)*f                          #delta_f on a Bloch amplitude
dOmega = -sp.Rational(3, 2)*W/a                          #dOmega/da, Keplerian
#background + eta * perturbation, linearized
delta_r = D + eta*Dl_(xid)
delta_w = dOmega*D + eta*Dl_(yd)/a                       #(w_{j+1} - w_{j-1})/2
P = -nu*(lam/delta_r)*a*delta_w/delta_r
dP = sp.diff(P, eta).subs(eta, 0)                        #perturbation amplitude of P_j
At = sp.simplify(-Dl_(dP)/lam)
print('O4 shear stencil  A_t = -nu sin^2(phi)/Delta^2 (ydot + 3 Omega xi):',
      sp.simplify(At + nu*sp.sin(phi)**2/D**2*(yd + 3*W*xid)) == 0)
Pb = -nu1*(lam/delta_r)*(eta*Dl_(ud))/delta_r
Ar = sp.simplify(-Dl_(sp.diff(Pb, eta).subs(eta, 0))/lam)
print('   bulk stencil   A_r = -nu1 sin^2(phi)/Delta^2 u:', sp.simplify(Ar + nu1*sp.sin(phi)**2/D**2*ud) == 0)
print('   so the lattice is O2-O3 with K = sin^2(phi)/Delta^2, which vanishes at the zigzag phi = pi')

#--- O5: self-gravity strength at the zigzag, x = X(pi)/kappa^2 = R/2 (P6's X, R = 2 pi^2 G sigma/kappa^2 Delta) ---
Xpi = 4*G*sig/D*(sp.pi*sp.pi/2 - sp.pi**2/4)
R = 2*sp.pi**2*G*sig/(W**2*D)
print('O5 x(zigzag) = R/2:', sp.simplify(Xpi/W**2 - R/2) == 0)
print('   so with nu_b = nu_s the lattice is stable for R < 2/7, with nu_b = 2 nu_s for R < 4/5')
