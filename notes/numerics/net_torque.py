#net_torque.py
#
#Net torque of the ring's self-gravity, sum_i m r A_theta, relative to sum_i |m r A_theta|, for the
#wire kernel variants: lambda of the perturbed particle (epi_int_lite) or the perturber's own, and the
#wire's direction from the velocity (epi_int_lite) or from the shape (Lagrange slope).
#(notes/numerics/angular_momentum.md)
#
#Run from the repo root:  uv run python notes/numerics/net_torque.py
import sys, numpy as np, runpy, pickle
sys.path.insert(0, 'source')
from helper_fns import wrap_ring, advance, sidestep, lagrange_poly_fit, get_lambda, sort_particles, ring_gravity

def lag_deriv(x0, x1, x2, y0, y1, y2, x):
    d0 = (2*x - x1 - x2)/((x0 - x1)*(x0 - x2))
    d1 = (2*x - x0 - x2)/((x1 - x0)*(x1 - x2))
    d2 = (2*x - x0 - x1)/((x2 - x0)*(x2 - x1))
    return y0*d0 + y1*d1 + y2*d2

def neighbors(t, f, n):
    t1 = advance(t, -n); f1 = advance(f, -n)
    return sidestep(t1, 1), t1, sidestep(t1, -1), sidestep(f1, 1), f1, sidestep(f1, -1)

def gravity(lam, r, t, vr, vt, own_lambda, shape):
    Ar = 0; At = 0
    for n in range(1, r.shape[0]):
        t0, t1, t2, r0, r1, r2 = neighbors(t, r, -n)
        rp = lagrange_poly_fit(t0, t1, t2, r0, r1, r2, t)
        dr = rp - r
        L = lagrange_poly_fit(*neighbors(t, lam, -n)[:3], *neighbors(t, lam, -n)[3:], t) if own_lambda else lam
        A = 2*L/dr
        if shape:
            s = lag_deriv(t0, t1, t2, r0, r1, r2, t)/r          # s_r = dr'/(r dtheta)
            Ar += A; At -= A*s                                   # exact planar wire, (1, -s_r)
        else:
            vri = lagrange_poly_fit(*neighbors(t, vr, -n)[:3], *neighbors(t, vr, -n)[3:], t)
            vti = lagrange_poly_fit(*neighbors(t, vt, -n)[:3], *neighbors(t, vt, -n)[3:], t)
            vi = np.sqrt(vri**2 + vti**2)
            Ar += A*vti/vi; At -= A*vri/vi
    return Ar, At

d = 'benchmark/runs/lamfix/off/'
p = runpy.run_path(d + 'inputs.py'); p.update(pickle.load(open(d + 'output/monitor.pkl', 'rb'))['modified_params'])
for k in (0, -1):
    r, t, vr, vt = (np.load(d + f'output/{x}.npy')[k].copy() for x in ('r', 't', 'vr', 'vt'))
    r, t, vr, vt = sort_particles(r, t, vr, vt)
    lam = get_lambda(p['total_ring_mass'], r.shape[0], p['J2'], p['Rp'], r, t, vr, vt)
    w = [wrap_ring(x, longitude=(i == 1)) for i, x in enumerate((r, t, vr, vt, lam))]
    rw, tw, vrw, vtw, lw = w
    ref = ring_gravity(lw, 1.0, rw, tw, vrw, vtw, False)
    mine = gravity(lw, rw, tw, vrw, vtw, False, False)
    assert np.allclose(ref[1][:, 1:-1], mine[1][:, 1:-1], rtol=1e-12, atol=1e-30)
    mp = p['total_ring_mass']/r.size
    print(f"{'initial' if k == 0 else 'final'} state:")
    for own in (False, True):
        for shape in (False, True):
            Ar, At = gravity(lw, rw, tw, vrw, vtw, own, shape)
            tq = mp*rw[:, 1:-1]*At[:, 1:-1]
            print(f"   lambda {'own' if own else 'perturbed'}  direction {'shape   ' if shape else 'velocity'}"
                  f"  net/sum|torque| = {tq.sum()/np.abs(tq).sum():+.2e}")
