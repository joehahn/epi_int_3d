#vertical_analysis.py
#
#Diagnostics of epi_int_3d's vertical motion (fast/epi_int_3d.py), for the inclined test
#(tests/inclined/inclined.ipynb) and the regression suite (test_vertical.py):
#  - inclined_metrics: each particle's inclination should stay constant and its node should
#    regress at Omega_a - nu_a (notes/equations.md V1-V4, V8);
#  - exact_orbits: the same particles integrated exactly in the 3D J2 potential, an independent
#    check whose nodal rates differ from Omega_a - nu_a by O(e^2, i^2).

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO/'source'))
sys.path.insert(0, str(REPO/'fast'))
from helper_fns import Omega
from epi_int_3d import nu_a, restore_output, vertical_elements
from analysis import read_inputs


def load3d(test_folder, output_folder=None, modified_params=None):
    """A test's parameters, and its output with the elements a, e, wt, M, inclination and node.
    Output 0 is skipped: like the reference driver's, it holds barycentric velocities."""
    p = read_inputs(test_folder, modified_params)
    folder = Path(output_folder) if output_folder else Path(test_folder)/p['output_folder']
    r, t, z, vr, vt, vz, times, monitor = restore_output(str(folder))
    r, t, z, vr, vt, vz, times = (x[1:] for x in (r, t, z, vr, vt, vz, times))
    a, e, wt, M, inc, node = vertical_elements(p['J2'], p['Rp'], r, t, z, vr, vt, vz)
    return p, dict(r=r, t=t, z=z, vr=vr, vt=vt, vz=vz, times=times, monitor=monitor, a=a, e=e, wt=wt, M=M,
                   inc=inc, node=np.unwrap(node, axis=0))


def nodal_rate(p, a):
    """Omega_a - nu_a (V4)."""
    return Omega(p['J2'], p['Rp'], a) - nu_a(p['J2'], p['Rp'], a)


def inclined_metrics(p, o):
    """Largest change in each particle's inclination, and largest departure of its node from
    node_0 + (Omega_a - nu_a)(t - t_0), over the run."""
    dt = o['times'] - o['times'][0]
    rate = nodal_rate(p, o['a'][0])
    node_theory = o['node'][0] + rate*dt[:, None, None]
    return {'max_delta_inc_rel': (np.abs(o['inc'] - o['inc'][0])/o['inc'][0]).max(),
            'max_delta_node': np.abs(o['node'] - node_theory).max(),
            'max_delta_a': np.abs(o['a'] - o['a'][0]).max(),
            'regression_orbits': (2*np.pi/np.abs(rate)).min()/(2*np.pi)}


def _rhs(t, s, J2, Rp):
    """Exact equations of motion in the J2 potential, Cartesian, for many particles (GM = 1)."""
    x, y, z, vx, vy, vz = s.reshape(6, -1)
    rho2 = x*x + y*y + z*z
    f = 1.0/rho2**1.5
    q = 1.5*J2*Rp*Rp/rho2
    zz = 5*z*z/rho2
    return np.concatenate([vx, vy, vz, -f*x*(1 + q*(1 - zz)), -f*y*(1 + q*(1 - zz)), -f*z*(1 + q*(3 - zz))])


def exact_orbits(p, o, idx, rtol=1e-11):
    """Integrate particles idx (index arrays into a streamline-by-particle array) exactly from their
    state at the first output; return their nodes (V3) at the output times, and the code's."""
    k = (slice(None),) + idx
    r, t, z, vr, vt, vz = (o[c][0][idx] for c in ('r', 't', 'z', 'vr', 'vt', 'vz'))
    c, s = np.cos(t), np.sin(t)
    s0 = np.concatenate([r*c, r*s, z, vr*c - vt*s, vr*s + vt*c, vz])
    times = o['times'] - o['times'][0]
    sol = solve_ivp(_rhs, (0, times[-1]), s0, t_eval=times, args=(p['J2'], p['Rp']), rtol=rtol, atol=1e-13,
                    method='DOP853')
    x, y, z, vx, vy, vz = sol.y.reshape(6, len(r), -1)
    R, T = np.hypot(x, y), np.arctan2(y, x)                  #node by V3, as in the code's output
    VR, VT = (x*vx + y*vy)/R, (x*vy - y*vx)/R
    node = vertical_elements(p['J2'], p['Rp'], R, T, z, VR, VT, vz)[-1].T
    return np.unwrap(node, axis=0), o['node'][k]


def forced_vertical(p, o):
    """Each particle's z(t) for a zero-mass ring of circular test particles under fast/satellite.py's vertical
    Fourier satellite, from z = vz = 0 at t = 0 while the satellite's mass grows as 1 - exp(-t/time_grow): the
    exact solution B11 (notes/equations.md), and its steady forced part B6, z_forced = F_z/(nu^2 - omega_D^2) at the
    final mass. A circular particle is at theta = theta_0 + Omega_a t, so the forcing's phase chi_0 + omega_D t
    follows from its current longitude, whatever its index."""
    import satellite as satellite_module
    sat = satellite_module.Satellite(p['satellite'], p['J2'], p['Rp'])
    a, theta, time = o['r'], o['t'], o['times'][:, None, None]
    m, gam = sat.m, 1.0/sat.time_grow
    nu, w = nu_a(p['J2'], p['Rp'], a), m*(Omega(p['J2'], p['Rp'], a) - sat.Omega_p)
    lam, node = sat.longitudes(time)
    chi = m*theta - (m + 1)*lam + node                                     #chi_0 + omega_D t
    F = 1j*sat.mass_final*sat.A_V_per_mass*np.exp(1j*(chi - w*time))       #F_z = Re[F e^{i omega_D t}](1 - e^{-gam t})
    s = 1j*w - gam
    zp = lambda t: F*(np.exp(1j*w*t)/(nu**2 - w**2) - np.exp(s*t)/(s**2 + nu**2))
    vp = lambda t: F*(1j*w*np.exp(1j*w*t)/(nu**2 - w**2) - s*np.exp(s*t)/(s**2 + nu**2))
    #free part C e^{i nu t} + D e^{-i nu t} with z(0) = z'(0) = 0 (for free_over_forced only)
    C = -(zp(0) + vp(0)/(1j*nu))/2
    D = -(zp(0) - vp(0)/(1j*nu))/2
    #B11 itself, evaluated as the equivalent Duhamel integral z = int_0^t sin(nu(t - t'))/nu f(t') dt', which does not
    #lose precision as nu - omega_D -> 0 (where the closed form's terms cancel): f is a sum of c e^{q t'}, and
    #int_0^t e^{q t'} dt' = t phi1(q t), phi1(u) = (e^u - 1)/u
    def phi1(u):
        small = np.abs(u) < 1.0e-3
        us = np.where(small, 1.0, u)
        return np.where(small, 1 + u/2 + u**2/6 + u**3/24, np.expm1(us)/us)
    t = time
    z = 0.0
    for c, q in ((F/2, 1j*w), (np.conj(F)/2, -1j*w), (-F/2, s), (-np.conj(F)/2, np.conj(s))):
        z = z + c/(2j*nu)*(np.exp(1j*nu*t)*t*phi1((q - 1j*nu)*t) - np.exp(-1j*nu*t)*t*phi1((q + 1j*nu)*t))
    z = z.real
    forced = (F*np.exp(1j*w*time)/(nu**2 - w**2)).real
    return sat, dict(z=z, forced=forced, x=(a[0].mean(axis=1) - sat.r_V)/sat.r_V,
                     free_over_forced=(np.abs(C) + np.abs(D))*np.abs(nu**2 - w**2)/np.abs(F))
