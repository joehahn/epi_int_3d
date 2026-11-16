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


def paper_figure(name):
    """The paper's figure script paper/figures/<name>.py as a module, so a notebook draws its figure with the same
    code, figure(p, o), that writes the paper's PDF: the two are identical (style/testing.md T28)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location('paper_figure_' + name, REPO/'paper'/'figures'/(name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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
    exact solution B11 (notes/equations.md), and its steady forced part B6, z_forced = F_z/(nu^2 - omega_m^2) at the
    final mass. A circular particle is at theta = theta_0 + Omega_a t, so the forcing's phase chi_0 + omega_m t
    follows from its current longitude, whatever its index."""
    import satellite as satellite_module
    sat = satellite_module.Satellite(p['satellite'], p['J2'], p['Rp'])
    a, theta, time = o['r'], o['t'], o['times'][:, None, None]
    m, gam = sat.m, 1.0/sat.time_grow
    nu, w = nu_a(p['J2'], p['Rp'], a), m*(Omega(p['J2'], p['Rp'], a) - sat.Omega_p)
    lam, node = sat.longitudes(time)
    chi = m*theta - (m + 1)*lam + node                                     #chi_0 + omega_m t
    F = 1j*sat.mass_final*sat.A_V_per_mass*np.exp(1j*(chi - w*time))       #F_z = Re[F e^{i omega_m t}](1 - e^{-gam t})
    s = 1j*w - gam
    zp = lambda t: F*(np.exp(1j*w*t)/(nu**2 - w**2) - np.exp(s*t)/(s**2 + nu**2))
    vp = lambda t: F*(1j*w*np.exp(1j*w*t)/(nu**2 - w**2) - s*np.exp(s*t)/(s**2 + nu**2))
    #free part C e^{i nu t} + D e^{-i nu t} with z(0) = z'(0) = 0 (for free_over_forced only)
    C = -(zp(0) + vp(0)/(1j*nu))/2
    D = -(zp(0) - vp(0)/(1j*nu))/2
    #B11 itself, evaluated as the equivalent Duhamel integral z = int_0^t sin(nu(t - t'))/nu f(t') dt', which does not
    #lose precision as nu - omega_m -> 0 (where the closed form's terms cancel): f is a sum of c e^{q t'}, and
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


def ring_sigma(p):
    """The surface density of a ring of uniform sigma spanning 1 <= r <= 1 + radial_width."""
    return p['total_ring_mass']/(np.pi*((1 + p['radial_width'])**2 - 1))


def bending_wave(sat, sigma, x, G=1.0):
    """The steady forced bending wave h(x) of B8 (notes/equations.md), Z = Re[h e^{i chi}]:
    h = -i K F e^{-i L x^2/2} int_x^inf e^{i L s^2/2} ds, K = r_V/(2 pi G sigma), L = D_V K, F = i A_V, written with
    the Fresnel integrals: int_x^inf e^{i L s^2/2} ds = (pi/L)^{1/2} [(1 + i)/2 - C(u) - i S(u)], u = x (L/pi)^{1/2}."""
    from scipy.special import fresnel
    K = sat.r_V/(2*np.pi*G*sigma)
    L = sat.D_V*K
    F = 1j*sat.mass_final*sat.A_V_per_mass
    S, C = fresnel(x*np.sqrt(L/np.pi))
    return -1j*K*F*np.exp(-1j*L*x**2/2)*np.sqrt(np.pi/L)*((1 + 1j)/2 - C - 1j*S)


def simulated_h(p, o, k):
    """Each streamline's complex amplitude h at output k, from a least-squares fit of its particles'
    z = Re[h e^{i chi}] = Re(h) cos chi - Im(h) sin chi, and its fractional distance x from r_V."""
    import satellite as satellite_module
    sat = satellite_module.Satellite(p['satellite'], p['J2'], p['Rp'])
    lam, node = sat.longitudes(o['times'][k])
    chi = sat.m*o['t'][k] - (sat.m + 1)*lam + node
    h = np.empty(chi.shape[0], dtype=complex)
    for j in range(chi.shape[0]):
        A = np.column_stack([np.cos(chi[j]), -np.sin(chi[j])])
        (re, im), *_ = np.linalg.lstsq(A, o['z'][k, j], rcond=None)
        h[j] = re + 1j*im
    x = (o['r'][k].mean(axis=1) - sat.r_V)/sat.r_V
    return sat, x, h


def satellite_torque(p, o, k):
    """The torque that the full resonant potential Phi = A_V z sin(chi), whose -dPhi/dz is fast/satellite.py's A_z,
    exerts on each streamline at output k: mu r A_theta = -mu dPhi/dtheta = -mu m A_V z cos(chi), summed over its
    particles (mu = particle mass). fast/satellite.py applies only A_z, so the ring itself does not feel this torque."""
    import satellite as satellite_module
    sat = satellite_module.Satellite(p['satellite'], p['J2'], p['Rp'])
    lam, node = sat.longitudes(o['times'][k])
    chi = sat.m*o['t'][k] - (sat.m + 1)*lam + node
    mu = p['total_ring_mass']/p['number_of_particles']
    A_V = sat.mass(o['times'][k])*sat.A_V_per_mass
    return -(mu*sat.m*A_V*o['z'][k]*np.cos(chi)).sum(axis=1)


def wave_wavenumber(p, sat, x, G=1.0):
    """The local wavenumber k(x) of the free bending wave at x < 0 in the continuum, near r_V where
    nu^2 - omega_m^2 = D_V x (B5): D_V|x| = 2 pi G sigma k (B7, the paper's Eqns. 25 and 37)."""
    return sat.D_V*np.abs(np.asarray(x, dtype=float))/(2*np.pi*G*ring_sigma(p))


#tests/bending_waves at three numbers of streamlines, the same ring otherwise (tools/run_simulations.py makes them):
#output folder -> modified inputs; number_of_particles is set from inputs.py's number_of_streamlines before -m applies
RESOLUTION_RUNS = {
    'output_N300': {'number_of_streamlines': 300, 'number_of_particles': 300*32},
    'output': {},
    'output_N1200': {'number_of_streamlines': 1200, 'number_of_particles': 1200*32},
}


def resolution_runs(test_folder):
    """tests/bending_waves' runs at 300, 600 and 1200 streamlines, as a list of (p, o) from load3d, coarsest first."""
    runs = []
    for folder, modified in RESOLUTION_RUNS.items():
        runs.append(load3d(test_folder, Path(test_folder)/folder, dict(modified, output_folder=folder)))
    return runs


def behind_front(p, o, sat, x):
    """Where the wave is in place at the end of the run: behind its front, which travels inward from r_V at
    c_g = pi G sigma/nu_a (Eqn. 26), and outside the resonance zone, |x| > 2 (2 pi sigma/(r_V D_V))^{1/2}."""
    sigma = ring_sigma(p)
    x_front = -np.pi*sigma*o['times'][-1]/(nu_a(p['J2'], p['Rp'], sat.r_V)*sat.r_V)
    return (x < -2.0*np.sqrt(2*np.pi*sigma/(sat.r_V*sat.D_V))) & (x > x_front), x_front


def resonance_D(p, a):
    """D = nu_a^2 - omega_m^2, omega_m = m(Omega_a - Omega_p) (paper Eqn. 25), at semimajor axes a, under the test's
    vertical Fourier satellite. Returns sat and D."""
    import satellite as satellite_module
    sat = satellite_module.Satellite(p['satellite'], p['J2'], p['Rp'])
    a = np.asarray(a, dtype=float)
    return sat, nu_a(p['J2'], p['Rp'], a)**2 - (sat.m*(Omega(p['J2'], p['Rp'], a) - sat.Omega_p))**2


LATE_TAU = 6.0                          #tests/vertical_resonance's averages are over t > LATE_TAU growth times


def late_outputs(p, o):
    """The outputs that tests/vertical_resonance averages over, t > LATE_TAU growth times tau, once the satellite is
    grown (its mass within e^-6 = 0.25% of its final mass); the paper's figure, numbers and the notebook all use it."""
    return o['times'] > LATE_TAU*p['satellite']['time_grow']


def ring_angular_momentum(p, o):
    """The ring's total angular momentum about the z axis at each output, sum of mu r v_theta over its particles, in
    the code's coordinates (planetocentric positions, barycentric velocities, whose sum is conserved: Chambers 1999).
    The vertical Fourier satellite pulls only along z, so it exerts no torque about z, and only the ring's gravity
    and the integrator can change this."""
    mu = p['total_ring_mass']/p['number_of_particles']
    return (mu*o['r']*o['vt']).sum(axis=(1, 2))




def wavenumber_profile(x, w, r_V):
    """The local radial wavenumber k = d(arg w)/dr of a complex wave w(x) sampled at x, at the midpoints xm, from
    the phase change between neighbors (so its sign is kept: here k > 0 is a trailing wave, paper Section 4.1.1)."""
    xm = 0.5*(x[1:] + x[:-1])
    return xm, np.angle(w[1:]/w[:-1])/np.diff(x)/r_V





