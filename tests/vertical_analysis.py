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
    """The local wavenumber k(x) of the free bending wave at x < 0, from three dispersion relations, near r_V
    where nu^2 - omega_m^2 = D_V x (B5): the continuum's D_V|x| = 2 pi G sigma k (B7); a lattice of streamlines
    spaced Delta, D_V|x| = X(k Delta) = (4 G sigma/Delta)(pi k Delta/2 - (k Delta)^2/4) (the paper's zigzag_X); and
    that lattice cut off at the ring's edges, which lacks the pull of the streamlines beyond them,
    D_V|x| = X(k Delta) - 2 G sigma (1/d_in + 1/d_out), with d the distance to each edge plus Delta/2 (B13)."""
    from scipy.optimize import brentq
    sig = ring_sigma(p)
    W, N = p['radial_width'], p['number_of_streamlines']
    Delta = W/(N - 1)
    r = sat.r_V*(1 + np.asarray(x, dtype=float))
    d_in, d_out = r - 1 + Delta/2, 1 + W - r + Delta/2
    lattice = lambda k: (4*G*sig/Delta)*(np.pi*k*Delta/2 - (k*Delta)**2/4)
    k_cont = sat.D_V*np.abs(x)/(2*np.pi*G*sig)
    k_lat = np.array([brentq(lambda k: lattice(k) - sat.D_V*abs(xx), 0.0, np.pi/Delta) for xx in np.ravel(x)])
    k_edge = np.array([brentq(lambda k: lattice(k) - 2*G*sig*(1/di + 1/do) - sat.D_V*abs(xx), 0.0, np.pi/Delta)
                       for xx, di, do in zip(np.ravel(x), np.ravel(d_in), np.ravel(d_out))])
    return k_cont, k_lat.reshape(np.shape(x)), k_edge.reshape(np.shape(x))


def resonance_frequencies(p, o, x=None, late_orbits=1.0e4):
    """nu_a^2 - omega_m^2 for each streamline of a zero-mass ring under the vertical Fourier satellite: simulated, from
    its forced amplitude, A_V/Im(h) (B6: z = Re[h e^{i chi}], Im h = A_V/(nu_a^2 - omega_m^2)), with h the least-squares
    amplitude of simulated_h averaged over the outputs after late_orbits, which averages away the free oscillation;
    and exact, at the streamline's mean semimajor axis, omega_m = m(Omega_a - Omega_p). Returns sat, x, simulated, exact.
    With x given, only the exact values at r = r_V(1 + x), and simulated = None."""
    import satellite as satellite_module
    sat = satellite_module.Satellite(p['satellite'], p['J2'], p['Rp'])
    exact = lambda a: nu_a(p['J2'], p['Rp'], a)**2 - (sat.m*(Omega(p['J2'], p['Rp'], a) - sat.Omega_p))**2
    if x is not None:
        x = np.asarray(x, dtype=float)
        return sat, x, None, exact(sat.r_V*(1 + x))
    late = np.where(o['times']/(2*np.pi) > late_orbits)[0]
    h = np.mean([simulated_h(p, o, k)[2] for k in late], axis=0)
    a = o['a'][late].mean(axis=(0, 2))
    A_V = sat.mass_final*sat.A_V_per_mass
    return sat, (a - sat.r_V)/sat.r_V, A_V/h.imag, exact(a)


def ring_angular_momentum(p, o):
    """The ring's total angular momentum about the z axis at each output, sum of mu r v_theta over its particles, in
    the code's coordinates (planetocentric positions, barycentric velocities, whose sum is conserved: Chambers 1999).
    The vertical Fourier satellite pulls only along z, so it exerts no torque about z, and only the ring's gravity
    and the integrator can change this."""
    mu = p['total_ring_mass']/p['number_of_particles']
    return (mu*o['r']*o['vt']).sum(axis=(1, 2))
