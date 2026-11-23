#density_analysis.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 23 November 2026.
#
#Diagnostics of a satellite's inner Lindblad resonance (fast/satellite.py's 'horizontal' Fourier satellite), the
#in-plane counterparts of tests/vertical_analysis.py's bending-wave ones, for tests/lindblad_resonance (massless
#particles) and tests/density_waves (a self-gravitating ring's spiral density wave), the paper's Section 4.2.
#Theory: notes/equations.md L1-L9 (notes/sympy/density_wave.py).

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
for folder in ('source', 'fast', 'tests'):
    if str(REPO/folder) not in sys.path:
        sys.path.insert(0, str(REPO/folder))
from helper_fns import Omega, Kappa
from satellite import Satellite
from vertical_analysis import load3d, ring_sigma, paper_figure, late_outputs, wavenumber_profile

__all__ = ['load3d', 'ring_sigma', 'paper_figure', 'late_outputs', 'wavenumber_profile', 'satellite', 'resonance_D',
           'forced_eccentricity', 'density_wave', 'L_L', 'wave_amplitude', 'simulated_xi', 'satellite_torque',
           'torque_theory', 'wave_wavenumber', 'tau_001', 'resonance_zone', 'wave_front',
           'steady', 'X_STEADY']


def satellite(p):
    """The test's horizontal Fourier satellite."""
    return Satellite(p['satellite'], p['J2'], p['Rp'])


def resonance_D(p, a):
    """D = kappa_a^2 - omega_m^2, omega_m = m(Omega_a - Omega_p), at semimajor axes a (the paper's Section 4.2.1).
    Returns sat and D."""
    sat = satellite(p)
    a = np.asarray(a, dtype=float)
    return sat, Kappa(p['J2'], p['Rp'], a)**2 - (sat.m*(Omega(p['J2'], p['Rp'], a) - sat.Omega_p))**2


def L_L(sat):
    """L_L = -Psi/(r_L D_L) at the satellite's final mass (L7)."""
    return -sat.mass_final*sat.Psi_per_mass/(sat.r_L*sat.D_L)


def forced_eccentricity(sat, x):
    """A test particle's forced eccentricity near r_L, e_f = |L_L/x| = |Psi|/(r_L^2 D_L |x|) (L3, L9)."""
    return np.abs(L_L(sat)/np.asarray(x, dtype=float))


def density_wave(sat, sigma, x, G=1.0):
    """The steady forced density wave xi(x) of L7, r = a + Re[xi e^{i chi}], chi = m(theta - Omega_p t):
    xi = -(i L_L/eps) e^{i x^2/2 eps} int_{-inf}^x e^{-i s^2/2 eps} ds, eps = 2 pi G sigma/(r_L D_L), written with the
    Fresnel integrals: int_{-inf}^x e^{-i s^2/2 eps} ds = (pi eps)^{1/2} [(1 - i)/2 + C(u) - i S(u)], u = x/(pi eps)^{1/2}."""
    from scipy.special import fresnel
    eps = 2*np.pi*G*sigma/(sat.r_L*sat.D_L)
    u = np.asarray(x, dtype=float)/np.sqrt(np.pi*eps)
    S, C = fresnel(u)
    integral = np.sqrt(np.pi*eps)*((1 - 1j)/2 + C - 1j*S)
    return -(1j*L_L(sat)/eps)*np.exp(1j*np.asarray(x)**2/(2*eps))*integral


def wave_amplitude(p, sat, G=1.0):
    """The wave's radial amplitude |Z_w| = (|Psi|/r_L)(r_L/(G sigma D_L))^{1/2} (L7)."""
    return abs(sat.mass_final*sat.Psi_per_mass)/sat.r_L*np.sqrt(sat.r_L/(G*ring_sigma(p)*sat.D_L))


def simulated_xi(p, o, k):
    """Each streamline's complex radial displacement xi at output k, from a least-squares fit of its particles'
    r = r_0 + Re[xi e^{i chi}] = r_0 + Re(xi) cos chi - Im(xi) sin chi, chi = m(theta - lambda_s), and its fractional
    distance x = (r_0 - r_L)/r_L from resonance, r_0 its mean radius."""
    sat = satellite(p)
    lam = sat.lambda0 + sat.Omega_a*o['times'][k]
    chi = sat.m*(o['t'][k] - lam)
    xi = np.empty(chi.shape[0], dtype=complex)
    r0 = np.empty(chi.shape[0])
    for j in range(chi.shape[0]):
        A = np.column_stack([np.ones_like(chi[j]), np.cos(chi[j]), -np.sin(chi[j])])
        (r0[j], re, im), *_ = np.linalg.lstsq(A, o['r'][k, j], rcond=None)
        xi[j] = re + 1j*im
    return sat, (r0 - sat.r_L)/sat.r_L, xi


def satellite_torque(p, o, k):
    """The torque about z that the satellite exerts on each streamline at output k, sum of mu r A_theta over its
    particles, with fast/satellite.py's own A_theta: the torque the ring feels."""
    sat = satellite(p)
    mu = p['total_ring_mass']/p['number_of_particles']
    A_r, A_t = sat.horizontal_acceleration(o['r'][k], o['t'][k], o['times'][k])
    return (mu*o['r'][k]*A_t).sum(axis=1)


def torque_theory(p, sat, G=1.0):
    """Linear theory's torque on the ring, T_L = -pi^2 m sigma Psi^2/D_L (L8; Goldreich & Tremaine 1979)."""
    return -np.pi**2*sat.m*ring_sigma(p)*(sat.mass_final*sat.Psi_per_mass)**2/sat.D_L


def wave_wavenumber(p, sat, x, G=1.0):
    """The local wavenumber of the free density wave at x > 0 in the continuum, D_L x = 2 pi G sigma k (L5)."""
    return sat.D_L*np.abs(np.asarray(x, dtype=float))/(2*np.pi*G*ring_sigma(p))


X_STEADY = 0.01                     #tests/density_waves: the wave is steady over the resonance zone < x < X_STEADY by the end


def tau_001(p, sat):
    """tau_0.01, the time the density wave takes to travel from r_L to x = 0.01 at its group speed
    c_g = pi G sigma/kappa_a (L6), with kappa_a at r_L."""
    return 0.01*sat.r_L/(np.pi*ring_sigma(p)/Kappa(p['J2'], p['Rp'], sat.r_L))


def resonance_zone(p, sat):
    """The resonance zone's half width in x, 2 (2 pi G sigma/(r_L D_L))^{1/2} = 2 eps^{1/2}."""
    return 2.0*np.sqrt(2*np.pi*ring_sigma(p)/(sat.r_L*sat.D_L))


def wave_front(p, o, k):
    """The wave front's x at output k, measured: going outward from the resonance zone, the first streamline where the
    wave part of xi, |xi - L_L/x|, falls below half of linear theory's |Z_w|."""
    sat, x, xi = simulated_xi(p, o, k)
    Z_w = wave_amplitude(p, sat)
    xi_wave = np.abs(xi - L_L(sat)/x)
    outward = np.where(x > resonance_zone(p, sat))[0]
    below = outward[xi_wave[outward] < 0.5*Z_w]
    return x[below[0]] if len(below) else x[outward[-1]]


def steady(p, sat, x):
    """Where the wave is compared with linear theory at the run's end: outside the resonance zone, x < X_STEADY."""
    return (x > resonance_zone(p, sat)) & (x < X_STEADY)
