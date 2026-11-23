"""A "Fourier satellite" for epi_int_3d: only the resonant m-th term of a satellite's potential.

The satellite is set by the dict `satellite` in inputs.py (notes/phase2c/bending_wave_theory.md, section 4):
    satellite = {
        'resonance': 'vertical',   #'vertical' (inner vertical resonance) or 'horizontal' (inner Lindblad)
        'm': 2,                    #number of spiral arms, m >= 2 for 'vertical'
        'mass_final': 1.0e-9,      #satellite mass/planet mass once grown; < 0 turns the satellite off
        'time_grow': 1.0e3,        #mass = mass_final*(1 - exp(-time/time_grow)); <= 0: mass_final from time = 0
        'r': 2.0911,               #satellite's orbit radius a_s (circular)
        't': 0.0,                  #its initial mean longitude lambda_s
        'inclination': 1.0e-3,     #i_s, for 'vertical'
        'node': 0.0,               #its initial longitude of ascending node Omega_s, for 'vertical'
    }
The satellite's mean longitude advances at Omega_{a_s}, its node at Omega_{a_s} - nu_{a_s} (notes/equations.md B4).

'vertical': an inclined circular satellite pulls a ring particle at longitude theta with
    A_z = -A_V sin(m theta - (m+1) lambda_s + Omega_s),   A_V = G M_s i_s b_{3/2}^(m)(beta_V)/(2 a_s^2)   (B2)
and nothing in the plane; b is evaluated once, at the inner vertical resonance r_V, beta_V = r_V/a_s, where
m(Omega_a - Omega_p) = nu_a with m Omega_p = m Omega_{a_s} + nu_{a_s} (B4), solved with J2.

'horizontal': a circular satellite in the ring plane, the m-th term of its potential (L1),
    Phi = -(G M_s/a_s) b_{1/2}^(m)(beta) cos m(theta - lambda_s),   Omega_p = Omega_{a_s},
which pulls a ring particle in the plane with A_r = -dPhi/dr, A_theta = -(1/r) dPhi/dtheta, and nothing along z. As in
'vertical', b and its derivative are evaluated once, at the inner Lindblad resonance r_L, where m(Omega_a - Omega_p) =
kappa_a (L2), solved with J2: Phi's radial profile is its tangent line there, b(beta) ~= b_L + b'_L (beta - beta_L), so
A_r = (G M_s/a_s^2) b'_L cos m(theta - lambda_s) and the force stays the gradient of a potential. This is
epi_int_lite's lindblad_resonance branch, which evaluates b, b' at the ring's mean radius and b' by finite
difference; here b' is analytic. Psi = r dPhi/dr + (2 Omega_a/(Omega_a - Omega_p)) Phi at r_L sets the forced
response (L3), and D_L = r d(kappa^2 - omega_m^2)/dr at r_L (L4).
"""
import json

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


def Omega(J2, Rp, a):
    """Angular velocity of a circular orbit (helper_fns.Omega, GM = 1)."""
    return np.sqrt((1.0/a**3)*(1.0 + 1.5*J2*(Rp/a)**2))


def nu(J2, Rp, a):
    """Vertical frequency (notes/equations.md V1, GM = 1)."""
    return np.sqrt((1.0/a**3)*(1.0 + 4.5*J2*(Rp/a)**2))


def kappa(J2, Rp, a):
    """Epicyclic frequency (helper_fns.Kappa, GM = 1)."""
    return np.sqrt((1.0/a**3)*(1.0 - 1.5*J2*(Rp/a)**2))


def laplace(beta, m, s):
    """Laplace coefficient b_s^(m)(beta) = (2/pi) int_0^pi cos(m phi)/(1 + beta^2 - 2 beta cos phi)^s dphi."""
    f = lambda phi: np.cos(m*phi)/(1.0 + beta**2 - 2.0*beta*np.cos(phi))**s
    return (2.0/np.pi)*quad(f, 0.0, np.pi, limit=200, epsabs=0.0, epsrel=1e-13)[0]


def dlaplace(beta, m, s):
    """d b_s^(m)/d beta, analytically (notes/equations.md L1)."""
    f = lambda phi: -2.0*s*(beta - np.cos(phi))*np.cos(m*phi)/(1.0 + beta**2 - 2.0*beta*np.cos(phi))**(s + 1)
    return (2.0/np.pi)*quad(f, 0.0, np.pi, limit=200, epsabs=0.0, epsrel=1e-13)[0]


class Satellite:
    """The satellite's inputs, derived constants and its state at the current time."""

    def __init__(self, sat, J2, Rp):
        self.inputs = dict(sat)
        self.resonance = sat.get('resonance', 'vertical')
        if self.resonance not in ('vertical', 'horizontal'):
            raise ValueError("satellite resonance must be 'vertical' or 'horizontal', not '%s'" % self.resonance)
        self.m = int(sat['m'])
        if self.m < 2:
            raise ValueError('an inner vertical or Lindblad resonance needs m >= 2')
        self.mass_final = sat['mass_final']
        self.time_grow = sat['time_grow']
        self.a = sat['r']
        self.lambda0 = sat.get('t', 0.0)
        self.Omega_a = Omega(J2, Rp, self.a)
        self.nu_a = nu(J2, Rp, self.a)
        m = self.m
        h = 1.0e-6                                                                  #for D, numerically
        if self.resonance == 'vertical':
            self.inclination = sat['inclination']
            self.node0 = sat.get('node', 0.0)
            self.Omega_p = self.Omega_a + self.nu_a/m                               #pattern speed (B4)
            f = lambda r: m*(Omega(J2, Rp, r) - self.Omega_p) - nu(J2, Rp, r)
            self.r_V = brentq(f, 1.0e-3*self.a, self.a, xtol=1e-15, rtol=1e-15)   #inner vertical resonance
            self.beta_V = self.r_V/self.a
            self.b = laplace(self.beta_V, m, 1.5)
            self.A_V_per_mass = self.inclination*self.b/(2.0*self.a**2)           #A_V = mass*this (G = 1)
            hV = h*self.r_V
            d = lambda F: (F(J2, Rp, self.r_V + hV) - F(J2, Rp, self.r_V - hV))/(2*hV)
            self.D_V = 2.0*nu(J2, Rp, self.r_V)*self.r_V*(d(nu) - m*d(Omega))      #D_V (B5)
        else:
            self.inclination, self.node0 = 0.0, 0.0
            self.Omega_p = self.Omega_a                                             #pattern speed (L2)
            f = lambda r: m*(Omega(J2, Rp, r) - self.Omega_p) - kappa(J2, Rp, r)
            r_kep = self.a*((m - 1.0)/m)**(2.0/3.0)                                 #Keplerian r_L (L2), to bracket it
            self.r_L = brentq(f, 0.5*r_kep, self.a, xtol=1e-15, rtol=1e-15)       #inner Lindblad resonance
            self.beta_L = self.r_L/self.a
            self.b = laplace(self.beta_L, m, 0.5)
            self.db = dlaplace(self.beta_L, m, 0.5)
            #Phi = -(mass/a_s)[b_L + b'_L (beta - beta_L)] cos m(theta - lambda_s) near r_L
            self.Phi_per_mass = -self.b/self.a                                      #Phi's amplitude at r_L, per mass
            Om_L = Omega(J2, Rp, self.r_L)
            #Psi = r dPhi/dr + (2 Omega/(Omega - Omega_p)) Phi at r_L (L3), per mass
            self.Psi_per_mass = -(self.beta_L*self.db + 2.0*Om_L/(Om_L - self.Omega_p)*self.b)/self.a
            hL = h*self.r_L
            d = lambda F: (F(J2, Rp, self.r_L + hL) - F(J2, Rp, self.r_L - hL))/(2*hL)
            self.D_L = 2.0*kappa(J2, Rp, self.r_L)*self.r_L*(d(kappa) - m*d(Omega))  #D_L (L4)
        self.time = 0.0

    def mass(self, time):
        """The satellite's mass at time: grown as 1 - exp(-time/time_grow), or mass_final throughout if time_grow <= 0."""
        if self.time_grow <= 0:
            return self.mass_final
        return self.mass_final*(1.0 - np.exp(-time/self.time_grow))

    def longitudes(self, time):
        """Mean longitude lambda_s and node Omega_s at time."""
        return self.lambda0 + self.Omega_a*time, self.node0 + (self.Omega_a - self.nu_a)*time

    def vertical_acceleration(self, theta, time):
        """A_z on particles at longitudes theta (B2)."""
        lam, node = self.longitudes(time)
        A_V = self.mass(time)*self.A_V_per_mass
        return -A_V*np.sin(self.m*theta - (self.m + 1)*lam + node)

    def horizontal_acceleration(self, r, theta, time):
        """A_r, A_theta on particles at r, theta (L1): minus the gradient of
        Phi = -(mass/a_s)[b_L + b'_L (r/a_s - beta_L)] cos m(theta - lambda_s)."""
        lam = self.lambda0 + self.Omega_a*time
        mass = self.mass(time)
        chi = self.m*(theta - lam)
        A_r = (mass*self.db/self.a**2)*np.cos(chi)
        A_t = -(mass*self.m/(self.a*r))*(self.b + self.db*(r/self.a - self.beta_L))*np.sin(chi)
        return A_r, A_t

    def summary(self):
        """The inputs and derived constants, for the log and output/satellite.json."""
        out = dict(self.inputs)
        out.update(resonance=self.resonance, Omega_a=self.Omega_a, nu_a=self.nu_a, Omega_p=self.Omega_p)
        if self.resonance == 'vertical':
            out.update(r_V=self.r_V, beta_V=self.beta_V, b_32=self.b, A_V_per_mass=self.A_V_per_mass,
                       A_V_final=self.mass_final*self.A_V_per_mass, D_V=self.D_V)
        else:
            out.update(r_L=self.r_L, beta_L=self.beta_L, b_12=self.b, db_12=self.db, Psi_per_mass=self.Psi_per_mass,
                       Psi_final=self.mass_final*self.Psi_per_mass, D_L=self.D_L)
        return {k: float(v) if isinstance(v, np.floating) else v for k, v in out.items()}

    def save(self, path):
        with open(path, 'w') as fp:
            json.dump(self.summary(), fp, indent=1)


def get_satellite(g):
    """The Satellite set by inputs.py's dict `satellite`, or None when it is absent or mass_final < 0."""
    sat = g.get('satellite')
    if sat is None or sat.get('mass_final', -1.0) < 0:
        return None
    return Satellite(sat, g['J2'], g['Rp'])
