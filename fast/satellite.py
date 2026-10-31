"""A "Fourier satellite" for epi_int_3d: only the resonant m-th term of a satellite's potential.

The satellite is set by the dict `satellite` in inputs.py (notes/phase2c/bending_wave_theory.md, section 4):
    satellite = {
        'resonance': 'vertical',   #'vertical' (inner vertical resonance) or 'horizontal' (inner Lindblad)
        'm': 2,                    #number of spiral arms, m >= 2 for 'vertical'
        'mass_final': 1.0e-9,      #satellite mass/planet mass once grown; < 0 turns the satellite off
        'time_grow': 1.0e3,        #mass = mass_final*(1 - exp(-time/time_grow))
        'r': 2.0911,               #satellite's orbit radius a_s (circular)
        't': 0.0,                  #its initial mean longitude lambda_s
        'inclination': 1.0e-3,     #i_s, for 'vertical'
        'node': 0.0,               #its initial longitude of ascending node Omega_s
    }
The satellite's mean longitude advances at Omega_{a_s}, its node at Omega_{a_s} - nu_{a_s} (notes/equations.md B4).

'vertical': an inclined circular satellite pulls a ring particle at longitude theta with
    A_z = -A_V sin(m theta - (m+1) lambda_s + Omega_s),   A_V = G M_s i_s b_{3/2}^(m)(beta_V)/(2 a_s^2)   (B2)
and nothing in the plane; b is evaluated once, at the inner vertical resonance r_V, beta_V = r_V/a_s, where
m(Omega_a - Omega_p) = nu_a with m Omega_p = m Omega_{a_s} + nu_{a_s} (B4), solved with J2. 'horizontal' (the
lindblad_resonance branch of epi_int_lite) is not ported yet.
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


def laplace(beta, m, s):
    """Laplace coefficient b_s^(m)(beta) = (2/pi) int_0^pi cos(m phi)/(1 + beta^2 - 2 beta cos phi)^s dphi."""
    f = lambda phi: np.cos(m*phi)/(1.0 + beta**2 - 2.0*beta*np.cos(phi))**s
    return (2.0/np.pi)*quad(f, 0.0, np.pi, limit=200, epsabs=0.0, epsrel=1e-13)[0]


class Satellite:
    """The satellite's inputs, derived constants and its state at the current time."""

    def __init__(self, sat, J2, Rp):
        self.inputs = dict(sat)
        self.resonance = sat.get('resonance', 'vertical')
        if self.resonance != 'vertical':
            raise NotImplementedError("satellite resonance '%s' is not ported yet: use 'vertical'" % self.resonance)
        self.m = int(sat['m'])
        if self.m < 2:
            raise ValueError('a vertical resonance needs m >= 2')
        self.mass_final = sat['mass_final']
        self.time_grow = sat['time_grow']
        self.a = sat['r']
        self.lambda0 = sat.get('t', 0.0)
        self.inclination = sat['inclination']
        self.node0 = sat.get('node', 0.0)
        self.Omega_a = Omega(J2, Rp, self.a)
        self.nu_a = nu(J2, Rp, self.a)
        m = self.m
        self.Omega_p = self.Omega_a + self.nu_a/m                                   #pattern speed (B4)
        f = lambda r: m*(Omega(J2, Rp, r) - self.Omega_p) - nu(J2, Rp, r)
        self.r_V = brentq(f, 1.0e-3*self.a, self.a, xtol=1e-15, rtol=1e-15)     #inner vertical resonance
        self.beta_V = self.r_V/self.a
        self.b = laplace(self.beta_V, m, 1.5)
        self.A_V_per_mass = self.inclination*self.b/(2.0*self.a**2)               #A_V = mass*this (G = 1)
        h = 1.0e-6*self.r_V                                                         #D_V (B5), numerically
        d = lambda F: (F(J2, Rp, self.r_V + h) - F(J2, Rp, self.r_V - h))/(2*h)
        self.D_V = 2.0*nu(J2, Rp, self.r_V)*self.r_V*(d(nu) - m*d(Omega))
        self.time = 0.0

    def mass(self, time):
        return self.mass_final*(1.0 - np.exp(-time/self.time_grow))

    def longitudes(self, time):
        """Mean longitude lambda_s and node Omega_s at time."""
        return self.lambda0 + self.Omega_a*time, self.node0 + (self.Omega_a - self.nu_a)*time

    def vertical_acceleration(self, theta, time):
        """A_z on particles at longitudes theta (B2)."""
        lam, node = self.longitudes(time)
        A_V = self.mass(time)*self.A_V_per_mass
        return -A_V*np.sin(self.m*theta - (self.m + 1)*lam + node)

    def summary(self):
        """The inputs and derived constants, for the log and output/satellite.json."""
        out = dict(self.inputs)
        out.update(resonance=self.resonance, Omega_a=self.Omega_a, nu_a=self.nu_a, Omega_p=self.Omega_p,
                   r_V=self.r_V, beta_V=self.beta_V, b_32=self.b, A_V_per_mass=self.A_V_per_mass,
                   A_V_final=self.mass_final*self.A_V_per_mass, D_V=self.D_V)
        return out

    def save(self, path):
        with open(path, 'w') as fp:
            json.dump(self.summary(), fp, indent=1)


def get_satellite(g):
    """The Satellite set by inputs.py's dict `satellite`, or None when it is absent or mass_final < 0."""
    sat = g.get('satellite')
    if sat is None or sat.get('mass_final', -1.0) < 0:
        return None
    return Satellite(sat, g['J2'], g['Rp'])
