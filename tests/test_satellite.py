#test_satellite.py
#
#Unit tests of fast/satellite.py's 'horizontal' Fourier satellite (notes/equations.md L1-L4): at its inner Lindblad
#resonance r_L its pull A_r, A_theta is the m-th Fourier term of a point-mass satellite's exact in-plane pull (computed
#here by a discrete Fourier transform in longitude), r_L satisfies m(Omega_a - Omega_p) = kappa_a, and its pull is the
#gradient of a potential (the torque it exerts is -dPhi/dtheta).
#
#Run:  uv run pytest tests/test_satellite.py

import sys

import numpy as np
import pytest

from analysis import REPO

sys.path.insert(0, str(REPO/'fast'))
from satellite import Satellite, Omega, kappa

J2, RP = 0.01, 0.5


def horizontal(m, a_s, mass=1.0):
    return Satellite({'resonance': 'horizontal', 'm': m, 'mass_final': mass, 'time_grow': -1.0, 'r': a_s}, J2, RP)


@pytest.mark.parametrize('m, a_s', [(2, 1.6), (3, 1.35), (4, 1.25)])
def test_point_mass_fourier_term(m, a_s):
    """At r_L, the m-th cos/sin coefficients of the exact pull of a unit-mass satellite at (a_s, 0) match A_r, A_theta."""
    sat = horizontal(m, a_s)
    theta = np.linspace(-np.pi, np.pi, 4096, endpoint=False)
    r = sat.r_L
    X, Y = r*np.cos(theta) - a_s, r*np.sin(theta)
    d3 = (X*X + Y*Y)**1.5
    Fx, Fy = -X/d3, -Y/d3
    F_r = Fx*np.cos(theta) + Fy*np.sin(theta)
    F_t = -Fx*np.sin(theta) + Fy*np.cos(theta)
    A_r, A_t = sat.horizontal_acceleration(r, theta, 0.0)
    coef = lambda f, g: 2.0*(f*g).mean()
    assert coef(A_r, np.cos(m*theta)) == pytest.approx(coef(F_r, np.cos(m*theta)), rel=1e-10)
    assert coef(A_t, np.sin(m*theta)) == pytest.approx(coef(F_t, np.sin(m*theta)), rel=1e-10)
    #and nothing at other harmonics
    assert abs(coef(A_r, np.cos((m + 1)*theta))) < 1e-12


def test_resonance_and_pattern_speed():
    sat = horizontal(2, 1.6)
    assert sat.Omega_p == Omega(J2, RP, 1.6)
    assert sat.m*(Omega(J2, RP, sat.r_L) - sat.Omega_p) == pytest.approx(kappa(J2, RP, sat.r_L), rel=1e-13)
    assert sat.r_L == pytest.approx(1.6*(1/2)**(2/3), rel=0.01)                    #Keplerian r_L, to O(J2)
    assert sat.D_L == pytest.approx(3.0*Omega(J2, RP, sat.r_L)**2, rel=0.02)       #Keplerian 3(m-1) Omega_L^2


def test_conservative():
    """A_r = -dPhi/dr and A_theta = -(1/r) dPhi/dtheta of one potential: the curl r A_theta's r-derivative equals
    A_r's theta-derivative, checked by finite differences off r_L."""
    sat = horizontal(2, 1.6, mass=1e-6)
    r, th, h = sat.r_L*1.01, 0.3, 1e-6
    rAt = lambda r, th: r*sat.horizontal_acceleration(r, th, 0.0)[1]
    Ar = lambda r, th: sat.horizontal_acceleration(r, th, 0.0)[0]
    d_rAt_dr = (rAt(r + h, th) - rAt(r - h, th))/(2*h)
    d_Ar_dth = (Ar(r, th + h) - Ar(r, th - h))/(2*h)
    assert d_rAt_dr == pytest.approx(d_Ar_dth, rel=1e-6)
