#overstability_analysis.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 27 November 2026.
#
#Diagnostics of tests/overstability: one axisymmetric epicyclic wave of lattice wavenumber phi seeded in a viscous,
#self-gravitating ring (tools/patches.py seed_mode), which grows or decays at the rate of notes/equations.md O7,
#gamma = -(K/2)[nu_b + (nu_eff - nu_b - 4 x nu_s/3)/(1 - x)], nu_eff the density wave's (O6, beta = 0), here on the
#streamline lattice: K = sin^2(phi)/Delta^2 (O4) and x = X(phi)/kappa^2 (P7), the paper's Appendix B.
#
#Kappa falls outward, so the wave winds up: its wavenumber rises by about 1.5 kappa t/r. The growth rate is measured
#locally, from each streamline's complex amplitude z_j = xi_j - i v_r,j/omega (xi_j and v_r,j its mean radius and radial
#velocity less those of the unseeded control run, which removes the initial ring's own small departures from
#equilibrium; omega the wave's lattice frequency), which for one wave is
#|z_j| = its amplitude whatever its wavenumber, and compared with gamma at the wave's measured wavenumber, the phase
#change of z_j from one streamline to the next.

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
for folder in ('source', 'fast', 'tests'):
    if str(REPO/folder) not in sys.path:
        sys.path.insert(0, str(REPO/folder))
from helper_fns import Omega, Kappa
from vertical_analysis import load3d, ring_sigma, paper_figure

__all__ = ['load3d', 'paper_figure', 'NU', 'PHIS', 'MIXES', 'OVERSTABILITY_RUNS', 'run_name', 'overstability_run',
           'spacing', 'lattice_K', 'lattice_x', 'growth_per_K', 'growth_rate', 'wave', 'measured_growth', 'WINDOW']

NU = 1.0e-11                                     #the shear (or bulk) viscosity of the 800-streamline runs
PHIS = {1: 0.25*np.pi, 2: 0.5*np.pi, 3: 0.75*np.pi}
#viscosity mixes: (nu_s, nu_b, ring gravity on, streamlines, outputs every so many timesteps, decaying?)
MIXES = {'S': (NU, -1.0, True, 800, 25),           #shear only: grows
         'S0': (NU, -1.0, False, 800, 25),         #shear only, no ring gravity: x = 0
         'C': (0.0, NU, True, 800, 25),            #bulk only: decays at K nu_b/2
         'A': (NU, NU, True, 800, 125),            #nu_b = nu_s: decays slowly, 500 orbits
         'D': (NU, 2*NU, True, 800, 25),           #nu_b = 2 nu_s: decays
         'AH': (NU/4, NU/4, True, 1600, 125)}      #nu_b = nu_s, 1600 streamlines: x up to 0.18, so phi >= pi/2 grows
WINDOW = (0.2, 0.8)                              #the streamlines measured, as fractions of the ring's width


def run_name(mix, i):
    return 'output' if (mix, i) == ('S', 2) else f'output_{mix}{i}'


def _mods(mix, i):
    nu_s, nu_b, gravity, n, steps = MIXES[mix]
    m = {'shear_viscosity': nu_s, 'bulk_viscosity': nu_b, 'seed_phi': PHIS[i], 'timesteps_per_output': steps}
    if not gravity:
        m['G_ring'] = -1.0
    if n != 800:
        m.update(number_of_streamlines=n, number_of_particles=32*n)
    if nu_b > 0.0:
        m['seed_amplitude'] = 0.003 if mix == 'AH' else 0.01     #decaying waves start larger, to end above the noise
    return m


OVERSTABILITY_RUNS = {run_name(mix, i): _mods(mix, i) for mix in MIXES for i in PHIS if run_name(mix, i) != 'output'}
#each mix's control run, unseeded: its streamlines' small departures from their circles (the initial ring's departure
#from equilibrium, which the viscosity also amplifies or damps) are subtracted from the seeded runs'
OVERSTABILITY_RUNS.update({f'output_{mix}0': dict(_mods(mix, 2), seed_amplitude=0.0) for mix in MIXES})


def overstability_run(test_folder, mix, i):
    """(p, o) of a run (i = 0: the mix's unseeded control), or None if its output folder is missing."""
    name = f'output_{mix}0' if i == 0 else run_name(mix, i)
    folder = Path(test_folder)/name
    if not folder.exists():
        return None
    return load3d(str(test_folder), str(folder), None if name == 'output' else OVERSTABILITY_RUNS[name])


def spacing(p):
    return p['radial_width']/(p['number_of_streamlines'] - 1)


def lattice_K(p, phi):
    """K = sin^2(phi)/Delta^2, the code's viscous stencils' k^2 (O4)."""
    return np.sin(phi)**2/spacing(p)**2


def lattice_x(p, phi, G=1.0):
    """x = X(phi)/kappa^2 at the ring's middle, X = (4 G sigma/Delta)(pi phi/2 - phi^2/4) the lattice's self-gravity (P7);
    0 without ring gravity."""
    if p['G_ring'] <= 0.0:
        return 0.0
    kappa = Kappa(p['J2'], p['Rp'], 1.0 + 0.5*p['radial_width'])
    X = 4.0*G*ring_sigma(p)/spacing(p)*(np.pi*phi/2.0 - phi**2/4.0)
    return X/kappa**2


def growth_per_K(p, x):
    """O7 per unit K, beta = 0: gamma/K = -(1/2)[nu_b + (nu_eff - nu_b - 4 x nu_s/3)/(1 - x)],
    nu_eff - nu_b = [10/3 - 4 (Omega/kappa)^2] nu_s, at the ring's middle."""
    nu_s = max(p['shear_viscosity'], 0.0)
    nu_b = max(p['bulk_viscosity'], 0.0)
    r = 1.0 + 0.5*p['radial_width']
    ratio = Omega(p['J2'], p['Rp'], r)/Kappa(p['J2'], p['Rp'], r)
    shear = (10.0/3.0 - 4.0*ratio**2)*nu_s
    return -0.5*(nu_b + (shear - 4.0*x*nu_s/3.0)/(1.0 - x))


def growth_rate(p, phi):
    """O7 on the lattice: K growth_per_K(x), K = lattice_K(phi), x = lattice_x(phi)."""
    return lattice_K(p, phi)*growth_per_K(p, lattice_x(p, phi))


def wave(p, o, o0):
    """The wave over the measured streamlines (WINDOW), the seeded run o less its unseeded control o0: times t,
    amplitude A(t) (the median |z_j|, in units of the seeded xi_0) and wavenumber phi(t) (the slope of a straight line
    fitted to the phase of z_j against j)."""
    n = p['number_of_streamlines']
    xi = o['r'].mean(axis=2) - o0['r'].mean(axis=2)
    vr = o['vr'].mean(axis=2) - o0['vr'].mean(axis=2)
    phi0 = p['seed_phi']
    kappa = Kappa(p['J2'], p['Rp'], 1.0 + 0.5*p['radial_width'])
    omega = kappa*np.sqrt(1.0 - lattice_x(p, phi0))
    z = (xi - 1j*vr/omega)[:, int(WINDOW[0]*n):int(WINDOW[1]*n)]
    A = np.median(np.abs(z), axis=1)/(p['seed_amplitude']*spacing(p))
    j = np.arange(z.shape[1])
    phi = np.array([abs(np.polyfit(j, np.unwrap(np.angle(zt)), 1)[0]) for zt in z])
    return o['times'], A, phi


def measured_growth(p, o, o0):
    """The fitted growth rate (ln A linear in t), the wave's mean wavenumber, O7 averaged over the run as the wave
    winds up, K and x at the mean wavenumber, and the rms residual of ln A."""
    t, A, phi = wave(p, o, o0)
    c = np.polyfit(t, np.log(A), 1)
    phi_mean = phi.mean()
    return {'gamma': c[0], 'phi': phi_mean, 'theory': np.mean([growth_rate(p, f) for f in phi]),
            'K': lattice_K(p, phi_mean), 'x': lattice_x(p, phi_mean),
            'residual': np.std(np.log(A) - np.polyval(c, t)), 'growth': A[-1]/A[0],
            'crossed': o['monitor'].get('streamline_crossing_timestep')}
