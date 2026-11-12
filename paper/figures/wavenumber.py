#wavenumber.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 12 November 2026.
#
#the paper's wavenumber figure (Section 4.1.3), after Fig. 5 of Hahn (2007): the simulated bending wave's
#dimensionless wavenumber |k| r_V versus fractional distance from resonance x at the end of tests/bending_waves,
#against the continuum dispersion relation of Shu (1984), D|x| = 2 pi G sigma |k| (notes/equations.md B7), and that of
#the finite lattice of streamlines with the ring's edges, B13. The simulated k is the radial gradient of the phase of
#the wave part of each streamline's complex amplitude h, h - L_V/x (tests/bending_waves/bending_waves.ipynb Section 4).
#Needs that test's output: uv run python tools/run_simulations.py tests/bending_waves
#
#usage, from the repo's top folder:
#   uv run python paper/figures/wavenumber.py     (writes wavenumber.pdf and wavenumber.png)

import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

#paths
figures_folder = os.path.dirname(os.path.abspath(__file__))
repo_folder = os.path.dirname(os.path.dirname(figures_folder))
sys.path.insert(0, os.path.join(repo_folder, 'source'))
sys.path.insert(0, os.path.join(repo_folder, 'fast'))
sys.path.insert(0, os.path.join(repo_folder, 'tests'))
from vertical_analysis import load3d, simulated_h, ring_sigma, wave_wavenumber, nu_a

#read the run, and the wave's local wavenumber from the phase of its wave part, at the end
p, o = load3d(os.path.join(repo_folder, 'tests', 'bending_waves'))
sat, x, h = simulated_h(p, o, -1)
sigma = ring_sigma(p)
A_V = sat.mass_final*sat.A_V_per_mass
forced = 1j*A_V/(sat.D_V*x)                                         #L_V/x
x_mid = 0.5*(x[1:] + x[:-1])
k_sim = np.angle((h - forced)[1:]/(h - forced)[:-1])/np.diff(x)/sat.r_V

#where the wave is in place at the end: behind its front, which travels inward at c_g = pi G sigma/nu
nu_V = nu_a(p['J2'], p['Rp'], sat.r_V)
x_front = -np.pi*sigma*o['times'][-1]/(nu_V*sat.r_V)
behind = (x_mid < -2.0*np.sqrt(2*np.pi*sigma/(sat.r_V*sat.D_V))) & (x_mid > x_front)

#theory
x_theory = np.linspace(x.min(), -1.0e-5, 400)
k_continuum, k_lattice, k_edges = wave_wavenumber(p, sat, x_theory)

#plot prep
sns.set(font_scale=1.5, font='DejaVu Sans')
color_seq = plt.rcParams['axes.prop_cycle'].by_key()['color']
f, ax = plt.subplots(figsize=(10, 6))
ax.plot(x_mid[behind]/1.0e-3, np.abs(k_sim[behind])*sat.r_V/1.0e4, marker='o', markersize=6, linestyle='none',
    alpha=0.7, color=color_seq[0], label='simulation')
ax.plot(x_theory/1.0e-3, k_continuum*sat.r_V/1.0e4, linestyle='dashed', linewidth=2, color=color_seq[0],
    label='continuum')
ax.plot(x_theory/1.0e-3, k_edges*sat.r_V/1.0e4, linestyle='dashed', linewidth=2, color=color_seq[1],
    label='lattice + edges')
ax.axvline(x_front/1.0e-3, color='gray', linewidth=1, linestyle='dotted')

#axes
ax.set_xlim(x.min()/1.0e-3, 0.0)
ax.set_ylim(0.0, 6.0)
ax.set_xlabel(r'fractional distance from resonance $x$    ($10^{-3}$)')
ax.set_ylabel(r'wavenumber $|k|r_V$    ($10^4$)')
ax.legend(loc='upper left')

#save
for extension in ['pdf', 'png']:
    path = os.path.join(figures_folder, 'wavenumber.' + extension)
    plt.savefig(path, bbox_inches='tight', dpi=130)
    print('wrote', os.path.relpath(path, repo_folder))
