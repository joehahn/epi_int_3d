#bending_wave.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 10 November 2026.
#
#the paper's bending-wave figure (Section 4.3): each streamline's inclination, averaged over its particles,
#versus its fractional distance from resonance x, at six times as the wave propagates inward, and linear
#theory's inclination |h|/r, with h the trailing solution of Shu (1984)'s Eqn. (88) (notes/equations.md B8).
#Drawn from tests/bending_waves/bending_waves.ipynb Section 1.
#Needs that test's output: uv run python tools/run_simulations.py tests/bending_waves
#
#usage, from the repo's top folder:
#   uv run python paper/figures/bending_wave.py     (writes bending_wave.pdf and bending_wave.png)

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
from vertical_analysis import load3d, simulated_h, bending_wave, ring_sigma, nu_a

#read the run
p, o = load3d(os.path.join(repo_folder, 'tests', 'bending_waves'))
sat, x, h = simulated_h(p, o, -1)
sigma = ring_sigma(p)
nu_V = nu_a(p['J2'], p['Rp'], sat.r_V)

#the wave's crossing time from r_V to the inner edge, at the group velocity pi G sigma/nu (B10)
c_g = np.pi*sigma/nu_V
crossing_time = (sat.r_V - 1.0)/c_g
print('crossing_time (orbits) =', np.round(crossing_time/(2.0*np.pi), 1))

#plot prep
sns.set(font_scale=1.5, font='DejaVu Sans')
color_seq = plt.rcParams['axes.prop_cycle'].by_key()['color']
f, ax = plt.subplots(figsize=(10, 6))

#each streamline's inclination at six times, in units of the crossing time
time_fractions = [0.05, 0.15, 0.25, 0.4, 0.6, 0.85]
for idx, time_fraction in enumerate(time_fractions):
    tidx = np.argmin(np.abs(o['times'] - time_fraction*crossing_time))
    x_now = ((o['a'][tidx] - sat.r_V)/sat.r_V).mean(axis=1)
    inc_now = o['inc'][tidx].mean(axis=1)
    linewidth = 1
    if (idx == len(time_fractions) - 1):
        linewidth = 2
    label = r'$t/\tau_c$ = ' + str(np.round(o['times'][tidx]/crossing_time, 2))
    ax.plot(x_now, inc_now/1.0e-6, linestyle='-', linewidth=linewidth, color=color_seq[idx], label=label)

#linear theory, dashed in the color of the latest time
xs = np.linspace(x.min(), x.max(), 20000)
inc_theory = np.abs(bending_wave(sat, sigma, xs))/(sat.r_V*(1.0 + xs))
ax.plot(xs, inc_theory/1.0e-6, linestyle='dashed', linewidth=2, color=color_seq[len(time_fractions) - 1],
    label='linear theory')

#axes
A_V = sat.mass_final*sat.A_V_per_mass
Z_w = A_V*np.sqrt(sat.r_V/(sigma*sat.D_V))
ax.set_xlim(x.min(), x.max())
ax.set_ylim(0.0, 1.9*Z_w/sat.r_V/1.0e-6)
ax.set_xlabel(r'fractional distance from resonance $x$')
ax.set_ylabel(r'inclination $i$    ($10^{-6}$)')
ax.legend(loc='upper right', ncol=2)

#save
for extension in ['pdf', 'png']:
    path = os.path.join(figures_folder, 'bending_wave.' + extension)
    plt.savefig(path, bbox_inches='tight', dpi=130)
    print('wrote', os.path.relpath(path, repo_folder))
