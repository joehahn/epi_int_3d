#torque.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 12 November 2026.
#
#the paper's torque figure (Section 4.1.3): the torque that the satellite's resonant potential A_V z sin(chi) exerts on
#the simulated ring of tests/bending_waves, -sum mu m A_V z cos(chi) over its particles, versus time, in units of the
#torque of linear theory, T = -pi^2 m sigma r_V^2 A_V^2/D (Shu 1984; notes/equations.md B9), which is also the
#angular momentum luminosity of the wave it launches (tests/bending_waves/bending_waves.ipynb Section 5).
#fast/satellite.py applies only the vertical pull -dPhi/dz, so this is the torque on the simulated wave, not one the
#ring feels.
#Needs that test's output: uv run python tools/run_simulations.py tests/bending_waves
#
#usage, from the repo's top folder:
#   uv run python paper/figures/torque.py     (writes torque.pdf and torque.png)

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
from vertical_analysis import load3d, ring_sigma, satellite_torque, nu_a
from satellite import Satellite

#read the run, and the torque at each output
p, o = load3d(os.path.join(repo_folder, 'tests', 'bending_waves'))
sat = Satellite(p['satellite'], p['J2'], p['Rp'])
sigma = ring_sigma(p)
A_V = sat.mass_final*sat.A_V_per_mass
torque_theory = -np.pi**2*sat.m*sigma*sat.r_V**2*A_V**2/sat.D_V
torque = np.array([satellite_torque(p, o, k).sum() for k in range(len(o['times']))])

#time in units of the wave's crossing time from r_V to the inner edge, at c_g = pi G sigma/nu
crossing_time = (sat.r_V - 1.0)/(np.pi*sigma/nu_a(p['J2'], p['Rp'], sat.r_V))
half = len(torque)//2
print('torque/theory over the last half of the run: mean', np.round(torque[half:].mean()/torque_theory, 4),
    ' range', np.round(torque[half:].min()/torque_theory, 4), np.round(torque[half:].max()/torque_theory, 4))

#plot
sns.set(font_scale=1.5, font='DejaVu Sans')
color_seq = plt.rcParams['axes.prop_cycle'].by_key()['color']
f, ax = plt.subplots(figsize=(10, 6))
ax.plot(o['times']/crossing_time, torque/torque_theory, marker='o', markersize=7, linestyle='-', linewidth=2,
    color=color_seq[0], label='simulation')
ax.plot(o['times']/crossing_time, o['times']*0 + 1.0, linestyle='dashed', linewidth=2, color=color_seq[0],
    label='linear theory')
ax.set_xlim(0.0, o['times'][-1]/crossing_time*1.02)
ax.set_ylim(0.0, 1.3)
ax.set_xlabel(r'time $t/\tau_c$')
ax.set_ylabel(r'torque $T/T_V$')
ax.legend(loc='lower right')

#save
for extension in ['pdf', 'png']:
    path = os.path.join(figures_folder, 'torque.' + extension)
    plt.savefig(path, bbox_inches='tight', dpi=130)
    print('wrote', os.path.relpath(path, repo_folder))
