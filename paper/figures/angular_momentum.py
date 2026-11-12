#angular_momentum.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 12 November 2026.
#
#the paper's angular-momentum figure (Section 4.1.3): the fractional change in the total angular momentum L_z of the
#self-gravitating ring of tests/bending_waves over the run. The satellite pulls only along z and so exerts no torque
#about z, so L_z changes only by the ring's own gravity, if it does not conserve angular momentum, and by the
#integrator's round-off (tests/vertical_analysis.py ring_angular_momentum).
#Needs that test's output: uv run python tools/run_simulations.py tests/bending_waves
#
#usage, from the repo's top folder:
#   uv run python paper/figures/angular_momentum.py     (writes angular_momentum.pdf and .png)

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
from vertical_analysis import load3d, ring_angular_momentum

#read the run, and L_z at each output
p, o = load3d(os.path.join(repo_folder, 'tests', 'bending_waves'))
L = ring_angular_momentum(p, o)
dL = (L - L[0])/L[0]
orbits = o['times']/(2.0*np.pi)
print('final dL/L =', dL[-1], ' largest |dL/L| =', np.abs(dL).max())

#plot
sns.set(font_scale=1.5, font='DejaVu Sans')
color_seq = plt.rcParams['axes.prop_cycle'].by_key()['color']
f, ax = plt.subplots(figsize=(10, 6))
ax.plot(orbits/1.0e3, dL/1.0e-13, linestyle='-', linewidth=2, color=color_seq[0])
ax.set_xlim(0.0, orbits[-1]/1.0e3)
ax.set_xlabel(r'time $t$    ($10^3$ orbits)')
ax.set_ylabel(r'$\Delta L_z/L_z$    ($10^{-13}$)')

#save
for extension in ['pdf', 'png']:
    path = os.path.join(figures_folder, 'angular_momentum.' + extension)
    plt.savefig(path, bbox_inches='tight', dpi=130)
    print('wrote', os.path.relpath(path, repo_folder))
