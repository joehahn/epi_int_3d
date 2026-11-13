#resonance_location.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 12 November 2026.
#
#the paper's resonance-location figure (Section 4.1.2): for the zero-mass ring of tests/vertical_resonance, each
#streamline's nu_a^2 - omega_m^2 inferred from its forced vertical amplitude, versus its fractional distance from
#resonance x, and the exact nu_a^2 - omega_m^2 at its semimajor axis, which vanishes at r_V. A particle forced by
#A_z = -A_V sin(chi) oscillates as z = Re[h e^{i chi}] with Im(h) = A_V/(nu_a^2 - omega_m^2) (notes/equations.md B6),
#so nu_a^2 - omega_m^2 = A_V/Im(h), with h each streamline's least-squares amplitude averaged over the last 10^4
#orbits, which averages away what free oscillation remains.
#Needs that test's output: uv run python tools/run_simulations.py tests/vertical_resonance
#
#usage, from the repo's top folder:
#   uv run python paper/figures/resonance_location.py     (writes resonance_location.pdf and .png)

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
from vertical_analysis import load3d, resonance_frequencies

#read the run, and each streamline's nu_a^2 - omega_m^2, simulated and exact, in units of D_V
p, o = load3d(os.path.join(repo_folder, 'tests', 'vertical_resonance'))
sat, x, simulated, exact = resonance_frequencies(p, o)
off = np.abs(x) > 1.0e-6                     #the streamline at r_V has no steady forced amplitude
print('largest |simulated/exact - 1| =', np.round(np.abs(simulated[off]/exact[off] - 1.0).max(), 4))

#plot prep
sns.set(font_scale=1.5, font='DejaVu Sans')
color_seq = plt.rcParams['axes.prop_cycle'].by_key()['color']
f, ax = plt.subplots(figsize=(10, 6))

#simulation, exact, and the linear approximation D_V x
ax.plot(x[off]/1.0e-2, simulated[off]/sat.D_V/1.0e-2, marker='o', markersize=6, linestyle='none', alpha=0.7,
    color=color_seq[0], label='simulation')
a_s = np.linspace(x.min(), x.max(), 400)
_, _, _, exact_curve = resonance_frequencies(p, o, x=a_s)
ax.plot(a_s/1.0e-2, exact_curve/sat.D_V/1.0e-2, linestyle='dashed', linewidth=2, color=color_seq[0], zorder=3,
    label='theory')
ax.plot(a_s/1.0e-2, a_s/1.0e-2, linestyle='dotted', linewidth=2, color=color_seq[1], label=r'${\cal D}x$')
ax.axhline(0.0, color='gray', linewidth=1)
ax.axvline(0.0, color='gray', linewidth=1)

#axes
ax.set_xlim(x.min()/1.0e-2*1.05, x.max()/1.0e-2*1.05)
ax.set_xlabel(r'fractional distance from resonance $x$    ($10^{-2}$)')
ax.set_ylabel(r'$D/{\cal D}$    ($10^{-2}$)')
ax.legend(loc='upper left')

#save
for extension in ['pdf', 'png']:
    path = os.path.join(figures_folder, 'resonance_location.' + extension)
    plt.savefig(path, bbox_inches='tight', dpi=130)
    print('wrote', os.path.relpath(path, repo_folder))
