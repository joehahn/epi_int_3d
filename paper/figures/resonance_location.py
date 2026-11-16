#resonance_location.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 12 November 2026; D from its definition, 16 November 2026.
#
#The paper's resonance-location figure (Section 4.1.3): for the zero-mass ring of tests/vertical_resonance,
#D = nu_a^2 - omega_m^2 (paper Eqn. 25) at each streamline's semimajor axis, averaged over its particles and the orbits
#after 10^4, versus its fractional distance from resonance x, with D at all x as a curve. The two streamlines whose
#inclinations, averaged the same way, are largest are marked: they straddle D = 0 (style/testing.md T25: D is computed
#from its definition, not inferred from the forced amplitude).
#Needs that test's output: uv run python tools/run_simulations.py tests/vertical_resonance
#
#usage, from the repo's top folder:
#   uv run python paper/figures/resonance_location.py     (writes resonance_location.pdf and .png)

import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO/'source'))
sys.path.insert(0, str(REPO/'fast'))
sys.path.insert(0, str(REPO/'tests'))
from vertical_analysis import load3d, resonance_D
from helper_fns import Omega

#each streamline's semimajor axis and inclination, averaged over its particles and the orbits after 10^4
p, o = load3d(REPO/'tests'/'vertical_resonance')
orbits = o['times']/(2*np.pi)
late = orbits > 1.0e4
a_avg = o['a'][late].mean(axis=(0, 2))
i_avg = o['inc'][late].mean(axis=(0, 2))
sat, D = resonance_D(p, a_avg)
x = (a_avg - sat.r_V)/sat.r_V
xs = np.linspace(x.min(), x.max(), 4001)
D_curve = resonance_D(p, sat.r_V*(1 + xs))[1]
Omega_V2 = Omega(p['J2'], p['Rp'], sat.r_V)**2
peak = np.argsort(i_avg)[-2:]                       #the two streamlines with the largest inclinations
j = np.where(np.diff(np.sign(D)) != 0)[0][0]
x_zero = x[j] - D[j]*(x[j + 1] - x[j])/(D[j + 1] - D[j])
print('D = 0 at x =', x_zero, 'between the streamlines at x =', x[j], x[j + 1])
print('largest inclinations at x =', np.sort(x[peak]))

def panel(ax, ms):
    ax.axhline(0.0, color='gray', lw=0.6)
    ax.axvline(0.0, color='gray', lw=0.6)
    ax.plot(xs, D_curve/Omega_V2, '-', color='C0', lw=1, label=r'$\nu_a^2 - \omega_m^2$')
    ax.plot(x, D/Omega_V2, 'k.', ms=ms, label='streamlines')
    ax.plot(x[peak], D[peak]/Omega_V2, 'o', mfc='none', mec='C3', ms=ms + 3, label=r'largest $\langle i\rangle$')

fig, ax = plt.subplots(figsize=(4.6, 3.6))
panel(ax, 4)
ax.set_xlim(-0.021, 0.021)
ax.set_xlabel(r'distance from resonance $x$')
ax.set_ylabel(r'$D/\Omega_a^2$ at $r_V$')
ax.legend(fontsize=8, loc='upper left', frameon=False)
#inset: the streamlines either side of r_V
w = 1.6e-3
ins = ax.inset_axes([0.6, 0.12, 0.37, 0.36])
panel(ins, 4)
ins.set_xlim(-w, w)
k = np.abs(xs) <= w
ins.set_ylim(1.1*D_curve[k].min()/Omega_V2, 1.1*D_curve[k].max()/Omega_V2)
ins.tick_params(labelsize=6)
ins.set_xticks([-1.0e-3, 0.0, 1.0e-3], ['-0.001', '0', '0.001'])
ins.set_yticks([-0.004, 0.0, 0.004], ['-0.004', '0', '0.004'])
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(HERE/f'resonance_location.{ext}', dpi=130)
    print('wrote', HERE/f'resonance_location.{ext}')
