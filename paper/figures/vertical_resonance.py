#vertical_resonance.py
#
#The paper's vertical-resonance figure (Section 4.2, "Massless particles"): each streamline's
#inclination, averaged over its particles and over the last 10^4 orbits of tests/vertical_resonance, versus its
#average fractional distance from resonance, and the forced inclination of Shu (1984), notes/equations.md B12.
#Needs that test's output: uv run python tools/run_simulations.py tests/vertical_resonance
#
#Run:  uv run python paper/figures/vertical_resonance.py   (writes vertical_resonance .pdf and .png)

import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO/'tests'))
sys.path.insert(0, str(REPO/'fast'))
from vertical_analysis import load3d
import satellite as satellite_module

p, o = load3d(REPO/'tests'/'vertical_resonance')
sat = satellite_module.Satellite(p['satellite'], p['J2'], p['Rp'])
orbits = o['times']/(2*np.pi)
late = orbits > 1.0e4
i_avg = o['inc'][late].mean(axis=(0, 2))
x_avg = ((o['a'][late] - sat.r_V)/sat.r_V).mean(axis=(0, 2))
xs = np.linspace(-0.0205, 0.0205, 4000)                                 #even: skips x = 0
i_f = sat.mass_final*sat.inclination*sat.beta_V**2*sat.b/(6*(sat.m - 1)*np.abs(xs))

fig, ax = plt.subplots(figsize=(4.6, 3.6))
ax.plot(xs, i_f, '-', color='C0', lw=1, label='forced inclination')
ax.plot(x_avg, i_avg, 'k.', ms=4, label='simulation')
ax.set_yscale('log')
ax.set_ylim(1e-8, 1e-5)
ax.set_xlim(-0.021, 0.021)
ax.set_xlabel(r'distance from resonance $\langle x\rangle$')
ax.set_ylabel(r'inclination $\langle i\rangle$')
ax.legend(fontsize=9, loc='upper left', frameon=False)
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(HERE/f'vertical_resonance.{ext}', dpi=130)
    print('wrote', HERE/f'vertical_resonance.{ext}')
ok = np.abs(x_avg) > 1e-6
ratio = i_avg[ok]/(sat.mass_final*sat.inclination*sat.beta_V**2*sat.b/(6*(sat.m - 1)*np.abs(x_avg[ok])))
print(f'simulation/B12: {ratio.min():.3f} to {ratio.max():.3f}')
