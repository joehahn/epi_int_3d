#vertical_resonance.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com; one figure with the notebook's, 16 November 2026.
#
#The paper's vertical-resonance figure (Section 4.1.3), which tests/vertical_resonance/vertical_resonance.ipynb draws
#with this same function: for the zero-mass ring of tests/vertical_resonance, with each streamline's semimajor axis a
#and inclination i averaged over its particles and the orbits after 10^4,
#  (a) D = nu_a^2 - omega_m^2 (Eqn. 25) at each streamline's a, in units of Omega_a^2 at r_V, versus x = (a - r_V)/r_V,
#      the curve D at all x, and the two streamlines with the largest <i> circled;
#  (b) <i> versus x, and the forced inclination i_f of Eqn. (38).
#Needs that test's output: uv run python tools/run_simulations.py tests/vertical_resonance
#
#usage, from the repo's top folder:
#   uv run python paper/figures/vertical_resonance.py     (writes vertical_resonance.pdf and .png)

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
for folder in ('source', 'fast', 'tests'):
    if str(REPO/folder) not in sys.path:
        sys.path.insert(0, str(REPO/folder))
from vertical_analysis import resonance_D
from helper_fns import Omega
sys.path.insert(0, str(HERE))
from figstyle import style

LATE_ORBITS = 1.0e4                     #averages over the orbits after this


def figure(p, o):
    """The figure, from the run's parameters p and output o (load3d)."""
    orbits = o['times']/(2*np.pi)
    late = orbits > LATE_ORBITS
    a_avg = o['a'][late].mean(axis=(0, 2))
    i_avg = o['inc'][late].mean(axis=(0, 2))
    sat, D = resonance_D(p, a_avg)
    x = (a_avg - sat.r_V)/sat.r_V
    xs = np.linspace(x.min(), x.max(), 4000)                    #an even number of points skips x = 0
    D_curve = resonance_D(p, sat.r_V*(1 + xs))[1]
    Omega_V2 = Omega(p['J2'], p['Rp'], sat.r_V)**2
    i_f = sat.mass_final*sat.inclination*sat.beta_V**2*sat.b/(6*(sat.m - 1)*np.abs(xs))
    peak = np.argsort(i_avg)[-2:]                               #the two streamlines with the largest <i>
    with style(font_scale=1.0):
        f, (ax, ax2) = plt.subplots(2, 1, figsize=(6, 7.5))
        ax.axhline(0.0, color='gray', linewidth=1)
        ax.axvline(0.0, color='gray', linewidth=1)
        ax.plot(xs, D_curve/Omega_V2, linewidth=1, color='black', label=r'$\nu_a^2 - \omega_m^2$')
        ax.plot(x, D/Omega_V2, marker='o', markersize=2.6, linestyle='none', color='gray', label='streamlines',
            zorder=3)
        ax.plot(x[peak], D[peak]/Omega_V2, marker='o', markersize=6.6, linestyle='none', markerfacecolor='none',
            markeredgecolor='black', label=r'largest $\langle i\rangle$', zorder=4)
        ax.set_xlabel(r'fractional distance from resonance $x$')
        ax.set_ylabel(r'frequency distance from resonance $D$' + '\n' + r'($\Omega_a^2$ at $r_V$)')
        ax.legend(loc='upper left', fontsize=9)
        ax2.axvline(0.0, color='gray', linewidth=1)
        ax2.plot(x, i_avg, marker='o', markersize=2.6, linestyle='none', color='gray', label='simulation', zorder=3)
        ax2.plot(xs, i_f, linewidth=1, color='black', label=r'$i_f$, Eqn. (38)')
        ax2.set_yscale('log')
        ax2.set_ylim(1.0e-8, 1.0e-5)
        ax2.set_xlabel(r'fractional distance from resonance $x$')
        ax2.set_ylabel(r'inclination $\langle i\rangle$')
        ax2.legend(loc='upper left', fontsize=9)
        for a_, label in ((ax, '(a)'), (ax2, '(b)')):
            a_.set_xticks([-0.02, -0.01, 0.0, 0.01, 0.02])
            a_.text(0.97, 0.93, label, transform=a_.transAxes, ha='right', va='top', fontsize=11)
        f.tight_layout(h_pad=1.5)
    return f


if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    from vertical_analysis import load3d
    p, o = load3d(REPO/'tests'/'vertical_resonance')
    f = figure(p, o)
    for ext in ('pdf', 'png'):
        f.savefig(HERE/f'vertical_resonance.{ext}', dpi=130)
        print('wrote', HERE/f'vertical_resonance.{ext}')
