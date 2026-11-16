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
#Drawn by figure(p, o), which tests/bending_waves/bending_waves.ipynb also calls, so the two are identical.
#
#usage, from the repo's top folder:
#   uv run python paper/figures/bending_wave.py     (writes bending_wave.pdf and bending_wave.png)

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
for folder in ('source', 'fast', 'tests'):
    if str(REPO/folder) not in sys.path:
        sys.path.insert(0, str(REPO/folder))
sys.path.insert(0, str(HERE))
from figstyle import style
from vertical_analysis import simulated_h, bending_wave, ring_sigma, nu_a


def figure(p, o):
    """The figure, from the run's parameters p and output o (load3d)."""
    sat, x, h = simulated_h(p, o, -1)
    sigma = ring_sigma(p)
    #the wave's crossing time from r_V to the inner edge, at the group velocity pi G sigma/nu (B10)
    nu_V = nu_a(p['J2'], p['Rp'], sat.r_V)
    crossing_time = (sat.r_V - 1.0)/(np.pi*sigma/nu_V)
    with style(font_scale=1.5) as color_seq:
        f, ax = plt.subplots(figsize=(10, 6))
        #each streamline's inclination at six times, in units of the crossing time
        time_fractions = [0.05, 0.15, 0.25, 0.4, 0.6, 0.85]
        for idx, time_fraction in enumerate(time_fractions):
            tidx = np.argmin(np.abs(o['times'] - time_fraction*crossing_time))
            x_now = ((o['a'][tidx] - sat.r_V)/sat.r_V).mean(axis=1)
            inc_now = o['inc'][tidx].mean(axis=1)
            linewidth = 2 if (idx == len(time_fractions) - 1) else 1
            label = r'$t/\tau_c$ = ' + str(np.round(o['times'][tidx]/crossing_time, 2))
            ax.plot(x_now, inc_now/1.0e-6, linestyle='-', linewidth=linewidth, color=color_seq[idx], label=label)
        #linear theory, dashed in the color of the latest time
        xs = np.linspace(x.min(), x.max(), 20000)
        inc_theory = np.abs(bending_wave(sat, sigma, xs))/(sat.r_V*(1.0 + xs))
        ax.plot(xs, inc_theory/1.0e-6, linestyle='dashed', linewidth=2, color=color_seq[len(time_fractions) - 1],
            label='linear theory')
        A_V = sat.mass_final*sat.A_V_per_mass
        Z_w = A_V*np.sqrt(sat.r_V/(sigma*sat.D_V))
        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(0.0, 1.9*Z_w/sat.r_V/1.0e-6)
        ax.set_xlabel(r'fractional distance from resonance $x$')
        ax.set_ylabel(r'inclination $i$    ($10^{-6}$)')
        ax.legend(loc='upper right', ncol=2)
    return f


if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    from vertical_analysis import load3d
    p, o = load3d(REPO/'tests'/'bending_waves')
    f = figure(p, o)
    for ext in ('pdf', 'png'):
        f.savefig(HERE/f'bending_wave.{ext}', bbox_inches='tight', dpi=130)
        print('wrote', HERE/f'bending_wave.{ext}')
