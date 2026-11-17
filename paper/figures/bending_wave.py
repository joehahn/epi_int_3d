#bending_wave.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 10 November 2026.
#
#the paper's bending-wave figure (Section 4.1.4): each streamline's inclination, averaged over its particles,
#versus its fractional distance from resonance x over -0.011 < x < 0.002, at six times, in units of the time
#tau_0.01 that the wave takes to reach x = -0.01, as the wave propagates inward, and linear theory's inclination
#|h|/r, with h the trailing solution of Shu (1984)'s Eqn. (88) (notes/equations.md B8).
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
from vertical_analysis import simulated_h, bending_wave, ring_sigma, tau_001


def figure(p, o):
    """The figure, from the run's parameters p and output o (load3d)."""
    sat, x, h = simulated_h(p, o, -1)
    sigma = ring_sigma(p)
    #the time tau_0.01 the wave takes to travel from r_V to x = -0.01, at the group velocity pi G sigma/nu (B10)
    tau = tau_001(p, sat)
    x_min, x_max = -0.011, 0.002
    with style(font_scale=1.5) as color_seq:
        f, ax = plt.subplots(figsize=(10, 6))
        #each streamline's inclination at six times, in units of tau_0.01
        time_fractions = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5]
        for idx, time_fraction in enumerate(time_fractions):
            tidx = np.argmin(np.abs(o['times'] - time_fraction*tau))
            x_now = ((o['a'][tidx] - sat.r_V)/sat.r_V).mean(axis=1)
            inc_now = o['inc'][tidx].mean(axis=1)
            linewidth = 2 if (idx == len(time_fractions) - 1) else 1
            label = r'$t/\tau_{0.01}$ = ' + str(np.round(o['times'][tidx]/tau, 2))
            ax.plot(x_now, inc_now/1.0e-6, linestyle='-', linewidth=linewidth, color=color_seq[idx], label=label)
        #linear theory, dashed in the color of the latest time
        xs = np.linspace(x_min, x_max, 20000)
        inc_theory = np.abs(bending_wave(sat, sigma, xs))/(sat.r_V*(1.0 + xs))
        ax.plot(xs, inc_theory/1.0e-6, linestyle='dashed', linewidth=2, color=color_seq[len(time_fractions) - 1],
            label='linear theory')
        A_V = sat.mass_final*sat.A_V_per_mass
        Z_w = A_V*np.sqrt(sat.r_V/(sigma*sat.D_V))
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(0.0, 2.0*Z_w/sat.r_V/1.0e-6)
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
