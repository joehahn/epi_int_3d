#wavenumber.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 12 November 2026.
#
#the paper's wavenumber figure (Section 4.1.4): the simulated bending wave's dimensionless wavenumber |k| r_V versus
#fractional distance from resonance x at the end of tests/bending_waves, in the same ring simulated with 400, 800 and
#1600 streamlines (resolution_runs), against the continuum dispersion relation of Shu (1984), D|x| = 2 pi G sigma |k|
#(Eqn. 37; notes/equations.md B7). The simulated k is the radial gradient of the phase of the wave part of each
#streamline's complex amplitude h, h - L_V/x: k = d arg(h - L_V/x)/dr, the change in that phase between neighboring
#streamlines divided by their separation (B21).
#Needs that test's outputs: uv run python tools/run_simulations.py tests/bending_waves
#
#Drawn by figure(runs), which tests/bending_waves/bending_waves.ipynb also calls, so the two are identical.
#
#usage, from the repo's top folder:
#   uv run python paper/figures/wavenumber.py     (writes wavenumber.pdf and wavenumber.png)

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
from vertical_analysis import simulated_h, wave_wavenumber, wavenumber_profile, behind_front


def figure(runs):
    """The figure, from the runs at each number of streamlines, [(p, o), ...] (resolution_runs), coarsest first."""
    with style(font_scale=1.5) as color_seq:
        f, ax = plt.subplots(figsize=(10, 6))
        for (p, o), color in zip(runs, color_seq):
            sat, x, h = simulated_h(p, o, -1)
            A_V = sat.mass_final*sat.A_V_per_mass
            x_mid, k_sim = wavenumber_profile(x, h - 1j*A_V/(sat.D_V*x), sat.r_V)    #the wave part, h - L_V/x
            behind, x_front = behind_front(p, o, sat, x_mid)
            ax.plot(x_mid[behind], np.abs(k_sim[behind])*sat.r_V/1.0e4, marker='o', markersize=4, linestyle='none',
                alpha=0.7, color=color, label=f"{p['number_of_streamlines']} streamlines")
            ax.axvline(x_front, color=color, linewidth=1.5, linestyle='dotted')
        x_theory = np.linspace(x.min(), 0.0, 400)
        ax.plot(x_theory, wave_wavenumber(p, sat, x_theory)*sat.r_V/1.0e4, linestyle='dashed', linewidth=2,
            color='black', label='Eqn. (37)')
        ax.set_xlim(x.min(), 0.0)
        ax.set_ylim(0.0, 4.5)
        ax.set_xlabel(r'fractional distance from resonance $x$')
        ax.set_ylabel(r'wavenumber $|k|r_V$    ($10^4$)')
        ax.legend(loc='upper right')
    return f


if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    from vertical_analysis import resolution_runs
    f = figure(resolution_runs(REPO/'tests'/'bending_waves'))
    for ext in ('pdf', 'png'):
        f.savefig(HERE/f'wavenumber.{ext}', bbox_inches='tight', dpi=130)
        print('wrote', HERE/f'wavenumber.{ext}')
