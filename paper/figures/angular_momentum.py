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
#Drawn by figure(p, o), which tests/bending_waves/bending_waves.ipynb also calls, so the two are identical.
#
#usage, from the repo's top folder:
#   uv run python paper/figures/angular_momentum.py     (writes angular_momentum.pdf and .png)

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
from vertical_analysis import ring_angular_momentum


def figure(p, o):
    """The figure, from the run's parameters p and output o (load3d)."""
    L = ring_angular_momentum(p, o)
    dL = (L - L[0])/L[0]
    orbits = o['times']/(2.0*np.pi)
    with style(font_scale=1.5) as color_seq:
        f, ax = plt.subplots(figsize=(10, 6))
        ax.plot(orbits/1.0e3, dL/1.0e-13, linestyle='-', linewidth=2, color=color_seq[0])
        ax.set_xlim(0.0, orbits[-1]/1.0e3)
        ax.set_xlabel(r'time $t$    ($10^3$ orbits)')
        ax.set_ylabel(r'$\Delta L_z/L_z$    ($10^{-13}$)')
    return f


if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    from vertical_analysis import load3d
    p, o = load3d(REPO/'tests'/'bending_waves')
    f = figure(p, o)
    for ext in ('pdf', 'png'):
        f.savefig(HERE/f'angular_momentum.{ext}', bbox_inches='tight', dpi=130)
        print('wrote', HERE/f'angular_momentum.{ext}')
