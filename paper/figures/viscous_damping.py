#viscous_damping.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 27 November 2026.
#
#the paper's viscous density-wave figure (Section 4.2.5): tests/density_waves' ring run with viscosity, runs A
#(nu_s = nu_b), C (bulk viscosity only) and D (nu_b = 2 nu_s). Each streamline's wave amplitude |xi - L_L/x| over the
#inviscid run's at the end of the runs, versus x, over the stretch where the inviscid wave is steady, against linear
#theory's exp[-(x/x_vis)^3] for constant kinematic viscosity (dashed; notes/equations.md O6).
#Needs tests/density_waves' viscous runs: uv run python tools/run_simulations.py tests/density_waves
#
#Drawn by figure(runs), which tests/density_waves/density_waves.ipynb also calls, so the two are identical.
#
#usage, from the repo's top folder:
#   uv run python paper/figures/viscous_damping.py     (writes viscous_damping.pdf and viscous_damping.png)

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
from density_analysis import satellite, damping_ratio, damping_length, resonance_zone, X_STEADY

LABELS = {'A': r'A: $\nu_b = \nu_s$', 'C': r'C: $\nu_s = 0$', 'D': r'D: $\nu_b = 2\nu_s$'}


def figure(runs):
    """The figure, from runs = {name: ((p, o), (p0, o0))}, each viscous run and the inviscid run (viscous_run)."""
    with style(font_scale=1.5) as color_seq:
        f, ax = plt.subplots(figsize=(10, 6))
        x_theory = np.linspace(0.0, X_STEADY, 300)
        for i, (name, ((p, o), (p0, o0))) in enumerate(runs.items()):
            sat = satellite(p)
            x, ratio = damping_ratio(p, o, p0, o0)
            use = (x > resonance_zone(p, sat)) & (x < X_STEADY)
            ax.plot(x[use], ratio[use], marker='o', markersize=4, linestyle='none', alpha=0.7, color=color_seq[i],
                label=LABELS[name])
            ax.plot(x_theory, np.exp(-(x_theory/damping_length(p, sat))**3), linestyle='dashed', linewidth=2,
                color=color_seq[i])
        ax.set_xlim(0.0, X_STEADY)
        ax.set_ylim(0.0, 1.05)
        ax.set_xlabel(r'fractional distance from resonance $x$')
        ax.set_ylabel(r'wave amplitude, viscous/inviscid')
        ax.legend(loc='upper right')
        f.tight_layout()
    return f


if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    from density_analysis import viscous_run
    f = figure({name: viscous_run(REPO/'tests'/'density_waves', name) for name in 'ACD'})
    for ext in ('pdf', 'png'):
        f.savefig(HERE/f'viscous_damping.{ext}', bbox_inches='tight', dpi=130)
        print('wrote', HERE/f'viscous_damping.{ext}')
