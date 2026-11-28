#overstability.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 27 November 2026.
#
#the paper's viscous-overstability figure (Appendix B): tests/overstability's seeded axisymmetric epicyclic waves. Each
#run's measured growth rate gamma in units of K nu (K = sin^2 phi/Delta^2, nu = nu_s, or nu_b for bulk viscosity only),
#versus the self-gravity x = X(phi)/kappa^2 at the wave's measured wavenumber, against notes/equations.md O7 (lines):
#gamma/K nu_s = -(1/2)[nu_b/nu_s + (nu_eff - nu_b - 4 x nu_s/3)/(nu_s (1 - x))]. Above zero the wave grows (overstable).
#Needs tests/overstability's runs: uv run python tools/run_simulations.py tests/overstability
#
#Drawn by figure(runs), which tests/overstability/overstability.ipynb also calls, so the two are identical.
#
#usage, from the repo's top folder:
#   uv run python paper/figures/overstability.py     (writes overstability.pdf and overstability.png)

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
from overstability_analysis import MIXES, measured_growth, growth_per_K

LABELS = {'S': r'$\nu_b = 0$', 'S0': r'$\nu_b = 0$, no ring gravity', 'C': r'$\nu_s = 0$',
          'A': r'$\nu_b = \nu_s$', 'AH': r'$\nu_b = \nu_s$, 1600 streamlines', 'D': r'$\nu_b = 2\nu_s$'}
COLORS = {'S': 0, 'S0': 0, 'C': 1, 'A': 2, 'AH': 2, 'D': 3}
MARKERS = {'S': 'o', 'S0': 's', 'C': 'o', 'A': 'o', 'AH': '^', 'D': 'o'}


def unit(p):
    """The viscosity that gamma is scaled by: nu_s, or nu_b when there is no shear viscosity."""
    return p['shear_viscosity'] if p['shear_viscosity'] > 0.0 else p['bulk_viscosity']


def figure(runs):
    """The figure, from runs = {(mix, i): ((p, o), (p0, o0))}, each seeded run and its mix's unseeded control
    (overstability_run(..., mix, 0))."""
    with style(font_scale=2.0) as color_seq:
        f, ax = plt.subplots(figsize=(13, 6))
        x_line = np.linspace(0.0, 0.2, 200)
        drawn = set()
        for (mix, i), ((p, o), (p0, o0)) in runs.items():
            m = measured_growth(p, o, o0)
            ax.plot(m['x'], m['gamma']/(m['K']*unit(p)), marker=MARKERS[mix], markersize=9, linestyle='none',
                    color=color_seq[COLORS[mix]], label=LABELS[mix] if mix not in drawn else None)
            drawn.add(mix)
        #O7 against x: gamma/(K nu) does not depend on K, so one line per viscosity mix (with ring gravity)
        for mix in ('S', 'C', 'A', 'D'):
            p = next(run[0][0] for (m_, i), run in runs.items() if m_ == mix)
            ax.plot(x_line, growth_per_K(p, x_line)/unit(p), linestyle='dashed', linewidth=2,
                    color=color_seq[COLORS[mix]])
        ax.axhline(0.0, color='gray', linewidth=1)
        ax.set_xlim(-0.005, 0.2)
        ax.set_xlabel(r'self-gravity $x = X(\phi)/\kappa^2$')
        ax.set_ylabel(r'growth rate $\gamma/K\nu$')
        ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1.0), fontsize='small')
        f.tight_layout()
    return f


if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    from overstability_analysis import overstability_run, PHIS
    folder = REPO/'tests'/'overstability'
    runs = {(mix, i): (overstability_run(folder, mix, i), overstability_run(folder, mix, 0)) for mix in MIXES for i in PHIS}
    f = figure({k: v for k, v in runs.items() if None not in v})
    for ext in ('pdf', 'png'):
        f.savefig(HERE/f'overstability.{ext}', bbox_inches='tight', dpi=130)
        print('wrote', HERE/f'overstability.{ext}')
