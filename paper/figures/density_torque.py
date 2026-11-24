#density_torque.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 23 November 2026.
#
#the paper's second density-wave figure (Section 4.2.4), the counterpart of wavenumber.py and torque.py for
#tests/density_waves: (a) the simulated wave's wavenumber |k| r_L versus x at the end of the run, k = d arg(xi - L_L/x)/dr
#between neighboring streamlines (the wave part of each streamline's radial displacement xi, as B21), from the resonance
#zone to the measured front, against the dispersion relation D = 2 pi G sigma |k| (Eqn. 40; notes/equations.md L5);
#(b) the satellite's torque on the ring, sum of mu r A_theta over its particles, versus time, in units of linear
#theory's T_L = -pi^2 m sigma Psi^2/D (Eqn. 46; L8).
#Needs tests/density_waves' output: uv run python tools/run_simulations.py tests/density_waves
#
#Drawn by figure(p, o), which tests/density_waves/density_waves.ipynb also calls, so the two are identical.
#
#usage, from the repo's top folder:
#   uv run python paper/figures/density_torque.py     (writes density_torque.pdf and density_torque.png)

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
from density_analysis import (simulated_xi, L_L, wave_wavenumber, wavenumber_profile, resonance_zone, wave_front,
    satellite_torque, torque_theory, tau_001)


def figure(p, o):
    """The figure, from the run's parameters p and output o (load3d)."""
    sat, x, xi = simulated_xi(p, o, -1)
    x_mid, k_sim = wavenumber_profile(x, xi - L_L(sat)/x, sat.r_L)            #the wave part, xi - L_L/x
    x_front = wave_front(p, o, -1)
    behind = (x_mid > resonance_zone(p, sat)) & (x_mid < x_front)
    tau = tau_001(p, sat)
    torque = np.array([satellite_torque(p, o, k).sum() for k in range(len(o['times']))])/torque_theory(p, sat)
    with style(font_scale=1.5) as color_seq:
        f, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
        #(a) the wavenumber behind the front, and the dispersion relation
        ax1.plot(x_mid[behind], np.abs(k_sim[behind])*sat.r_L/1.0e4, marker='o', markersize=4, linestyle='none',
            alpha=0.7, color=color_seq[0], label='simulation')
        x_theory = np.linspace(0.0, x.max(), 400)
        ax1.plot(x_theory, wave_wavenumber(p, sat, x_theory)*sat.r_L/1.0e4, linestyle='dashed', linewidth=2,
            color='black', label='theory')
        ax1.axvline(x_front, color=color_seq[0], linewidth=1.5, linestyle='dotted')
        ax1.set_xlim(0.0, x.max())
        ax1.set_ylim(0.0, 4.5)                                         #as wavenumber.py
        ax1.set_xlabel(r'fractional distance from resonance $x$')
        ax1.set_ylabel(r'wavenumber $|k|r_L$    ($10^4$)')
        ax1.legend(loc='upper left')
        ax1.text(0.97, 0.06, '(a)', transform=ax1.transAxes, ha='right', va='bottom')
        #(b) the torque, in units of T_L
        ax2.plot(o['times']/tau, torque, marker='o', markersize=7, linestyle='-', linewidth=2, color=color_seq[0],
            label='simulation')
        ax2.plot(o['times']/tau, o['times']*0 + 1.0, linestyle='dashed', linewidth=2, color='black', label='theory')
        ax2.set_xlim(0.0, o['times'][-1]/tau*1.02)
        ax2.set_ylim(0.0, 1.3)
        ax2.set_xlabel(r'time $t/\tau_{0.01}$')
        ax2.set_ylabel(r'torque $T/T_L$')
        ax2.legend(loc='lower right')
        ax2.text(0.03, 0.06, '(b)', transform=ax2.transAxes, ha='left', va='bottom')
        f.tight_layout()
    return f


if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    from density_analysis import load3d
    p, o = load3d(REPO/'tests'/'density_waves')
    f = figure(p, o)
    for ext in ('pdf', 'png'):
        f.savefig(HERE/f'density_torque.{ext}', bbox_inches='tight', dpi=130)
        print('wrote', HERE/f'density_torque.{ext}')
