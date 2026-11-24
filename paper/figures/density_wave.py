#density_wave.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 23 November 2026.
#
#the paper's density-wave figure (Section 4.2.4), the counterpart of bending_wave.py: each streamline's radial amplitude
#|xi|/r, xi from a least-squares fit of its particles' radii to r_0 + Re[xi e^{i m(theta - Omega_p t)}] (simulated_xi),
#versus its fractional distance from resonance x over -0.002 < x < 0.011 (bending_wave.py's range,
#mirrored), at six times, in units of the time tau_0.01 that the
#wave takes to reach x = 0.01, as the wave propagates outward, and linear theory's |xi|/r, with xi the trailing solution
#of Eqn. (45) (notes/equations.md L7). The fit keeps only the m-armed part of each streamline's shape: its
#osculating eccentricity also holds a free part, up to 3.5e-6 at the ring's edges, where the initial circular orbits
#are out of balance with the ring's own radial gravity (the same in tests/bending_waves' ring at its first output).
#Needs tests/density_waves' output: uv run python tools/run_simulations.py tests/density_waves
#
#Drawn by figure(p, o), which tests/density_waves/density_waves.ipynb also calls, so the two are identical.
#
#usage, from the repo's top folder:
#   uv run python paper/figures/density_wave.py     (writes density_wave.pdf and density_wave.png)

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
from density_analysis import satellite, simulated_xi, density_wave, ring_sigma, tau_001, wave_amplitude


def figure(p, o):
    """The figure, from the run's parameters p and output o (load3d)."""
    sat = satellite(p)
    sigma = ring_sigma(p)
    tau = tau_001(p, sat)
    x_min, x_max = -0.002, 0.011                                       #Fig. 3's range, mirrored
    with style(font_scale=1.5) as color_seq:
        f, ax = plt.subplots(figsize=(10, 6))
        #each streamline's |xi|/r at six times, in units of tau_0.01
        time_fractions = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5]
        for idx, time_fraction in enumerate(time_fractions):
            tidx = np.argmin(np.abs(o['times'] - time_fraction*tau))
            sat_unused, x_now, xi_now = simulated_xi(p, o, tidx)
            linewidth = 2 if (idx == len(time_fractions) - 1) else 1
            label = r'$t/\tau_{0.01}$ = ' + str(np.round(o['times'][tidx]/tau, 2))
            ax.plot(x_now, np.abs(xi_now)/(sat.r_L*(1.0 + x_now))/1.0e-7, linestyle='-', linewidth=linewidth,
                color=color_seq[idx], label=label)
        #linear theory, dashed in the color of the latest time
        xs = np.linspace(x_min, x_max, 20000)
        e_theory = np.abs(density_wave(sat, sigma, xs))/(sat.r_L*(1.0 + xs))
        ax.plot(xs, e_theory/1.0e-7, linestyle='dashed', linewidth=2, color=color_seq[len(time_fractions) - 1],
            label='linear theory')
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(0.0, 2.0*wave_amplitude(p, sat)/sat.r_L/1.0e-7)
        ax.set_xlabel(r'fractional distance from resonance $x$')
        ax.set_ylabel(r'radial amplitude $|\xi|/r$    ($10^{-7}$)')
        ax.legend(loc='upper left', ncol=2)
    return f


if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    from density_analysis import load3d
    p, o = load3d(REPO/'tests'/'density_waves')
    f = figure(p, o)
    for ext in ('pdf', 'png'):
        f.savefig(HERE/f'density_wave.{ext}', bbox_inches='tight', dpi=130)
        print('wrote', HERE/f'density_wave.{ext}')
