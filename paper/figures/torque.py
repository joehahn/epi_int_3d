#torque.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 12 November 2026.
#
#the paper's torque figure (Section 4.1.3): the torque that the satellite's resonant potential A_V z sin(chi) exerts on
#the simulated ring of tests/bending_waves, -sum mu m A_V z cos(chi) over its particles, versus time, in units of the
#torque of linear theory, T = -pi^2 m sigma r_V^2 A_V^2/D (Shu 1984; notes/equations.md B9), which is also the
#angular momentum luminosity of the wave it launches (tests/bending_waves/bending_waves.ipynb Section 5).
#fast/satellite.py applies only the vertical pull -dPhi/dz, so this is the torque on the simulated wave, not one the
#ring feels.
#Needs that test's output: uv run python tools/run_simulations.py tests/bending_waves
#
#Drawn by figure(p, o), which tests/bending_waves/bending_waves.ipynb also calls, so the two are identical.
#
#usage, from the repo's top folder:
#   uv run python paper/figures/torque.py     (writes torque.pdf and torque.png)

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
from vertical_analysis import ring_sigma, satellite_torque, nu_a
from satellite import Satellite


def figure(p, o):
    """The figure, from the run's parameters p and output o (load3d)."""
    #the torque at each output, in units of linear theory's
    sat = Satellite(p['satellite'], p['J2'], p['Rp'])
    sigma = ring_sigma(p)
    A_V = sat.mass_final*sat.A_V_per_mass
    torque_theory = -np.pi**2*sat.m*sigma*sat.r_V**2*A_V**2/sat.D_V
    torque = np.array([satellite_torque(p, o, k).sum() for k in range(len(o['times']))])
    #time in units of the wave's crossing time from r_V to the inner edge, at c_g = pi G sigma/nu
    crossing_time = (sat.r_V - 1.0)/(np.pi*sigma/nu_a(p['J2'], p['Rp'], sat.r_V))
    with style(font_scale=1.5) as color_seq:
        f, ax = plt.subplots(figsize=(10, 6))
        ax.plot(o['times']/crossing_time, torque/torque_theory, marker='o', markersize=7, linestyle='-', linewidth=2,
            color=color_seq[0], label='simulation')
        ax.plot(o['times']/crossing_time, o['times']*0 + 1.0, linestyle='dashed', linewidth=2, color=color_seq[0],
            label='linear theory')
        ax.set_xlim(0.0, o['times'][-1]/crossing_time*1.02)
        ax.set_ylim(0.0, 1.3)
        ax.set_xlabel(r'time $t/\tau_c$')
        ax.set_ylabel(r'torque $T/T_V$')
        ax.legend(loc='lower right')
    return f


if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    from vertical_analysis import load3d
    p, o = load3d(REPO/'tests'/'bending_waves')
    f = figure(p, o)
    for ext in ('pdf', 'png'):
        f.savefig(HERE/f'torque.{ext}', bbox_inches='tight', dpi=130)
        print('wrote', HERE/f'torque.{ext}')
