#!/usr/bin/env python3
"""Compare the Python 3 port of epi_int_lite with the Python 2 results in its test notebooks.

Reads the full-length runs in tests/<test>/output (made by running source/epi_int_lite.py in
each test folder), redraws key figures of each test notebook, and places each beside the
notebook's original figure (tests/<test>/reference/cellNN.png). Writes the figures to
notes/phase0/figs/ and prints the port's diagnostics.

Run:  uv run python tools/compare_port.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO/'tests'))
from analysis import (load, unperturbed_metrics, spreading_metrics, spreading_profile,
    spreading_angular_momentum, libration_metrics)
from libration import orbit_deltas

FIGS = REPO/'notes'/'phase0'/'figs'


def side_by_side(name, reference_png, plot_port, title):
    """Notebook figure on the left, the port's version on the right."""
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(14, 4.6),
        gridspec_kw={'width_ratios': [1.15, 1]})
    ax0.imshow(mpimg.imread(reference_png))
    ax0.set_axis_off()
    ax0.set_title('epi_int_lite notebook (Python 2)', fontsize=11)
    plot_port(ax1)
    ax1.set_title('this port (Python 3.14)', fontsize=11)
    fig.suptitle(title, fontsize=13)
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS/f'{name}.png', dpi=110)
    plt.close(fig)


def unperturbed():
    folder = REPO/'tests'/'unperturbed'
    p, o = load(folder)
    ref = folder/'reference'
    times, e, a = o['times'], o['e'], o['a']
    h = o['r']*o['vt']
    e0 = e[0, :, 0]
    order = e0.argsort()

    def de_vs_time(ax):
        for s in range(p['number_of_streamlines']):
            for k in range(p['particles_per_streamline']):
                ax.plot(times/(2*np.pi), e[:, s, k] - e[0, s, k], lw=0.5, alpha=0.5)
        ax.set(xscale='log', yscale='log', xlim=(100, 1.2e5), ylim=(2e-16, 2e-11),
               xlabel='time t    (orbits)', ylabel=r'$|\Delta e|$')

    def delta_vs_e0(x, ylim, ylabel):
        def plot(ax):
            dmax = np.abs(x[-1] - x[0]).max(axis=1)
            ax.plot(e0[order], dmax[order], marker='o', ls='-')
            ax.set(xscale='log', yscale='log', xlim=(5e-9, 3e-2), ylim=ylim,
                   xlabel='initial eccentricity', ylabel=ylabel)
        return plot

    side_by_side('unperturbed_de_vs_time', ref/'cell03.png', de_vs_time,
                 'Unperturbed orbits: e variations vs time')
    side_by_side('unperturbed_da', ref/'cell05.png',
                 delta_vs_e0(a, (2e-16, 2e-11), r'$|\Delta a|$'), 'Unperturbed orbits: a variations')
    side_by_side('unperturbed_dh', ref/'cell06.png',
                 delta_vs_e0(h, (2e-16, 2e-11), r'$|\Delta h|$'), 'Unperturbed orbits: h variations')
    return unperturbed_metrics(p, o)


def spreading():
    folder = REPO/'tests'/'spreading'
    p, o = load(folder)
    ref = folder/'reference'
    tau_nu = p['radial_width']**2/(12*p['shear_viscosity'])
    tv = o['times']/tau_nu
    a_inner = o['a'][:, 0, :].mean(axis=1)
    a_outer = o['a'][:, -1, :].mean(axis=1)

    def width(ax):
        da = a_outer - a_inner
        ax.plot(tv, da/da[0], lw=2)
        ax.set(xlabel=r'time $t/\tau_\nu$', ylabel=r'radial width $\Delta a/\Delta a_0$')

    def dL(ax):
        ax.plot(tv, spreading_angular_momentum(p, o)/1e-7, lw=2)
        ax.set(xlabel=r'time $t/\tau_\nu$', ylabel=r'$\Delta L/L_0$    $(10^{-7})$')

    def sigma(ax):
        for k, t_visc in enumerate([1.0, 2.5, 5.0, 10.0, 20.0, 50.0]):
            t_actual, r1, sd, sd_exp = spreading_profile(p, o, t_visc)
            c = f'C{k}'
            ax.plot(r1, sd_exp, lw=2, alpha=0.7, color=c, label=rf'$t/\tau_\nu$={t_visc:g}')
            ax.plot(r1, sd, marker='o', ms=3, ls='-.', lw=1, color=c)
        ax.set(xlabel=r'radial distance r', ylabel=r'surface density $\sigma/\bar{\sigma}$')
        ax.legend(loc='center right', fontsize=8)

    side_by_side('spreading_width', ref/'cell06.png', width, 'Viscous spreading: ringlet width')
    side_by_side('spreading_angular_momentum', ref/'cell10.png', dL,
                 'Viscous spreading: angular momentum conservation')
    side_by_side('spreading_sigma', ref/'cell11.png', sigma,
                 'Viscous spreading: surface density vs analytic solution (curves)')
    return spreading_metrics(p, o)


def libration():
    folder = REPO/'tests'/'libration'/'narrow_ring'/'ring_mass_1.0e-9'
    p, o = load(folder)
    ref = folder/'reference'
    times = o['times']
    m = libration_metrics(p, o)
    a_inner, a_outer, a_mean, da, e_inner, e_outer, e_mean, de, wt_inner, wt_outer, dwt = \
        orbit_deltas(times, o['r'], o['t'], o['vr'], o['vt'], o['a'], o['e'], o['wt'],
            p['J2'], p['Rp'])
    x = times/(2*np.pi)/1e4

    def dwt_plot(ax):
        ax.plot(x, dwt, marker='o', ms=4, ls='none')
        ax.plot(x, dwt, lw=2, alpha=0.5)
        ax.plot(x, 0*x + m['dwt_eq'], lw=2, color='green')
        ax.set(xlabel=r'time t    ($10^4$ orbits)', ylabel=r'$\Delta\tilde{\omega}$    (radians)')

    def phase_space(ax):
        Ix = (de - m['de_eq'])*a_mean/da
        Iy = e_mean*a_mean*(dwt - m['dwt_eq'])/da
        ax.plot(Ix, Iy, marker='o', ms=2, ls='-', alpha=0.5)
        ax.set(xlabel=r'$I_x = (\Delta e - \Delta e_{eq})a/\Delta a$',
               ylabel=r'$I_y = ea\Delta\tilde{\omega}/\Delta a$')
        ax.set_aspect('equal', adjustable='datalim')

    side_by_side('libration_dwt', ref/'cell09.png', dwt_plot,
                 r'Libration: $\Delta\tilde{\omega}$ at periapse')
    side_by_side('libration_phase_space', ref/'cell17.png', phase_space,
                 'Libration: integral-of-motion level curve')
    return m


if __name__ == '__main__':
    for name, fn in [('unperturbed', unperturbed), ('spreading', spreading),
                     ('libration', libration)]:
        print(f'\n{name}')
        for k, v in fn().items():
            print(f'    {k:24s} {v:.10g}')
    print(f'\nfigures written to {FIGS.relative_to(REPO)}/')
