#animate.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 23 November 2026.
#
#Animate tests/density_waves, the counterpart of tests/bending_waves/animate.py. Each streamline's radial displacement
#is first stripped of its mean radius and its m = 1 part, a least-squares fit of r = r_0 + c cos(theta) + s sin(theta):
#the m = 1 part is the free eccentricity, up to 3.5e-6 at the ring's edges, where the initial circular orbits are out
#of balance with the ring's own radial gravity, about 10 times the wave's amplitude |Z_w|. What is left, delta, holds
#the satellite's m = 2 wave. Seen in the frame that rotates with the satellite's pattern speed Omega_p, where the wave
#is stationary. Writes
#   density_wave.gif: the ring unrolled, longitude in the pattern's frame across and x = (r - r_L)/r_L up, grayscale,
#       the surface density's departure from its mean along each streamline, sigma/sigma_0 - 1 = -d delta/dr
#       (streamline spacing; white dense, black sparse), in percent. Each streamline's delta is interpolated along it
#       with a periodic cubic spline in longitude onto a fixed longitude grid, and d delta/dr is the centered
#       difference between neighboring streamlines. Not face-on, as tests/bending_waves' is: the wave winds about 60
#       times across the ring, which a face-on frame can't resolve. Nor is there a lines version: with about 4 of
#       every 3rd streamline per wavelength, the curves either cross or barely move.
#It shows every output and is committed. It is also written as docs/density_waves/density_wave.html, with play, pause
#and a frame slider, served by GitHub Pages at https://joehahn.github.io/epi_int_3d/density_waves/.
#
#Run from tests/density_waves:  uv run python animate.py [output_folder]

import os
import sys

REPO = os.path.abspath(os.path.dirname(__file__))
while not os.path.exists(os.path.join(REPO, 'pyproject.toml')):
    REPO = os.path.dirname(REPO)
for d in ('source', 'fast', 'tests'):
    sys.path.insert(0, os.path.join(REPO, d))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.interpolate import CubicSpline
from density_analysis import load3d, satellite, tau_001

HERE = os.path.dirname(os.path.abspath(__file__))
folder = sys.argv[1] if len(sys.argv) > 1 else None
p, o = load3d(HERE, output_folder=None if folder is None else os.path.join(HERE, folder))
sat = satellite(p)
tau = tau_001(p, sat)                                                   #the wave's travel time from r_L to x = 0.01

phi = np.linspace(-np.pi, np.pi, 361)                                   #longitude in the pattern's frame


def displacement(k):
    """Each particle's radial displacement delta at output k, less its streamline's mean radius and m = 1 part, its
    longitude in the pattern's frame, and each streamline's mean radius."""
    th = o['t'][k]
    delta = np.empty_like(o['r'][k])
    r0 = np.empty(th.shape[0])
    for j in range(th.shape[0]):
        A = np.column_stack([np.ones_like(th[j]), np.cos(th[j]), np.sin(th[j])])
        c, *_ = np.linalg.lstsq(A, o['r'][k, j], rcond=None)
        delta[j] = o['r'][k, j] - A@c
        r0[j] = c[0]
    rot = (th - sat.Omega_p*o['times'][k] + np.pi) % (2*np.pi) - np.pi
    return delta, rot, r0


def density(k):
    """sigma/sigma_0 - 1 = -d delta/dr on the (streamline, phi) grid at output k, and each streamline's mean radius."""
    delta, rot, r0 = displacement(k)
    dk = np.empty((rot.shape[0], phi.size))
    for j in range(rot.shape[0]):
        i = np.argsort(rot[j])
        th, d = rot[j, i], delta[j, i]
        th, d = np.append(th, th[0] + 2*np.pi), np.append(d, d[0])     #close the loop for the periodic spline
        dk[j] = CubicSpline(th, d, bc_type='periodic')((phi - th[0]) % (2*np.pi) + th[0])
    return -np.gradient(dk, r0, axis=0), r0


sk, rk = density(0)
x_edges = (np.concatenate([[rk[0] - (rk[1] - rk[0])/2], (rk[1:] + rk[:-1])/2, [rk[-1] + (rk[-1] - rk[-2])/2]]) - sat.r_L)/sat.r_L
phi_edges = np.concatenate([[phi[0] - (phi[1] - phi[0])/2], (phi[1:] + phi[:-1])/2, [phi[-1] + (phi[1] - phi[0])/2]])
s_max = 1.2                                                             #the gray scale's limit, percent
fig, ax = plt.subplots(figsize=(8, 11))
mesh = ax.pcolormesh(phi_edges/np.pi, x_edges, 100*sk, cmap='gray', vmin=-s_max, vmax=s_max, shading='flat')
ax.axhline(0, color='r', lw=0.8, ls='--')
ax.set_xlabel(r'longitude in the frame rotating at $\Omega_p$, $/\pi$')
ax.set_ylabel(r'$x = (r - r_L)/r_L$')
fig.colorbar(mesh, ax=ax, shrink=0.5, pad=0.02, label=r'surface density $\sigma/\sigma_0 - 1$  (%)')
title = ax.set_title('')
ax.text(1.03, 0.0, f"{p['number_of_streamlines']} streamlines\n{p['particles_per_streamline']} particles/streamline\n"
        f"dt = {p['dt']:g}\n\nm = {sat.m} satellite\nred: $r_L$\nm = 1 part removed", transform=ax.transAxes,
        fontsize=9, va='bottom')


def frame(k):
    sk, _ = density(k)
    mesh.set_array((100*sk).ravel())
    title.set_text(f'$t = {o["times"][k]/tau:.2f}\\,\\tau_{{0.01}}$  ({o["times"][k]/(2*np.pi):.0f} orbits)')
    return mesh, title


anim = animation.FuncAnimation(fig, frame, frames=range(len(o['times'])), blit=False)
out = os.path.join(HERE, 'density_wave.gif')
anim.save(out, writer=animation.PillowWriter(fps=5), dpi=80)
print('wrote', out, f'({os.path.getsize(out)/1e6:.1f} MB)')
SITE = os.path.join(REPO, 'docs', 'density_waves')
os.makedirs(SITE, exist_ok=True)
html = os.path.join(SITE, 'density_wave.html')
anim.save(html, writer=animation.HTMLWriter(fps=5, embed_frames=True, default_mode='loop'), dpi=80)
print('wrote', html)

