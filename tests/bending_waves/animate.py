#animate.py
#
#Animate tests/bending_waves: the ring's height z above or below the plane (grayscale, white above, black below,
#in units of the wave amplitude |Z_w| of B8), one frame per output, seen face-on in the frame that rotates with the
#satellite's pattern speed Omega_p, where the m = 2 spiral bending wave is stationary. The ring's width is stretched
#radially so that the wave is visible. Each streamline's particles are interpolated along the streamline with a
#periodic cubic spline in longitude onto a fixed longitude grid (the streamlines stay circular to ~1e-8 of their
#spacing, so each is plotted at its mean radius). Writes bending_wave.gif, which bending_waves.ipynb shows (committed),
#and figs/bending_wave_lines.gif (git-ignored):
#every 3rd streamline's particles joined by lines, z versus longitude in the pattern's frame, each streamline offset
#by its distance from resonance x, so the height is shown directly, without interpolation. Its z scale is the largest
#that keeps adjacent curves apart in every frame.
#
#Run from tests/bending_waves:  uv run python animate.py [output_folder]

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
from vertical_analysis import load3d, ring_sigma, nu_a
import satellite as satellite_module

HERE = os.path.dirname(os.path.abspath(__file__))
folder = sys.argv[1] if len(sys.argv) > 1 else None
p, o = load3d(HERE, output_folder=None if folder is None else os.path.join(HERE, folder))
sat = satellite_module.Satellite(p['satellite'], p['J2'], p['Rp'])
sigma = ring_sigma(p)
A_V = sat.mass_final*sat.A_V_per_mass
Z_w = A_V*np.sqrt(sat.r_V/(sigma*sat.D_V))                              #the wave's amplitude (B8)
tau_c = (sat.r_V - 1.0)/(np.pi*sigma/nu_a(p['J2'], p['Rp'], sat.r_V))   #the wave's crossing time, resonance to inner edge
W = p['radial_width']
inner = 0.35                                                            #the stretched ring spans inner <= rho <= 1
stretch = lambda r: inner + (1 - inner)*(r - 1.0)/W

phi = np.linspace(-np.pi, np.pi, 361)                                   #longitude in the pattern's frame


def height(k):
    """z on the (streamline, phi) grid at output k, and each streamline's mean radius."""
    rot = (o['t'][k] - sat.Omega_p*o['times'][k] + np.pi) % (2*np.pi) - np.pi
    zk = np.empty((rot.shape[0], phi.size))
    for j in range(rot.shape[0]):
        i = np.argsort(rot[j])
        th, z = rot[j, i], o['z'][k, j, i]
        th, z = np.append(th, th[0] + 2*np.pi), np.append(z, z[0])     #close the loop for the periodic spline
        zk[j] = CubicSpline(th, z, bc_type='periodic')((phi - th[0]) % (2*np.pi) + th[0])
    return zk, o['r'][k].mean(axis=1)


zk, rk = height(0)
rho_edges = stretch(np.concatenate([[rk[0] - (rk[1] - rk[0])/2], (rk[1:] + rk[:-1])/2, [rk[-1] + (rk[-1] - rk[-2])/2]]))
phi_edges = np.concatenate([[phi[0] - (phi[1] - phi[0])/2], (phi[1:] + phi[:-1])/2, [phi[-1] + (phi[1] - phi[0])/2]])
P, R = np.meshgrid(phi_edges, rho_edges)
fig, ax = plt.subplots(figsize=(7.5, 6.5))
mesh = ax.pcolormesh(R*np.cos(P), R*np.sin(P), zk/Z_w, cmap='gray', vmin=-1.5, vmax=1.5, shading='flat')
circle = np.linspace(0, 2*np.pi, 400)
ax.plot(stretch(sat.r_V)*np.cos(circle), stretch(sat.r_V)*np.sin(circle), 'r--', lw=0.8)
ax.set_aspect('equal')
ax.axis('off')
fig.colorbar(mesh, ax=ax, shrink=0.7, label=r'height $z/|Z_w|$')
title = ax.set_title('')
ax.text(-1.25, -1.05, f"{p['number_of_streamlines']} streamlines\n{p['particles_per_streamline']} particles/streamline\n"
        f"dt = {p['dt']:g}\n\nm = {sat.m} satellite, red: $r_V$\nring width stretched\n"
        f"frame rotates at $\\Omega_p$", fontsize=9, va='bottom')


def frame(k):
    zk, _ = height(k)
    mesh.set_array((zk/Z_w).ravel())
    title.set_text(f'$t = {o["times"][k]/tau_c:.2f}\\,\\tau_c$  ({o["times"][k]/(2*np.pi):.0f} orbits)')
    return mesh, title


os.makedirs(os.path.join(HERE, 'figs'), exist_ok=True)
anim = animation.FuncAnimation(fig, frame, frames=len(o['times']), blit=False)
out = os.path.join(HERE, 'bending_wave.gif')
anim.save(out, writer=animation.PillowWriter(fps=8), dpi=80)
print('wrote', out, f'({os.path.getsize(out)/1e6:.1f} MB)')


#the line version: every 3rd streamline, z(phi) offset by x
every = 3
picks = np.arange(0, p['number_of_streamlines'], every)
x0 = (o['r'][0][picks].mean(axis=1) - sat.r_V)/sat.r_V
gap = np.median(np.diff(x0))
#the largest height difference between adjacent shown streamlines at the same longitude, over all frames (on the
#spline grid), fills 90% of the gap between them, so the curves never cross
dz = max(np.abs(np.diff(height(k)[0][picks], axis=0)).max() for k in range(len(o['times'])))
scale = 0.9*gap*Z_w/dz                                                  #offset per unit z/|Z_w|
fig2, ax2 = plt.subplots(figsize=(8, 8))
lines = [ax2.plot([], [], '-', lw=0.6, color='k')[0] for _ in picks]
ax2.axhline(0, color='r', lw=0.8, ls='--')
ax2.set_xlim(-np.pi, np.pi)
ax2.set_ylim(x0.min() - gap, x0.max() + gap)
ax2.set_xlabel(r'longitude in the frame rotating at $\Omega_p$')
ax2.set_ylabel(r'$x = (r - r_V)/r_V$, plus each streamline height $z$')
ax2.text(1.01, 0.0, f"{p['number_of_streamlines']} streamlines,\nevery {every}rd shown\n{p['particles_per_streamline']} particles/streamline\n"
         f"dt = {p['dt']:g}\n\nred: $r_V$\nbar: $z = 2|Z_w|$", transform=ax2.transAxes, fontsize=9, va='bottom')
ax2.plot([np.pi*1.04]*2, [x0.max() - scale, x0.max() + scale], 'k-', lw=2, clip_on=False)
title2 = ax2.set_title('')
fig2.subplots_adjust(right=0.78)


def frame2(k):
    rot = (o['t'][k] - sat.Omega_p*o['times'][k] + np.pi) % (2*np.pi) - np.pi
    xk = (o['r'][k].mean(axis=1) - sat.r_V)/sat.r_V
    for line, j in zip(lines, picks):
        i = np.argsort(rot[j])
        th, z = rot[j, i], o['z'][k, j, i]
        line.set_data(np.append(th, th[0] + 2*np.pi), xk[j] + scale*np.append(z, z[0])/Z_w)
    title2.set_text(f'$t = {o["times"][k]/tau_c:.2f}\\,\\tau_c$  ({o["times"][k]/(2*np.pi):.0f} orbits)')
    return lines + [title2]


anim2 = animation.FuncAnimation(fig2, frame2, frames=len(o['times']), blit=False)
out2 = os.path.join(HERE, 'figs', 'bending_wave_lines.gif')
anim2.save(out2, writer=animation.PillowWriter(fps=8), dpi=80)
print('wrote', out2, f'({os.path.getsize(out2)/1e6:.1f} MB)')
