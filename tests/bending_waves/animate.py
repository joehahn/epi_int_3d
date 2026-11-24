#animate.py
#
#Animate tests/bending_waves: the ring's height z above or below the plane (grayscale, white above, black below,
#in units of the wave amplitude |Z_w| of B8), in the frame that rotates with the satellite's pattern speed Omega_p,
#where the m = 2 spiral bending wave is stationary, seen face-on (left, the ring's width stretched radially) and
#unrolled (right, longitude across and x = (r - r_V)/r_V up), as tests/density_waves/animate.py. Each streamline's
#particles are interpolated along the streamline with a periodic cubic spline in longitude onto a fixed longitude grid
#(the streamlines stay circular to ~1e-8 of their spacing, so each is plotted at its mean radius). The face-on panel
#is drawn 3 times finer than its pixels and averaged down to them, so the tightly wound wave doesn't alias into a
#moire pattern. Writes bending_wave.gif, and bending_wave_lines.gif:
#every 3rd streamline's particles joined by lines, z versus longitude in the pattern's frame, each streamline offset
#by its distance from resonance x, so the height is shown directly, without interpolation. Its z scale draws |Z_w| as
#two curve spacings, so curves overlap where neighbours are out of phase. Both show every 2nd output and are committed; bending_waves.ipynb
#links to them. Each is also written as docs/bending_waves/<name>.html, with play, pause and a frame slider, served by
#GitHub Pages at https://joehahn.github.io/epi_int_3d/bending_waves/ (GitHub's file view shows such a page's source).
#
#Run from tests/bending_waves:  uv run python animate.py [output_folder]
#   LAST_FRAME=<file.png> draws only the face-on and unrolled view's last frame, to that file

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
from scipy.ndimage import map_coordinates
from vertical_analysis import load3d, ring_sigma, nu_a
import satellite as satellite_module
from helper_fns import wrap_ring

HERE = os.path.dirname(os.path.abspath(__file__))
folder = sys.argv[1] if len(sys.argv) > 1 else None
p, o = load3d(HERE, output_folder=None if folder is None else os.path.join(HERE, folder))
sat = satellite_module.Satellite(p['satellite'], p['J2'], p['Rp'])
sigma = ring_sigma(p)
A_V = sat.mass_final*sat.A_V_per_mass
Z_w = A_V*np.sqrt(sat.r_V/(sigma*sat.D_V))                              #the wave's amplitude (B8)
tau_c = (sat.r_V - 1.0)/(np.pi*sigma/nu_a(p['J2'], p['Rp'], sat.r_V))   #the wave's crossing time, resonance to inner edge
inner = 0.1                                                             #face-on, the stretched ring spans inner <= rho <= 1

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
r_edges = np.concatenate([[rk[0] - (rk[1] - rk[0])/2], (rk[1:] + rk[:-1])/2, [rk[-1] + (rk[-1] - rk[-2])/2]])
x_edges = (r_edges - sat.r_V)/sat.r_V
phi_edges = np.concatenate([[phi[0] - (phi[1] - phi[0])/2], (phi[1:] + phi[:-1])/2, [phi[-1] + (phi[1] - phi[0])/2]])
z_max = 1.5                                                             #the gray scale's limit, in |Z_w|
DPI = 80
fig, (ax1, ax) = plt.subplots(1, 2, figsize=(17, 9), dpi=DPI, gridspec_kw=dict(width_ratios=[1.15, 1]))
fig.subplots_adjust(left=0.01, right=0.97, top=0.93, bottom=0.07, wspace=0.12)
#left: face-on, the ring's width stretched. Drawn as an image SUPER times finer than the panel's pixels, then averaged
#down to them, so a wave finer than a pixel fades to gray instead of aliasing into a moire pattern
SUPER = 3
lim = 1.01
n_pix = int(round(ax1.get_window_extent().height))                     #the panel's height in pixels, at DPI
u = np.linspace(-lim, lim, SUPER*n_pix + 1)
u = (u[1:] + u[:-1])/2                                                  #the subpixels' centers
X, Y = np.meshgrid(u, u[::-1])
rho = np.hypot(X, Y)
r_sub = r_edges[0] + (rho - inner)/(1 - inner)*(r_edges[-1] - r_edges[0])
outside = (rho < inner) | (rho > 1)
coords = np.array([np.interp(r_sub, rk, np.arange(rk.size)),           #each subpixel's place on the (streamline, phi) grid
                   (np.arctan2(Y, X) - phi[0])/(phi[1] - phi[0])])


def face_on(zk):
    """The face-on image of z/|Z_w|: zk sampled at each subpixel, averaged over SUPER x SUPER."""
    img = map_coordinates(zk/Z_w, coords, order=1, mode='nearest')
    img[outside] = np.nan
    return img.reshape(n_pix, SUPER, n_pix, SUPER).mean(axis=(1, 3))      #a pixel the ring's edge crosses stays blank


mesh1 = ax1.imshow(face_on(zk), cmap='gray', vmin=-z_max, vmax=z_max, extent=(-lim, lim, -lim, lim),
    interpolation='none')
circle = np.linspace(0, 2*np.pi, 400)
stretch_rV = inner + (1 - inner)*(sat.r_V - r_edges[0])/(r_edges[-1] - r_edges[0])
ax1.plot(stretch_rV*np.cos(circle), stretch_rV*np.sin(circle), 'r--', lw=0.8)
ax1.set_aspect('equal')
ax1.set_xlim(-lim, lim)
ax1.set_ylim(-lim, lim)
ax1.axis('off')
fig.text(0.905, 0.07, f"{p['number_of_streamlines']} streamlines\n{p['particles_per_streamline']} particles/streamline\n"
         f"dt = {p['dt']:g}\n\nm = {sat.m} satellite, red: $r_V$\nring width stretched\n"
         f"frame rotates at $\\Omega_p$", fontsize=9, va='bottom')
#right: unrolled
mesh = ax.pcolormesh(phi_edges/np.pi, x_edges, zk/Z_w, cmap='gray', vmin=-z_max, vmax=z_max, shading='flat')
ax.axhline(0, color='r', lw=0.8, ls='--')
ax.set_xlabel(r'longitude in the frame rotating at $\Omega_p$, $/\pi$')
ax.set_ylabel(r'$x = (r - r_V)/r_V$')
fig.colorbar(mesh, ax=ax, shrink=0.6, pad=0.02, label=r'height $z/|Z_w|$')
title = fig.suptitle('', fontsize=14)


def frame(k):
    zk, _ = height(k)
    mesh.set_array((zk/Z_w).ravel())
    mesh1.set_data(face_on(zk))
    title.set_text(f'$t = {o["times"][k]/tau_c:.2f}\\,\\tau_c$  ({o["times"][k]/(2*np.pi):.0f} orbits)')
    return mesh, mesh1, title


if os.environ.get('LAST_FRAME'):                                         #a quick look at the last frame
    frame(len(o['times']) - 1)
    fig.savefig(os.environ['LAST_FRAME'], dpi=DPI)
    sys.exit()
anim = animation.FuncAnimation(fig, frame, frames=range(1, len(o['times']), 2), blit=False)
out = os.path.join(HERE, 'bending_wave.gif')
anim.save(out, writer=animation.PillowWriter(fps=5), dpi=DPI)
print('wrote', out, f'({os.path.getsize(out)/1e6:.1f} MB)')
SITE = os.path.join(REPO, 'docs', 'bending_waves')
os.makedirs(SITE, exist_ok=True)
html = os.path.join(SITE, 'bending_wave.html')
anim.save(html, writer=animation.HTMLWriter(fps=5, embed_frames=True, default_mode='loop'), dpi=DPI)
print('wrote', html)


#the line version: every 3rd streamline, z(phi) offset by x
every = 3
picks = np.arange(0, p['number_of_streamlines'], every)
x0 = (o['r'][0][picks].mean(axis=1) - sat.r_V)/sat.r_V
gap = np.median(np.diff(x0))
scale = 2.0*gap                                                         #|Z_w| is drawn as two curve spacings
fig2, ax2 = plt.subplots(figsize=(8, 11))
lines = [ax2.plot([], [], '-', lw=0.5, color='k')[0] for _ in picks]
ax2.axhline(0, color='r', lw=0.8, ls='--')
ax2.set_xlim(-1, 1)
ax2.set_ylim(x0.min() - gap, x0.max() + gap)
ax2.set_xlabel(r'longitude in the frame rotating at $\Omega_p$, $/\pi$')
ax2.set_ylabel(r'$x = (r - r_V)/r_V$, plus each streamline height $z$')
ax2.text(1.01, 0.0, f"{p['number_of_streamlines']} streamlines,\nevery {every}rd shown\n{p['particles_per_streamline']} particles/streamline\n"
         f"dt = {p['dt']:g}\n\nred: $r_V$\nbar: $z = 2|Z_w|$", transform=ax2.transAxes, fontsize=9, va='bottom')
ax2.plot([1.04]*2, [x0.max() - scale, x0.max() + scale], 'k-', lw=2, clip_on=False)
title2 = ax2.set_title('')
fig2.subplots_adjust(right=0.78)


def frame2(k):
    rot = (o['t'][k] - sat.Omega_p*o['times'][k] + np.pi) % (2*np.pi) - np.pi
    xk = (o['r'][k].mean(axis=1) - sat.r_V)/sat.r_V
    for line, j in zip(lines, picks):
        i = np.argsort(rot[j])
        th = wrap_ring(rot[j, i][None, :], longitude=True)[0]           #the end particles repeated across the edges,
        z = wrap_ring(o['z'][k, j, i][None, :])[0]                      #so each curve spans -pi to pi
        line.set_data(th/np.pi, xk[j] + scale*z/Z_w)
    title2.set_text(f'$t = {o["times"][k]/tau_c:.2f}\\,\\tau_c$  ({o["times"][k]/(2*np.pi):.0f} orbits)')
    return lines + [title2]


anim2 = animation.FuncAnimation(fig2, frame2, frames=range(1, len(o['times']), 2), blit=False)
out2 = os.path.join(HERE, 'bending_wave_lines.gif')
anim2.save(out2, writer=animation.PillowWriter(fps=5), dpi=80)
print('wrote', out2, f'({os.path.getsize(out2)/1e6:.1f} MB)')
html2 = os.path.join(SITE, 'bending_wave_lines.html')
anim2.save(html2, writer=animation.HTMLWriter(fps=5, embed_frames=True, default_mode='loop'), dpi=80)
print('wrote', html2)
