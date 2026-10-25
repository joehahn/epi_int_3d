#wire_geometry.py
#
#No longer in the paper (its Section 3.1 now uses the vector form D = t x (d x t), 2026-10-25); kept
#for notes/phase2/streamline_gravity.md.
#hahn13 Fig. 3 moved and tilted into 3D: the gravity that a perturbing streamline a', treated as a
#straight wire, exerts on a particle (equation G12 of notes/sympy/streamline_gravity.py). Seen from above, the wire raised a height d_z above the particle
#but still horizontal is the planar geometry of the paper's Appendix B (wire_planar.py), with the
#horizontal perpendicular Delta_h = d_r cos(phi). Then:
#(a) the vertical plane through that horizontal perpendicular: the full perpendicular Delta is the
#    hypotenuse of Delta_h and d_z, and A_g = 2 G lambda/Delta acts along it;
#(b) the theta-z plane, seen from outside the ring: tilting the wire by its vertical slope
#    s_z = dz'/(r dtheta) moves the nearest point along theta by -d_z s_z.
#
#Run:  uv run python paper/figures/wire_geometry.py   (writes wire_geometry.pdf and .png alongside)

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
AXIS = dict(arrowstyle='->', lw=1)
PULL = dict(arrowstyle='-|>', lw=2.0, mutation_scale=14, color='C3')
STEP = dict(arrowstyle='-|>', lw=1.5, mutation_scale=12, color='C0')


def frame(ax, xlabel, ylabel, title, xsize=13):
    ax.set_aspect('equal')
    ax.axis('off')
    ax.annotate('', xy=(3.3, 0), xytext=(-1.0, 0), arrowprops=AXIS)
    ax.text(3.38, 0, xlabel, va='center', fontsize=xsize)
    ax.annotate('', xy=(0, 2.4), xytext=(0, -1.5), arrowprops=AXIS)
    ax.text(0, 2.5, ylabel, ha='center', fontsize=13)
    ax.set_xlim(-1.3, 3.9)
    ax.set_ylim(-1.7, 2.8)
    ax.set_title(title, fontsize=11)


def particle(ax):
    ax.plot(0, 0, 'ko', ms=9, zorder=5)
    ax.text(-0.15, 0.15, 'particle', ha='right', fontsize=10)


fig, (ax2, ax3) = plt.subplots(1, 2, figsize=(7.2, 3.8))

#(a) the vertical plane through the horizontal perpendicular Delta_h; the wire runs into the page
frame(ax2, r'along $\Delta_h$', '$z$', '(a) the wire raised by $d_z$', xsize=10)
W = np.array([2.0, 1.3])
ax2.plot([W[0], W[0]], [0, W[1]], 'k:', lw=1)
ax2.plot([0, W[0]], [0, W[1]], 'k-', lw=1)
ax2.text(W[0]/2, -0.35, r'$\Delta_h$', ha='center', fontsize=12)
ax2.text(W[0] + 0.1, W[1]/2, r'$d_z$', va='center', fontsize=12)
ax2.text(W[0]*0.72 - 0.2, W[1]*0.72 + 0.12, r'$\Delta$', fontsize=12)
u = W/np.hypot(*W)
ax2.annotate('', xy=0.95*u, xytext=(0, 0), arrowprops=PULL)
ax2.text(*(0.95*u + np.array([-0.35, 0.12])), r'$A_g$', color='C3', fontsize=11)
ax2.plot(*W, 'o', ms=13, mfc='white', mec='k', mew=1.5)
ax2.plot(*W, 'kx', ms=8, mew=1.5)
ax2.text(W[0] - 0.15, W[1] + 0.3, "$a'$", fontsize=12)
particle(ax2)

#(b) theta to the right, z up, seen from outside the ring: the raised wire (dotted), then tilted by s_z
frame(ax3, r'$\theta$', '$z$', '(b) the wire tilted by $s_z$')
dz, chi = 1.1, np.radians(25)
ax3.plot([-1.1, 3.1], [dz, dz], ':', color='0.5', lw=1.2)
T = np.array([np.cos(chi), np.sin(chi)])
C = np.array([0.0, dz])
s = np.array([-1.0, 2.5])
ax3.plot(C[0] + s*T[0], C[1] + s*T[1], 'k-', lw=1.6)
ax3.text(*(C + 2.55*T + np.array([0.05, 0.0])), "$a'$", fontsize=12)
F = C - (C @ T)*T
ax3.plot([0, F[0]], [0, F[1]], 'k-', lw=1)
ax3.plot([0, 0], [0, dz], 'k:', lw=1)
ax3.text(0.08, 0.45*dz, '$d_z$', fontsize=12)
ax3.annotate('', xy=(F[0], dz + 0.5), xytext=(0, dz + 0.5), arrowprops=dict(arrowstyle='->', lw=1))
ax3.plot([F[0], F[0]], [F[1], dz + 0.6], 'k:', lw=0.8)
ax3.text(F[0] - 0.1, dz + 0.65, r'$d_zs_z$', ha='center', fontsize=11)
V = 1.3*T
ax3.annotate('', xy=C + V, xytext=C, arrowprops=STEP)
ax3.plot([C[0] + V[0], C[0] + V[0]], [dz, dz + V[1]], ':', color='C0', lw=1)
ax3.text(C[0] + V[0]/2, dz - 0.3, r'$r\,d\theta$', ha='center', color='C0', fontsize=11)
ax3.text(C[0] + V[0] + 0.08, dz + V[1]/2 - 0.05, "$dz'$", color='C0', fontsize=11)
particle(ax3)

fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(HERE/f'wire_geometry.{ext}', dpi=130)
    print('wrote', HERE/f'wire_geometry.{ext}')
