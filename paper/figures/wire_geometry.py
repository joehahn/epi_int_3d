#wire_geometry.py
#
#hahn13 Fig. 3 moved and tilted into 3D: the gravity that a perturbing streamline a', treated as a
#straight wire, exerts on a particle (the paper's Section 3 and Figure 1; equation G12 of
#notes/sympy/streamline_gravity.py):
#(a) hahn13 Fig. 3, seen from above: the wire, raised a height d_z above the particle but still
#    horizontal, crosses the particle's longitude a distance d_r away, tilted from theta-hat by phi,
#    tan(phi) = s_r = dr'/(r dtheta); the horizontal perpendicular from the particle to it is
#    Delta_h = d_r cos(phi), at the angle phi to r-hat;
#(b) the vertical plane through that horizontal perpendicular: the full perpendicular Delta is the
#    hypotenuse of Delta_h and d_z, and A_g = 2 G lambda/Delta acts along it;
#(c) the theta-z plane, seen from outside the ring: tilting the wire by its vertical slope
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


fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10.5, 3.8))

#(a) hahn13 Fig. 3: r to the right, theta up
frame(ax1, '$r$', r'$\theta$', '(a) seen from above')
dx, phi = 2.0, np.radians(25)
T = np.array([np.sin(phi), np.cos(phi)])
s = np.array([-1.5, 2.2])
ax1.plot(dx + s*T[0], s*T[1], 'k-', lw=1.6)
for q in (-1.4, 1.95):
    ax1.plot(dx + q*T[0], q*T[1], 'ko', ms=6)
ax1.text(dx + 2.25*T[0], 2.25*T[1], "$a'$", fontsize=12)
F = np.array([dx, 0]) - (np.array([dx, 0]) @ T)*T
ax1.plot([0, F[0]], [0, F[1]], 'k-', lw=1)
ax1.text(*(0.7*F + np.array([0.05, -0.42])), r'$\Delta_h$', fontsize=12)
u = F/np.hypot(*F)
ax1.add_patch(matplotlib.patches.Arc((0, 0), 1.9, 1.9, theta1=np.degrees(np.arctan2(u[1], u[0])), theta2=0, lw=1))
ax1.text(1.0, -0.3, r'$\varphi$', fontsize=11)
ax1.text(dx/2, 0.12, r'$d_r$', ha='center', fontsize=12)
V = 1.45*T
ax1.annotate('', xy=(dx + V[0], V[1]), xytext=(dx, 0), arrowprops=STEP)
ax1.plot([dx, dx], [0, V[1]], ':', color='C0', lw=1)
ax1.plot([dx, dx + V[0]], [V[1], V[1]], ':', color='C0', lw=1)
ax1.text(dx - 0.08, V[1]/2, r'$r\,d\theta$', ha='right', color='C0', fontsize=11)
ax1.text(dx + V[0]/2, V[1] + 0.1, "$dr'$", ha='center', color='C0', fontsize=11)
ax1.add_patch(matplotlib.patches.Arc((dx, 0), 1.1, 1.1, theta1=90 - np.degrees(phi), theta2=90, lw=1))
ax1.text(dx + 0.06, 0.62, r'$\varphi$', fontsize=10)
particle(ax1)

#(b) the vertical plane through the horizontal perpendicular Delta_h; the wire runs into the page
frame(ax2, r'along $\Delta_h$', '$z$', '(b) the wire raised by $d_z$', xsize=10)
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

#(c) theta to the right, z up, seen from outside the ring: the raised wire (dotted), then tilted by s_z
frame(ax3, r'$\theta$', '$z$', '(c) the wire tilted by $s_z$')
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
