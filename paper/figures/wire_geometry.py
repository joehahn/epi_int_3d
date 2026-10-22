#wire_geometry.py
#
#The 3D analog of hahn13 Fig. 3: the gravity that a perturbing streamline, treated as a straight wire,
#exerts on a particle that lies off the streamline's plane, drawn as three projections that each show
#one angle (the paper's Section 3 and Figure 1; notes/phase2/streamline_gravity.md, equation G7).
#
#Run:  uv run python paper/figures/wire_geometry.py   (writes wire_geometry.pdf and .png alongside)

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
PULL = dict(arrowstyle='-|>', lw=1.8, mutation_scale=14, color='C3')
AXIS = dict(arrowstyle='->', lw=1)


def axes(ax, xlabel, ylabel, title):
    ax.set_aspect('equal')
    ax.axis('off')
    ax.annotate('', xy=(3.0, -0.9), xytext=(-1.2, -0.9), arrowprops=AXIS)
    ax.text(3.08, -0.9, xlabel, va='center', fontsize=13)
    ax.annotate('', xy=(-1.0, 2.6), xytext=(-1.0, -1.1), arrowprops=AXIS)
    ax.text(-1.0, 2.7, ylabel, ha='center', fontsize=13)
    ax.set_xlim(-1.4, 3.4)
    ax.set_ylim(-1.3, 3.1)
    ax.set_title(title, fontsize=11)


def tilted_wire(ax, angle, name, sym):
    """A wire through (0, 1.8) tilted by angle; the perpendicular from the particle at the origin, the
    pull along it, and the angle between the perpendicular and the vertical axis of the panel."""
    C = np.array([0.0, 1.8])
    T = np.array([np.cos(angle), np.sin(angle)])
    s = np.array([-1.3, 2.4])
    ax.plot(C[0] + s*T[0], C[1] + s*T[1], 'k-', lw=1.6)
    ax.text(*(C + 2.45*T + np.array([0.05, 0.05])), name, fontsize=12)
    ax.plot([0, 0], [0, C[1]], 'k:', lw=1)
    F = C - (C @ T)*T                            #foot of the perpendicular
    ax.plot([0, F[0]], [0, F[1]], 'k--', lw=1)
    u = F/np.hypot(*F)
    ax.annotate('', xy=1.0*u, xytext=(0, 0), arrowprops=PULL)
    ax.text(*(0.6*u + np.array([-0.5, -0.05])), r'$A_g$', color='C3', fontsize=12)
    ax.add_patch(matplotlib.patches.Arc((0, 0), 1.6, 1.6, theta1=90, theta2=np.degrees(np.arctan2(u[1], u[0])), lw=1))
    ax.text(-0.18, 0.95, sym, fontsize=12)
    ax.plot([C[0], C[0] + 1.2], [C[1], C[1]], 'k:', lw=1)
    ax.add_patch(matplotlib.patches.Arc(C, 1.9, 1.9, theta1=0, theta2=np.degrees(angle), lw=1))
    ax.text(C[0] + 1.02, C[1] + 0.12, sym, fontsize=12)
    ax.plot(0, 0, 'ko', ms=9)
    ax.text(0.1, -0.45, 'particle', ha='center', fontsize=10)


fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 3.6))

#(a) the particle's meridional plane: the wire runs along theta-hat, into the page, through (r', z')
axes(ax1, 'r', 'z', '(a) meridional plane')
W = np.array([2.0, 1.3])
ax1.plot([0, W[0]], [0, 0], 'k:', lw=1)
ax1.plot([W[0], W[0]], [0, W[1]], 'k:', lw=1)
ax1.plot([0, W[0]], [0, W[1]], 'k--', lw=1)
ax1.text(W[0]/2, -0.35, r'$\Delta_r$', ha='center', fontsize=12)
ax1.text(W[0] + 0.1, W[1]/2, r'$\Delta_z$', va='center', fontsize=12)
ax1.text(W[0]*0.55 - 0.1, W[1]*0.55 + 0.15, r'$\delta$', fontsize=12)
u = W/np.hypot(*W)
ax1.annotate('', xy=1.0*u, xytext=(0, 0), arrowprops=PULL)
ax1.text(*(0.55*u + np.array([-0.45, 0.2])), r'$A_g$', color='C3', fontsize=12)
ax1.add_patch(matplotlib.patches.Arc((0, 0), 1.5, 1.5, theta1=0, theta2=np.degrees(np.arctan2(W[1], W[0])), lw=1))
ax1.text(0.8, 0.12, r'$\psi$', fontsize=12)
ax1.plot(*W, 'o', ms=13, mfc='white', mec='k', mew=1.5)
ax1.plot(*W, 'kx', ms=8, mew=1.5)
ax1.text(W[0] - 0.2, W[1] + 0.35, "$a'$", fontsize=12)
ax1.plot(0, 0, 'ko', ms=9)
ax1.text(0.1, -0.45, 'particle', ha='center', fontsize=10)

#(b) seen from above the ring plane: hahn13 Fig. 3, the wire tilted by phi
axes(ax2, r'$\theta$', 'r', '(b) ring plane')
tilted_wire(ax2, np.radians(22), "$a'$", r'$\varphi$')
ax2.text(0.08, 0.6*1.8 + 0.25, r'$\Delta_r$', fontsize=12)

#(c) seen from outside the ring: the wire tilted by chi
axes(ax3, r'$\theta$', 'z', r'(c) from outside, along $-\hat r$')
tilted_wire(ax3, np.radians(22), "$a'$", r'$\chi$')
ax3.text(0.08, 0.6*1.8 + 0.25, r'$\Delta_z$', fontsize=12)

fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(HERE/f'wire_geometry.{ext}', dpi=130)
    print('wrote', HERE/f'wire_geometry.{ext}')
