#wire_geometry.py
#
#The 3D analog of hahn13 Fig. 3: the gravity that a perturbing streamline, treated as a straight wire,
#exerts on a particle that lies off the streamline's plane, in three steps, one per panel (the paper's
#Section 3 and Figure 1; notes/phase2/streamline_gravity.md, equation G7):
#(a) the r-z plane at the particle's longitude: the pull A_g toward the streamline's spot, split into
#    an in-plane part A_g cos psi and a vertical part A_g sin psi;
#(b) the ring plane: the wire's tilt phi splits the in-plane part into A_r and A_theta, as in hahn13;
#(c) the theta-z plane: the wire's tilt chi tips the vertical part, adding a little to A_theta.
#A_theta = A_theta,phi + A_theta,chi.
#
#Run:  uv run python paper/figures/wire_geometry.py   (writes wire_geometry.pdf and .png alongside)

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
PULL = dict(arrowstyle='-|>', lw=2.0, mutation_scale=14, color='C3')
PART = dict(arrowstyle='-|>', lw=1.3, mutation_scale=11, color='C0')
AXIS = dict(arrowstyle='->', lw=1)
L = 1.25                                         #length of the drawn pull


def axes(ax, xlabel, ylabel, title):
    ax.set_aspect('equal')
    ax.axis('off')
    ax.annotate('', xy=(3.0, -0.9), xytext=(-1.6, -0.9), arrowprops=AXIS)
    ax.text(3.08, -0.9, xlabel, va='center', fontsize=13)
    ax.annotate('', xy=(-1.4, 2.6), xytext=(-1.4, -1.1), arrowprops=AXIS)
    ax.text(-1.4, 2.7, ylabel, ha='center', fontsize=13)
    ax.set_xlim(-1.8, 3.4)
    ax.set_ylim(-1.3, 3.1)
    ax.set_title(title, fontsize=11)


def arrow(ax, tip, style, label, offset, **kw):
    ax.annotate('', xy=tip, xytext=(0, 0), arrowprops=style)
    ax.text(tip[0] + offset[0], tip[1] + offset[1], label, color=style['color'], fontsize=11, **kw)


def tilted_wire(ax, angle, sym, pull, parts):
    """A wire through (0, 1.8) tilted by angle, the pull along the perpendicular from the particle, and
    that pull's components along the panel's two axes."""
    C = np.array([0.0, 1.8])
    T = np.array([np.cos(angle), np.sin(angle)])
    s = np.array([-1.5, 2.4])
    ax.plot(C[0] + s*T[0], C[1] + s*T[1], 'k-', lw=1.6)
    ax.text(*(C + 2.45*T + np.array([0.05, 0.05])), "$a'$", fontsize=12)
    F = C - (C @ T)*T                            #foot of the perpendicular
    ax.plot([0, F[0]], [0, F[1]], 'k--', lw=1)
    ax.plot([0, 0], [L, C[1]], 'k:', lw=1)
    ax.plot([C[0], C[0] + 1.2], [C[1], C[1]], 'k:', lw=1)
    ax.add_patch(matplotlib.patches.Arc(C, 1.9, 1.9, theta1=0, theta2=np.degrees(angle), lw=1))
    ax.text(C[0] + 1.02, C[1] + 0.12, sym, fontsize=12)
    u = F/np.hypot(*F)
    ax.add_patch(matplotlib.patches.Arc((0, 0), 1.0, 1.0, theta1=90, theta2=np.degrees(np.arctan2(u[1], u[0])), lw=1))
    ax.text(-0.17, 0.58, sym, fontsize=11)
    tip = L*u
    arrow(ax, (0, tip[1]), PART, parts[0], (0.08, -0.25))
    arrow(ax, (tip[0], 0), PART, parts[1], (-0.25, -0.38), ha='center')
    ax.plot([tip[0], tip[0], 0], [0, tip[1], tip[1]], ':', color='C0', lw=0.8)
    arrow(ax, tip, PULL, pull, (-0.1, -0.75), ha='right')
    ax.plot(0, 0, 'ko', ms=9, zorder=5)
    ax.text(0.15, -0.5, 'particle', ha='center', fontsize=10)


fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10.5, 3.7))

#(a) the r-z plane at the particle's longitude; the wire runs along theta-hat, into the page
axes(ax1, 'r', 'z', "(a) r–z plane at the particle's longitude")
W = np.array([2.0, 1.3])
ax1.plot([W[0], W[0]], [-0.15, W[1]], 'k:', lw=1)
ax1.plot([0, W[0]], [-0.15, -0.15], 'k:', lw=1)
ax1.plot([0, W[0]], [0, W[1]], 'k--', lw=1)
ax1.text(W[0]/2, -0.48, r'$\Delta_r$', ha='center', fontsize=12)
ax1.text(W[0] + 0.1, W[1]/2, r'$\Delta_z$', va='center', fontsize=12)
ax1.text(W[0]*0.78 - 0.15, W[1]*0.78 + 0.12, r'$\delta$', fontsize=12)
u = W/np.hypot(*W)
tip = L*u
ax1.add_patch(matplotlib.patches.Arc((0, 0), 1.0, 1.0, theta1=0, theta2=np.degrees(np.arctan2(W[1], W[0])), lw=1))
ax1.text(0.55, 0.08, r'$\psi$', fontsize=11)
arrow(ax1, (tip[0], 0), PART, r'$A_g\cos\psi$', (-0.3, -0.75), ha='center')
arrow(ax1, (0, tip[1]), PART, r'$A_g\sin\psi$', (-0.15, 0.12), ha='right')
ax1.plot([tip[0], tip[0], 0], [0, tip[1], tip[1]], ':', color='C0', lw=0.8)
arrow(ax1, tip, PULL, r'$A_g$', (-0.3, 0.12))
ax1.plot(*W, 'o', ms=13, mfc='white', mec='k', mew=1.5)
ax1.plot(*W, 'kx', ms=8, mew=1.5)
ax1.text(W[0] - 0.55, W[1] + 0.3, "$a'$ at $(r', z')$", fontsize=11)
ax1.plot(0, 0, 'ko', ms=9, zorder=5)
ax1.text(-0.3, -0.5, 'particle', ha='center', fontsize=10)

#(b) the ring plane, as hahn13 Fig. 3: the in-plane part, split by the wire's tilt phi
axes(ax2, r'$\theta$', 'r', r'(b) ring plane ($r$–$\theta$)')
tilted_wire(ax2, np.radians(22), r'$\varphi$', r'$A_g\cos\psi$', (r'$A_r$', r'$A_{\theta,\varphi}$'))

#(c) the theta-z plane, seen from outside the ring: the vertical part, tipped by the wire's tilt chi
axes(ax3, r'$\theta$', 'z', r'(c) $\theta$–$z$ plane, from outside')
tilted_wire(ax3, np.radians(22), r'$\chi$', r'$A_g\sin\psi$', (r'$A_z$', r'$A_{\theta,\chi}$'))

fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(HERE/f'wire_geometry.{ext}', dpi=130)
    print('wrote', HERE/f'wire_geometry.{ext}')
