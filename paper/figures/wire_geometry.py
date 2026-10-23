#wire_geometry.py
#
#hahn13 Fig. 3 extended to 3D: the gravity that a perturbing streamline a', treated as a straight wire,
#exerts on a particle, in hahn13's layout and notation (the paper's Section 3 and Figure 1;
#notes/phase2/streamline_gravity.md, equations G10-G11):
#(a) the r-theta plane, as hahn13 Fig. 3 but with the wire along the streamline's shape: tilted by
#    phi from theta-hat, tan(phi) = s_r = dr'/(r dtheta), so the perpendicular Delta from the particle
#    makes the angle phi with r-hat, and A_g has components A_g cos(phi) and -A_g sin(phi);
#(b) the z-theta plane, the same construction with z in place of r: the slope s_z = dz'/(r dtheta);
#(c) the r-z plane at the particle's longitude: the separation d = (d_r, 0, d_z) to the crossing point.
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
VEL = dict(arrowstyle='-|>', lw=1.5, mutation_scale=12, color='C0')


def hahn13_panel(ax, x, xp, vx, angle, title):
    """hahn13 Fig. 3: the particle at the origin, the coordinate x (r or z) to the right and theta up.
    The wire crosses the particle's longitude a distance d_x away and is tilted from theta-hat by
    the streamline's slope (dx'/r dtheta); the perpendicular Delta from the particle makes the same angle with x."""
    ax.set_aspect('equal')
    ax.axis('off')
    ax.annotate('', xy=(3.3, 0), xytext=(-1.0, 0), arrowprops=AXIS)
    ax.text(3.38, 0, f'${x}$', va='center', fontsize=13)
    ax.annotate('', xy=(0, 2.4), xytext=(0, -1.5), arrowprops=AXIS)
    ax.text(0, 2.5, r'$\theta$', ha='center', fontsize=13)
    dx = 2.0
    T = np.array([np.sin(angle), np.cos(angle)])            #the wire's direction, along v'
    s = np.array([-1.5, 2.2])
    ax.plot(dx + s*T[0], s*T[1], 'k-', lw=1.6)
    for q in (-1.4, 1.95):
        ax.plot(dx + q*T[0], q*T[1], 'ko', ms=6)
    ax.text(dx + 2.25*T[0], 2.25*T[1], "$a'$", fontsize=12)
    F = np.array([dx, 0]) - (np.array([dx, 0]) @ T)*T        #foot of the perpendicular
    ax.plot([0, F[0]], [0, F[1]], 'k-', lw=1)
    ax.text(*(0.7*F + np.array([0.05, -0.4])), r'$\Delta$', fontsize=12)
    u = F/np.hypot(*F)
    ax.annotate('', xy=0.9*u, xytext=(0, 0), arrowprops=PULL)
    ax.text(*(0.9*u + np.array([-0.55, -0.4])), r'$A_g$', color='C3', fontsize=11)
    if x == 'r':                                             #hahn13's angle phi; panel (b) shows the slope only
        ax.add_patch(matplotlib.patches.Arc((0, 0), 1.9, 1.9, theta1=np.degrees(np.arctan2(u[1], u[0])), theta2=0, lw=1))
        ax.text(1.0, -0.3, r'$\varphi$', fontsize=11)
    ax.text(dx/2, 0.12, rf'$d_{x}$', ha='center', fontsize=12)
    #a step along the streamline from the crossing: r dtheta along theta, dr' (or dz') outward
    V = 1.45*T
    ax.annotate('', xy=(dx + V[0], V[1]), xytext=(dx, 0), arrowprops=VEL)
    ax.plot([dx, dx], [0, V[1]], ':', color='C0', lw=1)
    ax.plot([dx, dx + V[0]], [V[1], V[1]], ':', color='C0', lw=1)
    ax.text(dx - 0.08, V[1]/2, r"$r\,d\theta$", ha='right', color='C0', fontsize=11)
    ax.text(dx + V[0]/2, V[1] + 0.1, f"$d{vx}'$", ha='center', color='C0', fontsize=11)
    if x == 'r':
        ax.add_patch(matplotlib.patches.Arc((dx, 0), 1.1, 1.1, theta1=90 - np.degrees(angle), theta2=90, lw=1))
        ax.text(dx + 0.06, 0.62, r'$\varphi$', fontsize=10)
    ax.plot(0, 0, 'ko', ms=9, zorder=5)
    ax.text(-0.15, 0.15, 'particle', ha='right', fontsize=10)
    ax.set_xlim(-1.3, 3.7)
    ax.set_ylim(-1.7, 2.8)
    ax.set_title(title, fontsize=11)


fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10.5, 3.8))
hahn13_panel(ax1, 'r', "r'", 'r', np.radians(25), r'(a) $r$–$\theta$ plane')
hahn13_panel(ax2, 'z', "z'", 'z', np.radians(25), r'(b) $z$–$\theta$ plane')

#(c) the r-z plane at the particle's longitude; the wire crosses it at (r', z'), into the page
ax3.set_aspect('equal')
ax3.axis('off')
ax3.annotate('', xy=(3.3, 0), xytext=(-1.0, 0), arrowprops=AXIS)
ax3.text(3.38, 0, '$r$', va='center', fontsize=13)
ax3.annotate('', xy=(0, 2.4), xytext=(0, -1.5), arrowprops=AXIS)
ax3.text(0, 2.5, '$z$', ha='center', fontsize=13)
W = np.array([2.0, 1.3])
ax3.plot([W[0], W[0]], [0, W[1]], 'k:', lw=1)
ax3.annotate('', xy=W, xytext=(0, 0), arrowprops=dict(arrowstyle='-|>', lw=1.2, mutation_scale=12, color='k', shrinkB=7))
ax3.text(W[0]/2, -0.35, r'$d_r$', ha='center', fontsize=12)
ax3.text(W[0] + 0.1, W[1]/2, r'$d_z$', va='center', fontsize=12)
ax3.text(W[0]*0.5 - 0.2, W[1]*0.5 + 0.12, r'$\mathbf{d}$', fontsize=12)
ax3.plot(*W, 'o', ms=13, mfc='white', mec='k', mew=1.5)
ax3.plot(*W, 'kx', ms=8, mew=1.5)
ax3.text(W[0] - 0.45, W[1] + 0.3, "$a'$ at $(r', z')$", fontsize=11)
ax3.plot(0, 0, 'ko', ms=9, zorder=5)
ax3.text(-0.15, 0.15, 'particle', ha='right', fontsize=10)
ax3.set_xlim(-1.3, 3.7)
ax3.set_ylim(-1.7, 2.8)
ax3.set_title(r"(c) $r$–$z$ plane at the particle's $\theta$", fontsize=11)

fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(HERE/f'wire_geometry.{ext}', dpi=130)
    print('wrote', HERE/f'wire_geometry.{ext}')
