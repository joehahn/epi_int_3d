#wire_planar.py
#
#hahn13 Fig. 3, redrawn: the gravity that a perturbing streamline a', treated as a straight wire,
#exerts on a particle in the ring plane (the paper's Appendix B; equations G3, G6, G14 of
#notes/sympy/streamline_gravity.py). The wire crosses the particle's longitude a distance d_r away,
#along the streamline's shape, tilted from theta-hat by phi, tan(phi) = s_r = dr'/(r dtheta); the
#perpendicular from the particle to it, Delta = d_r cos(phi), makes the angle phi with r-hat, and the
#wire pulls along it with A_g = 2 G lambda'/Delta. The streamline's velocity v' (dashed) is tilted by
#less than phi when its shape rotates (here as for an m = 2 pattern, v_r'/v_theta' = s_r/2).
#
#Run:  uv run python paper/figures/wire_planar.py   (writes wire_planar.pdf and .png alongside)

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
AXIS = dict(arrowstyle='->', lw=1)
PULL = dict(arrowstyle='-|>', lw=2.0, mutation_scale=14, color='C3')
STEP = dict(arrowstyle='-|>', lw=1.5, mutation_scale=12, color='C0')

fig, ax = plt.subplots(figsize=(4.6, 4.2))
ax.set_aspect('equal')
ax.axis('off')
#r to the right, theta up
ax.annotate('', xy=(3.3, 0), xytext=(-1.0, 0), arrowprops=AXIS)
ax.text(3.38, 0, '$r$', va='center', fontsize=13)
ax.annotate('', xy=(0, 2.4), xytext=(0, -1.5), arrowprops=AXIS)
ax.text(0, 2.5, r'$\theta$', ha='center', fontsize=13)
ax.set_xlim(-1.3, 3.7)
ax.set_ylim(-1.7, 2.8)

dx, phi = 2.0, np.radians(25)
T = np.array([np.sin(phi), np.cos(phi)])
s = np.array([-1.5, 2.2])
ax.plot(dx + s*T[0], s*T[1], 'k-', lw=1.6)
for q in (-1.4, 1.95):
    ax.plot(dx + q*T[0], q*T[1], 'ko', ms=6)
ax.text(dx + 2.25*T[0], 2.25*T[1], "$a'$", fontsize=12)
#the crossing point, and the radial gap d_r to it
ax.plot(dx, 0, 'o', ms=7, mfc='white', mec='k', mew=1.3, zorder=4)
ax.text(dx/2, 0.12, r'$d_r$', ha='center', fontsize=12)
#the perpendicular Delta = d_r cos(phi) to the near point, and the pull along it
F = np.array([dx, 0]) - (np.array([dx, 0]) @ T)*T
ax.plot([0, F[0]], [0, F[1]], 'k-', lw=1)
ax.plot(*F, 'ks', ms=4)
ax.text(0.45, -1.1, r'$\Delta = d_r\cos\varphi$', fontsize=11)
u = F/np.hypot(*F)
ax.annotate('', xy=0.85*u, xytext=(0, 0), arrowprops=PULL)
ax.text(0.2, -0.62, r'$A_g$', color='C3', fontsize=11)
ax.add_patch(matplotlib.patches.Arc((0, 0), 2.3, 2.3, theta1=np.degrees(np.arctan2(u[1], u[0])), theta2=0, lw=1))
ax.text(1.2, -0.22, r'$\varphi$', fontsize=11)
#the shape's slope at the crossing point: a step r dtheta moves the streamline dr'
V = 1.45*T
ax.annotate('', xy=(dx + V[0], V[1]), xytext=(dx, 0), arrowprops=STEP)
ax.plot([dx, dx], [0, V[1]], ':', color='C0', lw=1)
ax.plot([dx, dx + V[0]], [V[1], V[1]], ':', color='C0', lw=1)
ax.text(dx - 0.08, V[1]/2, r'$r\,d\theta$', ha='right', color='C0', fontsize=11)
ax.text(dx + V[0]/2, V[1] + 0.1, "$dr'$", ha='center', color='C0', fontsize=11)
ax.add_patch(matplotlib.patches.Arc((dx, 0), 1.1, 1.1, theta1=90 - np.degrees(phi), theta2=90, lw=1))
ax.text(dx + 0.06, 0.62, r'$\varphi$', fontsize=10)
#the velocity, tilted by less when the shape rotates
phv = np.arctan(0.5*np.tan(phi))
Vv = 1.0*np.array([np.sin(phv), np.cos(phv)])
ax.annotate('', xy=(dx + Vv[0], Vv[1]), xytext=(dx, 0),
            arrowprops=dict(arrowstyle='-|>', lw=1.3, ls='--', mutation_scale=11, color='C2'))
ax.text(dx + Vv[0] - 0.12, Vv[1] + 0.06, "$\\mathbf{v}'$", color='C2', fontsize=12)
ax.plot(0, 0, 'ko', ms=9, zorder=5)
ax.text(-0.15, 0.15, 'particle', ha='right', fontsize=10)

fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(HERE/f'wire_planar.{ext}', dpi=130)
    print('wrote', HERE/f'wire_planar.{ext}')
