#wire_geometry.py
#
#The 3D analog of hahn13 Fig. 3: the geometry of the gravity that a perturbing streamline, treated
#as a straight wire, exerts on a particle that lies off the streamline's plane
#(the paper's Section 3 and Figure 1; notes/phase2/streamline_gravity.md, equations G1-G3).
#
#Run:  uv run python paper/figures/wire_geometry.py   (writes wire_geometry.pdf and .png alongside)

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ARROW = dict(arrowstyle='-|>', lw=1.6, mutation_scale=14)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
for ax in (ax1, ax2):
    ax.set_aspect('equal')
    ax.axis('off')

#(a) the particle's meridional plane (fixed longitude theta), seen along theta
P = np.array([0.0, 0.0])                     #particle
W = np.array([2.4, 1.3])                     #where the perturbing streamline crosses this plane
ax1.annotate('', xy=(3.6, -1.2), xytext=(-2.9, -1.2), arrowprops=dict(arrowstyle='->', lw=1))
ax1.text(3.65, -1.2, 'r', va='center', fontsize=13)
ax1.annotate('', xy=(-2.6, 2.3), xytext=(-2.6, -1.5), arrowprops=dict(arrowstyle='->', lw=1))
ax1.text(-2.6, 2.4, 'z', ha='center', fontsize=13)
ax1.plot([-2.9, 3.4], [-1.2, -1.2], alpha=0)
ax1.text(-2.45, -1.5, 'midplane z = 0 lies below; planet $M_p$ far to the left', fontsize=8.5, color='0.4')
ax1.plot(*P, 'ko', ms=9)
ax1.text(P[0] - 0.15, P[1] - 0.35, 'particle (r, z)', ha='center', fontsize=10)
ax1.plot(*W, 'o', ms=13, mfc='white', mec='k', mew=1.5)
ax1.plot(*W, 'kx', ms=8, mew=1.5)
ax1.text(W[0] + 0.25, W[1] + 0.2, "streamline $a'$ at $(r', z')$,\nheading into the page", fontsize=10)
ax1.plot([P[0], W[0]], [P[1], P[1]], 'k--', lw=1)
ax1.plot([W[0], W[0]], [P[1], W[1]], 'k--', lw=1)
ax1.plot([P[0], W[0]], [P[1], W[1]], 'k-', lw=1)
ax1.text(W[0]/2, -0.3, r"$\Delta_r = r' - r$", ha='center', fontsize=11)
ax1.text(W[0] + 0.1, W[1]/2, r"$\Delta_z = z' - z$", va='center', fontsize=11)
ax1.text(W[0]*0.68 + 0.05, W[1]*0.68 - 0.32, r'$\delta$', fontsize=12)
u = (W - P)/np.hypot(*(W - P))
ax1.annotate('', xy=P + 1.35*u, xytext=P, arrowprops=dict(color='C3', **ARROW))
ax1.text(*(P + 1.35*u + np.array([-0.75, 0.15])), r'$A_g$', color='C3', fontsize=12)
psi = np.degrees(np.arctan2(*(W - P)[::-1]))
ax1.add_patch(matplotlib.patches.Arc(P, 1.4, 1.4, theta1=0, theta2=psi, lw=1))
ax1.text(0.78, 0.17, r'$\psi$', fontsize=12)
ax1.set_xlim(-3.0, 5.6)
ax1.set_ylim(-1.8, 3.4)
ax1.set_title('(a) meridional plane through the particle', fontsize=11)

#(b) the same seen from outside the ring, looking inward along -r: theta to the right, z up.
#The streamline crosses the particle's longitude a height dz above it, tilted by chi.
ax2.annotate('', xy=(3.3, -1.2), xytext=(-3.0, -1.2), arrowprops=dict(arrowstyle='->', lw=1))
ax2.text(3.35, -1.2, r'$\theta$', va='center', fontsize=13)
ax2.annotate('', xy=(-2.7, 2.3), xytext=(-2.7, -1.5), arrowprops=dict(arrowstyle='->', lw=1))
ax2.text(-2.7, 2.4, 'z', ha='center', fontsize=13)
P = np.array([0.0, 0.0])
dz = 1.9
chi = np.radians(25)
T = np.array([np.cos(chi), np.sin(chi)])     #the wire's direction, slope v_z'/v_theta'
C = np.array([0.0, dz])                      #crossing point, at the particle's longitude
s = np.linspace(-2.6, 2.4, 2)
ax2.plot(C[0] + s*T[0], C[1] + s*T[1], 'k-', lw=1.5)
for q in (-1.9, 0.5, 1.8):
    ax2.plot(*(C + q*T), 'ko', ms=7)
ax2.text(*(C + 2.45*T + np.array([0.05, -0.05])), "$a'$", fontsize=12)
ax2.plot(*P, 'ko', ms=9)
ax2.text(P[0], P[1] - 0.4, 'particle', ha='center', fontsize=10)
F = C + ((P - C) @ T)*T                      #foot of the perpendicular from the particle
ax2.plot([P[0], C[0]], [P[1], C[1]], 'k--', lw=1)
ax2.text(C[0] + 0.1, dz*0.35, r'$\Delta_z$', ha='left', va='center', fontsize=11)
ax2.plot([P[0], F[0]], [P[1], F[1]], 'k-', lw=1)
ax2.text(*(P + 0.8*(F - P) + np.array([-0.5, 0.0])), r'$\Delta_\perp$', fontsize=11)
u = (F - P)/np.hypot(*(F - P))
ax2.annotate('', xy=P + 0.7*u, xytext=P, arrowprops=dict(color='C3', **ARROW))
ax2.text(*(P + 0.35*u + np.array([-0.5, -0.1])), r'$A_g$', color='C3', fontsize=12)
#angle chi between the vertical and the perpendicular, at the particle, and the wire's tilt
ax2.add_patch(matplotlib.patches.Arc(P, 1.5, 1.5, theta1=90, theta2=np.degrees(np.arctan2(u[1], u[0])), lw=1))
ax2.text(-0.12, 0.95, r'$\chi$', fontsize=12)
ax2.plot([C[0] + 0.0, C[0] + 1.4], [dz, dz], 'k:', lw=1)
ax2.add_patch(matplotlib.patches.Arc(C, 2.2, 2.2, theta1=0, theta2=np.degrees(chi), lw=1))
ax2.text(C[0] + 1.18, dz + 0.17, r'$\chi$', fontsize=12)
ax2.text(*(C + 1.6*T + np.array([0.1, 0.25])), r"$v'$", fontsize=11)
ax2.set_xlim(-3.1, 4.3)
ax2.set_ylim(-1.8, 3.4)
ax2.set_title(r'(b) seen from outside the ring, along $-\hat r$', fontsize=11)

fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(HERE/f'wire_geometry.{ext}', dpi=130)
    print('wrote', HERE/f'wire_geometry.{ext}')
