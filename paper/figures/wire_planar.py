#wire_planar.py
#
#hahn13 Fig. 3, redrawn, in two figures (equations G3, G6, G14, G16 of notes/sympy/streamline_gravity.py):
#wire_planar: the paper's Figure 1 (Section 3.1). The perturbing streamline a' crosses the particle's
#    longitude at d0 = d_r r-hat (the crossing point, x), tilted from theta-hat by phi; the
#    near point (black square) is d1 further along it, and Delta = d0 + d1 is perpendicular to the
#    streamline, with length d_r cos(phi) and at the angle phi to r-hat; the wire pulls along it.
#    r and r' come from the planet, far off to the left.
#    Vectors are drawn heavier than the streamline.
#wire_velocity: the paper's Appendix B figure. At the crossing point the streamline's shape has the slope
#    tan(phi) = s_r = dr'/(r dtheta), while its velocity v' (dashed) is tilted less when the shape
#    rotates (here as for an m = 2 pattern, v_r'/v_theta' = s_r/2).
#
#Run:  uv run python paper/figures/wire_planar.py   (writes wire_planar and wire_velocity .pdf and .png)

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
AXIS = dict(arrowstyle='->', lw=1)
VEC = dict(arrowstyle='-|>', lw=1.8, mutation_scale=14, color='k', shrinkA=0, shrinkB=0)
PULL = dict(arrowstyle='-|>', lw=2.4, mutation_scale=16, color='C3', shrinkA=0, shrinkB=0)
STEP = dict(arrowstyle='-|>', lw=1.5, mutation_scale=12, color='C0')
dx, phi = 2.0, np.radians(25)
T = np.array([np.sin(phi), np.cos(phi)])              #the streamline's direction, tilted phi from theta-hat


def save(fig, name):
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(HERE/f'{name}.{ext}', dpi=130)
        print('wrote', HERE/f'{name}.{ext}')


#Figure 1: d0, d1 and Delta
fig, ax = plt.subplots(figsize=(4.6, 4.0))
ax.set_aspect('equal')
ax.axis('off')
ax.annotate('', xy=(3.3, 0), xytext=(0, 0), arrowprops=AXIS)
ax.text(3.38, 0, r'$\hat{r}$', va='center', fontsize=13)
ax.annotate('', xy=(0, 2.2), xytext=(0, -1.5), arrowprops=AXIS)
ax.text(0, 2.3, r'$\hat{\theta}$', ha='center', fontsize=13)
ax.set_xlim(-1.35, 3.7)
ax.set_ylim(-1.7, 2.6)
X = np.array([dx, 0.0])                               #the crossing point
F = X - (X @ T)*T                                     #the near point
s = np.array([-1.5, 2.2])
ax.plot(dx + s*T[0], s*T[1], 'k-', lw=1.0)
for q in (-1.4, 0.85, 1.95):                          #the three particles of a' nearest in longitude
    ax.plot(dx + q*T[0], q*T[1], 'ko', ms=6)
ax.text(dx + 2.25*T[0], 2.25*T[1], "$a'$", fontsize=12)
ax.annotate('', xy=X, xytext=(0, 0), arrowprops=dict(VEC, shrinkB=5))
ax.text(dx/2, 0.12, r'$\mathbf{d}_0$', ha='center', fontsize=10)
ax.annotate('', xy=F, xytext=X, arrowprops=dict(VEC, shrinkA=5, shrinkB=6))
ax.text(*((X + F)/2 + 0.15*np.array([T[1], -T[0]])), r'$\mathbf{d}_1$', fontsize=10, va='center')
ax.annotate('', xy=F, xytext=(0, 0), arrowprops=dict(VEC, shrinkB=6))
u = F/np.hypot(*F)
ax.text(*(0.5*F + 0.14*np.array([u[1], -u[0]])), r'$\boldsymbol{\Delta}$', fontsize=10, ha='center', va='center')
ax.add_patch(matplotlib.patches.Arc((0, 0), 2.3, 2.3, theta1=np.degrees(np.arctan2(u[1], u[0])), theta2=0, lw=1))
am = np.arctan2(u[1], u[0])/2                         #midway along the arc
ax.text(1.27*np.cos(am), 1.27*np.sin(am), r'$\varphi$', fontsize=11, ha='center', va='center')
ax.plot(*X, 'kx', ms=9, mew=2, zorder=5)
ax.text(dx + 0.25, 0.1, 'crossing\npoint', fontsize=9, va='bottom')
ax.plot(*F, 'ks', ms=9, zorder=5)
ax.text(F[0] + 0.12, F[1] - 0.05, 'near\npoint', fontsize=9, va='top')
ax.plot(0, 0, 'ko', ms=9, zorder=5)
#r and r', from the planet far off to the left, to the particle and the near point
RV = dict(arrowstyle='-|>', lw=1.2, mutation_scale=12, color='k', shrinkA=0, shrinkB=6)
ax.annotate('', xy=(0, 0), xytext=(-1.25, 0), arrowprops=RV)
ax.text(-0.95, 0.07, r'$\mathbf{r}$', fontsize=10, va='bottom')
ax.annotate('', xy=F, xytext=(-1.25, F[1]), arrowprops=RV)
ax.text(-0.95, F[1] + 0.07, r"$\mathbf{r}'$", fontsize=10, va='bottom')
#the near point's longitude offset, r dtheta' along theta-hat (negative here)
ax.annotate('', xy=(-0.45, F[1]), xytext=(-0.45, 0), arrowprops=dict(arrowstyle='<->', lw=0.8, shrinkA=0, shrinkB=0))
ax.text(-0.5, F[1]/2, r"$r\,\delta\theta'$", fontsize=10, ha='right', va='center')
ax.text(-0.12, 0.3, 'particle', ha='right', fontsize=10)
save(fig, 'wire_planar')

#Appendix B: the shape's slope and the velocity at the crossing point
fig, ax = plt.subplots(figsize=(2.6, 3.0))
ax.set_aspect('equal')
ax.axis('off')
ax.set_xlim(-0.6, 1.4)
ax.set_ylim(-0.5, 2.0)
s = np.array([-0.45, 1.85])
ax.plot(s*T[0], s*T[1], 'k-', lw=1.6)
ax.text(1.88*T[0] + 0.05, 1.88*T[1], "$a'$", fontsize=12)
ax.plot([0, 0], [-0.4, 1.8], ':', color='0.5', lw=1)
ax.text(-0.05, 1.8, r'$\hat{\theta}$', ha='right', color='0.4', fontsize=11)
V = 1.45*T
ax.annotate('', xy=V, xytext=(0, 0), arrowprops=STEP)
ax.plot([0, V[0]], [V[1], V[1]], ':', color='C0', lw=1)
ax.text(-0.06, V[1]/2, r'$r\,d\theta$', ha='right', color='C0', fontsize=11)
ax.text(V[0]/2, V[1] + 0.07, "$dr'$", ha='center', color='C0', fontsize=11)
ax.add_patch(matplotlib.patches.Arc((0, 0), 1.0, 1.0, theta1=90 - np.degrees(phi), theta2=90, lw=1))
ax.text(0.07, 0.56, r'$\varphi$', fontsize=10)
phv = np.arctan(0.5*np.tan(phi))
Vv = 1.0*np.array([np.sin(phv), np.cos(phv)])
ax.annotate('', xy=Vv, xytext=(0, 0), arrowprops=dict(arrowstyle='-|>', lw=1.3, ls='--', mutation_scale=11, color='C2'))
ax.text(Vv[0] + 0.03, Vv[1] - 0.1, "$\\mathbf{v}'$", color='C2', fontsize=12)
ax.plot(0, 0, 'kx', ms=9, mew=2, zorder=5)
ax.text(0.1, -0.12, 'crossing point', fontsize=9, va='top')
save(fig, 'wire_velocity')
