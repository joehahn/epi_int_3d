#!/usr/bin/env python3
"""Figure for the paper's appendix on the streamline method's resolution limit (zigzag.pdf).

(a) The zigzag mode: alternate streamlines displaced by +-xi, so the gaps alternate Delta +- 2 xi
    while r_{j+1} - r_{j-1} = 2 Delta stays fixed.
(b) Axisymmetric restoring and destabilizing terms versus phi = k Delta for a uniform lattice of
    streamlines, in units of kappa^2, for the nominal ringlet at 1000 streamlines
    (G sigma/kappa^2 Delta = 0.048): discrete self-gravity X (Eq. P7) versus the continuum
    2 pi G sigma k, and kappa^2 plus the discrete pressure term (Eq. P6) for Q = 1, 2.

Run:  uv run python paper/figures/zigzag.py   (writes paper/figures/zigzag.pdf)
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

g = 0.048                                   #G sigma/(kappa^2 Delta), nominal ringlet, 1000 streamlines
phi = np.linspace(0, np.pi, 400)
X = 4*g*(np.pi*phi/2 - phi**2/4)            #discrete self-gravity
X_cont = 2*np.pi*g*phi                      #continuum 2 pi G sigma k
def pressure(Q):                            #c = Q pi G sigma/kappa, so (c/Delta)^2/kappa^2 = (Q pi)^2 g^2
    return (Q*np.pi*g)**2*np.sin(phi)**2
def pressure_cont(Q):
    return (Q*np.pi*g)**2*phi**2

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={'width_ratios': [1, 1.6]})

#(a) zigzag schematic: radius increases to the right, each streamline a vertical segment
xi = 0.25
pos = lambda j: j + xi*(-1)**j                      #even streamlines move out, odd ones in
for j in range(-2, 3):
    ax0.plot([j, j], [0, 1], color='0.7', linestyle='--', linewidth=1)
    ax0.plot([pos(j), pos(j)], [0, 1], color='C0', linewidth=2)
    ax0.annotate('', xy=(pos(j), 0.9), xytext=(j, 0.9),
                 arrowprops=dict(arrowstyle='->', color='C3', linewidth=1.2, shrinkA=0, shrinkB=0))
    ax0.text(j, 1.03, rf'$j{j:+d}$' if j else r'$j$', ha='center', va='bottom', fontsize=7)
for j, y in ((-1, 0.62), (0, 0.42)):
    xa, xb = pos(j), pos(j + 1)
    ax0.annotate('', xy=(xb, y), xytext=(xa, y), arrowprops=dict(arrowstyle='<->', linewidth=0.8,
                 shrinkA=0, shrinkB=0))
    ax0.text((xa + xb)/2, y + 0.03, r'$\Delta+2\xi$' if j % 2 else r'$\Delta-2\xi$', ha='center',
             va='bottom', fontsize=8)
xa, xb = pos(-1), pos(1)
ax0.annotate('', xy=(xb, 0.15), xytext=(xa, 0.15), arrowprops=dict(arrowstyle='<->', linewidth=0.8,
             shrinkA=0, shrinkB=0))
ax0.text(0, 0.19, r'$r_{j+1}-r_{j-1}=2\Delta$', ha='center', va='bottom', fontsize=8)
ax0.set_xlim(-2.7, 2.7)
ax0.set_ylim(-0.05, 1.18)
ax0.set_xlabel(r'radius $r$ $\rightarrow$')
ax0.set_yticks([]); ax0.set_xticks([])
ax0.set_title('(a) zigzag mode', fontsize=9)

#(b) terms of omega^2 = kappa^2 + omega_p^2 - X
x = phi/np.pi
ax1.plot(x, X, color='C3', linewidth=2, label=r'self-gravity $X$')
ax1.plot(x, X_cont, color='C3', linestyle='--', linewidth=1, label=r'continuum $2\pi G\sigma k$')
for Q, ls in ((1, '-'), (2, '-.')):
    ax1.plot(x, 1 + pressure(Q), color='C0', linestyle=ls, linewidth=1.5,
             label=rf'$\kappa^2+\omega_p^2$, $Q={Q}$')
    ax1.plot(x, 1 + pressure_cont(Q), color='C0', linestyle=':', linewidth=1)
ax1.axhline(1, color='0.5', linewidth=0.8)
ax1.set_xlim(0, 1)
ax1.set_ylim(0, 1.4)
ax1.set_xlabel(r'$\phi/\pi = k\Delta/\pi$')
ax1.set_ylabel(r'units of $\kappa^2$')
ax1.legend(fontsize=7, loc='upper left')
ax1.set_title(r'(b) $G\sigma/\kappa^2\Delta = 0.048$', fontsize=9)

fig.tight_layout()
fig.savefig(Path(__file__).with_suffix('.pdf'))
fig.savefig(Path(__file__).with_suffix('.png'), dpi=120)
