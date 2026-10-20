#animate.py
#
#Animate runs of the nominal ringlet up to (and through) a streamline crossing: each streamline's
#eccentricity, periapse relative to streamline 1's (just exterior to the innermost) and raggedness, and
#q of each adjacent pair, versus distance from the middle streamline, one frame per output. The panels
#are those of tests/nominal/nominal.ipynb's animation, plus raggedness, which shows the short-wavelength
#zigzag that precedes a crossing. Writes <run>/figs/periapse.html (with play/pause) and periapse.gif.
#
#Run from explore/<study>:  uv run python ../animate.py run1 run2 ...

import os
import sys

REPO = os.path.abspath(os.path.dirname(__file__))
while not os.path.exists(os.path.join(REPO, 'pyproject.toml')):
    REPO = os.path.dirname(REPO)
sys.path.insert(0, os.path.join(REPO, 'source'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import animation
import re
from helper_fns import restore_output, Kappa


def shape(r, t):
    """Each streamline's mean radius a, eccentricity e and periapse w from a fit to r = a[1 - e cos(theta - w)]."""
    a, c1, s1 = np.array([np.linalg.lstsq(np.array([np.ones_like(tk), np.cos(tk), np.sin(tk)]).T, rk, rcond=None)[0]
        for tk, rk in zip(t, r)]).T
    return a, np.hypot(c1, s1)/a, np.unwrap(np.arctan2(-s1, -c1))


def raggedness(r, t, mmax=3):
    """Each streamline's rms residual of r(theta) about a Fourier fit through m = mmax."""
    res = []
    for tk, rk in zip(t, r):
        X = np.array([np.ones_like(tk)] + [f(m*tk) for m in range(1, mmax + 1) for f in (np.cos, np.sin)]).T
        res += [np.sqrt(((rk - X@np.linalg.lstsq(X, rk, rcond=None)[0])**2).mean())]
    return np.array(res)


def animate(folder):
    cwd = os.getcwd()
    os.chdir(folder)
    g = {}
    exec(open('inputs.py').read(), g)
    r, t, vr, vt, times, monitor = restore_output(g['output_folder'])
    g.update(monitor['modified_params'])                    #the run's -m overrides
    Q = g['Q_ring']
    taper = re.search(r'patches\.apply\(taper=([0-9.]+)', open('inputs.py').read())
    cross = monitor['streamline_crossing_timestep']
    steps_per_orbit = 2*np.pi/g['dt']
    orbits = times/(2*np.pi)
    N = np.isfinite(r).all(axis=(1, 2)).sum()
    A, E, W = np.array([shape(r[k], t[k]) for k in range(N)]).transpose(1, 0, 2)
    spacing = np.median(np.diff(A[0]))
    EV = E*np.exp(1j*W)
    q_pair = 0.5*(A[:, 1:] + A[:, :-1])*np.abs(np.diff(EV, axis=1))/np.diff(A, axis=1)
    R = np.array([raggedness(r[k], t[k]) for k in range(N)])/spacing
    Nr = r.shape[1]
    mid = Nr//2
    x = (A - A[:, mid:mid + 1])/(A[0, -1] - A[0, 0])
    xq = 0.5*(x[:, 1:] + x[:, :-1])
    e4 = E
    w_in = ((W - W[:, 1:2] + np.pi) % (2*np.pi) - np.pi)/np.pi
    nu_s, nu_b = g['shear_viscosity'], g['bulk_viscosity']
    tau_nu = (A[0, -1] - A[0, 0])**2/(12*nu_s)/(2*np.pi) if nu_s > 0 else None
    #resolution parameter lambda_crit/2 Delta in the middle of the ringlet, lambda_crit = 4 pi^2 G sigma/kappa^2
    a_mid = A[0, mid]
    sigma = g['total_ring_mass']/(Nr*2*np.pi*a_mid*spacing)
    res = 2*np.pi**2*max(g['G_ring'], 0.0)*sigma/(Kappa(g['J2'], g['Rp'], a_mid)**2*spacing)
    io = g['initial_orbits']
    lines = [f'{Nr} streamlines', f"{r.shape[2]} particles/streamline", f"dt = {g['dt']:g}", '',
             f"width {g['radial_width']:g}"] + ([f'edge taper {taper.group(1)}'] if taper else []) + [
             f"mass {g['total_ring_mass']:g}", f"ring gravity {'on' if g['G_ring'] > 0 else 'off'}", '',
             f"orbits: {io['shape']}"]
    if io['shape'] == 'eccentric':
        lines += [f"  e = {io['e']:g}", f"  a de/da = {io['e_prime']:g}", f"  a e dw/da = {io['w_prime']:g}"]
    lines += ['', f"J2 = {g['J2']:g}", f"Rp = {g['Rp']:g}", '',
              f'shear viscosity {nu_s:g}', f'bulk viscosity {nu_b:g}', 'pressure off' if Q < 0 else f'Q_ring = {Q:g}', '',
              rf'$\lambda_{{crit}}/2\Delta = {res:.2f}$', '  (mid-ringlet)']
    settings = '\n'.join(lines)

    plt.rcParams.update({'font.size': 15, 'animation.embed_limit': 500})     #MB; the default 20 drops frames
    f, axs = plt.subplots(2, 2, figsize=(19, 11))
    f.subplots_adjust(left=0.07, right=0.76, bottom=0.08, top=0.9, wspace=0.3, hspace=0.3)
    f.text(0.79, 0.88, settings, ha='left', va='top', fontsize=15, color='0.2', linespacing=1.4)
    panels = [(axs[0, 0], x, e4, r'eccentricity $e$', 'log'),
              (axs[0, 1], x, w_in, r'$(\tilde{\omega} - \tilde{\omega}_1)/\pi$', 'linear'),
              (axs[1, 0], xq, q_pair, r'$q$ of adjacent pair (cross at 1)', 'log'),
              (axs[1, 1], x, R, 'raggedness/spacing', 'log')]
    artists = []
    for ax, X, Y, label, scale in panels:
        line, = ax.plot(X[0], Y[0], color='0.5', linewidth=1)
        dots = ax.scatter(X[0], Y[0], color='C0', s=10, zorder=3)
        ax.set_xlim(x.min(), x.max())
        ax.set_yscale(scale)
        if scale == 'log':
            ax.set_ylim(0.5*np.percentile(Y[Y > 0], 1), 2*Y.max())     #from the quiet start to the blow-up
        ax.set_xlabel(r'$(a - a_{mid})/\Delta a_0$')
        ax.set_ylabel(label)
        artists += [(line, dots, X, Y)]
    axs[0, 1].set_ylim(-1, 1)
    axs[1, 0].axhline(1, color='C3', linewidth=0.8)
    crossing = 'no crossing' if cross is None else f'crossing at step {cross} ({cross/steps_per_orbit:.0f} orbits)'
    title = f.suptitle('', x=0.42, y=0.965, fontsize=18)

    def update(k):
        for line, dots, X, Y in artists:
            line.set_data(X[k], Y[k])
            dots.set_offsets(np.c_[X[k], Y[k]])
        ex = int(np.floor(np.log10(orbits[k]))) if orbits[k] > 0 else 0
        when = rf'$t = {orbits[k]/10**ex:.2f}\times10^{{{ex}}}$ orbits' + (rf' $= {orbits[k]/tau_nu:.4f}\,\tau_\nu$' if tau_nu else '')
        title.set_text(f'{when}, step {int(times[k]/g["dt"] + 0.5)},  {crossing}')
        return [a for line, dots, X, Y in artists for a in (line, dots)] + [title]

    os.makedirs('figs', exist_ok=True)
    anim = animation.FuncAnimation(f, update, frames=range(N))
    anim.save('figs/periapse.gif', writer='pillow', fps=10, dpi=60)
    anim.save('figs/periapse.html', writer=animation.HTMLWriter(fps=10, embed_frames=True, default_mode='loop'), dpi=60)
    plt.close(f)
    print(f'{folder}: {N} frames, figs/periapse.html')
    os.chdir(cwd)


if __name__ == '__main__':
    for folder in sys.argv[1:]:
        animate(folder)
