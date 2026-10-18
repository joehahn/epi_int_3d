#!/usr/bin/env python3
"""Check that the benchmark ringlets stayed well behaved: smooth r(theta), velocities and orbit
elements, rather than the ragged particle distribution of an unstable ring.

For each benchmark run in benchmark/runs/ (the n2 run of each case), this compares the initial
and final states:
- along each streamline, r, vr and vt are fit by a Fourier series in theta through m = 3 (an
  eccentric streamline is m = 1 plus small harmonics); the rms residual measures raggedness, in
  units of the streamline spacing (r) or of n times the spacing (vr, vt);
- the scatter of each particle's a and e about its streamline's mean;
- the streamlines' e and periapse longitude (from r(theta)) versus a, and how smooth they are: the largest
  second difference between adjacent streamlines over the typical first difference.
It also plots e(a), periapse(a) and r(theta) residuals to notes/phase0b/figs/.

Run:  uv run python benchmark/check_smoothness.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BENCH = Path(__file__).resolve().parent
REPO = BENCH.parent
sys.path.insert(0, str(REPO/'source'))
sys.path.insert(0, str(BENCH))
from helper_fns import restore_output, coords2elem, unwrap_angle
from benchmark import CASES, SPACING, case_name

J2, Rp = 0.01, 0.5
FIGS = REPO/'notes'/'phase0b'/'figs'


def fourier_residual(t, f, mmax=3):
    """rms residual of a least-squares Fourier fit of f(t) through mmax, per streamline."""
    res = []
    for ti, fi in zip(t, f):
        cols = [np.ones_like(ti)]
        for m in range(1, mmax + 1):
            cols += [np.cos(m*ti), np.sin(m*ti)]
        X = np.array(cols).T
        coef, *_ = np.linalg.lstsq(X, fi, rcond=None)
        res.append(np.sqrt(np.mean((fi - X@coef)**2)))
    return np.array(res)


def shape(r, t):
    """Each streamline's mean radius, eccentricity and periapse longitude from the m = 1 Fourier
    term of r(theta) = a[1 - e cos(theta - w)]. This is geometric; averaging the particles'
    osculating periapses instead jumps whenever a particle's mean anomaly wraps past +-pi."""
    a0, c1, s1 = [], [], []
    for ti, ri in zip(t, r):
        X = np.array([np.ones_like(ti), np.cos(ti), np.sin(ti)]).T
        coef, *_ = np.linalg.lstsq(X, ri, rcond=None)
        a0.append(coef[0]); c1.append(coef[1]); s1.append(coef[2])
    a0, c1, s1 = np.array(a0), np.array(c1), np.array(s1)
    return a0, np.hypot(c1, s1)/a0, np.arctan2(-s1, -c1)


def roughness(y):
    """Largest |second difference| over the median |first difference| across streamlines."""
    d1 = np.abs(np.diff(y))
    d2 = np.abs(np.diff(y, 2))
    return d2.max()/np.median(d1) if len(y) > 2 else 0.0


def metrics(r, t, vr, vt):
    a, e, wt, M = coords2elem(J2, Rp, r, t, vr, vt)
    a_mean, e_mean, w = shape(r, t)
    w = unwrap_angle(w)
    speed = a_mean.mean()**-0.5*SPACING              #n times the streamline spacing, with a ~ 1
    return {'r_resid': fourier_residual(t, r).max()/SPACING,
            'vr_resid': fourier_residual(t, vr).max()/speed,
            'vt_resid': fourier_residual(t, vt).max()/speed,
            'a_scatter': (a - a.mean(axis=1)[:, None]).std(axis=1).max()/SPACING,
            'e_scatter': (e - e.mean(axis=1)[:, None]).std(axis=1).max(),
            'e_rough': roughness(e_mean), 'w_rough': roughness(w),
            'nan': bool(np.isnan(r).any() or np.isnan(vr).any())}, (a_mean, e_mean, w, a, e)


def main():
    FIGS.mkdir(parents=True, exist_ok=True)
    keys = ['r_resid', 'vr_resid', 'vt_resid', 'a_scatter', 'e_scatter', 'e_rough', 'w_rough']
    print(f"{'case':>13} {'state':>7} " + ' '.join(f'{k:>10}' for k in keys) + '   nan')
    for Nr, Nt, n1, n2 in CASES:
        folder = BENCH/'runs'/f'{case_name(Nr, Nt)}_{n2}'/'output'
        r, t, vr, vt, times, monitor = restore_output(str(folder))
        curves = []
        for k, label in ((0, 'initial'), (-1, 'final')):
            m, curve = metrics(r[k], t[k], vr[k], vt[k])
            curves.append(curve)
            print(f'{case_name(Nr, Nt):>13} {label:>7} ' +
                  ' '.join(f'{m[key]:10.3g}' for key in keys) + f"   {m['nan']}")
        if Nt == 20 and curves[0][1].mean() > 1e-6:     #e(a) and periapse(a) mean nothing for circular orbits
            plot(Nr, Nt, n2, r, t, curves)
    print(f'figures written to {FIGS.relative_to(REPO)}/')


def plot(Nr, Nt, n2, r, t, curves):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    for (a_mean, e_mean, w, a, e), label, style in zip(curves, ('initial', f'after {n2} steps'),
                                                        ('o', '.')):
        x = (a_mean - a_mean[0])/SPACING
        axes[0].plot(x, e_mean, style, ms=3, label=label)
        axes[1].plot(x, w - w.mean(), style, ms=3, label=label)
    axes[0].set_xlabel('(a - a_inner)/spacing'); axes[0].set_ylabel('streamline e')
    axes[1].set_xlabel('(a - a_inner)/spacing'); axes[1].set_ylabel('periapse - mean (rad)')
    axes[0].legend(); axes[1].legend()
    #r(theta) relative to each streamline's Fourier fit, for a few streamlines
    for i in np.linspace(0, Nr - 1, 5).astype(int):
        ti, ri = t[-1][i], r[-1][i]
        cols = [np.ones_like(ti)] + [f(m*ti) for m in (1, 2, 3) for f in (np.cos, np.sin)]
        X = np.array(cols).T
        coef, *_ = np.linalg.lstsq(X, ri, rcond=None)
        o = ti.argsort()
        axes[2].plot(ti[o], (ri - X@coef)[o]/SPACING, '.-', label=f'streamline {i}')
    axes[2].set_xlabel('longitude (rad)'); axes[2].set_ylabel('r - Fourier fit (spacings)')
    axes[2].legend(fontsize=8)
    fig.suptitle(f'benchmark ringlet, Nr = {Nr}, Nt = {Nt}')
    fig.tight_layout()
    fig.savefig(FIGS/f'smoothness_{case_name(Nr, Nt)}.png', dpi=100)
    plt.close(fig)


if __name__ == '__main__':
    main()
