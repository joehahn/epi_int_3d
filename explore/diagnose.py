#!/usr/bin/env python3
"""Is a ringlet run well behaved? Time series of its shape and smoothness, from output/.

For each output: each streamline's mean radius a, eccentricity e and periapse w from the m = 1
Fourier term of r(theta); the ringlet's e' = a de/da and w' = a e dw/da (straight-line fits
across the streamlines) and q, the largest nonlinearity a |d(e exp(i w))|/da of any adjacent
pair of streamlines (they touch when q -> 1); the smallest gap between adjacent
streamlines over the initial spacing; the width a_outer - a_inner over its initial value; and the raggedness, the largest rms residual of a
streamline's r(theta) about a Fourier fit through m = 3, in units of the initial spacing.

Run:  uv run python explore/diagnose.py explore/ladder/nr0010 [more folders...]
Writes <folder>/diagnose.png and prints a one-line summary per folder.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO/'source'))
sys.path.insert(0, str(REPO/'benchmark'))
from helper_fns import restore_output, unwrap_angle, wrap_ring, interpolate_fn
from check_smoothness import fourier_residual, shape


def series(folder):
    r, t, vr, vt, times, monitor = restore_output(str(Path(folder)/'output'))
    Nr = r.shape[1]
    a0, _, _ = shape(r[0], t[0])
    spacing = (a0[-1] - a0[0])/(Nr - 1)
    out = {k: [] for k in ('orbits', 'e_prime', 'w_prime', 'q', 'e_mean', 'gap', 'ragged', 'width')}
    snaps = []
    for k in range(len(times)):
        if np.isnan(r[k]).any():
            break
        a, e, w = shape(r[k], t[k])
        w = np.unwrap(w)                   #continuous across streamlines (adjacent ones differ by << pi)
        w = w - w[len(w)//2]               #relative to the middle streamline
        abar, ebar = a.mean(), e.mean()
        e_prime = abar*np.polyfit(a, e, 1)[0]
        w_prime = abar*ebar*np.polyfit(a, w, 1)[0]
        rw = wrap_ring(r[k], longitude=False)
        tw = wrap_ring(t[k], longitude=True)
        gap = (rw - interpolate_fn(tw, rw, -1, interpolate=True))[1:, 1:-1].min()/spacing
        out['orbits'].append(times[k]/(2*np.pi))
        out['e_prime'].append(e_prime)
        out['w_prime'].append(w_prime)
        ev = e*np.exp(1j*w)               #eccentricity vectors; adjacent streamlines touch when q_pair -> 1
        out['q'].append((0.5*(a[1:] + a[:-1])*np.abs(np.diff(ev))/np.diff(a)).max())
        out['e_mean'].append(ebar)
        out['width'].append((a[-1] - a[0])/(a0[-1] - a0[0]))
        out['gap'].append(gap)
        out['ragged'].append(fourier_residual(t[k], r[k]).max()/spacing)
        snaps.append((a, e, w))
    return {k: np.array(v) for k, v in out.items()}, snaps, monitor, spacing, Nr


def main():
    for folder in sys.argv[1:]:
        s, snaps, monitor, spacing, Nr = series(folder)
        fig, ax = plt.subplots(2, 3, figsize=(16, 8))
        ax[0, 0].plot(s['orbits'], s['e_prime'], label="e'")
        ax[0, 0].plot(s['orbits'], s['w_prime'], label="w'")
        ax[0, 0].plot(s['orbits'], s['q'], label='q (max over pairs)')
        ax[0, 0].set_ylabel("e', w', q"); ax[0, 0].legend()
        ax[0, 1].plot(s['orbits'], s['gap'])
        ax[0, 1].set_ylabel('min gap / initial spacing')
        ax[0, 2].semilogy(s['orbits'], s['ragged'])
        ax[0, 2].set_ylabel('raggedness (spacings)')
        for k in sorted({0, len(snaps)//2, len(snaps) - 1}):
            a, e, w = snaps[k]
            x = (a - a.mean())/spacing
            ax[1, 0].plot(x, e, '.-', ms=3, label=f"{s['orbits'][k]:.0f} orbits")
            ax[1, 1].plot(x, w - w.mean(), '.-', ms=3, label=f"{s['orbits'][k]:.0f} orbits")
        ax[1, 0].set_ylabel('e'); ax[1, 1].set_ylabel('periapse - mean (rad)')
        for a_ in (ax[1, 0], ax[1, 1]):
            a_.set_xlabel('(a - mean a)/spacing'); a_.legend(fontsize=8)
        ax[1, 2].plot(s['orbits'], s['width'])
        ax[1, 2].set_ylabel('width / initial width')
        for a_ in (*ax[0], ax[1, 2]):
            a_.set_xlabel('orbits')
        stop = monitor['streamline_crossing_timestep'] or monitor['nan_timestep']
        title = (f"{folder}: Nr = {Nr}, {s['orbits'][-1]:.0f} orbits" +
                 (f', stopped at timestep {stop}' if stop else ''))
        fig.suptitle(title)
        fig.tight_layout()
        fig.savefig(Path(folder)/'diagnose.png', dpi=90)
        plt.close(fig)
        print(f"{folder}: Nr {Nr}, {s['orbits'][-1]:.0f} orbits, q {s['q'][0]:.3f} -> "
              f"{s['q'][-1]:.3f} (max {s['q'].max():.3f}), min gap {s['gap'].min():.3f}, "
              f"max raggedness {s['ragged'].max():.2g} spacings, width x{s['width'][-1]:.3f}, "
              f"crossing {monitor['streamline_crossing_timestep']}, nan {monitor['nan_timestep']}")


if __name__ == '__main__':
    main()
