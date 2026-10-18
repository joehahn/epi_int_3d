#!/usr/bin/env python3
"""Changes to the reference implementation (source/) that are applied at run time, leaving source/
unchanged. Used by explore/run.py and by tests/nominal/inputs.py.

- The streamline-crossing check wraps the ring in longitude first, because the original misfires
  near theta = +-pi for eccentric, closely spaced streamlines and stops the run (issue 1 in
  notes/epi_int_lite_issues.md). A run still stops on a real crossing or a NaN. Always applied.
- taper=f tapers the ringlet's edges: streamlines still have equal masses, but are placed at
  equal-mass quantiles of a surface density that is flat across the middle and falls as cos^2 to
  zero over a fraction f of the width at each edge, so the initial sigma has no step at the edges.
- fix_edge='own' gives the outer streamline its own one-sided pressure or viscous flux at its own
  longitudes, the mirror image of the inner edge (issue 4); fix_edge='hahn13' instead uses hahn13's
  p(r - Delta), the inner neighbor's flux at the outer particles' longitudes.

Usage, before source/epi_int_lite.py runs `from helper_fns import *`:
    import patches; patches.apply(taper=0.2)
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent/'source'))
import helper_fns
from helper_fns import wrap_ring, interpolate_fn, df_dr


def monitor_streamlines(monitor, r, t, timestep):
    """As helper_fns.monitor_streamlines, but with r and t wrapped in longitude first."""
    monitor['current_timestep'] = timestep
    if np.isnan(r).any() and monitor['nan_timestep'] is None:
        print('nan coordinate at timestep = ', timestep)
        monitor['nan_timestep'] = timestep
    rw = wrap_ring(r, longitude=False)
    tw = wrap_ring(t, longitude=True)
    dr = (rw - interpolate_fn(tw, rw, -1, interpolate=True))[1:, 1:-1]
    if (dr < 0).any() and monitor['streamline_crossing_timestep'] is None:
        print('steamlines cross at timestep = ', timestep)
        monitor['streamline_crossing_timestep'] = timestep
    return monitor


def tapered_spacing(n, f):
    """n points in [0, 1] at equal-mass quantiles (k + 1/2)/n of the tapered profile."""
    x = np.linspace(0.0, 1.0, 100001)
    edge = np.minimum(x, 1.0 - x)/f
    sigma = np.where(edge < 1.0, np.sin(0.5*np.pi*np.minimum(edge, 1.0))**2, 1.0)
    cdf = np.concatenate(([0.0], np.cumsum(0.5*(sigma[1:] + sigma[:-1])*np.diff(x))))
    q = (np.arange(n) + 0.5)/n
    return np.interp(q*cdf[-1], cdf, x)


def tapered(initialize, f):
    """initialize_streamline() with its uniform streamline spacing replaced by tapered_spacing()."""
    def wrapped(number_of_streamlines, particles_per_streamline, radial_width, *args, **kwargs):
        linspace = np.linspace
        def spacing(start, stop, num=50, endpoint=True, **kw):
            if num == number_of_streamlines and endpoint and start == 1.0:
                return start + (stop - start)*tapered_spacing(num, f)
            return linspace(start, stop, num=num, endpoint=endpoint, **kw)
        np.linspace = spacing
        try:
            return initialize(number_of_streamlines, particles_per_streamline, radial_width, *args, **kwargs)
        finally:
            np.linspace = linspace
    return wrapped


def A_P_hahn13(lambda_, sd, P, t, delta_P, delta_r):
    """As helper_fns.A_P, but the outer edge uses the inner neighbor's P at the outer particles' longitudes."""
    A = -df_dr(delta_P, delta_r)/sd
    A[0] = -P[0]/lambda_[0]
    A[-1] = interpolate_fn(t, P, -1)[-1]/lambda_[-1]
    return A


def A_P_own(lambda_, sd, P, t, delta_P, delta_r):
    """As helper_fns.A_P, but the outer edge uses its own (one-sided) P at its own longitudes,
    the mirror image of the inner edge's -P[0]/lambda."""
    A = -df_dr(delta_P, delta_r)/sd
    A[0] = -P[0]/lambda_[0]
    A[-1] = P[-1]/lambda_[-1]
    return A


def apply(taper=None, fix_edge=None):
    """Patch helper_fns: the wrapped crossing check always, and optionally an edge taper and edge fix."""
    helper_fns.monitor_streamlines = monitor_streamlines
    if fix_edge is not None:
        print('outer-edge flux fixed (issue 4):', fix_edge)
        helper_fns.A_P = {'own': A_P_own, 'hahn13': A_P_hahn13}[fix_edge]
    if taper is not None and not getattr(helper_fns.initialize_streamline, 'tapered', False):
        print('edge taper fraction =', taper)
        helper_fns.initialize_streamline = tapered(helper_fns.initialize_streamline, taper)
        helper_fns.initialize_streamline.tapered = True
