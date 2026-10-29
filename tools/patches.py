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
- fix_viscosity=True computes the viscous accelerations whenever either viscosity is on (issue 6):
  the reference skips bulk viscosity when shear_viscosity = 0 and pressure is off, and skips shear
  viscosity's radial stress when bulk_viscosity = 0. Runs with both viscosities on, or both off,
  are unchanged bit for bit.
- fix_seam=True computes pressure and viscosity on the ring wrapped with two ghost columns per side
  instead of one (issue 10). delta_f's three-point stencil takes a particle's longitude neighbours by
  rolling the row around without the 2 pi offset, so it is wrong in the outermost column, and each
  nested delta_f (P from delta_r, then delta_P, and A_P's edge interpolation of P) spreads that one
  column inward: with one ghost column the particles beside theta = +-pi got wrong viscous forces.
  Gravity still sees one ghost column, as its stencil is not nested. Particles away from the seam get
  the same accelerations bit for bit. With confine_inner_edge or confine_outer_edge, the edge torque
  is averaged over the real particles only, not also over the ghost columns as before.

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


def accelerations(lambda_, G_ring, shear_viscosity, bulk_viscosity, c, r, t, vr, vt, fast_gravity, confine_inner_edge, confine_outer_edge):
    """As helper_fns.accelerations, but sigma is computed when either viscosity is on, and the radial
    viscous stress (bulk viscosity plus shear viscosity's radial part) when either is on (issue 6)."""
    h = helper_fns
    rw = h.wrap_ring(r, longitude=False)
    tw = h.wrap_ring(t, longitude=True)
    vrw = h.wrap_ring(vr, longitude=False)
    vtw = h.wrap_ring(vt, longitude=False)
    vw = np.sqrt(vrw*vrw + vtw*vtw)
    lw = h.wrap_ring(lambda_, longitude=False)
    Ar = 0
    At = 0
    if (G_ring > 0.0):
        A = h.ring_gravity(lw, G_ring, rw, tw, vrw, vtw, fast_gravity)
        Ar += A[0]
        At += A[1]
    if ((c > 0.0) or (shear_viscosity > 0.0) or (bulk_viscosity > 0.0)):          #was: c or shear only
        delta_rw = h.delta_f(rw, tw)
        sdw = h.surface_density(lw, delta_rw)
        if (c > 0.0):
            A = h.ring_pressure(c, lw, sdw, rw, tw, vrw, vtw, vw, delta_rw)
            Ar += A[0]
            At += A[1]
        if (shear_viscosity > 0.0):
            At += h.ring_shear_viscosity(shear_viscosity, lw, sdw, rw, tw, vrw, vtw, vw, delta_rw)
        if ((shear_viscosity > 0.0) or (bulk_viscosity > 0.0)):                     #was: bulk only
            Ar += h.ring_bulk_viscosity(shear_viscosity, bulk_viscosity, lw, sdw, rw, tw, vrw, vtw, vw, delta_rw)
    if (confine_inner_edge or confine_outer_edge):
        At = h.edge_torques(rw, At, confine_inner_edge, confine_outer_edge)
    if (type(Ar) != int):
        Ar = Ar[:, 1:-1]
    if (type(At) != int):
        At = At[:, 1:-1]
    return Ar, At


def wrap(c, n, longitude=False):
    """helper_fns.wrap_ring with n ghost columns per side (n = 1 is wrap_ring, bit for bit)."""
    left = c[:, -n:].copy()
    right = c[:, :n].copy()
    if (longitude):
        twopi = 2*np.pi
        left -= twopi
        right += twopi
    return np.concatenate((left, c, right), axis=1)


def seam_accelerations(fix_viscosity):
    """accelerations (helper_fns' or, with fix_viscosity, this module's) with issue 10 fixed: gravity on
    the ring wrapped with one ghost column per side, pressure and viscosity with two."""
    def accelerations(lambda_, G_ring, shear_viscosity, bulk_viscosity, c, r, t, vr, vt, fast_gravity,
                      confine_inner_edge, confine_outer_edge):
        h = helper_fns
        rw = wrap(r, 2)
        tw = wrap(t, 2, longitude=True)
        vrw = wrap(vr, 2)
        vtw = wrap(vt, 2)
        lw = wrap(lambda_, 2)
        Ar = 0
        At = 0
        if (G_ring > 0.0):
            one = (slice(None), slice(1, -1))                 #the one-ghost wrap, as wrap_ring gives
            A = h.ring_gravity(np.ascontiguousarray(lw[one]), G_ring, np.ascontiguousarray(rw[one]),
                np.ascontiguousarray(tw[one]), np.ascontiguousarray(vrw[one]), np.ascontiguousarray(vtw[one]),
                fast_gravity)
            Ar += A[0][:, 1:-1]
            At += A[1][:, 1:-1]
        viscous = (shear_viscosity > 0.0) or (fix_viscosity and bulk_viscosity > 0.0)
        if ((c > 0.0) or viscous):
            vw = np.sqrt(vrw*vrw + vtw*vtw)
            delta_rw = h.delta_f(rw, tw)
            sdw = h.surface_density(lw, delta_rw)
            if (c > 0.0):
                A = h.ring_pressure(c, lw, sdw, rw, tw, vrw, vtw, vw, delta_rw)
                Ar += A[0][:, 2:-2]
                At += A[1][:, 2:-2]
            if (shear_viscosity > 0.0):
                At += h.ring_shear_viscosity(shear_viscosity, lw, sdw, rw, tw, vrw, vtw, vw, delta_rw)[:, 2:-2]
            if ((bulk_viscosity > 0.0) or (fix_viscosity and shear_viscosity > 0.0)):
                Ar += h.ring_bulk_viscosity(shear_viscosity, bulk_viscosity, lw, sdw, rw, tw, vrw, vtw, vw,
                    delta_rw)[:, 2:-2]
        if (confine_inner_edge or confine_outer_edge):
            At = h.edge_torques(r, At, confine_inner_edge, confine_outer_edge)
        return Ar, At
    accelerations.fix_seam = True
    return accelerations


def apply(taper=None, fix_edge=None, fix_viscosity=False, fix_seam=False):
    """Patch helper_fns: the wrapped crossing check always, and optionally an edge taper, edge fix,
    viscosity fix and seam fix."""
    helper_fns.monitor_streamlines = monitor_streamlines
    if fix_viscosity:
        print('viscosity gating fixed (issue 6)')
        helper_fns.accelerations = accelerations
    if fix_seam:
        apply_seam()
    if fix_edge is not None:
        print('outer-edge flux fixed (issue 4):', fix_edge)
        helper_fns.A_P = {'own': A_P_own, 'hahn13': A_P_hahn13}[fix_edge]
    if taper is not None and not getattr(helper_fns.initialize_streamline, 'tapered', False):
        print('edge taper fraction =', taper)
        helper_fns.initialize_streamline = tapered(helper_fns.initialize_streamline, taper)
        helper_fns.initialize_streamline.tapered = True


def apply_seam():
    """Fix issue 10 in helper_fns.accelerations, keeping issue 6's fix if it is applied; idempotent."""
    if getattr(helper_fns.accelerations, 'fix_seam', False):
        return
    print('seam fixed (issue 10)')
    helper_fns.accelerations = seam_accelerations(helper_fns.accelerations is accelerations)

