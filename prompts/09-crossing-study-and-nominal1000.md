# 09 — Why 333 streamlines cross, and a 1000-streamline nominal ringlet

(Assembled from a run of prompts on 2026-10-19; it is the final version of what was asked.)

Tabulate the nominal-ringlet runs (Nr, λ_crit/2Δ, orbits, slow or fast code, crossed or not), and
confirm that the reference code crosses at 333 streamlines at the same step as the fast code.

Then find what drives the 333-streamline crossing, one suspect at a time, before running anything
that has not been discussed:

- ring pressure, Q_ring = 0.5, 1, 2, 10;
- particles per streamline, Nt = 20, 30, 40, 50, 60, 80, 160, 320, extending the runs that survive;
- ring gravity off, viscosity off, shear viscosity only, bulk viscosity only, bulk viscosity doubled,
  circular orbits;
- if a code fix is needed, show it first, then apply it outside `source/` and rerun.

Review what the viscous overstability literature says about avoiding or delaying it, and redo the
stability calculation with the code's own finite-difference viscosity operators. Sympy-check the
results and check the new citations.

Put "fix the tangential zigzag" on the to-do list, without analyzing or fixing it yet. Instead choose
nominal-ringlet parameters for which viscous overstability is a non-issue and the tangential zigzag does
not appear at 1000 streamlines, with viscosity and pressure off, since the first bending-wave runs will
have them off (vertical accelerations for pressure and viscosity come later). Screen candidates, pick
the best, run it long enough for J2 to matter (5×10⁴ orbits), and animate it with all relevant settings
listed beside the plots. It replaces the old nominal ringlet; also update the paper's resolution appendix.
