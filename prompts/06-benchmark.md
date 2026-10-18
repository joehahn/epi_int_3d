# 06 — Benchmark the reference implementation

Start Phase 0b of `PLAN.md` with step 1: measure how fast the Python 3 port in `source/` is and
where its time goes, before any of it is made faster.

- Write a benchmark that times the reference implementation on a few problem sizes: Nr = 10,
  100 and 1000 streamlines, with Nt = 20 particles per streamline (the bending-wave target is
  Nr ≈ 400–800, Nt ≈ 20). Use a self-gravitating, viscous ringlet with the wide 10⁻⁹ libration
  ringlet's surface density and circular streamlines (`initial_orbits` `'circular'`; if that
  makes NaNs, lower the ring mass); `G_ring = 1`, `fast_gravity = False`, shear and bulk
  viscosity 1e-13, `Q_ring = -1` (pressure off), dt = 0.5. Keep the streamline spacing and surface density fixed
  as Nr grows: `radial_width` = 1e-4, 1e-3, 1e-2 and `total_ring_mass` = 1e-10, 1e-9, 1e-8.
  Run each size for a few hundred timesteps, enough to get a stable time per step, and check
  that the ringlet stays well behaved: no NaNs or crossing streamlines, and r(θ), velocities and
  orbit elements that vary smoothly.
- Profile those runs (cProfile, and line_profiler on the hot functions) and report which
  functions take the time and how each cost scales with Nr and Nt. Check `PLAN.md`'s guess that
  `ring_gravity` dominates and grows as Nr² × Nt.
- From the measured scaling, estimate the wall-clock time of one bending-wave run of 10⁵–10⁶
  timesteps with the current code.
- Save the benchmark so it can be rerun unchanged on each faster version, and record this
  laptop's results (CPU, cores, Python and NumPy versions) as the baseline for the speedup table
  in the README. Also save each size's final state, so faster versions can be checked against
  the reference to ~1e-12.
- Don't change `source/`. Write up the results in `notes/phase0b/` and show me them, with a
  recommendation for what to speed up first.
