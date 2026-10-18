# 07 — A nominal ringlet for benchmarking and testing

(Assembled from a run of prompts on 2026-10-18; it is the final version of what was asked.)

This project should have its own nominal ringlet for benchmarking and testing. Explore, starting from
epi_int_lite's nominal self-confining ringlet but letting it spread rather than self-confine:

- Initial `radial_width` = 0.01, viscosity chosen so that Δa spreads by about ×2 over the run, and
  `initial_orbits` eccentric with `e_prime` = 0 and a `w_prime` chosen so that the periapse rotates by π
  from the inner to the outer edge at t = 0, so that the ringlet's self-gravity is non-axisymmetric.
- Choose the mass, if possible, so that self-gravity dominates J2: enough to keep J2 from twisting the
  periapses so much that streamlines cross. If viscosity is too high the ringlet spreads a lot and then
  self-gravity can't stop J2 from twisting it further.
- 20–40 particles per streamline. Aim for 10⁴ orbits, shorter if the run times get excessive (lengthen
  it after the 3D code is parallelized).
- Ladder up from ~10 to 100 to 1000 streamlines, tweaking mass and viscosity so the ringlet stays
  coherent and stable. If it goes unstable, find out whether that is Toomre instability or numerical.
- Try tapered edges. Ring pressure is in hahn13: confirm its equations, confirm that the code
  encodes them correctly (z = 0), then try pressure on with Q = 0.5, 1.0, 1.2, 2.0. Check whether
  the outer-edge flux bug matters for the nominal ringlet. (Viscous spreading is already tested by
  `tests/spreading`; no further spreading tests are needed.)
- To check stability and coherence over time, plot a(t) per streamline as in
  `tests/spreading/spreading.ipynb` cell [9], plus e(a) and ϖ̃(a) versus time.
- Put the nominal ringlet in `tests/nominal/` with a notebook, try 300 streamlines, write it up, and
  commit. Don't change `source/`.
