# Tie-break: code

2026-11-08 · 2 disputed findings · verdicts `reviews/2026-11-08-code-verdicts.json` · cost $0.0426

## gpt-1: claude (judge gemini)

- **Finding (gpt):** The crossing check's 3-point Lagrange stencil on the inner streamline is chosen by particle index, not by longitude, so if neighboring streamlines' particles drift apart in longitude it extrapolates and can report a false crossing.
- **Claude (partly):** Right in principle, and the same index-based stencil is epi_int_lite's own design (source/helper_fns.py interpolate_fn), used also by delta_f and so by the viscous and pressure forces; the kernel must match the reference bit for bit. Measured in our runs, same-index particles on neighboring streamlines stay within 1.37 particle spacings (self_confining, 135,000 orbits) and 0.78 (bending_waves), so the stencil brackets or nearly brackets the target and the smooth streamline makes the mild extrapolation harmless. Not changing the kernel; logging the assumption as epi_int_lite issue 11 and adding a check of the offset to the regression tests.
- **Reason:** The Numba kernel in `fast/bookkeeping.py` is explicitly designed to match the reference implementation (`source/helper_fns.py`) bit-for-bit, which uses the same index-based stencil via `interpolate_fn`. Furthermore, because particles on each streamline are sorted by longitude at every velocity kick (`sort_streamlines`), the j-th particle on any streamline stays near the same longitude across streamlines, making Claude's decision to document the assumption and add a regression check rather than alter the reference-matching kernel appropriate.
- **Fix:** 

## gpt-2: claude (judge gemini)

- **Finding (gpt):** _adjust shifts by 2 pi only once, so an angle outside [-3 pi, 3 pi] is not brought into [-pi, pi]; e.g. dt = 4 pi leaves t near 2 pi, breaking the longitude convention that sorting and wrapping need.
- **Claude (disagree):** It is helper_fns.adjust_angle exactly (one shift each way), which the kernel must match bit for bit. Its argument cannot leave [-3 pi, 3 pi] in use: M from atan2 is in [-pi, pi], wt in [-pi, pi], Omega_a/kappa_a = 1 + O(J2 (Rp/a)^2), e << 1, and kappa dt is a small fraction of an orbit (timesteps are ~1/100 of an orbit or less), so the sum is below 2 pi + small. The dt = 4 pi example is outside the code's use.
- **Reason:** `_adjust` is an exact scalar translation of `helper_fns.adjust_angle` required for bit-for-bit equivalence with the reference code. In all call sites within `drift`, the inputs are combinations of bounded angles in `[-pi, pi]` and a single sub-orbital timestep `Kap * dt`, so the argument never exceeds `[-3*pi, 3*pi]` in valid physical simulations.
- **Fix:** 

