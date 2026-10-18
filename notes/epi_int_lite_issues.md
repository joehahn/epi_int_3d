# Issues found in epi_int_lite

Problems found in epi_int_lite (master, commit 5e0942a) while porting and extending it. `source/`
stays a straight port of epi_int_lite, so these are fixed in the faster code and epi_int_3d,
not there, unless noted. Changes made by the Python 3 port itself are in `port_changes.md`.

## 1. The streamline-crossing check misfires near θ = ±π

- **Found:** Phase 0b benchmark (`notes/phase0b/benchmark.md`), 2026-10-18.
- **What:** `monitor_streamlines()` (`source/helper_fns.py:423`) estimates each streamline's gap
  to its inner neighbor with `interpolate_fn()` on arrays that are not wrapped in longitude. For
  a streamline's first and last particle, `sidestep()` borrows a point from the other end of the
  array without adding ±2π to its longitude, so the 3-point Lagrange fit spans
  θ ≈ (2.8, 3.1, −3.1) and extrapolates wildly. The check then reports crossing streamlines
  that have not crossed, and `update_display()` stops a self-interacting run.
- **When it bites:** the error grows with eccentricity and with the longitude gap between
  particles, and is compared with the streamline spacing. In the benchmark ringlet (spacing
  10⁻⁵, e = 0.005, 10–40 particles per streamline) it fired after 7–143 steps in every case,
  while the true smallest gap stayed at 0.5–0.7 of the spacing. Bending-wave runs are in the
  same regime.
- **Effect on hahn25:** none found. None of the 21 self-interacting test runs (`tests/spreading`,
  the 20 libration runs) was flagged. The forces are unaffected, since `accelerations()` wraps the
  ring with `wrap_ring()` first.
- **Fix:** run the check on `wrap_ring()`-padded arrays, as `benchmark/benchmark.py`'s
  `min_separation()` does. To be made in the faster driver (Phase 0b), with tests that it does
  not misfire on the benchmark ringlet and still catches a ringlet whose streamlines do cross.
  Not yet fixed.

## 2. `tests/unperturbed/inputs.py` predates the current driver

- **Found:** Python 3 port (Phase 0), 2026-10-17.
- **What:** it defines `confine_edges`, which the driver no longer reads, and lacks
  `confine_inner_edge` and `confine_outer_edge`, which the driver requires, so the test stops
  with a `NameError`.
- **Effect on hahn25:** none; it is a test input only.
- **Fix:** both set to `False` (the old setting) in this repo's `tests/unperturbed/inputs.py`;
  see `port_changes.md`. Fixed.

## 3. Closely spaced eccentric streamlines go unstable when sparsely sampled

- **Found:** checking that the Phase 0b benchmark ringlet stayed well behaved, 2026-10-18.
- **What:** an eccentric ringlet with streamline spacing 10⁻⁵ (e = 0.005, a de/da = 0.336, the
  wide 10⁻⁹ libration ringlet's surface density, 10–100 streamlines) stays smooth for a few tens
  of orbits; then a kink localized in longitude appears in every streamline's r(θ), grows by
  ~20× in ~25 orbits, the streamlines cross and the code makes NaNs. Run with the reference
  helper functions, without the driver's (misfiring, issue 1) stop:

  | variant (Nr = 30, Nt = 20 unless noted) | unstable after |
  |---|---|
  | as above | ~100 orbits |
  | Nt = 40 | ~70 orbits |
  | Nt = 80 | ~50 orbits |
  | Nt = 241 (the libration tests' sampling) | stable for 480 orbits |
  | Nr = 10 | ~140 orbits |
  | Nr = 100 | ~90 orbits |
  | ring gravity off (viscosity only) | ~110 orbits |
  | viscosity off (gravity only) | ~200 orbits |
  | dt = 0.25 instead of 0.5 | ~100 orbits |
  | a de/da = 0 | ~130 orbits |
  | circular streamlines | stable for 380 orbits |

- **So far:** it needs eccentric streamlines and some ring force, but not self-gravity, and it is
  not the timestep. It is not monotonic in Nt. The cause (numerical, e.g. the 3-point Lagrange
  interpolation between sparsely sampled eccentric streamlines, or physical) is unknown.
- **Effect on hahn25:** none known; the libration tests (2 streamlines, Nt = 241, spacing
  5×10⁻⁴–10⁻³) are stable.
- **Effect on epi_int_3d:** bending-wave runs have closely spaced streamlines with Nt ≈ 20, and
  satellite forcing makes them eccentric, so this needs to be understood before Phase 5. The
  Phase 0b benchmark uses circular streamlines meanwhile. Open.

## 4. The outer-edge flux in `A_P()` is evaluated at the wrong longitudes

- **Found:** checking the ring-pressure code against hahn13, 2026-10-18.
- **What:** hahn13 (Eq. 23 and text) gives the outermost streamline the acceleration
  p(r − Δ)/λ, the pressure of its inner neighbor at the particle's longitude. `A_P()`
  (`source/helper_fns.py:189`) uses `interpolate_fn(t, P, 1)[-2]`: the outer streamline's own P,
  interpolated to the longitudes of the inner neighbor's particles, then applied column by column to
  the outer streamline's particles. Using the outer streamline's own (one-sided) P would match the
  inner edge, `A[0] = -P[0]/λ`, but the longitudes are off by up to ~1.6 particle spacings in an
  eccentric ringlet, which changed the outer-edge value by up to 70% in one test state.
- **Scope:** `A_P()` also gives the shear- and bulk-viscous accelerations, so this affects
  viscosity, used in hahn25, as well as pressure (unused there).
- **Test:** the tapered π-twist nominal ringlet (`explore/taper/nr0100_pi*`, 100 streamlines,
  10⁴ orbits, viscosity on, pressure off) with three outer-edge terms:

  | outer-edge term | max q | width | last gap, outer/inner | max raggedness |
  |---|---|---|---|---|
  | as coded | 0.387 | ×2.655 | 23.67/23.51 | 1.0×10⁻⁴ |
  | own P, own longitudes (`--fix-edge-own`) | 0.387 | ×2.655 | 23.67/23.51 | 3.4×10⁻⁵ |
  | hahn13's p(r − Δ) (`--fix-edge`) | 0.297 | ×2.792 | 34.89/23.51 | 3.3×10⁻⁵ |

- **So:** the code's choice of the outer streamline's own one-sided flux is right: it mirrors the
  inner edge and keeps the edges symmetric, while hahn13's literal p(r − Δ) does not. The only
  error is the longitudes, which adds a little raggedness (×3) and nothing else visible.
- **Effect on hahn25:** probably negligible (smooth streamlines, Nt = 241); not checked.
- **Fix:** `A[-1] = P[-1]/lambda_[-1]` in the faster code, with a test. Minor. Open.

## 5. Ring-pressure scheme: omitted term and two properties of the discretization

- **Found:** same check; sympy checks P1–P7 in `notes/sympy/ring_pressure.py`.
- **Omitted term:** hahn13 also has the pressure gradient along a streamline,
  A_p = −(∂p/∂θ)/(rσ). The code omits it; it is smaller than the radial term by ~qΔa/a.
- **Momentum:** with the edge terms, the net radial pressure force on the ring is
  (−p₀ + p₁ + p_{N−2} − p_{N−1})/2 ≠ 0 (P3). The same holds for the viscous fluxes.
- **Odd–even decoupling:** central differences do not see a grid-scale zigzag of the
  streamlines, so pressure exerts no restoring force on it (P6), and likewise for the viscous
  fluxes. This may be why viscously spreading ringlet edges develop alternating gaps (seen in
  `explore/nominal`, 1000 streamlines).
- **Status:** properties of the published scheme, recorded for epi_int_3d's design. Open.

Issues 4 and 5 are minor for the nominal simulation. Decision (2026-10-18): if the z ≠ 0 runs use
pressure, pressure gets a deep dive then (physics and code) before it is relied on.
