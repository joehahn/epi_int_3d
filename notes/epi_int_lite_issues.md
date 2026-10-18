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
