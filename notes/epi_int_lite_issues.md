# Issues found in epi_int_lite

Problems found in epi_int_lite (master, commit 5e0942a) while porting and extending it. `source/`
stays a straight port of epi_int_lite, so these are fixed in the faster code and epi_int_3d,
not there, unless noted. Changes made by the Python 3 port itself are in `port_changes.md`.

## Summary of the fixes

The fixes are applied at run time (`tools/patches.py`, called from a test's `inputs.py`) or in the
faster driver (`fast/epi_int_3d.py`). All but the crossing check are opt-in, so a default run gives
epi_int_lite's results bit for bit.

| issue | what is wrong | turn the fix on with | changes results? |
|---|---|---|---|
| 1 | the streamline-crossing check misfires near θ = ±π and stops runs | `patches.apply()` | no (only when a run stops) |
| 4 | the outer edge's viscous/pressure flux is taken at the wrong longitudes | `patches.apply(fix_edge='own')` | slightly |
| 6 | bulk viscosity is skipped when shear viscosity is 0, and shear's radial stress when bulk is 0 | `patches.apply(fix_viscosity=True)` | only with one viscosity off |
| 7 | each streamline's gravity uses the perturbed particle's linear density λ, not its own | input `fix_lambda = True` (`fast/epi_int_3d.py` only) | slightly |

Still open: closely spaced eccentric streamlines go unstable when sampled by too few particles
(issue 3), and properties of the pressure and viscosity scheme (issue 5). Issue 2 was an outdated
test input, fixed in this repo's copy.

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
  **Done** in `tools/patches.py` (`patches.apply()`, used by `tests/nominal/inputs.py` and
  `explore/run.py`); `fast/epi_int_3d.py` runs that wrapped check as a kernel. Runs whose inputs.py
  does not call `patches.apply()` still use the original check. Fixed.

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
- **Update (2026-10-19), the nominal ringlet at 333 streamlines** (`explore/twist0_n333_*`): every
  crossing checked (7 runs), except shear-only (issue 6), is a sawtooth along the streamlines at the particle
  scale, azimuthal wavenumber m ≈ Nt/2 (Nt = 20: m = 9; 40: m = 17–20; 320: m = 157), growing near
  (a − a_mid)/Δa₀ ≈ +0.2 from round-off. It needs neither eccentricity (circular streamlines crossed
  at ~216 orbits, m = 20) nor self-gravity (~168 orbits) nor viscosity (~849 orbits), and more bulk
  viscosity barely delays it. The code's axisymmetric modes (`notes/numerics/viscous_stability.py`,
  exact interpolation) are stable for all these cases, so the instability lives in the discretization
  along the streamlines, i.e. the 3-point Lagrange interpolation between neighbouring streamlines'
  particles. That is also why it depends on Nt. The shear-only crossings are different: a viscous
  overstability of the code's axisymmetric modes (growth 0.2–0.5 per orbit, wavelength ~4 spacings).

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
- **Fix:** `A[-1] = P[-1]/lambda_[-1]`: `tools/patches.py`, `apply(fix_edge='own')` (opt-in; the
  table's second row). Not yet in a kernel or a regression test. Minor. Fixed (opt-in).

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

## 6. Bulk viscosity is skipped when shear viscosity is zero, and vice versa for the radial stress

- **Found:** switching viscosity terms off separately in the 333-streamline nominal ringlet
  (`explore/twist0_n333_off`), 2026-10-19.
- **What:** `accelerations()` in `source/helper_fns.py` computes the viscous accelerations only
  `if ((c > 0.0) or (shear_viscosity > 0.0))`, so with pressure off and `shear_viscosity = 0` the
  bulk viscosity is silently skipped: the bulk-only run is bit-for-bit the run with no viscosity.
  Conversely the radial viscous stress, which carries shear viscosity's 4ν_s/3 and −2ν_s/3 parts as
  well as the bulk terms, is computed only `if (bulk_viscosity > 0.0)`, so with `bulk_viscosity = 0`
  shear viscosity exerts only its azimuthal stress.
- **Effect on hahn25:** none: its runs have both viscosities on.
- **Effect on epi_int_3d:** setting one viscosity to zero does not do what it says.
- **Fix:** `tools/patches.py`, `apply(fix_viscosity=True)` (opt-in): σ is computed when any of c, ν_s,
  ν_b is positive, and the radial stress when ν_s or ν_b is. Runs with both viscosities on are
  unchanged bit for bit (checked, nominal ringlet, 300 steps). To build into the faster code, with a
  test.
- **What it showed** (333-streamline nominal ringlet, Nt = 40, gravity on, to 10³ orbits; this is
  probably issue 3's instability): both viscosities, crossed at ~145 orbits; shear only, full stress,
  ~57 orbits; shear's azimuthal stress only (the unfixed shear-only run), ~23 orbits; bulk only, no
  crossing; no viscosity, ~849 orbits; viscosity but no gravity, ~168 orbits. So the shear stress
  drives the crossing, and the radial stresses (shear's radial part, and bulk) slow it.

## 7. Each perturbing streamline's gravity uses the perturbed particle's own λ

- **Found:** deriving the vertical streamline gravity (`notes/phase2/streamline_gravity.md`),
  2026-10-20.
- **What:** `ring_gravity()` computes A = 2Gλ/Δ_r with λ = `lambda_` of the perturbed particle,
  m₁Ω/(2πv) evaluated at that particle, for every perturbing streamline. The pull of streamline a′
  should use that streamline's own λ′ at the particle's longitude.
- **Effect:** λ′ − λ ≈ (a′ − a)∂λ/∂a, so each term is off by about 2G∂λ/∂a, whichever side a′ lies
  on. For circular, equal-mass streamlines λ = m₁/(2πa), so the fix adds −Σ_{a′} Gm₁/(πaa′) ≈
  −Gm_r/(πa²): a nearly uniform inward pull, about m_r/M_p of the planet's. That is the same order as
  the streamlines' curvature, which the wire approximation already neglects.
  - Nominal ringlet at 100 streamlines, 5000 orbits (`benchmark/runs/lamfix`, git-ignored): the
    mean precession rate changes by 1.4×10⁻⁷ of itself, the periapse twist (oscillating between
    −0.39 and −0.01 rad) by at most 9×10⁻⁶ rad, e by 6×10⁻⁵ of itself and a by 10⁻¹¹.
- **Fix:** in `fast/epi_int_3d.py`, input `fix_lambda = True` (opt-in, default False so the
  reference results stay bit for bit). The kernel `kernels.ring_gravity_lam_parallel` interpolates λ
  of the perturbing streamline to the particle's longitude with the same Lagrange weights as r′.
  Tested by `tests/test_kernels.py` (closed-form sums on a circular ringlet, including the
  correction above). With the fix off, all benchmark final states still match the reference. Fixed.
