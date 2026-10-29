# Phase 3, step 2: vertical streamline gravity

`fast/epi_int_3d.py` now has ring gravity in 3D: each perturbing streamline is a straight wire along
its shape, and it pulls the particle toward the wire's near point in r, θ and z. This is the paper's
Eqn. (23) (`notes/equations.md` G16, G17), the kernel `ring_gravity_3d` in `fast/kernels.py`.

## What the code does

- **Kernel.** At the particle's longitude, the same three-particle Lagrange polynomials that give the
  perturbing streamline's r′ also give z′ and the shape's slopes s_r = (dr′/dθ)/r′ and
  s_z = (dz′/dθ)/r′. With the crossing-point gaps d_r = r′ − r, d_z = z′ − z, N² = 1 + s_r² + s_z² and
  c = (d_r s_r + d_z s_z)/N², the pull is 2Gλ′(d_r − c s_r, −c, d_z − c s_z)/Δ² with
  Δ² = d_r² + d_z² − c(d_r s_r + d_z s_z). λ′ is the perturbing streamline's (issue 7) unless
  `fix_lambda = False`.
- **In the plane.** A term with d_z = s_z = 0 exactly uses the planar kernel's arithmetic,
  (2Gλ′/d_r)(1, −s_r, 0), which is Eqn. (23) at z = 0. So a ring at z = v_z = 0 gets the planar
  `nearpoint` and `shape_slope` kernel's result bit for bit, and stays at z = 0.
- **Driver.** `Ring.velocity_kick` wraps z in longitude like r and t, the `ring_gravity` drop-in
  that `helper_fns.accelerations` calls passes it to the 3D kernel, returns A_r and A_θ to
  `accelerations` and leaves A_z, which kicks v_z (planetocentric, like the in-plane kick). The
  reference's `source/` is untouched. During initialization, before z exists, the planar kernel runs.
- **Input.** `vertical_gravity` (default: `nearpoint and shape_slope`, i.e. on unless
  `epi_int_lite_gravity`); it is the exact wire along the shape, so it requires both and refuses
  otherwise. Off, ring gravity acts as if every particle were at z = 0, as before.

## Tests

- `tests/test_kernels.py`:
  - `test_3d_in_plane_is_planar`: at z = 0 the 3D kernel equals the planar `nearpoint` + `shape_slope`
    kernel bit for bit (both λ options, serial and parallel), and A_z = 0.
  - `test_3d_wire`: m = 2 streamlines also warped in z (slopes up to ~0.1), against the near point
    computed as a vector, Δ = d₀ − (d₀·t̂)t̂, with the analytic shape and slopes: agreement to 4×10⁻⁹
    (r), 3×10⁻⁷ (θ), 3×10⁻⁸ (z) of the largest pull, the three-point Lagrange fit's error.
  - `test_3d_tilted_ringlet`: a flat circular ringlet tilted by 0.1 rad is pulled within its own
    plane, |A·n|/|A| = 1.6×10⁻⁸, while A_z is 10% of |A|.
- The regression suite passes (`uv run pytest`, 9 tests).
- A z = 0 run with the default gravity (now through the 3D kernel), the libration ringlet
  `narrow_ring/ring_mass_1.0e-9`, is bit for bit the earlier run with the planar kernel
  (`explore/libration_gravity/.../new`): r, θ, v_r, v_θ, z, v_z identical.

## First physics: an inclined narrow ringlet

The same libration ringlet (two streamlines, Δa = 5×10⁻⁴, m_r = 10⁻⁹, e = 0.0025), tilted rigidly
by i = 0.002 (node 0), run for t = 5×10⁴ (about 8000 orbits); scratch runs, not kept:

| | node_outer − node_inner at t = 5000, 12500, 25000, 37500, 50000 | Δi range |
|---|---|---|
| `vertical_gravity = False` | +0.039, +0.098, +0.195, +0.293, +0.390 rad | ~10⁻¹¹ |
| `vertical_gravity = True` | −0.0002, +0.012, −0.0008, −0.012, +0.0015 rad | −8×10⁻⁸ to +5×10⁻⁵ |

Without vertical gravity the nodes spread at the differential nodal regression rate of J₂; with it,
self-gravity locks them: the node difference oscillates within ±0.012 rad, as the periapses of an
eccentric ringlet librate, and a small inclination difference Δi develops, the vertical analog of the
equilibrium Δe. The mean inclination stays 0.00200. This is the Phase 4 test "narrow inclined
self-gravitating ringlet: rigid nodal precession / libration", qualitatively; the quantitative version
(the libration frequency against a vertical analog of Borderies, Goldreich & Tremaine's) is next.
