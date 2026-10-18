# Phase 0b, step 1: benchmark of the reference implementation

Prompt: `prompts/06-benchmark.md`. Code: `benchmark/benchmark.py` (timing, cProfile, reference final
states) and `benchmark/line_profile.py` (line-by-line profile). Raw results:
`benchmark/results/reference_Apple-M5.json`.

## Setup

- Machine: Apple M5 laptop, 10 cores (4 performance + 6 efficiency), 16 GB, macOS 26; Python
  3.14.4, NumPy 2.5. The reference implementation is single-threaded.
- Ringlet (`benchmark/inputs.py`): the wide 10⁻⁹ libration ringlet, eccentric (e = 0.005,
  a de/da = 0.336, aligned periapses), with ring gravity (`G_ring = 1`, `fast_gravity = False`),
  shear and bulk viscosity 10⁻¹³, pressure off (`Q_ring = -1`), dt = 0.5. As Nr grows,
  `radial_width` = 10⁻⁵(Nr − 1) and `total_ring_mass` = 10⁻⁶ × `radial_width`, so the streamline
  spacing and surface density stay fixed.
- Each case runs `source/epi_int_lite.py` unchanged, as its own process, for n1 and then n2
  timesteps; time per step = (wall time difference)/(n2 − n1), so start-up, initialization and
  output cancel. 1 timestep = 1 velocity kick + 1 drift; 12.6 timesteps per orbit.
- No streamlines crossed: the smallest gap between adjacent streamlines stayed at 0.49–0.68 of
  the initial spacing (the column below).

## Time per step

| Nr | Nt | s/step | μs per (Nr² Nt) | min gap |
|---:|---:|---:|---:|---:|
| 10 | 20 | 0.00044 | 0.22 | 0.68 |
| 30 | 20 | 0.0013 | 0.072 | 0.51 |
| 100 | 20 | 0.0063 | 0.032 | 0.61 |
| 300 | 20 | 0.047 | 0.026 | 0.65 |
| 1000 | 20 | 0.43 | 0.021 | 0.65 |
| 100 | 10 | 0.0045 | 0.045 | 0.49 |
| 100 | 40 | 0.011 | 0.027 | 0.62 |

All seven cases fit, to ±20%,

    s/step ≈ (Nr − 1) × (33 μs + 0.017 μs × Nr (Nt + 2))

that is, a sum over the Nr − 1 other streamlines, each costing a fixed 33 μs of Python and
NumPy call overhead plus 0.017 μs per array element (each streamline has Nt + 2 particles once
wrapped in longitude). The cost becomes Nr² × Nt only when Nr Nt ≳ 2000; below that the fixed
overhead dominates, which is why time per step grows more slowly than Nt at Nr = 100.

## Where the time goes

**Ring gravity is 99.6% of the force calculation at Nr = 1000** (line profile, 20 steps); shear
and bulk viscosity are 0.2% each. `PLAN.md`'s guess is confirmed. Inside `ring_gravity`
(`source/helper_fns.py:166`), for each of the Nr − 1 shifts:

| line | share of `ring_gravity` |
|---|---:|
| `interpolate_fn(t, r, -shft)`, the neighbor's r at each particle's longitude | 30.5% |
| `interpolate_fn(t, vr, -shft)` | 29.6% |
| `interpolate_fn(t, vt, -shft)` | 29.6% |
| `sqrt`, divisions and sums | 10% |

Each `interpolate_fn` call shifts the arrays with `advance` and `sidestep` (6 `np.concatenate`
copies, 30% of its time) and then evaluates a 3-point Lagrange polynomial (70%), whose
weights are 2/3 of the work. The three calls per shift interpolate at the same longitudes, so
they recompute identical shifted longitudes and identical weights; only the interpolated
quantity differs. The whole sum makes about 90 full-size temporary arrays per shift.

At small Nr the same functions lead (cProfile, Nr = 10: `lagrange_poly_fit`, `sidestep`,
`advance`), but the drift, coordinate conversions and sorting take a larger share.

## A bending-wave run with the current code

From the fit, with Nt = 20 and dt = 0.5:

| Nr | s/step | 10⁵ steps (8,000 orbits) | 10⁶ steps (80,000 orbits) |
|---:|---:|---:|---:|
| 400 (adaptive spacing) | 0.073 | 2 h | 20 h |
| 800 (uniform spacing) | 0.26 | 7 h | 3 days |

A parameter study of tens of such runs is weeks of laptop time, so the speedup matters.

## The driver's crossing check misfires

`monitor_streamlines()` (`source/helper_fns.py:423`) reported crossing streamlines in every
case, after 7–143 steps, although none crossed. It interpolates each streamline's inner
neighbor with `interpolate_fn` on unwrapped arrays: for the first and last particle,
`sidestep` wraps the array around without adding ±2π to the longitude, so the Lagrange fit
spans θ ≈ (2.8, 3.1, −3.1) and extrapolates wildly (in one case a gap of −5.5×10⁻⁴, ~50
spacings, where the true gap was 0.66 spacings). The forces are unaffected: `accelerations()`
wraps the ring with `wrap_ring()` first.

The error grows with the eccentricity and the longitude gap between particles and is compared
with the streamline spacing, so it bites for closely spaced, eccentric, sparsely sampled
streamlines: the benchmark and the bending-wave runs, not the hahn25 tests. None of the 21
self-interacting test runs was flagged (the unperturbed run is flagged at step 1, but it isn't
self-interacting, so it doesn't stop). Because a flagged self-interacting run is terminated
(`update_display`), a long bending-wave run would stop early. `source/` is left as is; the
benchmark uses its own check on the wrapped arrays (`min_separation()`), and the faster code
should do the same.

## Reference states

The final states after n2 steps for Nr = 10, 100 and 1000 (Nt = 20) are in
`benchmark/reference/` (700 kB). `benchmark.py --driver <faster driver> --label <name> --check`
times other code with the same harness and compares its final states with these to relative
precision 10⁻¹². Rerunning the reference implementation reproduces them exactly.

## Recommendation: what to speed up first

1. **Ring gravity as a compiled Numba kernel** (Phase 0b, step 2). One loop over target
   particle × source streamline that finds the neighbor's three nearest particles in longitude
   (wrapping at ±π), computes the Lagrange weights once, and applies them to r, vr and vt: no
   temporaries, no Python overhead per shift, and a third of the arithmetic. This is 99.6% of the
   time at the target size, so it alone sets the speedup. To match the reference to 10⁻¹², keep
   the order of the sum over shifts and the same arithmetic per element, with no `fastmath`.
2. **Multicore** (step 3): `prange` over target streamlines in that kernel.
3. Everything else (viscosity, drift, conversions, sorting) is under 1% at Nr = 1000 and can
   follow later. At Nr ≈ 400–800 it becomes visible only after gravity is ~100× faster.
4. Fix the crossing check in the faster driver, as above.
