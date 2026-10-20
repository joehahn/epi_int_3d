# Phase 0b, steps 2–3: ring gravity as a Numba kernel

`fast/kernels.py` computes ring self-gravity, which takes 99.6% of the reference implementation's time
at Nr = 1000 (`benchmark.md`), as one compiled loop. `fast/epi_int_fast.py` is the faster driver: it
runs the reference driver, `source/epi_int_lite.py`, unchanged, with `helper_fns.ring_gravity` replaced
by the kernel. `source/` is not modified. Numba 0.68.0 (llvmlite 0.50.0), Python 3.14.4, NumPy 2.5.4.

## The kernel

For each target particle (i, j) of the longitude-wrapped arrays and each source streamline
s = i − shft (mod Nr), shft = 1, …, Nr − 1, it computes the 2nd-order Lagrange weights at t_ij from
particles j − 1, j, j + 1 of streamline s, applies them to r, v_r and v_θ, and adds 2Gλ_ij/Δr along
streamline s's unit normal: what `ring_gravity` does with `fast_gravity = False`, with no temporary
arrays and the weights computed once instead of three times.

**Same results as the reference, bit for bit.** Each particle's sum runs over shft in the reference's
order, each term uses the reference's operations in the reference's order, and there is no `fastmath`,
so no operation is reordered or fused. The loops are ordered i (prange), shft, j: the innermost loop
runs over contiguous, independent particles, and splitting the outer loop across threads leaves
every sum's order unchanged. Only the interior columns are computed; the wrapped edge columns, which
`accelerations()` drops, are left zero.

**Caching pitfall.** Compiling the serial and the parallel versions from one Python function with
`cache=True` made the parallel one load the serial one's cached machine code (no threading at all).
The kernels are therefore compiled at start-up, about 1 s.

## Checks

| check | result |
|---|---|
| kernel vs `ring_gravity`, nominal ringlet's first and last states | identical (`np.array_equal`) |
| regression suite, `EPI_INT_DRIVER=fast/epi_int_fast.py uv run pytest` | 3 passed |
| benchmark final states, Nr = 10, 100, 1000 (`--check`) | max relative difference 0.0 |
| `tests/nominal`, 126 000 steps, all 201 outputs | identical; 3.3 min versus 24.4 min |

## Speed (Apple M5: 4 performance + 6 efficiency cores)

Kernel alone, Nr = 1000, Nt = 40 (reference `ring_gravity` 0.75 s):

| | time | speedup |
|---|---:|---:|
| serial `@njit` | 124 ms | ×6 |
| parallel, 1 thread | 72 ms | ×10 |
| 2 threads | 36 ms | ×21 |
| 4 threads | 28 ms | ×27 |
| 10 threads | 15 ms | ×50 |

The parallel build runs faster than the serial one even on one thread, probably because the
outlined parallel loop lets LLVM vectorize the inner loop (not verified). Scaling is good to the
4 performance cores and weaker beyond.

Whole driver, `benchmark.py --driver fast/epi_int_fast.py --label numba --check`, Nt = 20, 10 threads
(`benchmark/results/numba_Apple-M5.json`):

| Nr | reference | fast | speedup |
|---:|---:|---:|---:|
| 10 | 0.47 ms | 0.14 ms* | ×3.4* |
| 30 | 1.28 ms | 0.47 ms | ×2.7 |
| 100 | 6.35 ms | 0.88 ms | ×7.2 |
| 300 | 51.5 ms | 2.62 ms | ×20 |
| 1000 | 447 ms | 14.1 ms | ×32 |

*Noise: measured over more steps, Nr = 10 takes 0.49 ms (reference) and 0.38 ms (fast); see
`driver.md`, which also has the next step, epi_int_3d's own driver.

The nominal ringlet (`--ringlet nominal`: eccentric, tapered edges, Nt = 40, width and mass fixed;
`benchmark/results/*_nominal_Apple-M5.json`), at streamline counts a factor ~3.3 apart:

| Nr | reference | fast | speedup | final state |
|---:|---:|---:|---:|---|
| 33 | 1.83 ms | 0.66 ms | ×2.8 | identical |
| 100 | 12.4 ms | 1.7 ms | ×7.4 | identical |
| 333 | 105 ms | 6.3 ms | ×16.5 | identical |

The 5×10⁴-orbit nominal run (629 000 steps) takes 16 min with the fast driver.

## Where the time goes now

cProfile, 400 steps at Nr = 1000, Nt = 40: gravity 6.4 s (16 ms/step), everything else ~13 ms/step:
Lagrange interpolation for the viscosity (`delta_f`, `lagrange_poly_fit`, `sidestep`) 1.3 s,
coordinate conversions (`xy2rt`, `rt2xy`, `coords2elem`) 1.7 s, `sort_particles` 0.4 s. At Nr = 100
the NumPy bookkeeping is nearly all of it. So the next gain, ×2 at Nr = 1000 and more at Nr = 100,
is step 4: viscosity (and pressure) as kernels, and cheaper sorting and coordinate conversions.

What this buys now: the narrow 1000-streamline ringlet (`nominal.md`) for 10⁴ orbits at Nt = 40 is
~126 000 × 29 ms ≈ 1 h, versus ~28 h with the reference.
