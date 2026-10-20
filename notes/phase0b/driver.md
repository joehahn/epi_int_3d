# Phase 0b, step 4: epi_int_3d's own driver

`fast/epi_int_3d.py` is epi_int_3d's driver, planar for now: Phase 3 adds the vertical motion to it.
Unlike `fast/epi_int_fast.py`, which runs the reference driver with only ring gravity swapped out, it
has its own main loop and bookkeeping. It reads the same `inputs.py` (including what
`tools/patches.py` changes) and gives the reference's results bit for bit. Numba 0.68.0, Python
3.14.4, NumPy 2.5.4, Apple M5 (4 performance + 6 efficiency cores).

## What it changes

- **The state lives in preallocated arrays, updated in place**, instead of new arrays for every
  operation.
- **The bookkeeping is Numba kernels** (`fast/bookkeeping.py`):
  - the polar ↔ Cartesian conversions of `coordinate_kick`, `mixed2planeto` and `planeto2mixed`;
  - the drift (`coords2elem`, `drift`, `elem2coords` fused into one kernel);
  - each streamline's linear density;
  - the longitude wrap and the wrapped crossing check;
  - sorting by insertion sort, which is nearly free because particles barely move in longitude
    per step.
- **The whole-ring sums stay NumPy's.** The centre-of-mass sums and means between the kernels are
  NumPy's own pairwise sums, as in the reference.
- **Ring gravity** is the existing kernel in `fast/kernels.py`.
- **Still the reference's NumPy code:** pressure and viscosity (`helper_fns.accelerations`), the
  initialization, the output kicks and storage. Whatever `inputs.py` patches stays patched:
  - the edge taper;
  - `fix_viscosity` and `fix_edge`;
  - the crossing check, which runs as a kernel when it is the patched (wrapped) one and as the
    reference's own code otherwise.

## Bit for bit

Every kernel repeats the reference's arithmetic in the reference's order (e.g. `a*a*Omg/r` is
`((a*a)*Omg)/r`), without `fastmath`. Two traps:

- **sincos.** When sin and cos of the same angle are taken in one loop, LLVM fuses them into one
  `sincos` call. Its sin differs from NumPy's (and libm's) sin in the last bit for 0.35% of angles.
  NumPy's cos, arctan2 and sqrt match Numba's, and so does sin when it is computed alone, so the
  kernels take sin and cos of an angle in separate loops.
- **Sums.** NumPy sums pairwise, Numba sequentially, so every sum is left to NumPy.

| check | result |
|---|---|
| regression suite, `EPI_INT_DRIVER=fast/epi_int_3d.py uv run pytest` | 3 passed |
| benchmark final states, Nr = 10, 100, 1000 and the former nominal ringlet at 33, 100, 333 (`--check`) | identical |
| all outputs versus `fast/epi_int_fast.py`, short runs (`benchmark/runs/compare.sh`): nominal at 100 and 1000 streamlines; viscous benchmark ringlet; `tests/spreading`; `tests/unperturbed` (no ring forces); pressure Q = 1 with both viscosities and both edge torques; `fix_viscosity` with shear viscosity only | identical |
| the former nominal ringlet at 333 streamlines, 3000 steps | identical; both cross at step 1825 |
| the full nominal run, 629 000 steps (5×10⁴ orbits), versus `tests/nominal/output` | identical: r, t, v_r, v_θ and times at all 501 outputs |

## Speed

`benchmark.py --driver fast/epi_int_3d.py --label driver --check`, 10 threads for gravity, time per
step (`benchmark/results/driver_*.json`; the previous step's numbers are from `numba_*.json`):

| ringlet | Nr × Nt | reference | gravity kernel only (`epi_int_fast.py`) | this driver | gain over previous |
|---|---|---:|---:|---:|---:|
| benchmark (viscous) | 10 × 20 | 0.49 ms | 0.38 ms | 0.30 ms | ×1.3 |
| | 100 × 20 | 6.35 ms | 0.88 ms | 0.73 ms | ×1.2 |
| | 300 × 20 | 51.5 ms | 2.62 ms | 1.95 ms | ×1.3 |
| | 1000 × 20 | 447 ms | 14.1 ms | 12.3 ms | ×1.15 |
| former nominal (viscous) | 33 × 40 | 1.83 ms | 0.66 ms | 0.49 ms | ×1.3 |
| | 100 × 40 | 12.4 ms | 1.7 ms | 1.10 ms | ×1.5 |
| | 333 × 40 | 105 ms | 6.3 ms | 4.3 ms | ×1.5 |

At Nr = 10 a step takes too little time for `benchmark.py`'s 1000-step difference (n2 − n1) to
resolve it: repeated runs gave 0.03–0.38 ms for the same driver. The Nr = 10 row is instead the
wall time of 21 000 steps minus that of 1000 steps (11 000 for the reference), three repeats each,
on an idle machine; they agree to ±5%. The earlier 0.14 ms for `epi_int_fast.py` in `numba.md` was
such noise.

The full nominal run took 177 min, sharing the machine with these benchmarks, so it is not a
timing; `epi_int_fast.py` took 169 min for the same run.

The gain is modest because at 1000 × 20 the gravity kernel is already 8.7 ms of the step. The
bookkeeping is now ~3.4 ms, almost all of it the ~6 arctan2 and ~11 sin/cos pairs that each
particle needs per step. Bit-for-bit results forbid skipping any of them.

## Threads for the bookkeeping: not yet

Each kernel can also be compiled with its loops spread over all cores (`bookkeeping.PARALLEL`); the
iterations are independent, so the results are the same. But Numba's only threading layer available here,
workqueue, costs ~60 µs to launch even an empty parallel loop. That is more than these kernels take
at 20 000 particles: in the benchmark, parallel bookkeeping made the step slower at 300 and 1000
streamlines (4.5 and 14.3 ms). Measured crossover: serial versus parallel per kernel

| particles | rt2xy | xy2rt | drift |
|---:|---|---|---|
| 6 000 | 35 vs 197 µs | 73 vs 268 µs | 103 vs 281 µs |
| 20 000 | 155 vs 223 µs | 384 vs 330 µs | 340 vs 344 µs |
| 80 000 | 900 vs 366 µs | 2261 vs 627 µs | 1474 vs 649 µs |

So the driver runs them serially below 50 000 particles (`PARALLEL_MIN_PARTICLES`), which covers
the planned bending-wave runs (8 000–16 000 particles). An OpenMP or TBB threading layer, with lower
launch costs, needs `intel-openmp` or `tbb`, which this project does not install.

## What is left

At 1000 × 20 the step is ~70% gravity and ~30% bookkeeping, so the remaining levers are:

1. **Gravity.** The algorithmic option of step 5 (near-field sum plus a smooth far field), which
   could not stay bit for bit.
2. **Viscosity and pressure as kernels.** These help viscous runs at small Nr, where the NumPy
   viscosity is now the largest cost. The first bending-wave runs have both off.
3. **Reuse sin/cos.** Two of the ~11 sin/cos pairs per step are of a longitude whose sin and cos
   were just computed; reusing them saves ~0.4 ms at 1000 × 20.
4. **A faster threading layer**, which would let the bookkeeping run on all cores.
