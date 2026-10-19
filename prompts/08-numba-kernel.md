# 08 — Ring gravity as a Numba kernel

(Assembled from a run of prompts on 2026-10-19; it is the final version of what was asked.)

Start Phase 0b step 2, the Numba gravity kernel (after recording the narrow-ringlet results and
correcting the "one aligned ellipse" wording). Explain how a Numba kernel speeds up execution: it
compiles the loop to machine code, removes the temporary arrays and the Python overhead per pass, and
can use all cores with `prange`.

- Rewrite `ring_gravity` as a compiled loop, without `fastmath`, outside `source/`, which stays the
  unchanged reference implementation.
- Check that its results match the reference on the benchmark ringlet and on `tests/nominal`.
- Measure the speedup at Nr = 10, 100 and 1000 for the README.
- Write it up and update the plan.

Then make the nominal ringlet the one with zero initial periapse twist, evolved for 5×10⁴ orbits
(~2 viscous timescales) if it can safely go that far. If it does, rerun short simulations of it with
the slow (reference) and fast codes at 33, 100 and 333 streamlines (factors of ~3.3 rather than 10),
make sure their outputs are equivalent, and update the README. Could the 333-streamline run go
5×10⁴ orbits?
