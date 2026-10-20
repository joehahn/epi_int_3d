# 10 — epi_int_3d's own driver, with the bookkeeping as Numba kernels

(Assembled from a run of prompts on 2026-10-19/20; it is the final version of what was asked.)

Start Phase 0b step 4, written as the skeleton of the 3D driver rather than a polished 2D one: a new
driver in `fast/` that no longer wraps `source/`, with the step-4 savings (keep particles sorted between
steps, convert coordinates once per step, viscosity interpolation in Numba where it pays), and state
arrays laid out so that z and v_z can be added. Check it bit for bit, or to ~1e-12, against the
reference on the benchmark ringlets and the regression suite, and add a speedup row to the README.
Afterwards: the vertical drift (Phase 2a) in that driver, with z ≡ 0 reproducing the 2D results.
