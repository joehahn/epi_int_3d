# code: google/gemini-4-argon

Input: `fast/bookkeeping.py` · 2026-11-08 · tokens in 14572, out 895 · cost $0.0381

---

# Review of `fast/bookkeeping.py`

**Summary on equations, units, numerics and bugs:**
1. **Equations:** Every equation in `fast/bookkeeping.py` matches `notes/equations.md` (`V1`, `V8`, `H1`–`H3`) and the standard epicyclic/Lagrange formulas exactly—every term, factor, sign, and index.
2. **Units and normalization:** Consistent with $G M_p = 1$ in `_omega_kappa` and `vertical_drift`, while `two_GM_Rp` is passed in explicitly for `H2`.
3. **Numerics and bugs:** One low-severity angle-wrapping edge case in `_adjust` when `dt` is not small. No race conditions in `prange` loops (including the parallel reduction in `crossed`), no aliasing bugs, and the separation of `sin` and `cos` loops properly prevents LLVM's `sincos` fusion.

---

### Findings

1. **Type:** numerics
   - **Severity:** low
   - **Where:** `fast/bookkeeping.py`, lines 31–37, `def _adjust(angle):`
   - **Finding:** `_adjust` subtracts or adds `TWOPI` only once, which wraps `angle` into $[-\pi, \pi]$ only if the input is already in $[-3\pi, 3\pi]$. In `drift` (line 125), the argument `(Omg/Kap)*(M[k] + 2.0*e_sin_M) + wt[k]` can exceed $3\pi$ if a caller passes a multi-orbit timestep `dt` (e.g., `t[k] = 3.0`, `wt[k] ≈ 3.0`, `M[k] = Mk + Kap*dt ≈ pi + 2*pi` for `dt = 2*pi/Kap` gives an argument $> 3\pi$, leaving `t[k] > pi` after `_adjust` and breaking `sort_streamlines`).
   - **Suggestion:** None if `dt` is always a small fraction of an orbital period and `helper_fns.adjust_angle` does a single shift that must be matched bit-for-bit; otherwise note the `dt` bound in a comment.

2. **Type:** style
   - **Severity:** low
   - **Where:** `fast/bookkeeping.py`, lines 1–24 and all functions (`_adjust`, `_omega_kappa`, `rt2xy`, `xy2rt`, `polar`, `shift`, `shift4`, `drift`, `shift1`, `vertical_drift`, `sort_streamlines`, `linear_density`, `kick`, `wrap`, `crossed`)
   - **Finding:** As noted in `style/code.md` ("`fast/` predates this guide and is written in Claude's own style; restyling it is on the to-do list in `PLAN.md`"), the file uses module and function docstrings instead of `#` header and function comments (breaching **C1**, **C6**), uses `lam` instead of `lambda_` in `linear_density` (**C10**), uses `k`, `i`, `j`, `s` as loop/streamline indices instead of `pidx`, `sidx` (**C13**), uses `TWOPI` instead of `twopi` (**C14**), uses integer `0` instead of `0.0` in `if rw[i, j] - (...) < 0:` at line 209 (**C15**), omits parentheses around `if` and `while` conditions in `_adjust`, `sort_streamlines`, and `crossed` (**C18**), and takes output/work arrays as mutated arguments with parameters after arrays rather than returning tuples with physical parameters first (**C19**, though required here for preallocated Numba kernels).
   - **Suggestion:** Address when the scheduled restyling of `fast/` in `PLAN.md` is carried out.
