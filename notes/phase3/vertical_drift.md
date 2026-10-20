# Phase 3, step 1: the vertical drift

`fast/epi_int_3d.py` now follows each particle's vertical motion, z and v_z, through the planet's
potential. This is the drift step of the plan's Phase 3 (theory: Phase 2a, `notes/equations.md`
V1–V8). There are no vertical ring forces yet (Phase 2b), and ring gravity still acts as if every
particle were at z = 0.

## What the code does

- **Drift.** The vertical part of the drift is a rotation of (z, v_z/ν_a) through ν_aΔt
  (`bookkeeping.vertical_drift`, equation V8), with ν_a from V1 at the semimajor axis that the
  in-plane drift has just computed. V8 is V2–V4 exactly: the mean longitude advances at Ω_a and the
  node regresses at Ω_a − ν_a, so the phase λ − Ω advances at ν_a. It needs no atan2, and z = v_z = 0
  stays exactly 0. The in-plane drift kernel is unchanged.
- **Barycentric bookkeeping.** z and v_z take part in it like x and y: the coordinate kick moves z by
  dt·m·mean(v_z); `mixed2planeto` and `planeto2mixed` shift v_z by the ring's mean momentum.
- **Sorting.** The insertion sort carries z and v_z along with r, θ, v_r and v_θ.
- **Initial conditions.** Optional inputs `inclination` (i) and `node` (default 0) give every particle
  z = a i sin(λ − Ω), v_z = a ν_a i cos(λ − Ω) (V2).
- **Output.** `z.npy` and `vz.npy` (planetocentric and sorted, like the other outputs) join the
  reference's files. `epi_int_3d.restore_output` reads them, with zeros for older output, and
  `epi_int_3d.vertical_elements` gives i and the node (V3).

The vertical work costs about 5% of a step at 1000 × 20 (13.0 versus 12.3 ms).

## z ≡ 0 reproduces the 2D code

None of the in-plane arithmetic uses z or v_z, so the in-plane results are the reference's bit for
bit, with or without inclination.

| check | result |
|---|---|
| regression suite, `EPI_INT_DRIVER=fast/epi_int_3d.py uv run pytest` | 3 passed |
| benchmark final states, Nr = 10, 100, 1000, and the former nominal ringlet at 33, 100, 333 (`--check`) | identical to the reference |
| all outputs versus `fast/epi_int_fast.py` (`benchmark/runs/compare.sh`): nominal ringlet at 100 streamlines, viscous benchmark ringlet, `tests/spreading`, `tests/unperturbed`, former nominal ringlet with pressure Q = 1 and edge confinement | identical; z and v_z exactly 0 |
| the same with an inclination: `tests/spreading` (i = 10⁻³), `tests/unperturbed` (i = 0.01), nominal ringlet at 100 streamlines (i = 10⁻³) | in-plane outputs identical |

## Inclined test particles (`tests/inclined/`)

tests/unperturbed's 51 eccentric streamlines (10⁻⁸ < e < 5×10⁻³), without ring forces, all at
i = 0.01, for 1000 orbits (`inclined.ipynb`; also `tests/test_vertical.py`, 200 orbits):

- **The code keeps its own equations exactly.** Inclinations stay constant to 2×10⁻¹³, and the nodes
  follow node₀ + (Ω_a − ν_a)t to 3×10⁻¹¹ rad, through 3.7 cycles of regression.
- **Exact integration.** The same particles were integrated exactly in the 3D J2 potential (DOP853),
  with their nodes measured by the same elements (V3). For e < 10⁻⁴, the exact nodes regress more
  slowly than the code's, by the fraction

  | i | (exact rate)/(code rate) − 1 | scatter |
  |---:|---:|---:|
  | 0.01 | −4.04×10⁻⁴ | 2×10⁻⁶ |
  | 0.003 | −3.64×10⁻⁵ | 5×10⁻⁷ |
  | 0.001 | −4.04×10⁻⁶ | 6×10⁻⁷ |

  That is −4.04 i², the first-order theory's O(i²) error. It has a simple origin: the code takes a
  from h_z = r v_θ = h cos i, so a_code = a cos²i. The nodal rate goes as a^(−7/2) and, in the exact
  problem, also as cos i, so exact/code = cos⁸i ≈ 1 − 4i².
- **Larger e.** For e ≳ 10⁻³ the difference grows roughly as e², with coefficient −3 to −4. A few
  particles are outliers, up to −7×10⁻⁴ at i = 0.01. That is the measurement: the exact orbits' V3
  nodes carry short-period terms of amplitude ≈ 2e rad, the theory's O(e·i) error in z (V6). The
  10-orbit output spacing aliases them into slow wobbles, which bias the fitted rates.

So the errors are the expected O(i², e², e·i) of the first-order theory. For bending waves,
i ≲ 10⁻³, so the nodal rate is good to about 10⁻⁵ of itself.

## Next

- The vertical self-gravity of streamlines (Phase 2b): derive and sympy-check it, then add it as a
  kernel. Until then an inclined self-gravitating ringlet regresses differentially, each streamline
  at its own Ω_a − ν_a.
- The Phase 4 test "narrow inclined self-gravitating ringlet" (rigid nodal precession) needs 2b.
