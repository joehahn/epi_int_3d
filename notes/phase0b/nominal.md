# epi_int_3d's nominal ringlet

A ringlet for benchmarking and testing this project's code (`tests/nominal/`), found by the parameter
search in `explore/` (2026-10-18). The goal was a ringlet that is eccentric, has a periapse twist that
makes its self-gravity non-axisymmetric, is massive enough that self-gravity rather than J2 controls
the twist, spreads viscously by about ×2, and stays smooth for 10⁴ orbits with as many streamlines as
possible, ideally 1000.

All runs use the reference implementation (`source/`) unchanged, through `explore/run.py`, which
applies `tools/patches.py` at run time: always the longitude-wrapped crossing check (issue 1 in
`epi_int_lite_issues.md`), and optionally an edge taper. Diagnostics are from `explore/diagnose.py`:
q = a|Δ(e e^{iϖ})|/Δa, the largest over adjacent streamline pairs (they touch when q → 1); the
smallest gap between adjacent streamlines over the initial spacing; the raggedness, the largest rms
residual of a streamline's r(θ) about a Fourier fit through m = 3, in spacings.

## The nominal ringlet

| parameter | value |
|---|---|
| streamlines × particles | 100 × 40 |
| width Δa₀ | 0.0125, edges tapered (below) |
| mass | 3×10⁻⁸ |
| e, a de/da | 10⁻⁴, 0 |
| periapse twist | π from inner to outer edge (`w_prime` = π e ā/Δa₀ = 0.02529) |
| ν_s = ν_b | 7×10⁻¹¹ |
| pressure, fast_gravity | off, off |
| run | 10⁴ orbits (126 000 steps), about 30 min on one M5 core |

Result (`tests/nominal/nominal.ipynb`): no crossing or NaN; max q 0.39; smallest gap 0.97; raggedness
≤ 1.2×10⁻⁴ spacings; width ×2.66; |ΔL/L₀| ≤ 1.3×10⁻⁵. The ringlet's body holds its twist: its
periapses librate about the middle streamline's rather than winding up. The outermost few streamlines
on each edge spread away from the body and their e and ϖ wander.

**Edge taper.** Streamlines keep equal masses but sit at equal-mass quantiles of a surface density that
is flat across the middle and falls as cos² to zero over 20% of the width at each edge. The width is
1.25× that of an untapered ringlet with the same middle spacing and surface density.

## How it was found

**Width 10⁻³, e = 10⁻², 241 particles per streamline** (`explore/ladder/`). At 10 streamlines,
mass 10⁻¹⁰ does not hold the twist (q 0.2 → 0.56), 10⁻⁸ does (it librates, q 0.13–0.28), 10⁻⁷ crosses at
once. More viscosity makes J2 win. At 100 streamlines 10⁻⁸ crosses at step 14 and 3×10⁻⁹ at 391
orbits; at 1000 streamlines 3×10⁻⁹ and 10⁻⁸ cross within 4 steps. A π twist is impossible at e = 10⁻²
(q ≈ 31 at this width), so e was lowered.

**Width 0.01, π twist, 40 particles per streamline** (`explore/nominal/`). At 10 streamlines and
e = 10⁻³, mass 10⁻⁷ held best (q 0.26–0.63) and 3×10⁻⁷ crossed at 6100 orbits; at e = 10⁻⁴,
3×10⁻⁸ held the twist (q ≤ 0.10). At 100 streamlines that ringlet (ν = 7×10⁻¹¹) was smooth for
10⁴ orbits (q ≤ 0.40, width ×2.92). At 1000 streamlines it crossed at step 302:

| 1000 streamlines, 400 steps | result |
|---|---|
| e = 10⁻⁴, π twist, 3×10⁻⁸ | crossed at step 302 |
| circular, 3×10⁻⁸ | crossed at step 308 |
| 241 particles per streamline | crossed at step 308 |
| mass 10⁻⁸ | smooth |
| ring gravity off | smooth |

The cause was the edges: viscous spreading moved the edge streamlines out by ~10 spacings with
alternating large and small gaps (odd–even decoupling, issue 5), and at 3×10⁻⁸ a compression front
then ran inward from the outer edge at ~2 streamlines per orbit, squeezing gaps to 0.42 of the spacing
before they crossed. The zigzag self-gravity scales as 1/gap², so a squeezed region collapses.
100 streamlines are spaced 10× wider, and their edges spread by only ~1 spacing.

**Tapered edges, width 0.0125** (`explore/taper/`):

| run | 100 streamlines, 10⁴ orbits | 1000 streamlines, ~32 orbits |
|---|---|---|
| circular | – | smooth |
| e = 10⁻⁴, no twist | smooth, q ≤ 0.35 | smooth |
| e = 10⁻⁴, π twist | **smooth, q ≤ 0.39** | crossed at step 361 |
| e = 10⁻³, no twist | crossed at ~8600 orbits | crossed at step 89 |
| e = 10⁻³, 10° twist | crossed at ~8400 orbits (q climbing) | crossed at step 55 |

The taper removes the edge collapse. At e = 10⁻³ J2 winds the ringlet up faster than 3×10⁻⁸ of
self-gravity can hold it. The 1000-streamline π-twist failure starts inside the ringlet (streamlines
460–860), not at the edges: the raggedness grows ×2 every ~2.5 orbits while the periapse span stays
at 2.95 rad. With 241 particles per streamline it grows at the same rate from a smaller start; with
ring gravity off it does not grow. So it is gravitational, not longitude sampling (unlike issue 3).

**Why 1000 streamlines fail.** A pressureless self-gravitating ringlet is a cold disk, unstable at all
wavelengths shorter than λ_crit = 4π²Gσ/κ²; only the streamline spacing Δ hides them. Here
λ_crit/2Δ ≈ 0.94 at 1000 streamlines and ≈ 0.09 at 100. A uniform-lattice analysis of the code's
forces (P6, P7 in `notes/equations.md`) says axisymmetric modes are still stable at 1000 streamlines,
because the discrete self-gravity at the grid scale is half the continuum value; the failure is
non-axisymmetric, so that analysis does not cover it. The self-gravity that holds the twist against
J2 scales as M/Δa², but λ_crit/2Δ ≈ πGM N_r/(aΔa²) scales the same way times N_r, so for a given
twist-holding strength, enough streamlines always resolve the instability.

**Ring pressure does not help** (`explore/pressure/`, Q_ring = 0.5, 1, 1.2, 2): the 1000-streamline
run crosses at steps 361, 366, 369, 384 (361 without pressure), and the 100-streamline run is
unchanged to three digits. At Q ≈ 1 the Toomre wavelength 2c²/Gσ ≈ Q²Δ is at the grid scale, where
the coded pressure is weak, and it exerts no force on the zigzag (P6).

**The outer-edge flux (issue 4)** changes nothing visible when its longitudes are corrected
(`explore/taper/nr0100_pi_fixedge_own`), but hahn13's literal p(r − Δ) does
(`nr0100_pi_fixedge`): max q 0.30, width ×2.79, asymmetric edges.

## Open

- 300 streamlines (`explore/taper/nr0300_pi`, λ_crit/2Δ ≈ 0.28) also go unstable, more slowly:
  crossed at step 2355 (~190 orbits); the raggedness was already 3.6×10⁻³ spacings at 150 orbits,
  versus ≤ 1.2×10⁻⁴ for 10⁴ orbits at 100 streamlines.
- So the nominal ringlet tops out between 100 and 300 streamlines, below the 400–800 that bending
  waves need; bending-wave runs have different surface densities, so this limit is specific to it.
- Whether the 1000-streamline failure is physical (Toomre/swing amplification in a cold, eccentric,
  twisted ringlet) or a property of the discrete scheme is not settled.
