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
≤ 1.2×10⁻⁴ spacings; width ×2.66; |ΔL/L₀| ≤ 1.3×10⁻⁵. The initial π twist does not persist:
within ~1000 orbits self-gravity straightens it, and the middle ~80% of the ringlet then stays nearly
aligned (e ≈ 0.4×10⁻⁴, periapses spread by ≲ 0.3π) instead of being wound up by J2. The tapered edges, the outer ~10% of streamlines on each side, carry larger eccentricities (up to
≳ 3×10⁻⁴) and librating periapses; the inner edge's periapse circulates once every ~4800 orbits. All the
large q is at the edges, mostly the outer one; the body stays below ~0.05. In the notebook's
animation (e, periapse and q versus a) that body is two sub-ringlets, the inner ~25% and the next ~65%,
precessing slowly relative to each other, with the outer edge decoupled. (Space-time diagrams,
surface-density maps and the optional animation in the notebook.) The 10⁴ orbits span 0.43 of the
viscous timescale τ_ν = Δa₀²/12ν_s (hahn25 Eq. 9), ≈ 2.3×10⁴ orbits.

**A 30% wider ringlet** (`explore/wide/`: width 0.01625 at the same surface density, mass 3.9×10⁻⁸, π
twist) also stays smooth for 10⁴ orbits (max q 0.41, smallest gap 0.98, raggedness ≤ 8×10⁻⁵, width
×2.19), but still forms only two sub-ringlets, not three, and spans only 0.25 τ_ν; so the nominal
width was kept.

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
non-axisymmetric, so that analysis does not cover it. A rough estimate that fits these runs: J2's differential precession across the
ringlet scales as J2 Δa, self-gravity's (for small q) as m_r Δϖ/Δa², so keeping J2 from winding the
ringlet up needs m_r ∝ J2 Δa³/Δϖ, and then λ_crit/2Δ ∝ J2 Δa N_r/Δϖ. A wider ringlet does not help;
a narrower one (with smaller e, to keep q and the interpolation error small) does.

**A 10× narrower ringlet** (`explore/narrow/`: width 0.00125, e = 10⁻⁵ so q starts at 0.025 as
before, π twist, ν = 7×10⁻¹³ for the same ×2 spreading, J2 still 0.01; the scaling gives mass ≈ 3×10⁻¹¹).
"Body twist" is the periapse span across the middle 80% of the streamlines:

| run | result | max q | body twist (rad), 0 → 4000 → 10⁴ orbits | body e at end |
|---|---|---|---|---|
| 100 streamlines, 10⁻¹¹ | smooth for 10⁴ orbits | 0.026 | 2.0 → 1.6 → 0.9 | 0.97×10⁻⁵ |
| 100 streamlines, 3×10⁻¹¹ | smooth for 10⁴ orbits | 0.026 | 2.0 → 1.5 → 0.5 | 0.87×10⁻⁵ |
| 100 streamlines, 10⁻¹⁰ | smooth for 10⁴ orbits | 0.11 | 2.0 → 1.05 → 1.6 (librates) | 0.81×10⁻⁵ |
| 1000 streamlines, 3×10⁻¹¹ | no crossing in 159 orbits | 0.027 | 2.0 throughout | 1.0×10⁻⁵ |

All 100-streamline runs widen ×2.65 with raggedness ≤ 4.5×10⁻⁵ spacings. At 1000 streamlines and
3×10⁻¹¹ (λ_crit/2Δ ≈ 0.09) the raggedness falls from 1.2×10⁻⁴ to 0.8×10⁻⁴ spacings over 159 orbits,
where the wide ringlet at 300 streamlines was already growing. So the narrow ringlet is the
1000-streamline candidate; its 10⁴-orbit run (~28 h with the reference code) waits for the faster
code. The catch: J2's differential precession scales with the width, so across this ringlet it alone
would twist it by only ~1 rad in 10⁴ orbits (versus ~10 rad for the nominal), and at 3×10⁻¹¹
self-gravity unwinds the twist slowly instead of aligning the body within ~1000 orbits. At 10⁻¹⁰
self-gravity controls the twist, but that is λ_crit/2Δ ≈ 0.31 at 1000 streamlines, where the wide
ringlet went unstable at 300; untested.

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
