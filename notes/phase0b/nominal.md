# epi_int_3d's nominal ringlet

A ringlet for testing this project's code and starting the bending-wave runs (`tests/nominal/`), found by the parameter
search in `explore/` (2026-10-18/19). The goal was a ringlet that is eccentric, has a periapse twist that
makes its self-gravity non-axisymmetric, is massive enough that self-gravity rather than J2 controls
the twist, spreads viscously by about ×2, and stays smooth for 10⁴ orbits with as many streamlines as
possible, ideally 1000.

All runs use the reference implementation (`source/`) unchanged, through `explore/run.py`, which
applies `tools/patches.py` at run time: always the longitude-wrapped crossing check (issue 1 in
`epi_int_lite_issues.md`), and optionally an edge taper. Diagnostics are from `explore/diagnose.py`:
q = a|Δ(e e^{iϖ})|/Δa, the largest over adjacent streamline pairs (they touch when q → 1); the
smallest gap between adjacent streamlines over the initial spacing; the raggedness, the largest rms
residual of a streamline's r(θ) about a Fourier fit through m = 3, in spacings.

## The nominal ringlet (since 2026-10-19, evening)

The 100-streamline viscous ringlet below could not be run at more than ~300 streamlines (it crosses
at 333, from the m ≈ Nt/2 tangential sawtooth along each streamline and the viscous terms; see issues 3
and 6 in `notes/epi_int_lite_issues.md`). The bending-wave runs need ~1000 streamlines and will start
with viscosity and pressure off, so the nominal ringlet is now a narrower, lighter ringlet that is well
resolved at 1000 streamlines (screened in `explore/nominal1000/`, all four candidates smooth to 10³
orbits; the chosen one is `explore/nominal1000/long_nt20_m3e-11`, whose output is the test's).

| parameter | value |
|---|---|
| streamlines × particles | 1000 × 20 |
| width Δa₀ | 0.00125, edges tapered (20% each side) |
| mass | 3×10⁻¹¹ (λ_crit/2Δ ≈ 0.09 mid-ringlet) |
| e, a de/da | 10⁻⁵, 0 |
| periapse twist | none (`w_prime` = 0) |
| ν_s, ν_b, pressure | off |
| J2, R_p, dt | 0.01, 0.5, 0.5 |
| run | 5×10⁴ orbits (629 000 steps, 500 outputs), 2.8 h with `fast/epi_int_fast.py` (the reference driver would take days) |

Result (`tests/nominal/nominal.ipynb`): no crossing or NaN; max q 0.20; smallest gap 0.58 of the
spacing (at ~3×10⁴ orbits, 0.83 at the end); |ΔL/L₀| ≤ 1.7×10⁻¹¹; width constant to 3×10⁻⁶. J2 alone
would twist the middle 80% of the ringlet by ~5 rad over the run (differential precession
≈ 0.08 rad per 10³ orbits); self-gravity holds the twist within 0.38 rad, the twist and an eccentricity
gradient (e from ~0.5 to ~1.8×10⁻⁵ across the ringlet) oscillating with a period of ~2–4×10⁴ orbits, the
ringlet staying one aligned body with no sub-ringlets. The raggedness grows from round-off (~10⁻⁹
spacings) by ~10⁷ over the first 2×10⁴ orbits, then levels off at ≤ 0.05 spacings (median ~5×10⁻³ at
the end). The largest q sits in narrow streaks that start near the taper boundaries, (a − ā₀)/Δa₀ ≈ ±0.3,
and run inward across the ringlet, short-wavelength waves launched where the surface density profile
bends. Neither is a threat over 5×10⁴ orbits, but both are noted for the bending-wave runs.

## The former nominal ringlet (100 streamlines, viscous; now the benchmark's)

Until the evening of 2026-10-19 the nominal ringlet was the one below; it is frozen in
`benchmark/nominal_inputs.py` so that the nominal benchmark results stay valid. From earlier that day it
had no initial periapse twist; the π-twisted ringlet found first (below) is kept as a variant.

| parameter | value |
|---|---|
| streamlines × particles | 100 × 40 |
| width Δa₀ | 0.0125, edges tapered (below) |
| mass | 3×10⁻⁸ |
| e, a de/da | 10⁻⁴, 0 |
| periapse twist | none (`w_prime` = 0) |
| ν_s = ν_b | 7×10⁻¹¹ |
| pressure, fast_gravity | off, off |
| run | 5×10⁴ orbits (629 000 steps, 500 outputs), 16 min with `fast/epi_int_fast.py` |

The 5×10⁴ orbits are 2.2 viscous timescales τ_ν = Δa₀²/12ν_s (hahn25 Eq. 9) ≈ 2.3×10⁴ orbits, with Δa₀
the initial distance between the edge streamlines.

Result (`tests/nominal/nominal.ipynb`): no crossing or NaN; max q 0.38 (at ~10⁴ orbits, lower after);
smallest gap 0.996 of the spacing; raggedness ≤ 4.3×10⁻⁴ spacings, creeping up over the last 10⁴
orbits; width ×5.35; |ΔL/L₀| ≤ 6.5×10⁻⁵. J2 starts to twist the ringlet, but self-gravity holds the body to a small, smooth
twist (~0.1–0.3π) with no sub-ringlets, while an eccentricity gradient builds up: e grows outward,
to ~6×10⁻⁴ at the outer edge and ~10⁻⁵ at the inner edge by 10⁴ orbits. The body then circularizes
(e ≈ 0 by ~3×10⁴ orbits, after which its periapses are meaningless), and only the outer few streamlines
stay eccentric (e ≈ 4–5×10⁻⁴), carrying all the large q. At 333 streamlines (`explore/twist0_n333`,
λ_crit/2Δ ≈ 0.31) it crosses at step 1825 (~145 orbits), like the π-twisted ringlet at 300, so the
limit is resolution, not the twist. The reference code (`explore/twist0_n333_ref`, 3.0 min versus
0.18 min with `fast/epi_int_fast.py`) crosses at the same step, with all outputs bit-for-bit identical. Ring pressure does not prevent it
(`explore/twist0_n333_Q`): Q_ring = 0.5, 1, 2, 10 cross at steps 1829, 1835, 1836, 2210 (~176 orbits).
The particles per streamline Nt matter, not monotonically (`explore/twist0_n333_Nt`, to 10³ orbits):
Nt = 20, 40, 320 cross at ~375, 145, 473 orbits, while Nt = 30, 50, 60, 80, 160 do not, their
raggedness levelling off at ~5×10⁻³, 10⁻³, 8×10⁻⁴, 3×10⁻⁴, 4×10⁻⁵ spacings. So 333 streamlines are
marginal, and whether a run crosses depends on Nt; 100 streamlines stay well clear.

**The π-twisted variant** (`w_prime` = π e ā/Δa₀ = 0.02529, 10⁴ orbits; `explore/taper/nr0100_pi`):
no crossing; max q 0.39; smallest gap 0.97; raggedness ≤ 1.2×10⁻⁴ spacings; width ×2.66;
|ΔL/L₀| ≤ 1.3×10⁻⁵. Within ~1000 orbits self-gravity straightens the twist, and the middle ~80% of the
ringlet then stays nearly aligned (e ≈ 0.4×10⁻⁴, periapses spread by ≲ 0.3π), as two sub-ringlets, the
inner ~25% and the next ~65%, precessing slowly relative to each other. The tapered edges carry larger
eccentricities (≳ 3×10⁻⁴) and librating periapses; the inner edge's periapse circulates once every
~4800 orbits. All the large q is at the edges. With zero twist (`explore/twist0`, 10⁴ orbits) the
same ringlet has no sub-ringlets: max q 0.35, raggedness ≤ 5.7×10⁻⁵, width ×2.65. So the sub-ringlets
are left over from the initial twist.

**A 30% wider π-twisted ringlet** (`explore/wide/`: width 0.01625 at the same surface density, mass
3.9×10⁻⁸) also stays smooth for 10⁴ orbits (max q 0.41, smallest gap 0.98, raggedness ≤ 8×10⁻⁵, width
×2.19), but still forms only two sub-ringlets, not three, and spans only 0.25 τ_ν.

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
