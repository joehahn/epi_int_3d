# Phase 2c: satellite vertical forcing and linear bending-wave theory

Derived in Goldreich & Tremaine notation and checked by `notes/sympy/bending_wave.py`
(B1–B10 in `notes/equations.md`; all pass, 2026-10-31). This is the theory for the Phase 4 test
"linear bending wave from a vertical resonance". Shu's papers (`papers/shu83.pdf`, `papers/shu84.pdf`) are
scans without a text layer, so this is derived from first principles. Read page by page afterwards (2026-10-31), it agrees with Shu (1984): his vertical forcing f_M (Eq. A21) is A_V, his Eq. 88, iε dh/dx + xh = L_V, is B8's equation (complex-conjugate convention, ε = 2πGσ/(r_V D_V), L_V = f_M/D_V), and the L_V/x term of his Eqs. 91–92 is B6's test-particle forced response near r_V. The B3 check Fourier-analyzes an
inclined satellite's exact force, and B8 and B9 confirm the forced solution independently.

Notation: the satellite has mass M_s and a circular orbit with radius a_s, inclination i_s (to the planet's
equator, the ring plane), longitude of ascending node Ω_s and mean longitude λ_s. Ω_a and ν_a are the angular
and vertical frequencies at radius a (V1). σ is the ring's surface density. x = (r − r_V)/r_V.

## 1. Forcing (B1–B4)

- **Vertical pull.** To first order in i_s, the satellite pulls a ring particle at (r, θ, z) vertically with
  F_z = GM_s(z_s − z)/D₀³, where z_s = a_s i_s sin(λ_s − Ω_s) and D₀² = r² + a_s² − 2ra_s cos(θ − λ_s) (B1).
- **The z_s term is the forcing.** The −z term is a small secular shift of ν (its k = 0 part) plus
  parametric, second-order terms. The Fourier satellite drops both.
- **Laplace expansion.** Expanding 1/D₀³ = (1/2a_s³) Σ_k b_{3/2}^{(k)}(β) cos k(θ − λ_s), with β = r/a_s,
  gives terms ∝ sin[(λ_s − Ω_s) + k(θ − λ_s)] (B2).
- **The resonant term is k = −m:**

  **F_z = −A_V sin(mθ − (m+1)λ_s + Ω_s),  A_V = GM_s i_s b_{3/2}^{(m)}(β)/(2a_s²)**

  The amplitude uses b_{3/2}^{(m)} itself, not its derivative.
- **No indirect term.** The indirect term is uniform in θ, so it contributes only to k = 0.
- **Numerical check (B3).** For m = 2, 3, 5, a 2D Fourier transform of the exact force of a satellite with
  i_s = 10⁻⁵ recovers A_V to 2×10⁻¹⁰–1.5×10⁻⁹.
- **Pattern speed.** λ_s advances at Ω_{a_s}, and the argument of latitude λ_s − Ω_s at ν_{a_s} (V8). Hence
  mΩ_p = mΩ_{a_s} + ν_{a_s} (B4).
- **Inner vertical resonance (IVR).** A ring particle feels the forcing at ω_D = m(Ω_a − Ω_p), so the IVR is
  at m(Ω_a − Ω_p) = ν_a, i.e. m(Ω_a − Ω_{a_s}) = ν_a + ν_{a_s}.
- **Keplerian limit.** The resonance condition becomes (m−1)Ω_V = (m+1)Ω_{a_s}, so
  r_V = a_s[(m−1)/(m+1)]^{2/3}. For example, m = 4 is a 5:3 resonance (Mimas 5:3).
- **m ≥ 2 only.** For m = 1 the IVR is where the ring's nodal regression matches 2Ω_{a_s} − Ω̇_s, far from
  any satellite, so there are no m = 1 bending waves of this kind.

## 2. Response (B5–B10)

- **Distance from resonance (B5).** Near r_V, ν² − ω_D² = D_V x, with D_V = 2ν_V r_V(dν/dr − m dΩ/dr) at r_V.
  - Keplerian: D_V = 3(m−1)Ω_V².
  - With J₂: D_V = (GM_p/r³)[3(m−1) + (3/2)(8m−15)J₂(R_p/r)²] + O(J₂²).
  - The code should evaluate D_V numerically from Ω_a and ν_a.
- **Test particles (B6).** Without self-gravity, the forced response is Z = F_z/(ν² − ω_D²) ≈ F_z/(D_V x),
  singular at r_V. This is the first check of the satellite code.
- **Test particles while the satellite grows (B11).** With mass ∝ 1 − e^{−t/τ} and z = v_z = 0 at t = 0, a
  particle's exact z(t) is the ramped forced response plus a free oscillation at ν. The free oscillation never
  decays, and its size relative to the forced response is ~1/(τ|ν − ω_D|). So particles settle into B6 where
  |x| ≫ 2ν_V/(D_V τ), and the zero-mass test compares each particle's z(t) with this exact solution everywhere.
- **Self-gravity (B7).** A sheet corrugated by Z = h cos kr pulls itself back with −2πGσ|k|Z. This gives the
  dispersion relation (ω_D² − ν²) = 2πGσ|k|.
  - Waves exist only inside r_V (x < 0). Bending waves at an IVR propagate inward, toward the planet and away
    from the satellite, the opposite direction to density waves at an inner Lindblad resonance.
  - Wavenumber: |k| = D_V|x|/(2πGσ). Wavelength: λ_b = 4π²Gσ/(D_V|x|).
- **Forced wave (B8).** For a trailing wave, |k|h = −i dh/dr. The forced equation becomes
  −i(2πGσ/r_V) dh/dx + D_V x h = F̂, where F_z = Re[F̂ e^{i(mθ − mΩ_p t)}] and F̂ = iA_V e^{iφ₀}. Its solution is

  **h(x) = −iKF̂ e^{−iLx²/2} ∫ₓ^∞ e^{iLs²/2} ds,  K = r_V/(2πGσ), L = D_V K**

  The integral is a Fresnel integral (scipy `fresnel`).
  - **Outside r_V (x ≫ L^{−1/2}):** h → F̂/(D_V x), the test-particle response.
  - **Inside r_V:** a trailing wave of constant amplitude |h| = A_V(r_V/(GσD_V))^{1/2} and local wavenumber
    k = −Lx/r_V, plus the local forced response.
  - **Resonance zone:** its width is |x| ~ L^{−1/2}.
- **Torque (B9).** The torque on the ring is T = −π²mσr_V²A_V²/D_V. It is negative (the ring loses angular
  momentum) and independent of the ring's self-gravity. This is a second, integral test: the ring's angular
  momentum should drop at this rate once the wave is established.
- **Group velocity (B10).** c_g = −πGσ/ω_D ≈ −πGσ/ν_V: inward, at constant speed. The front is at
  x_f ≈ −πGσ t/(ν_V r_V) after the forcing turns on. The run must stop before the front reaches the inner
  edge (Phase 2e), and the comparison uses only |x| < |x_f|.
- **Damping.** Not included. epi_int_3d's viscosity acts in the plane only (Phase 2d), so the simulated waves
  should be undamped, like this theory.

## 3. Scale of a test run (J₂ = 0.01, R_p = 0.5, r_V at 90% of the ring's width)

| m | width W | ring mass | a_s | L | resonance zone | crests in ring | shortest λ_b | streamlines for 10 per λ_b | front crosses in |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 0.004 | 10⁻⁹ | 2.09 | 1.2×10⁷ | 2.9×10⁻⁴ | 12 | 1.5×10⁻⁴ | ~270 | 4600 orbits |
| 2 | 0.01 | 10⁻⁸ | 2.10 | 3.0×10⁶ | 5.8×10⁻⁴ | 19 | 2.4×10⁻⁴ | ~420 | 2900 orbits |
| 3 | 0.004 | 10⁻⁹ | 1.59 | 2.4×10⁷ | 2.0×10⁻⁴ | 25 | 7.3×10⁻⁵ | ~550 | 4600 orbits |

The crest count scales as L x² ∝ W²/σ, and the streamlines needed as W·L·x ∝ W²/σ too. So a heavier ring
needs fewer streamlines and gives a faster front, at the cost of a wider resonance zone. Particles per
streamline can stay small (~10–20 for m = 2), since the pattern is m-armed.
- **Satellite inclination.** A_V sets |h|. Linear theory needs |kh| ≪ 1, which is largest at the inner edge.
  i_s, or M_s, should be chosen to keep max |kh| ≲ 0.1.

## 4. Inputs (epi_int_lite's `lindblad_resonance` branch, and what epi_int_3d adds)

On the branch, the satellite is a dict in `inputs.py`, the same file as the rest of the run (e.g.
`runs/resonance/inputs.py`):

```python
satellite = {'mass_final': 6.5969e-8, 'time_grow': 1.0e6, 'mass': 0.0, 'r': 1.58, 't': 0.0, 'm': 2}
```

- **Mass.** `mass_final < 0` turns the satellite off. The mass grows as mass_final(1 − e^{−t/time_grow}),
  so the wave starts smoothly.
- **Orbit.** `r` and `t` are the satellite's radius and longitude. The driver advances `t` at Ω(r) with J₂,
  and overwrites `mass` each step.
- **Laplace coefficients.** At setup, `lc` = b_{1/2}^{(m)} and `dlc` = its β-derivative are evaluated once,
  at the ring's mean radius. The derivative uses `scipy.misc.derivative`, which no longer exists in SciPy.
- **Forces.** In-plane only (the m-th Lindblad term, with the m = 1 indirect term), added in `accelerations`.

epi_int_3d uses the same dict with these additions (`fast/satellite.py`; `'vertical'` is implemented, `'horizontal'` and `'potential': 'direct'` not yet):

```python
#satellite: only the resonant m-th term of its potential (a "Fourier satellite"); mass_final < 0 turns it off
satellite = {
    'resonance': 'vertical',   #'vertical' (inner vertical resonance) or 'horizontal' (the branch's inner Lindblad)
    'm': 2,                    #number of spiral arms, m >= 2 for 'vertical'
    'mass_final': 1.0e-9,      #satellite mass/planet mass once grown
    'time_grow': 1.0e3,        #mass = mass_final*(1 - exp(-time/time_grow))
    'r': 2.0911,               #satellite's orbit radius a_s (circular)
    't': 0.0,                  #its initial mean longitude lambda_s
    'inclination': 1.0e-3,     #NEW: i_s, for 'vertical'
    'node': 0.0,               #NEW: its initial longitude of ascending node Omega_s
}
```

- **New keys.** `resonance`, `inclination` and `node`.
- **Satellite motion.** The driver advances λ_s at Ω_{a_s} and Ω_s at Ω_{a_s} − ν_{a_s} (nodal regression).
- **Vertical forcing.** It adds A_z = −A_V(t) sin(mθ − (m+1)λ_s + Ω_s) to the v_z kick, and nothing in the plane.
- **Laplace coefficient.** b_{3/2}^{(m)} is evaluated once, at r_V, which the driver solves for with J₂ and
  prints. That matches the theory, which needs b only at r_V.
- **Horizontal option.** `'horizontal'` ports the branch's force, with the analytic β-derivative of b_{1/2}^{(m)}.
- **Cross-check option.** An optional `'potential': 'direct'` replaces the Fourier term with the full point-mass
  satellite (all its terms) for the cross-check in PLAN.
- **Choosing `time_grow`.** The growing mass is the satellite's. A switch-on that is slow compared with the
  forcing's detuning, time_grow ≫ 1/|ν − ω_D| ≈ 2ν_V/(D_V|x|), leaves test particles in the purely forced
  state B6 instead of forced plus free. Near r_V, though, no growth is slow enough: the detuning there is
  set by self-gravity, and the free part becomes the wave that propagates away. So `time_grow` should also be
  short compared with the front's crossing time. Otherwise the wave's amplitude varies along it, following
  the mass at the retarded time t − r_V|x|/|c_g|. In the m = 2, W = 0.004 example, adiabatic at
  |x| ≳ 10⁻³ needs ≫ 700, and the crossing time is 2.9×10⁴, so time_grow ≈ 2–3×10³.
- **Run-time state.** The branch's run-time state (`mass`, `lc`, `dlc`) should not be an input: the driver
  keeps it in its own dict.
