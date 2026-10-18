# Equation verification tracker

Status of every derived equation: **derived → sympy-verified (CC) → hand-verified (JMH)**.
An equation may be used in the code once sympy-verified; all must be hand-verified before the
paper is submitted.

Notation: Ω_a, κ_a, ν_a are the orbital, epicyclic and vertical frequencies at semimajor axis a;
Ω is the longitude of ascending node. M_p and R_p are the planet's mass and radius (M is the
mean anomaly).

Sympy checks: `notes/sympy/vertical_motion.py` (V1–V7, H1–H3) and `notes/sympy/ring_pressure.py`
(P1–P8); run with `uv run python notes/sympy/<file>.py`.

| ID | Equation | Where used | Derived | Sympy | Hand | Notes |
|---|---|---|---|---|---|---|
| V1 | ν_a² = (GM_p/a³)[1 + (9/2)J₂(R_p/a)²] (vertical frequency, ∂²Φ/∂z² at z=0) | drift | ✅ | ✅ | ☐ | satisfies κ_a² + ν_a² = 2Ω_a² |
| V2 | z = a i sin(λ − Ω), v_z = a ν_a i cos(λ − Ω), with λ = ω̃ + (Ω_a/κ_a)M | drift (elements → coords) | ✅ | ✅ | ☐ | first-order vertical analog of hahn25 Eq. A9 |
| V3 | i = √[(z/a)² + (v_z/aν_a)²], λ − Ω = atan2(z, v_z/ν_a) | drift (coords → elements) | ✅ | ✅ | ☐ | exact inverse of V2 |
| V4 | ΔΩ = (Ω_a − ν_a)Δt | drift | ✅ | ✅ | ☐ | nodal regression ≈ −(3/2)J₂(R_p/a)²Ω_a |
| V5 | Φ(r,z) = Φ(r,0) + ½ν²(r)z² + O(z⁴), ν²(r) = (GM_p/r³)[1 + (9/2)J₂(R_p/r)²] | derivation of the vertical EOM | ✅ | ✅ | ☐ | Φ is even in z |
| V6 | z̈ + ν_a²z ≃ (3GM_p/a³)[1 + (15/2)J₂(R_p/a)²] δ z, with r = a(1 + δ) | accuracy (O(e·i) residual) | ✅ | ✅ | ☐ | equals −a(dν²/da)δz |
| V7 | ∇²Φ = r⁻¹∂(r ∂Φ/∂r)/∂r + ∂²Φ/∂z² = 0, and κ_a² = 3Ω_a² + ∂²Φ/∂r², hence κ_a² + ν_a² = 2Ω_a² | frequency relation | ✅ | ✅ | ☐ | Laplace's equation outside the planet |
| H1 | in-plane solution (hahn25 Eq. A9) conserves h = r v_θ = a²Ω_a exactly | paper Appendix A | published | ✅ | n/a | transcription check |
| H2 | a = R_p[C + √(C² − 3J₂/2)], C = h²/2GM_pR_p | paper Appendix A, step 1 | published | ✅ | n/a | transcription check |
| H3 | e sin M = (v_r/aκ_a)(a/r)³, e cos M = ½[(a/r)² − 1] | paper Appendix A, step 3 | published | ✅ | n/a | transcription check |
| P1 | A_p = [p(r−Δ) − p(r)]/λ ≃ −(1/σ)∂p/∂r, σ = λ/Δ (hahn13 Eq. 23) | `A_P` | published | ✅ | n/a | transcription check |
| P2 | σ_j = 2λ/(r_{j+1} − r_{j−1}), A_j = −(∂p/∂r)_j/σ_j = −(p_{j+1} − p_{j−1})/2λ = −(1/σ)∂p/∂r + O(Δ²) (hahn13 Eqs. 24–25) | `surface_density`, `delta_f`, `A_P` | published | ✅ | ☐ | code matches the paper in the interior |
| P3 | net radial pressure force Σλ_jA_j = (−p₀ + p₁ + p_{N−2} − p_{N−1})/2 with edges A₀ = −p₀/λ, A_{N−1} = p_{N−2}/λ | `A_P` | ✅ | ✅ | ☐ | ≠ 0: Newton's 3rd law fails at O(Δ ∂p/∂r) per edge |
| P4 | A_r = A v_θ/v, A_θ = −A v_r/v is normal to the streamline | `ring_pressure` | published | ✅ | n/a | hahn13 uses v_θ ≈ v |
| P5 | ω² = κ² − 2πGσk + c²k² is stable iff Q = cκ/πGσ ≥ 1; code sets c = QπGσ/Ω | `initialize_streamline` | standard | ✅ | ☐ | Ω for κ is exact to O(J₂); σ is the ring-wide mean |
| P6 | coded pressure on a radial mode ξ_j ∝ cos(jkΔ): ω_p² = (c/Δ)² sin²(kΔ) | stability of the scheme | ✅ | ✅ | ☐ | → c²k² as Δ → 0; = 0 for the zigzag kΔ = π (odd–even decoupling) |
| P7 | coded gravity on that mode: X = (4Gσ/Δ)(πφ/2 − φ²/4), φ = kΔ; max π²Gσ/Δ at φ = π | stability of the scheme | ✅ | ✅ (P7 sum numerically) | ☐ | half the continuum 2πGσk at k = π/Δ; so ω² = κ² − X + ω_p² |
| P8 | 𝓡 = λ_crit/2Δ = 2π²Gσ/κ²Δ ≃ πGm_rN_r/κ²aΔa²; cold lattice axisymmetrically stable iff 𝓡 < 2 | paper App. B (Eqs. B3–B6) | ✅ | ✅ | ☐ | P6, P7 are Eqs. B4, B3 |
