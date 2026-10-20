# Equation verification tracker

Status of every derived equation: **derived → sympy-verified (CC) → hand-verified (JMH)**.
An equation may be used in the code once sympy-verified; all must be hand-verified before the
paper is submitted.

Notation: Ω_a, κ_a, ν_a are the orbital, epicyclic and vertical frequencies at semimajor axis a;
Ω is the longitude of ascending node. M_p and R_p are the planet's mass and radius (M is the
mean anomaly).

Sympy checks: `notes/sympy/vertical_motion.py` (V1–V7, H1–H3) and `notes/sympy/ring_pressure.py`
(P1–P8) and `notes/sympy/viscous_overstability.py` (O1–O5); run with `uv run python notes/sympy/<file>.py`.

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
| O1 | cold viscous self-gravitating sheet, κ = Ω, constant ν (β = 0): s³ + (ν_b + 7ν_s/3)k²s² + [κ² − 2πGσk + (ν_b + 4ν_s/3)ν_sk⁴]s + ν_sk²(3Ω² − 2πGσk) = 0 | viscous overstability (paper App. B, planned) | ✅ | ✅ | ☐ | Eulerian; perturbation of the background shear stress ∝ (1 + β) |
| O2 | the same in Lagrangian (Hill) form: ξ̈ − 2Ωẏ − 3Ω²ξ = gξ − (ν_b + 4ν_s/3)Kξ̇, ÿ + 2Ωξ̇ = −ν_sK(ẏ + 3Ωξ), with K = k², g = 2πGσk | as O1 | ✅ | ✅ | ☐ | equals O1 |
| O3 | growth rate to first order in ν: γ = K[ν_s(4x + 2) − 3ν_b(1 − x)]/[6(1 − x)], x = g/κ²; stable iff ν_b/ν_s > (4x + 2)/3(1 − x) | as O1 | ✅ | ✅ | ☐ | marginal x = 1/7 for ν_b = ν_s, 2/5 for ν_b = 2ν_s; ν_b/ν_s > 2/3 without gravity. With β kept and no gravity, onset at β_c = (ν_b/ν_s − 2/3)/3, exactly the criterion of Lehmann et al. (2017, Sec. 2) and Schmit & Tscharnuter (1995) (check O3b) |
| O4 | the code's viscous stencils on a lattice mode ∝ e^{ijφ}: A_θ = −ν_s(sin²φ/Δ²)(ẏ + 3Ωξ), A_r = −(ν_b + 4ν_s/3)(sin²φ/Δ²)ξ̇ | stability of the scheme | ✅ | ✅ | ☐ | so O2–O3 hold with K = sin²φ/Δ², = 0 at the zigzag φ = π; v_r/r term (O(Δ/a)) dropped; numerically confirmed by `notes/numerics/viscous_stability.py` (0.21 vs 0.23 per orbit) |
| O5 | at the zigzag x = X(π)/κ² = 𝓡/2 (P7, P8); so the lattice is stable for 𝓡 < 2/7 with ν_b = ν_s, 𝓡 < 4/5 with ν_b = 2ν_s | as O1 | ✅ | ✅ | ☐ | lattice modes with sin²φ > 0 feel x < 𝓡/2 |
