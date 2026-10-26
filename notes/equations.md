# Equation verification tracker

Status of every derived equation: **derived → sympy-verified (CC) → hand-verified (JMH)**.
An equation may be used in the code once sympy-verified; all must be hand-verified before the
paper is submitted.

Notation: Ω_a, κ_a, ν_a are the orbital, epicyclic and vertical frequencies at semimajor axis a;
Ω is the longitude of ascending node. M_p and R_p are the planet's mass and radius (M is the
mean anomaly).

Sympy checks: `notes/sympy/vertical_motion.py` (V1–V8, H1–H3) and `notes/sympy/ring_pressure.py`
(P1–P8) and `notes/sympy/viscous_overstability.py` (O1–O5) and `notes/sympy/streamline_gravity.py` (G1–G17); run with `uv run python notes/sympy/<file>.py`.

| ID | Equation | Where used | Derived | Sympy | Hand | Notes |
|---|---|---|---|---|---|---|
| V1 | ν_a² = (GM_p/a³)[1 + (9/2)J₂(R_p/a)²] (vertical frequency, ∂²Φ/∂z² at z=0) | drift | ✅ | ✅ | ☐ | satisfies κ_a² + ν_a² = 2Ω_a² |
| V2 | z = a i sin(λ − Ω), v_z = a ν_a i cos(λ − Ω), with λ = ω̃ + (Ω_a/κ_a)M | drift (elements → coords) | ✅ | ✅ | ☐ | first-order vertical analog of hahn25 Eq. A9 |
| V3 | i = √[(z/a)² + (v_z/aν_a)²], λ − Ω = atan2(z, v_z/ν_a) | drift (coords → elements) | ✅ | ✅ | ☐ | exact inverse of V2 |
| V4 | ΔΩ = (Ω_a − ν_a)Δt | drift | ✅ | ✅ | ☐ | nodal regression ≈ −(3/2)J₂(R_p/a)²Ω_a |
| V5 | Φ(r,z) = Φ(r,0) + ½ν²(r)z² + O(z⁴), ν²(r) = (GM_p/r³)[1 + (9/2)J₂(R_p/r)²] | derivation of the vertical EOM | ✅ | ✅ | ☐ | Φ is even in z |
| V6 | z̈ + ν_a²z ≃ (3GM_p/a³)[1 + (15/2)J₂(R_p/a)²] δ z, with r = a(1 + δ) | accuracy (O(e·i) residual) | ✅ | ✅ | ☐ | equals −a(dν²/da)δz |
| V8 | drift of V2–V4 as a rotation: z′ = z cos(ν_aΔt) + (v_z/ν_a) sin(ν_aΔt), v_z′ = v_z cos(ν_aΔt) − ν_a z sin(ν_aΔt) | drift (`bookkeeping.vertical_drift`) | ✅ | ✅ | ☐ | λ − Ω advances by Ω_aΔt − (Ω_a − ν_a)Δt = ν_aΔt; equal to V2–V4 exactly; z = v_z = 0 stays 0 |
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
| G1 | straight wire: A_g = 2Gλ D/\|D\|², D = d − (d·T̂)T̂, d = (Δ_r, 0, Δ_z), T̂ = (ε, 1, η)/√(1 + ε² + η²), ε = v_r′/v_θ′, η = v_z′/v_θ′ | vertical streamline gravity (planned kernel) | ✅ | ✅ | ☐ | 3D hahn13 Fig. 3 (`notes/phase2/streamline_gravity.md`; paper Sec. 3); integral of the point-mass pull along the wire |
| G2 | to first order in ε, η, any Δ_z/Δ_r: A_r = 2GλΔ_r/δ², A_z = 2GλΔ_z/δ², A_θ = −2Gλ(εΔ_r + ηΔ_z)/δ², δ² = Δ_r² + Δ_z² | as G1 | ✅ | ✅ | ☐ | A_θ has no O(slope²) term; A_r, A_z do (printed by the check) |
| G3 | planar limit Δ_z = η = 0: exact wire A = (2Gλ/Δ_r)(1, −ε, 0); hahn13 Eq. 20's (2Gλ/Δ_r)(cos φ, −sin φ) is that /√(1 + ε²) | consistency with `kernels.ring_gravity` | ✅ | ✅ | ☐ | agree to O(e²) |
| G4 | continuum sheet Z = Z₀ cos kr, \|Δ_z\| ≪ \|Δ_r\|: ΣA_z = 2GσZ₀∫(cos kx − 1)/x² dx = −2πGσ\|k\|Z | check of G2 against bending-wave theory | ✅ | ✅ | ☐ | the self-gravity term of (ω − mΩ)² = ν² + 2πGσ\|k\| |
| G5 | lattice Z_j = Z₀ cos jφ: A_z = (4GσZ₀/Δ)Σ(cos nφ − 1)/n² = −X(φ)Z₀, X of P7 | stability of the scheme, vertical | ✅ | ✅ | ☐ | ω² = ν² + X; → 2πGσk as φ → 0 |
| G6 | pattern rotating at Ω_p: v_r′/v_θ′ = (1 − Ω_p/Ω)(1/r)∂r′/∂θ, same for z | choice of the wire's direction | ✅ | ✅ | ☐ | velocity ≠ shape for m ≥ 2; the velocity tangent is part of the wire model and is kept |
| G7 | hahn13 form + tilts: A_r = A_g cos ψ cos φ, A_θ = −A_g(cos ψ sin φ + sin ψ sin χ), A_z = A_g sin ψ; A_g = 2Gλ/δ, tan ψ = Δ_z/Δ_r, sin φ = v_r′/v′, sin χ = v_z′/√(v_θ′² + v_z′²); coded as A = 2Gλ/(Δ_r(1 + q²)), A_r = A cos φ, A_θ = −A(sin φ + q sin χ), A_z = Aq | vertical streamline gravity (planned kernel); paper Sec. 3.1 | ✅ | ✅ | ☐ | = G1 + O(e′², i′², e′i′); hahn13 Eq. 20 exactly at Δ_z = 0 |
| G8 | vector form of G1: D = Δ − (Δ·v̂′)v̂′ = v̂′ × (Δ × v̂′), A_g = 2Gλ(Δ − (Δ·v̂′)v̂′)/(δ² − (Δ·v̂′)²); in the ring plane G7 = this × cos φ | exact result that G7 approximates; paper Sec. 3.1 | ✅ | ✅ | ☐ | BAC − CAB identity; kernel keeps G7 for bit-for-bit agreement with epi_int_lite |
| G9 | hahn13 notation: Δ = perpendicular to the wire, A_g = 2Gλ/Δ along it, A_g,i = A_g Δ_i/Δ; to first order Δ_r = r′ − r, Δ_z = z′ − z, Δ_θ = −(Δ_r v_r′ + Δ_z v_z′)/v_θ′ (from Δ·v′ = 0), Δ = √(Δ_r² + Δ_z²); planar Δ = Δ_r cos φ (cos φ, −sin φ, 0) exactly | paper Sec. 3.1, Fig. 1 | ✅ | ✅ | ☐ | same content as G7; epi_int_lite uses Δ_r for hahn13's Δ |
| G10 | wire along the shape: t̂′ ∝ (s_r, 1, s_z), s_r = ∂r′/r∂θ, s_z = ∂z′/r∂θ; Δ_θ = −Δ_r s_r − Δ_z s_z; planar Δ = Δ_r cos φ (cos φ, −sin φ, 0), tan φ = s_r; at an m-armed resonance 1 − Ω_p/Ω = ±1/m | paper Sec. 3.1, 3.3, Fig. 1 | ✅ | ✅ | ☐ | replaces the velocity tangent (G9) so m ≥ 2 waves are not off by ~m |
| G11 | exact straight wire, local frame: d = (r′ − r, 0, z′ − z), t = (s_r, 1, s_z), c = (d_r s_r + d_z s_z)/(1 + s_r² + s_z²), Δ = d − c t, Δ² = d_r² + d_z² − c(d_r s_r + d_z s_z), A = 2GλΔ/Δ²; planar A = (2Gλ/d_r)(1, −s_r, 0) | paper Sec. 3.1 (Eqns. 14–19: t̂, D = t̂×(d×t̂), A = 2Gλ′D/D², planar limit), planned kernel | ✅ | ✅ | ☐ | only approximations: straight infinite wire, \|r′ − r\|, \|z′ − z\| ≪ r |
| G12 | hahn13 Fig. 3 raised by d_z (exact): Δ = (Δ_h cos φ, −Δ_h sin φ, d_z), Δ_h = d_r cos φ, Δ² = d_r² cos²φ + d_z², tan φ = s_r; tilted by s_z: Δ_θ −= d_z s_z; A_r = 2Gλ d_r cos²φ/Δ², A_θ = −2Gλ(d_r sin φ cos φ + d_z s_z)/Δ², A_z = 2Gλ d_z/Δ² | paper Sec. 3.1 (Eqns. 16–17), Fig. 1; planned kernel | ✅ | ✅ | ☐ | = G11 up to O(s_r s_z, s_z²); replaced in the paper by the exact vector form G8/G11 (2026-10-25) |
| G13 | planned kernel, distance at the crossing point: Δ² = d_r² + d_z², A_r = 2Gλ′d_r cos φ/Δ², A_θ = −2Gλ′(d_r sin φ + d_z s_z)/Δ², A_z = 2Gλ′d_z/Δ² | paper Sec. 3.3 (Eqn. 20); Phase 3 step 2 kernel | ✅ | ✅ | ☐ | = G12 + O(slope²); planar = epi_int_lite (2Gλ/d_r)(cos φ, −sin φ); φ from v′ by default, s_z from the shape |
| G14 | 2D kernel options: shape slope s_r = (∂r′/∂θ)/r′ from the Lagrange weights' derivatives dl₀ = (dx₁ + dx₂)/(dx₁₀dx₂₀), dl₁ = −(dx₀ + dx₂)/(dx₁₀dx₂₁), dl₂ = (dx₀ + dx₁)/(dx₂₀dx₂₁); `nearpoint`: A = (2Gλ/d_r)(1, −tan φ), the exact planar wire (G3); otherwise (2Gλ/d_r)(cos φ, −sin φ); tan φ = s_r (`shape_slope`) or v_r′/v_θ′ | `fast/kernels.py` (`nearpoint`, `shape_slope`); epi_int_lite issue 9 | ✅ | ✅ | ☐ | r′ vs r in s_r differs at O(d_r/r), which the wire model drops |
| G15 | accuracy for a curved streamline (curvature 1/R): near point at u = −p + κpq, so D = d − (d·t̂)t̂ is rotated by (d·t̂)/R and its length off by (d·t̂)²/2RD; a circular ring pulls with (2Gλ/D)[1 ± (D/2R)(ln(8R/D) − 2)] (+ outside), by quadrature | notes only (removed from paper Sec. 3.1, 2026-10-25) | ✅ | ✅ | ☐ | both O(D/R): the wire approximation's own error |
| G16 | paper's solve: Δ = d₀ + d₁(δθ′) to the nearest point, d₁ ≃ N r δθ′ t̂; Δ·t̂ = 0 at δθ′ = −c/r, c = (d₀·t̂)/N = (d_r s_r + d_z s_z)/N²; Δ = d₀ − (d₀·t̂)t̂ = t̂×(d₀×t̂) = (d_r − c s_r, −c, d_z − c s_z), Δ² = d_r² + d_z² − c(d_r s_r + d_z s_z) | paper Sec. 3.1 (Eqns. 14–18) | ✅ | ✅ | ☐ | = G11 |
| G17 | the recipe: Lagrange value and slope weights w_k = x_i x_j/((x_k - x_i)(x_k - x_j)), w'_k = -(x_i + x_j)/((x_k - x_i)(x_k - x_j)) with x_k = theta_k - theta, exact for a quadratic; A_g = 2G lambda' (d_r - c s_r, -c, d_z - c s_z)/Delta^2 | paper Sec. 3.1, Eqns. 20-21 | ✅ | ✅ | ☐ | the code's planar kernel uses the same weights (G14); the z terms are not yet coded |
