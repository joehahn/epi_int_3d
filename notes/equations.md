# Equation verification tracker

Status of every derived equation: **derived → sympy-verified (CC) → hand-verified (JMH)**.
An equation may be used in the code once sympy-verified; all must be hand-verified before the
paper is submitted.

Notation: Ω_a, κ_a, ν_a are the orbital, epicyclic and vertical frequencies at semimajor axis a;
Ω is the longitude of ascending node. M_p and R_p are the planet's mass and radius (M is the
mean anomaly).

Sympy checks: `notes/sympy/vertical_motion.py` (V1–V6, H1–H3); run with
`uv run python notes/sympy/vertical_motion.py`.

| ID | Equation | Where used | Derived | Sympy | Hand | Notes |
|---|---|---|---|---|---|---|
| V1 | ν_a² = (GM_p/a³)[1 + (9/2)J₂(R_p/a)²] (vertical frequency, ∂²Φ/∂z² at z=0) | drift | ✅ | ✅ | ☐ | satisfies κ_a² + ν_a² = 2Ω_a² |
| V2 | z = a i sin(λ − Ω), v_z = a ν_a i cos(λ − Ω), with λ = ω̃ + (Ω_a/κ_a)M | drift (elements → coords) | ✅ | ✅ | ☐ | first-order vertical analog of hahn25 Eq. A9 |
| V3 | i = √[(z/a)² + (v_z/aν_a)²], λ − Ω = atan2(z, v_z/ν_a) | drift (coords → elements) | ✅ | ✅ | ☐ | exact inverse of V2 |
| V4 | ΔΩ = (Ω_a − ν_a)Δt | drift | ✅ | ✅ | ☐ | nodal regression ≈ −(3/2)J₂(R_p/a)²Ω_a |
| V5 | Φ(r,z) = Φ(r,0) + ½ν²(r)z² + O(z⁴), ν²(r) = (GM_p/r³)[1 + (9/2)J₂(R_p/r)²] | derivation of the vertical EOM | ✅ | ✅ | ☐ | Φ is even in z |
| V6 | z̈ + ν_a²z ≃ (3GM_p/a³)[1 + (15/2)J₂(R_p/a)²] δ z, with r = a(1 + δ) | accuracy (O(e·i) residual) | ✅ | ✅ | ☐ | equals −a(dν²/da)δz |
| H1 | in-plane solution (hahn25 Eq. A9) conserves h = r v_θ = a²Ω_a exactly | paper Appendix A | published | ✅ | n/a | transcription check |
| H2 | a = R_p[C + √(C² − 3J₂/2)], C = h²/2GM_pR_p | paper Appendix A, step 1 | published | ✅ | n/a | transcription check |
| H3 | e sin M = (v_r/aκ_a)(a/r)³, e cos M = ½[(a/r)² − 1] | paper Appendix A, step 3 | published | ✅ | n/a | transcription check |
