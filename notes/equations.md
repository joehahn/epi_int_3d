# Equation verification tracker

Status of every derived equation: **derived → sympy-verified (CC) → hand-verified (JMH)**.
An equation may be used in the code once sympy-verified; all must be hand-verified before the
paper is submitted.

Notation: Ω_a, κ_a, ν_a are the orbital, epicyclic and vertical frequencies at semimajor axis a;
Ω is the longitude of ascending node. M_p and R_p are the planet's mass and radius (M is the
mean anomaly).

Sympy checks: `notes/sympy/vertical_motion.py` (V1–V4); run with
`uv run python notes/sympy/vertical_motion.py`.

| ID | Equation | Where used | Derived | Sympy | Hand | Notes |
|---|---|---|---|---|---|---|
| V1 | ν_a² = (GM_p/a³)[1 + (9/2)J₂(R_p/a)²] (vertical frequency, ∂²Φ/∂z² at z=0) | drift | ✅ | ✅ | ☐ | satisfies κ_a² + ν_a² = 2Ω_a² |
| V2 | z = a i sin(λ − Ω), v_z = a ν_a i cos(λ − Ω), with λ = ω̃ + (Ω_a/κ_a)M | drift (elements → coords) | ✅ | ✅ | ☐ | first-order vertical analog of hahn25 Eq. A9 |
| V3 | i = √[(z/a)² + (v_z/aν_a)²], λ − Ω = atan2(z, v_z/ν_a) | drift (coords → elements) | ✅ | ✅ | ☐ | exact inverse of V2 |
| V4 | ΔΩ = (Ω_a − ν_a)Δt | drift | ✅ | ✅ | ☐ | nodal regression ≈ −(3/2)J₂(R_p/a)²Ω_a |
