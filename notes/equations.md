# Equation verification tracker

Status of every derived equation: **derived → sympy-verified (CC) → hand-verified (JMH)**.
An equation may be used in the code once sympy-verified; all must be hand-verified before the
paper is submitted.

Sympy checks: `notes/sympy/vertical_motion.py` (V1–V4); run with
`uv run python notes/sympy/vertical_motion.py`.

| ID | Equation | Where used | Derived | Sympy | Hand | Notes |
|---|---|---|---|---|---|---|
| V1 | ν²(a) = (GM/a³)[1 + (9/2)J₂(R/a)²] (vertical frequency, ∂²Φ/∂z² at z=0) | drift | ✅ | ✅ | ☐ | satisfies κ² + ν² = 2Ω² |
| V2 | z = a i sin(λ − Ω̂), v_z = a ν i cos(λ − Ω̂), with λ = ω̃ + (Ω/κ)M | drift (elements → coords) | ✅ | ✅ | ☐ | first-order vertical analog of hahn25 Eq. A9 |
| V3 | i = √[(z/a)² + (v_z/aν)²], λ − Ω̂ = atan2(z, v_z/ν) | drift (coords → elements) | ✅ | ✅ | ☐ | exact inverse of V2 |
| V4 | ΔΩ̂ = (Ω − ν)Δt | drift | ✅ | ✅ | ☐ | nodal regression ≈ −(3/2)J₂(R/a)²Ω |
