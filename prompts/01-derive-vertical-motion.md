# 01 — Derive the vertical motion

`papers/hahn25.pdf` is my paper about the epi_int_lite N-body code. That code assumes all
rings and particles are coplanar, i.e. z = 0 always. epi_int_lite's physics is also detailed
in `papers/hahn13.pdf`; the most important part is Section 2.2, Eqs. 9, which describe the
particles' in-plane motions r, θ and their horizontal velocities v_r, v_θ.

I want to extend that solution, and the N-body code, into the vertical direction so that
particles and rings can have vertical displacement z ≠ 0 and vertical velocity v_z ≠ 0.
As a first task, derive expressions for z and v_z that are analogous to Eqs. 9 in hahn13.
