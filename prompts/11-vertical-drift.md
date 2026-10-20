# 11 — The vertical drift in epi_int_3d's driver

(Assembled from prompts on 2026-10-20; it is the final version of what was asked.)

Do the next step of the plan: add the vertical motion (Phase 2a) to `fast/epi_int_3d.py`, starting
with the drift. With z ≡ 0 the driver must reproduce the 2D results exactly. Test that inclined orbits
keep their inclination and that their nodes regress at the J2 rate, and compare with an exact
integration in the 3D J2 potential. No animations unless a result needs one.
