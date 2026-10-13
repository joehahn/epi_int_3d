# 02 — Make the vertical equations first-order and exactly reversible

epi_int_lite uses first-order equations, per Appendix A of `papers/hahn25.pdf`. Eq. A9 gives
positions and velocities (coordinates) as functions of the orbit elements a, e, ω̃, M, and
steps 1–4 invert those to give the elements as functions of the coordinates. The two are
exactly reversible: converting coordinates into elements and back involves no approximation
and introduces no error.

Can the same be achieved when extending epi_int_lite to first-order equations for z and v_z?
Derive them, and test that the round trip is exact and that the errors in the motion are of
the expected order.
