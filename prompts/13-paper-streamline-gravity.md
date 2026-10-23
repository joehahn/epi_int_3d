# 13 — Putting the streamline-gravity derivation in the paper

(2026-10-21)

put derivation and figure in paper

Follow-up: the hahn13 derivation of Ag and corresponding figure seems much simpler than the 3d
derivation; can the new derivation be revised so that it is as simple as hahn13 + small modification
for the out of plane motions? and ditto for figure? to same level of accuracy ie errors of order
e'^2 and i'^2.

Follow-up: (after the user's own version, A_g = 2Gλ'/D toward the spot, A_r = A_g cos φ,
A_θ = −A_g sin φ, A_z = A_g sin Φ, was corrected by the cos Φ factor and the χ term) how can we make
the above more transparent in the paper and plots? and what is meridional plane?

Follow-up: suppose vec(Ag) = Ag hat(n) is the force on the particle at vec(r) due to the streamline at
vec(r') with velocity vec(v'), and vec(Delta) = vec(r') - vec(r). Is there a simple vector identity
that derives the force direction from r, r', v'? (Answer: D = Delta - (Delta.v)v = v x (Delta x v);
then: yes, add it to Section 3.1 as the exact result that the angle form approximates.)

Follow-up: Delta_r = separation between r' and r when both are projected to the z=0 plane? wait, I
thought Delta_r = Delta seen in Fig 3 of hahn13, which doesn't share longitude theta. Can the paper use
the notation and derivation and Fig 3 [of hahn13] as much as possible when extending to 3d?

Follow-up: should [the derivation of the nearest point, s = -d.v, from the crossing point and v'] be
added to paper?

Follow-up: is there a simpler way to derive section 3.1? (Answer: Delta_r, Delta_z are the
crossing-point separations to first order, and Delta.v' = 0 gives Delta_theta; then: yes.)

Follow-up: (on keeping the v_z' term) angular momentum transfer between corrugated streamlines is
probably important. Don't make any approximations that will be problematic for m >= 2 spiral waves.

Follow-up: wire approximation doesn't care about spiral wave theory. Assume infinitely long straight
wires having linear mass density lambda, and follow that assumption to its logical conclusion; ignore
spiral wave theory. The other reasonable approximation is that the perturbed particle is close to the
perturber, |r'-r| << r and |z'-z| << r. Then derive in the simplest/easiest way possible the
accelerations Ar, At, Az consistent with those approximations.

Follow-up: the 2013 paper approximates A_g,r = A_g and A_g,theta = -A_g v_r'/v_theta' (its Eqns. 20),
but epi_int_lite didn't make that approximation, correct? So in the new paper assume Eqn. (14), and
we'll deal with the inconsistency between the code and the 2013 paper later.

Follow-up: don't keep referring to the Lagrange polynomial for interpolation. Maybe early on say that
r' and v' at the site in question are determined by Lagrange polynomial, then thereafter write as if
r' and v' are known exactly.

Follow-up: I wonder if we can get an even easier derivation. What if we start with Fig. 3 of hahn13?
To adapt it to the z != 0 case, we have to move Fig. 3 up by amount z, which doesn't change that
figure; then also move the perturbing streamline upwards by amount z' - z; and then tilt the
streamline's velocity vector vertically to account for vz'. Can those transformations be used to make
the derivation conceptually simpler? (Then: yes.)
