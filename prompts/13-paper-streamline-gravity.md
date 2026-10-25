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

Follow-up: will doing gravity more correctly improve the code's angular momentum conservation (~1e-5
or -6)? [Measured: no; the wire forces are not reciprocal.] Yes, log it. What is the fast/easy/most
correct way forward: keep the code as is, i.e. evaluate forces at the crossing point not the near
point, but fix the lambda' issue? Revise the paper accordingly; I'll keep reviewing, then we'll start
coding the kernel. Push.

Follow-up: I want to confirm my understanding of hahn13 and epi_int_lite. epi_int_lite uses Eqn 19 for
the magnitude of A_g and Fig. 3 for A_g,r = A_g cos phi, A_g,theta = -A_g sin phi, sin phi = v_r'/v', at
the perturbing streamline's near point; it doesn't use the approximations (22a,b), correct? Yet the code
seems to treat the middle particle of the perturbing streamline as the near point? [No: Eqns 20 are the
approximations; 22 is the own-streamline term; the code evaluates at the crossing point.] Is it
effectively setting phi = 0? [Only in the distance.] So epi_int_lite slides the near point to the
crossing point without changing the direction of A_g, i.e. a force with Delta = r' - r? How hard would it
be to compute A_g,r and A_g,theta correctly at the near point, inferred from the streamline's shape
rather than velocity? A quick remedy is Delta = d_r cos phi; the next fix is phi from the shape. I'd like
to get the 2D version of the code working the way we want, then deal with the z != 0 case.

Follow-up: for the paper, document these 3 fixes in an appendix (the lambda' fix, the shape inferred
from the streamline's shape rather than from v', and the Delta fix), and include hahn13 Fig. 3 in that
appendix so the discussion can refer to it; then the earlier parts of the paper generalize that appendix
to z != 0. We won't backtest this solution against epi_int_lite; we'll instead backtest by comparing
simulations of linear density and spiral bending waves to theory. [Agreed: redraw the figure (it is
(c) AAS); keep the bit-for-bit regression as a software check.] Yes, make the fixes the default, then
draft the appendix, push.

Follow-up: when a new user uses epi_int_3d, do they need to set a flag to run as described in the main
part of the paper, or is a flag set to run the epi_int_lite equivalent? [No flag; but run_simulations.py
gave the planar tests epi_int_lite's gravity.] Yes, make that change (epi_int_3d's own tests run with
its default gravity) and push.

Follow-up: Section 3 still seems too complex. Can we not set vec(D) = vec(r') - vec(r), let hat(t) be the
unit vector tangent to streamline a' at the near point so that vec(D) . hat(t) = 0, and use conventional
vector geometry to solve for hat(t) and then vec(D)? Don't assume the streamline is a straight wire; the
Lagrange polynomials can be used to compute hat(t), so it can be regarded as known, correct? From BAC-CAB
and vec(D) . hat(t) = 0 I get hat(t) x (r x r') = -(hat(t) . r) vec(D), please confirm. [Confirmed, but
r' there must be the unknown near point, and t . r ~ r s_r is nearly 0; with the origin at the crossing
point it becomes D = t x (d x t).] Rewrite Section 3.1 that way; note the accuracy in the paper.

Follow-up: [quoting "The perturbing streamline's position, shape and velocity at the particle's longitude
are found by fitting Lagrange polynomials ..."] means that the perturbing streamline's r' and z' can be
regarded as functions of delta theta', the longitude relative to the perturbed particle's, correct?
Should the paper introduce that relative longitude, perhaps via the symbol delta theta'? Yes, make those
edits.

Follow-up: merge the relevant parts of Appendix B with Section 3 to derive the 3D pull, then let a
revised Appendix B speak to how that solution reduces to epi_int_lite's z = 0 solution after accounting
for the discrepancies. Yes, do the restructure.

Follow-up: shouldn't D be Delta, to be consistent with past papers? [Yes; changed.]

Follow-up: should Eqn. (13) refer to vector Delta, then explain how the Lagrange polynomial is used to
solve for vec(Delta), which points to the near point? Yes, make that change.

Follow-up: the sentence "The perturbing streamline's position and linear density at longitudes
theta' = theta + delta theta' ..." should be revised to indicate that delta theta' gives the near point's
longitude relative to theta.

Follow-up: I'm not sure what "The pull of a displaced streamline" means; maybe change to "The perturbing
streamline's geometry"? Does "|r' - r| << r and |z' - z| << r" need to say |delta theta'| << 1? And what
is "it" in "Near there it lies at d + (s_r, 1, s_z) r delta theta', so its unit tangent is"?

Follow-up: can "as shown below, " be dropped? Change "The streamline crosses the particle's longitude,
delta theta' = 0, at the separation d = (d_r, 0, d_z)" to "The streamline crosses the particle's
longitude at delta theta' = 0, and that crossing point lies at d = (d_r, 0, d_z)", and "A point on the
streamline at a small offset delta theta'" to "other points on that locally straight streamline can be
parameterized by the small longitudinal offset delta theta'".

Follow-up: should "d = (d_r, 0, d_z) from the particle, where d_r = r'(0) - r and d_z = z'(0) - z"
indicate that d points from the particle to the open circle in Fig. 1? Would the discussion be simpler
if we set d = d0 + d1, where d0 = (d_r, 0, d_z) and d1 = [r'(dtheta') - r'(0), r dtheta',
z'(dtheta') - z'(0)]? Yes, rewrite it that way.

Follow-up: "Let d(delta theta') be the vector from the particle to the point of the streamline at the
longitude offset delta theta'" should be "the vector from the particle to the nearest point on streamline
a', which is the black symbol in Fig. 1", and make that symbol more prominent; then d = d0 + d1, where
d0 points from the particle to where the streamline crosses the particle's longitude. [That vector is
Delta itself, so it is written Delta = d0 + d1(delta theta').]

Follow-up: should Eqn. 14 be Delta = Delta_0 + Delta_1(dtheta')? [Recommended keeping d0, d1.] Should
d0 and d1 be shown as vectors in Fig. 1? [Yes.]

Follow-up: in Fig. 1 change d0 = d_r hat(r) to just vec(d0), move "crossing point" above the r axis,
move d1 to the other side of its vector, and change purple to black.

Follow-up: drop the d1 arrow and put d1's arrowhead on the line in question; use thick lines to indicate
vectors, so the line downward of the near point should be less thick. Does Fig. 1 need the v' vector and
the dr', r dtheta triangle? If not, drop them or move them to the appendix. [Moved to a new Appendix B
figure.]

Follow-up: indicate the Delta vector, drop "Delta = d_r cos phi", drop the red vector, make the vectors
less heavy but heavier than the non-vector lines. Should "nearest point" be "near point"? Move "crossing
point" a smidge right.

Follow-up: make the vector labels (Delta etc.) less heavy. Should the theta and r axis labels have hats?
[Yes: they are the local frame's directions.]

Follow-up: move Delta closer to its line and phi midway along its arc; move "nearest point" up a smidge.

Follow-up: move phi slightly closer to its arc, and push.

Follow-up: should there be a 3rd particle shown in a'? And should the crossing point be a non-circle so it
doesn't get confused with a particle? [Yes to both: three particles, crossing point marked x.]

Follow-up: change "points from there to the nearest point" to "points from the crossing to the nearest
point". Should the approx-equal in d1 = [...] ~ (s_r, 1, s_z) r dtheta' be an equality, with approx-equal
used in Eqns. 16 instead? [Recommended keeping it as is.]
