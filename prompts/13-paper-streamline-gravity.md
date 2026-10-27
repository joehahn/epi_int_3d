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

Follow-up: should we instead say Delta = (r'(theta'), theta', z'(theta')), d0 = (r'(0), 0, z'(0) - z), so
Delta = d0 + [(r'(theta'), theta', z'(theta')) - d0] = d0 + d1, where d1 = [r'(dtheta') - r'(0),
r dtheta', z'(dtheta') - z'(0)] ~ (s_r, 1, s_z) r dtheta'? [Yes, with the components measured from the
particle: Delta = [r'(dtheta') - r, r dtheta', z'(dtheta') - z].]

Follow-up: should that discussion instead start with Delta = r1 - r = d0 + d1, where d0 = (d_r, 0, d_z)?
[Yes, as Delta = r' - r ~ [r' - r, r dtheta', z' - z] = d0 + d1 in the local frame.]

Follow-up: should the ends of the r and r1 vectors be shown in Fig. 1? [Labelled: particle at r, nearest
point at r'.] r'

Follow-up: should "The streamline's nearest point lies at a longitude theta' = theta + dtheta' ... as
smooth functions of the longitude offset" be in Section 3.1? [Yes; moved.]

Follow-up: see the desktop photo of a sketch of how I think the vectors should be represented in the
paper (r and r'(dtheta') from off the figure to the particle and the nearest point; r'(dtheta') =
(r'(dtheta'), dtheta', z'(dtheta')); d0 = r'(dtheta' = 0) - r), what do you think? [Agreed, with
r'(dtheta') at cylindrical coordinates (r', theta + dtheta', z').] Yes, make both changes.

Follow-up: delete "that are treated as known below." Should "r'(dtheta') the point of streamline a' at
the cylindrical coordinates (r'(dtheta'), theta + dtheta', z'(dtheta'))" be "r'(dtheta') =
(r'(dtheta'), theta + dtheta', z'(dtheta')) represent the nearpoint in cylindrical coordinates, which is
parameterized by dtheta'"? Note that I used nearpoint rather than nearest point; if that works then
change the paper and figure accordingly. [Used "near point", defined once as the point nearest the
particle.]

Follow-up: should angle dtheta' be indicated in Fig. 1? [As the offset r dtheta' along theta-hat between
the r and r' arrows; the angle itself is subtended at the planet, off the figure.]

Follow-up: should vec(r') = (r', theta', z') refer to points on a', and should that be introduced near
"The perturbing streamline's near point lies at a longitude theta' = theta + dtheta' ..."? [Yes.]

Follow-up: 'The points of a' near the particle are then r'(dtheta') = (r', theta', z')' should be
introduced early in that paragraph.

Follow-up: should "r'(dtheta'), z'(dtheta') and lambda'(dtheta')" be "... and linear density
lambda'(dtheta')", to avoid confusion with a longitude? [Yes.]

Follow-up: change "dtheta' is now the near point's longitude offset" to "dtheta' is the near point's
longitude relative to the particle's longitude", and to "d0 = (d_r, 0, d_z) where d_r = r'(0) - r and
d_z = z'(0) - z".

Follow-up: can "The tangent is thus taken from the streamline's shape rather than its velocity, which is
tangent to the streamline only when the streamline's shape is stationary." be deleted? [Yes.]

Follow-up: delete "(Section 2.3.2 of J. M. Hahn & J. N. Spitale 2013; Eqn. 1 of J. M. Hahn et al. 2025)"
and say that those codes treat the gravity that a perturbing streamline exerts on a particle as a
wire, resulting in the acceleration of Eqn 13. [Done.]

Follow-up: delete "in another streamline" and change "Δ is the perpendicular from the particle to the
streamline's near point, the point of the streamline nearest the particle" to something like "Δ is the
distance from the particle to the streamline's near point, Fig 1". [Done, as "Δ is the separation of the
particle from the streamline's near point, its point nearest the particle (Fig. 1)": bold Δ is a vector,
and the near point is defined only here.]

Follow-up: should hahn13 and hahn25 be reviewed so this paper can be written in my style? [Yes:
extract a one-page style guide, apply it to new text, and propose any style pass on existing sections
separately.] Then: yes, write the style guide. [notes/writing_style.md.] Should the style guide be
pushed or kept local? [Committed to the repo, pushed when the author says "push".]

Follow-up: Equation (13) when that starts the sentence, Eqn otherwise; guide looks good. [Six
sentence-initial references changed to "Equation(s) (n)"; guide, CLAUDE.md and a memory updated.]

Follow-up: leave ("Eqn. 15") for now. Is the paper consistent with my style? [Partly: mechanics match,
prose is terser and more symbol-dense than the author's. Proposed: capitalize sentence-initial code
names, "Fig." to "Figure", "neighbouring" to "neighboring", the past-tense "We tested"; then a style
pass on Sec. 3.1 (physical reasons first, fewer colon chains, "acceleration" for "pull").] Yes. [Done;
Sec. 3.1 pass for review.]

Follow-up: change to "Epi_int and epi_int_lite assume that the gravitational acceleration of a particle
by perturbing streamline a' is that due to a straight and infinitely long wire:" and to "Δ is the
separation of the particle to the streamline's nearest point, designated the near point in Fig. 1."
Also start sentences with Figure but abbreviate Fig elsewhere. [Done; no figure reference starts a
sentence, so all 13 are now "Fig."; style guide and memory updated.]

Follow-up: delete "Its magnitude is Ag = 2Gλ'/Δ (Eqn. 19 of Hahn & Spitale 2013)." and change "This
section finds Δ in three dimensions, and Appendix B relates its planar limit to epi_int and
epi_int_lite." to "The following derives Δ in three dimensions, and Appendix B relates that result to
the co-planar limit examined by epi_int and epi_int_lite." Change "A particle in a narrow ring is always
close to the nearby portions of the other streamlines" to say that all epi_int_3d simulations are of
narrow annuli in a broader ring, so all particles are radially close to the near parts of all
perturbing streamlines. [Done.]

Follow-up: change to "and Appendix B relates that result to the co-planar problem considered by epi_int
and epi_int_lite" and change to "so every particle is radially close to the nearby portions of all
perturbing streamlines. Which allows us to treat the acceleration of that due to a perturbing streamline
as that due to a straight wire". [Done, reading "the particle's acceleration due to a perturbing
streamline".]

Follow-up: change to "All epi_int_3d simulations span only a narrow annulus within a broader planetary
ring". Can "The task is then to find the point of that streamline that is nearest the particle, ... and
the following generalizes it to a streamline that is also displaced and tilted vertically." be deleted?
[Yes: the Sec. 3 opening defines the near point and says Δ is derived in 3D, and the planar-limit
paragraph cites HS13 Fig. 3. Deleted.]

Follow-up: should "All epi_int_3d simulations span only a narrow annulus ... as that due to a straight
wire." then say that we'll build upon the treatment of J. M. Hahn & J. N. Spitale (2013), which should
be referred to as Hahn & Spitale (2013)? Doesn't aastex manage the reference text for you? [Yes to the
sentence: "The following builds upon the treatment of Hahn & Spitale (2013)." AASTeX does manage it:
aasjournalv7.bst adds first initials to in-text citations on purpose (AAS house style since 2024, as in
hahn25). paper/aasjournalv7_jmh.bst is a copy without them.]

Follow-up: use ApJ reference style. [aasjournalv7_jmh.bst removed; back to aasjournalv7.bst, with
first initials in in-text citations.]

Follow-up: change to "Which allows us to treat the particle's acceleration due to a perturbing
streamline as that due to a straight wire, and the following 3D builds treatment is an extension of the
2D force calculations used in J. M. Hahn & J. N. Spitale (2013)." [Done, reading "the following 3D
treatment is an extension of".]

Follow-up: change to "relates that to the co-planar problem" and revise "In cylindrical coordinates the
particle is at r = (r, θ, z), and the points of streamline a' near it are at r'(δθ') = (r', θ', z'),
where θ' = θ + δθ' is offset from the particle's longitude by δθ'" to say that the perturbing
streamline's cylindrical coordinates r' = (r', θ', z') are to be regarded as functions of θ' = θ + δθ'
since streamline a' is sampled by its particles, which allows us to use Lagrange interpolation to
estimate that near point. [Done.]

Follow-up: and "while the perturbing streamline's coordinates r' = (r', θ', z') are to be regarded as
functions of its longitude θ' = θ + δθ'" relative to the particle's longitude θ. [Done.]

Follow-up: change to "Streamline a' is sampled by its particles, so, as in J. M. Hahn & J. N. Spitale
(2013), Lagrange polynomials are fit to the three particles in a' that are nearest in longitude", and
"r'(δθ'), z'(δθ') and the linear density λ'(δθ') as smooth functions of δθ', which allows the
streamline's near point to be estimated below" to "which allows r'(δθ'), z'(δθ') and the linear density
λ'(δθ') to be regarded as smooth functions of δθ'". [Done.]

Follow-up: delete "J. M. Hahn & J. N. Spitale (2013)," from that sentence, change "This splits the trip
to the near point into two legs" to "across two legs", and to "travels from the particle to". Should we
insert "r'(0) points at the crossing point" into the discussion? [Yes: "The first leg, d0 = r'(0) - r,
travels from the particle to the crossing point r'(0), where the streamline crosses the particle's
longitude (the cross in Fig. 1)".]

Follow-up: change to "the crossing point r'(0) which is where" and to "then runs along the streamline
from the crossing point". Should the paragraph state that r'(0) = (r'(0), 0, z'(0))? [Yes, but as
r'(0) = (r'(0), θ, z'(0)) in cylindrical coordinates: its longitude is the particle's θ, not 0, which is
why d0 has no θ-component.]

Follow-up: ∂δθ' looks awkward. [Eqn. 16 now s_r = (1/r) dr'/dθ', s_z = (1/r) dz'/dθ' at θ' = θ: r' and z'
are functions of one variable, and d/dθ' = d/dδθ' since θ is fixed.]

Follow-up: should r'(δθ') be written as r'(θ') where θ' = θ + δθ'? [Yes: it matches Eqn. 16's d/dθ', and
the crossing point reads r'(θ) rather than the misleading r'(0). The near point is r'(θ + δθ'); δθ'
stays as the offset that is solved for.]

Follow-up: change to "to the streamline's nearest point, which is labeled the near point in Fig. 1." Can
"The following derives Δ in three dimensions, and Appendix B relates that to the co-planar problem
considered by epi_int and epi_int_lite." be deleted? [Yes, but it was the main text's only pointer to
Appendix B, so that pointer moved to just after Eqns. 19, the planar limit: "Appendix B relates these
to the in-plane accelerations computed by epi_int and epi_int_lite."]

Follow-up: in "where δθ' is the near point's longitude relative to the particle's longitude" is δθ'
measured in the z=0 plane or in the r, r' plane? [The z = 0 plane: a difference of cylindrical
longitudes about the spin axis.] Maybe that sentence should say "as measured in the z=0 plane". [Done.]

Follow-up: should "(the cross in Fig. 1)" say "(the cross in Fig. 1 where theta'=theta)"? [No: the
sentence already gives r'(θ) = (r'(θ), θ, z'(θ)); offered instead to label the figure.] Yes. [Fig. 1
labels now read "crossing point θ' = θ" and "near point θ' = θ + δθ'".]

Follow-up: can "Equations (13) and (14–19) are exact for a straight wire. But the streamline is curved,
... and smaller still for nearer streamlines." be deleted? [Yes: a side note on the wire approximation's
~1e-4 accuracy that the derivation doesn't need. Deleted; G15 stays in the sympy checks, its
equations.md row now marked notes only.]

Follow-up: epi_int_3d provides particles' r, theta coordinates, plus Lagrange polynomial fits to
particles along streamlines. Add a paragraph to the end of Section 3.1 explaining how those are used to
compute the factors needed to compute the r, theta, z components of Ag. [Added, with Eqns. 20 (Lagrange
value and slope weights) and 21 (the three components); sympy check G17.]

Follow-up: can that paragraph do what hahn13 does and skip the discussion of how Lagrange weights are
derived, and merely state that Lagrange polynomials are used to interpolate as needed? [Yes: the weights
equation is gone; the paragraph says the polynomials fit to the three bracketing particles give r'(θ),
z'(θ), λ'(θ) by interpolation and s_r, s_z from their derivatives. The weights stay in sympy G17.]

Follow-up: should λ' be evaluated at the crossing point? [Yes, it is fine: the near point is |c| ≪ Δ
along the streamline from it, inside the ~Δ stretch that dominates a wire's pull; the fractional
difference is ~(Δ/r)(e m)², below the wire approximation's O(Δ/R) error. Offered to sharpen the paper's
sentence:] no. [Left as is.]

revise "Epi_int and epi_int_lite assume that the gravitational acceleration of a particle by perturbing streamline a′ is that due to a straight and infinitely long wire:" to say that epi_int_3d as well as its predecessor codes. and what if Eqn 13 used hat(Delta)/Delta instead?

[Now "Epi_int_3d, like its predecessors epi_int and epi_int_lite, assumes...". Eqn. 13 is now A_g = (2Gλ′/Δ) Δ̂, which is the same thing (Δ̂/Δ = **Δ**/Δ²) but shows the magnitude 2Gλ′/Δ directly. The "where" clause now defines Δ as the distance to the near point and Δ̂ = **Δ**/Δ as the unit vector toward it. Eqns. 17–20 are unchanged. Also fixed a mid-sentence "Equation (13)" to "Eqn. (13)".]

should delta be unbold (scalar) and fig 1 and should delta-hat be added ?

[No. **Δ** stays bold in Eqn. 14 and Fig. 1: the figure draws the vector **Δ** = **d**₀ + **d**₁, and the plain Δ is its length. I advised against adding a Δ̂ arrow, since it would lie on top of the **Δ** arrow, and suggested changing the caption instead.]

yes, update the caption

[Fig. 1's caption now reads "The wire's acceleration **A**_g = (2Gλ′/Δ)**Δ̂** points along it (Eqns. 19)", which matches Eqn. 13.]

should "All epi_int_3d simulations span only a narrow annulus within a broader planetary ring, so every particle is radially close to the nearby portions of all perturbing streamlines. Which allows us to treat the particle's acceleration due to a perturbing streamline as that due to a straight wire," be above in section 3?

[Yes, since it is the physical reason for the wire approximation. Section 3 now opens with it: "All epi_int_3d simulations span ... perturbing streamlines. Which is why epi_int_3d, like its predecessors epi_int and epi_int_lite, assumes that ... a straight and infinitely long wire:" Section 3.1 now starts "The following 3D treatment is an extension of the 2D force calculations used in HS13." and goes straight on to the coordinates.]

make "while the perturbing streamline's coordinates r′ = (r′, θ′, z′) are to be regarded as functions of its longitude θ′ = θ + δθ′ relative to the particle's longitude θ." its own sentence and refer to r' = (r'(theta'), theta', z'(theta'))

[Now "In cylindrical coordinates the particle is at **r** = (r, θ, z). The perturbing streamline's coordinates **r**′ = (r′(θ′), θ′, z′(θ′)) are to be regarded as functions of its longitude θ′ = θ + δθ′ relative to the particle's longitude θ."]

how about changing "The perturbing streamline's coordinates r′ = (r′(θ′), θ′, z′(θ′)) are to be regarded as functions of its longitude θ′ = θ + δθ′ relative to the particle's longitude θ." to something like "while The perturbing streamline's coordinates r′, θ′, z′ are to be regarded as functions of its longitude θ′ ie r' = (r′(θ′), θ′, z′(θ′))" the later sentenc says theta'= θ + δθ′ relative to the particle's longitude θ.

[Now "The perturbing streamline's coordinates are to be regarded as functions of its longitude θ′, i.e., **r**′ = (r′(θ′), θ′, z′(θ′))." With δθ′ no longer introduced there, the next sentence's |δθ′| ≪ 1 became |θ′ − θ| ≪ 1, since Eqn. 14 is where δθ′ is defined.]

change to "But in the following the perturbing streamline's coordinates are to be regarded as functions of its longitude θ′, i.e., r′ = (r′(θ′), θ′, z′(θ′))." and introduce "Streamline a′ " in the paragraph above

[Done as worded. Section 3's opening paragraph now starts by introducing streamlines and their labels, following hahn25 Sec. 2.1: "Like its predecessors epi_int and epi_int_lite, epi_int_3d assigns its particles to streamlines, each a closed eccentric path around the planet whose particles initially share a common semimajor axis, and a streamline is labeled by that semimajor axis." The wire sentence now reads "Which is why all three codes assume that the gravitational acceleration of a particle by perturbing streamline a′, i.e., the streamline whose semimajor axis is a′, is that due to ...".]

change "the gravitational acceleration of a particle by perturbing streamline a′, i.e., the streamline whose semimajor axis is a′, " to something like "the gravitational acceleration of a particle by those particles inhhaniting the perturbing streamline at semimajor axis a′

[Now "...assume that the gravitational acceleration of a particle by the particles inhabiting the perturbing streamline at semimajor axis a′ is that due to a straight and infinitely long wire:"]

move "Like its predecessors epi_int and epi_int_lite," to end of its sentence. and change "epi_int_3d assigns its particles to streamlines, each a closed eccentric path around the planet whose particles initially share a common semimajor axis" to epi_int_3d assigns its particles to streamlines, each a closed eccentric path around the planet whose particles share a common semimajor axis a

[Now "Epi_int_3d assigns its particles to streamlines, each a closed eccentric path around the planet whose particles share a common semimajor axis a, and a streamline is labeled by that semimajor axis, like its predecessors epi_int and epi_int_lite."]

some papers use typwriter fonts when referring to codes. what is standard practice for ApJ papers?

[No firm ApJ rule. I said roman was safer and matched hahn25 and hahn13.]

i thought hahn13 refered to typewriter(epi_int)

[Right. The fonts embedded in the PDFs show hahn13 sets every code name in typewriter (CMTT10: epi_int ×16, SYMBA, MERCURY), while hahn25 uses roman, so ApJ has accepted both.]

do the macro

[Added \newcommand{\code}[1]{\texttt{#1}} to the preamble, and all 35 code names in the text, titles and captions are now \code{...}. Changing the macro to {#1} brings back roman. notes/writing_style.md records the rule.]

change to "Epi_int_3d assigns particles to streamlines that are a closed eccentric path around the planet whose particles share a common semimajor axis a"

[Done as worded.]

can we delete "and a streamline is labeled by that semimajor axis, like its predecessors epi_int and epi_int_lite"? and "Which is why all three codes assume " to "Which is why epi_int_3d and its predecessors

[Both done. The predecessors are named there instead, "Which is why epi_int_3d and its predecessors epi_int and epi_int_lite assume ...", because they are no longer named earlier in Section 3.]

change "predecessors epi_int and epi_int_lite assume " to "predecessors (epi_int and epi_int_lite) assume "

[Done.]

change to something like "the ring matter that is represented by the particles in the perturbing streamline a' will have coordinates that are coordinates r' theta' z' that are to be regarded as functions of longitude theta' ie r′ = (r′(θ′), θ′, z′(θ′)).

[Now "But in the following, the ring matter that is represented by the particles in perturbing streamline a′ has coordinates r′, θ′ and z′ that are to be regarded as functions of longitude θ′, i.e., **r**′ = (r′(θ′), θ′, z′(θ′))."]

that first paragraph in section 3.1 seems slightly repetititve

[Tightened. It said "the following" twice, "regarded as functions of θ′" twice, and introduced streamline a′'s particles twice. It now reads "In cylindrical coordinates the particle is at **r** = (r, θ, z), and the ring matter that is represented by the particles in perturbing streamline a′ is at **r**′ = (r′(θ′), θ′, z′(θ′)), with its coordinates regarded as functions of longitude θ′. These functions, and the streamline's linear density λ′(θ′), are Lagrange polynomials fit to the three particles in a′ that are nearest in longitude. And because ..."]

split

[Now "...the particle is at **r** = (r, θ, z). The ring matter that is represented by the particles in perturbing streamline a′ is at **r**′ = (r′(θ′), θ′, z′(θ′)), with its coordinates regarded as functions of longitude θ′."]

i preferred introducing r' = (r', θ', z') then saying they are functions of longitude θ′ then writing r′ = (r′(θ′), θ′, z′(θ′)),

[Now "The ring matter that is represented by the particles in perturbing streamline a′ is at **r**′ = (r′, θ′, z′), and these coordinates are to be regarded as functions of longitude θ′, i.e., **r**′ = (r′(θ′), θ′, z′(θ′))."]

maybe "streamline a′ is at r′ = (r′, θ′, z′) traces a continuouse wire of matter so these coordinates are to be regarded as functions of longitude θ′, i.e., r′ = (r′(θ′), θ′, z′(θ′))." and then later say that Lagrange polynomials fit to the particles allow the code (or use) to estimate where the wire is

[Now "The ring matter that is represented by the particles in perturbing streamline a′ is at **r**′ = (r′, θ′, z′) and traces a continuous wire of matter, so these coordinates are to be regarded as functions of longitude θ′, i.e., **r**′ = (r′(θ′), θ′, z′(θ′)). Lagrange polynomials fit to the three particles in a′ that are nearest in longitude then provide estimates of where that wire is, r′(θ′) and z′(θ′), as well as its linear density λ′(θ′)."]

r′ = (r′, θ′, z′) that traces a continuous wire . and to something like Lagrange polynomials fit to the three particles in a′ that are nearest in longitude to the parturbed particle at vec(r) are used to estimates of where that wire is,

[Now "...is at **r**′ = (r′, θ′, z′) that traces a continuous wire, so these coordinates are to be regarded as functions of longitude θ′, i.e., **r**′ = (r′(θ′), θ′, z′(θ′)). Lagrange polynomials fit to the three particles in a′ that are nearest in longitude to the perturbed particle at **r** are used to estimate where that wire is, r′(θ′) and z′(θ′), as well as its linear density λ′(θ′)."]

can "Lagrange polynomials fit to the three particles in a′ that are nearest in longitude to the perturbed particle at r are used to estimate where that wire is, r′(θ′) and z′(θ′), as well as its linear density λ′(θ′)." be changed to someting like "Lagrange polynomials will be fit to the particles in streamline a' which means that r′(θ′) and z′(θ′), as well as its linear density λ′(θ′) can be regard as known functions in the following

[Now "Lagrange polynomials will be fit to the particles in streamline a′, which means that r′(θ′), z′(θ′) and the streamline's linear density λ′(θ′) can be regarded as known functions in the following." The recipe paragraph at the end of Sec. 3.1 still says which three particles are fit.]

change "The vector from the particle to the streamline's near point (the black square in Fig. 1) is then" to "... can be written" and "where δθ′ is the near point's longitude relative to the particle's longitude, as measured in the z = 0 plane" to the same without the comma

[Both done.]

change to "∆= r′(θ′)− r = d0 + d1, " and "where δθ′ is the near point's longitude" to "where theta'=theta + δθ′ is the near point's longitude" and add |δθ′|<<1

[Eqn. 14 is now **Δ** = **r**′(θ′) − **r** = **d**₀ + **d**₁, followed by "where θ′ = θ + δθ′ is the near point's longitude as measured in the z = 0 plane, so δθ′ is its longitude relative to the particle's, and |δθ′| ≪ 1." The wording keeps "relative to the particle's" attached to δθ′, since θ′ itself is not a relative longitude.]

and to "and δθ′ is its longitude relative to the particle's with |δθ′| ≪ 1"

[Now "where θ′ = θ + δθ′ is the near point's longitude as measured in the z = 0 plane, and δθ′ is its longitude relative to the particle's with |δθ′| ≪ 1."]

change "The second leg, d1 = r′(θ + δθ′)− r′(θ)," to "The second leg, d1 = r′(θ′)− r′(θ) - d_0 = r'(theta')-r'(theta)'

[The middle expression as typed, r′(θ′) − r′(θ) − **d**₀, isn't equal to the last one, so I used **d**₁ = **Δ** − **d**₀ = **r**′(θ′) − **r**′(θ), which is presumably what was meant.]

d1 in cylindrical coordinates is vec(d1)=(r′(θ′)− r′(θ), θ′ - θ, z′(θ′)− z′(θ)), correct?

[No. In the particle's local frame the middle component is a length, r′ sin δθ′ ≃ r δθ′, not the angle θ′ − θ. Coordinate differences aren't vector components.]

so d1 in cyclindrical coordinat can be written vec(d1)=(r′(θ′), θ', z′(θ′)) - (r′(θ), θ, z′(θ)) correct?

[Yes, if each triple is read as a point (position vector), but not subtracted entry by entry, since the two points have different r̂ and θ̂.]

fig 1 also says that d0=delta*cos(varphi)-d1*sin(varphi) and delta*sin(varphi)=d1*cos(varphi)=r δθ′. does that simplify the math or discussion any

[Only for the plane. That's the right triangle in Fig. 1, and with signed d₁ = N r δθ′ the second relation is Δ sin φ = −d₁ cos φ = −r δθ′. In 3D there are two slopes and no single φ, so Eqn. 17's projection is still needed. Not added.]

does reader need to be reminded that so d1 ≃ [r′(θ + δθ′)− r′(θ), r δθ′, z′(θ + δθ′)− z′(θ)] is using cartesian coordinates?

[A brief reminder helps, since the paragraph also writes cylindrical coordinates as triples.]

do what's best

[Now "...so in the local frame **d**₁ ≃ [r′(θ′) − r′(θ), r δθ′, z′(θ′) − z′(θ)] ≃ ...", which matches the **d**₀ wording and Eqn. 14's θ′.]

can we instead say d1 ≃ (r′(θ′)− r′(θ)*r-hat +r δθ′*theta-hat + (z′(θ′)− z′(θ))z-hat because the curvilinear coordine system can be treatd as cartesian when all of the displacements are small?

[Yes. Now "**d**₁ ≃ [r′(θ′) − r′(θ)]r̂ + r δθ′ θ̂ + [z′(θ′) − z′(θ)]ẑ ≃ (s_r r̂ + θ̂ + s_z ẑ) r δθ′ = N r δθ′ t̂". The explicit unit vectors make "in the local frame" unnecessary. The reason the frame can be treated as Cartesian (small displacements) is already given at the end of Sec. 3.1's first paragraph, so it isn't repeated.]

the hatted unit vectors should be bold

[Done throughout the paper: r̂, θ̂, ẑ are now bold, as t̂ and Δ̂ already were. Fig. 1 and the Appendix B figure (wire_planar.py) were regenerated with bold axis labels, and notes/writing_style.md records the rule.]

should the r in Eqn 16 be r'?

[Both are fine to the derivation's order, since they differ at O(|d_r|/r). r′ = r′(θ), the crossing point's radius, gives the streamline's true slope, makes tan φ = s_r exact in Fig. 1, and matches the planar kernel (shape_slope divides by the interpolated r′). I recommended switching.]

yes

[Eqn. 16 is now s_r = (1/r′) dr′/dθ′ and s_z = (1/r′) dz′/dθ′ at θ′ = θ. The θ̂ leg of **d**₁ is r′(θ)δθ′, so **d**₁ ≃ N r′(θ)δθ′ t̂, the near point is at δθ′ = −c/r′(θ), and |c| ≪ r′. Fig. 1 and its caption now label r′δθ′. Eqns. 17–20 are unchanged. Sympy G10 (comments) and G16 were updated and all checks pass; the equations.md rows G10 and G16 were updated too.]

i think t hat should be written in terms or rhat etc intead of a tripple

[Done for t̂ in Eqn. 15 and for every other component triple in Sec. 3.1: **d**₀ = d_r **r̂** + d_z **ẑ**, t̂ = (s_r **r̂** + **θ̂** + s_z **ẑ**)/N, Eqn. 18's **Δ** = (d_r − c s_r)**r̂** − c **θ̂** + (d_z − c s_z)**ẑ**, and the planar **Δ** = d_r cos φ (cos φ **r̂** − sin φ **θ̂**). The cylindrical-coordinate triples for points, **r** = (r, θ, z) and **r**′(θ) = (r′(θ), θ, z′(θ)), stay as they are.]

change to "which are obtained from the derivatives of the Lagrange polynomials." instead of stating that t-hat is tangent to streamline should paper prove it?

[Wording changed. No proof needed: the d₁ line already shows it. d₁ is the chord from the crossing point to a nearby point on the streamline, and to first order in δθ′ it is N r′δθ′ t̂, so t̂ is the chord's direction as δθ′ → 0, which is the tangent.]

should text instead say that fig 1 shows that t-hat is the streamline's tangent

[I advised against it: Fig. 1 doesn't draw t̂ and can't show its s_z part. It does show **d**₁ along the streamline.]

what part of derivation in the text is what establishes that t̂ is the tangent?

[The **d**₁ line. But it only stated the result: the Taylor step r′(θ′) − r′(θ) ≃ (dr′/dθ′)δθ′ and the chord-to-tangent limit were both implicit, and the slopes were defined after they were used.]

yes add it

[The **d**₁ line now goes on "...to first order in δθ′, since r′(θ′) − r′(θ) ≃ (dr′/dθ′)δθ′ and likewise for z′, where the streamline's slopes are [Eqn. 15, s_r, s_z] at θ′ = θ, which are obtained from the derivatives of the Lagrange polynomials. As δθ′ → 0 the chord **d**₁ points along the streamline, so [Eqn. 16, t̂] is the streamline's unit tangent." The slopes and t̂ swap equation numbers, so the slopes are now defined before t̂ uses them.]

t̂ is the perturbing streamline's unit tangent:

[Now "As δθ′ → 0 the chord **d**₁ points along the streamline, so t̂ is the perturbing streamline's unit tangent: [Eqn. 16]".]

change to where the streamline's radial and vertical slopes are. should "s_r = (1/r′) dr′/dθ′" be written with vertical bar and theta'=theta subscript? and change "at θ′ = θ, which are obtained from the derivatives of the Lagrange polynomials." to "which are evaluated at the perturbed particle's longitide with θ′ = θ, which are obtained from the derivatives of the Lagrange polynomials.

[Yes to the bar: it shows that r′ as well as the derivative is evaluated at θ′ = θ. Eqn. 15 is now s_r = (1/r′)(dr′/dθ′)|_{θ′=θ} and likewise for s_z, introduced by "where the streamline's radial and vertical slopes are". It's followed by "which are evaluated at the perturbed particle's longitude θ′ = θ and are obtained from the derivatives of the Lagrange polynomials." I used "and are" to avoid two "which" clauses in a row.]

change "As δθ′ → 0 the chord d1 points along the streamline, so t̂ is the perturbing streamline's unit tangent:" to "note that as δθ′ → 0 the chord d1 points along the streamline, so t̂ is also the perturbing streamline's unit tangent:" and should the comma in eqn 16 be replaced with where?

[Sentence changed, starting "Note that ...". Yes to "where": Eqn. 16 is now just t̂ = (s_r **r̂** + **θ̂** + s_z **ẑ**)/N, followed by "where N² = 1 + s_r² + s_z²." in the text, as the style guide's "where" clauses do.]

change "Submitted to ApJ" to "to be Submitted to ApJ"

[AASTeX's \submitjournal{ApJ} always prints "Submitted to ApJ", so the preamble now sets the line directly, to "To be submitted to ApJ".]

d0 = dr r̂+dz ẑ,

[In the PDF the line was justified so tightly that TeX shrank the spaces around + and = in "**d**₀ = d_r **r̂** + d_z **ẑ**, where d_r = r′(θ) − r" to almost nothing. The preamble now limits how far inline-math spacing around +, − and = can shrink (\medmuskip, \thickmuskip), which fixes this line and others like it. No new overfull boxes.]

should d0 = dr r̂+dz ẑ, and d1=N r′(θ) δθ′ t̂ be displayed as numbered equations?

[I recommended both: **d**₁ is the derivation's main step and its inline chain broke across three lines, and **d**₀ defines d_r and d_z, which every later result uses.]

yes

[**d**₀ = d_r **r̂** + d_z **ẑ** is now Eqn. 15, followed by "where d_r = ... and d_z = ...". **d**₁ is the two-line Eqn. 16, ending in N r′(θ)δθ′ **t̂**, introduced by "so to first order in δθ′". The later equations are renumbered: slopes 17, t̂ 18, **Δ** 19–20, planar 21, components 22. The paper numbers quoted in notes/equations.md for G11, G16 and G17 were updated.]
