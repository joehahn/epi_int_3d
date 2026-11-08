# derive: openai/gpt-6-astra

Problem: `reviews/problems/vertical_motion.md` · 2026-11-08 · tokens in 3884, out 17322 · cost $0.9146

---

APPENDIX A
VERTICAL EPICYCLIC MOTION IN THE PLANET’S OBLATE POTENTIAL

A.1. Assumptions and ordering

During the drift, the particle moves solely under the planet’s gravity. We neglect ring gravity, collisions, satellites, and zonal harmonics beyond $J_2$. We use cylindrical coordinates $(r,\theta,z)$, so $r$ remains the distance from the rotation axis, not the spherical distance from the planet. The orbit is prograde, nearly circular, and nearly equatorial.

The supplied planar solution is retained without modification. Because the potential is symmetric about the equatorial plane, the vertical displacement changes the horizontal acceleration first at order $i^2$. Ignoring that feedback is therefore consistent with a calculation linear in inclination.

There is, however, an eccentricity-dependent term in the vertical acceleration at order $ei$. We retain it: this makes the vertical solution first order in eccentricity as well as inclination, and recovers the tilt of an eccentric Kepler orbit when $J_2=0$. The resulting vertical acceleration errors are $O(e^2i,i^3)$, while the horizontal acceleration errors remain $O(e^2,i^2)$. This is not a complete second-order expansion in the two small quantities.

We first derive the vertical force, then solve its periodically modulated oscillator equation. We subsequently construct an exactly invertible coordinate map and give the drift algorithm. Here $i$ and $\Omega$ will be epicyclic, or mean, inclination and nodal elements; their relation to an actual ascending crossing is stated below.

A.2. The vertical force

Let $\mu=GM_p$ and $s=(r^2+z^2)^{1/2}$. The planet’s monopole-plus-quadrupole potential is

$$
\Phi(r,z)
=
-\frac{\mu}{s}
\left[
1-J_2\left(\frac{R_p}{s}\right)^2
P_2\left(\frac{z}{s}\right)
\right],
\qquad
P_2(u)=\frac{3u^2-1}{2}.
\tag{A1}
$$

Because $|z|/r\ll1$, we can expand the potential about the equatorial plane. The identity used for each inverse power of $s$ is

$$
s^{-m}
=
r^{-m}
\left[
1-\frac{m}{2}\frac{z^2}{r^2}
+O\left(\frac{z^4}{r^4}\right)
\right].
\tag{A2}
$$

Substituting $P_2$ into Equation (A1) and then using Equation (A2) gives

$$
\begin{aligned}
\Phi
&=
-\frac{\mu}{s}
-\frac{\mu J_2R_p^2}{2s^3}
+\frac{3\mu J_2R_p^2z^2}{2s^5}
\\
&=
-\frac{\mu}{r}
-\frac{\mu J_2R_p^2}{2r^3}
+\frac{\mu z^2}{2r^3}
+\frac{9\mu J_2R_p^2z^2}{4r^5}
\\
&\quad
+O\left(
\frac{\mu z^4}{r^5}
\left[1+|J_2|\left(\frac{R_p}{r}\right)^2\right]
\right).
\end{aligned}
\tag{A3}
$$

The vertical equation of motion is $\ddot z=-\partial\Phi/\partial z$, so

$$
\ddot z+\nu^2(r)z
=
O\left(
\frac{\mu z^3}{r^5}
\left[1+|J_2|\left(\frac{R_p}{r}\right)^2\right]
\right),
\qquad
\nu^2(r)
=
\frac{\mu}{r^3}
\left[
1+\frac{9}{2}J_2\left(\frac{R_p}{r}\right)^2
\right].
\tag{A4}
$$

Thus the linear vertical restoring force is stronger than the monopole restoring force when $J_2>0$.

For completeness, define the equatorial potential $\Phi_0(r)=\Phi(r,0)$. Its derivatives are

$$
\Phi_0'(r)
=
\frac{\mu}{r^2}
+\frac{3\mu J_2R_p^2}{2r^4},
\qquad
\Phi_0''(r)
=
-\frac{2\mu}{r^3}
-\frac{6\mu J_2R_p^2}{r^5}.
\tag{A5}
$$

Axisymmetry conserves the axial specific angular momentum $h=r^2\dot\theta$. A circular equatorial orbit obeys $h^2/a^3=\Phi_0'(a)$. Expanding the radial equation about that orbit gives $\ddot{\delta r}=-[\Phi_0''(a)+3\Phi_0'(a)/a]\delta r$ to first order in $\delta r/a$. Consequently,

$$
\Omega_a^2=\frac{\Phi_0'(a)}{a},
\qquad
\kappa_a^2=\Phi_0''(a)+\frac{3\Phi_0'(a)}{a}.
\tag{A6}
$$

Introduce the dimensionless oblateness strength $q_a$ and the monopole frequency $n_a$,

$$
q_a=J_2\left(\frac{R_p}{a}\right)^2,
\qquad
n_a=\sqrt{\frac{\mu}{a^3}}.
\tag{A7}
$$

Equations (A4)–(A6) then yield

$$
\boxed{
\begin{aligned}
\Omega_a^2&=n_a^2\left(1+\frac{3}{2}q_a\right),\\
\kappa_a^2&=n_a^2\left(1-\frac{3}{2}q_a\right),\\
\nu_a^2&=n_a^2\left(1+\frac{9}{2}q_a\right).
\end{aligned}}
\tag{A8}
$$

Here $\nu_a$ is the vertical oscillation frequency of a circular, infinitesimally inclined orbit. No expansion in $q_a$ has been made. These frequencies also satisfy the exact identity $\kappa_a^2+\nu_a^2=2\Omega_a^2$, as required by Laplace’s equation outside an axisymmetric gravitating body.

We assume $\kappa_a^2>0$, so the guiding circular orbit is radially stable. For an oblate planet this requires $0\leq q_a<2/3$; ordinary planetary rings are well away from this upper boundary.

A.3. Eccentric modulation of the vertical oscillator

The planar epicycle supplied in the question, corresponding to the construction in Hahn et al. (2025, Appendix A), gives

$$
\begin{aligned}
r&=\frac{a}{\sqrt{1+2e\cos M}}
=a\left[1-e\cos M+O(e^2)\right],
\qquad
\dot M=\kappa_a,\\
\lambda&\equiv\frac{\Omega_a}{\kappa_a}M+\tilde\omega,
\qquad
\dot\lambda=\Omega_a,\\
\theta&=\lambda+2\frac{\Omega_a}{\kappa_a}e\sin M.
\end{aligned}
\tag{A9}
$$

The mean longitude $\lambda$ is useful because it advances uniformly, whereas $\theta$ does not.

Because the vertical restoring coefficient depends on $r$, the radial epicycle modulates it. Substituting Equation (A9) into Equation (A4) gives

$$
\begin{aligned}
\nu^2[r(M)]
&=
n_a^2(1+2e\cos M)^{3/2}
+\frac{9}{2}n_a^2q_a(1+2e\cos M)^{5/2}\\
&=
\nu_a^2+eC_a\cos M
+O\!\left(n_a^2(1+q_a)e^2\right),\\
C_a&\equiv n_a^2\left(3+\frac{45}{2}q_a\right).
\end{aligned}
\tag{A10}
$$

Since $z/a=O(i)$, the second term in Equation (A10) produces an acceleration of order $ei$. Retaining that term, while neglecting $e^2i$ and the cubic vertical force, gives

$$
\ddot z+\nu_a^2z
=
-eC_a\cos M\,z
+
O\!\left(an_a^2[e^2i+i^3]\right).
\tag{A11}
$$

Thus, even after linearization in $z$, the eccentric orbit has a periodically varying vertical restoring force. A constant-frequency sinusoid alone would omit an $O(ei)$ term.

Define a uniformly advancing vertical phase $\psi$ and a mean nodal longitude $\Omega$ by

$$
\dot\psi=\nu_a,
\qquad
\Omega\equiv\lambda-\psi,
\qquad
\dot\Omega=\Omega_a-\nu_a.
\tag{A12}
$$

The corrections to this vertical frequency and mean nodal rate begin at $O(e^2,i^2)$ and are not calculated here.

We now seek a solution of Equation (A11). Because its forcing is proportional to $\cos M\sin\psi$, differentiation couples that function to $\sin M\cos\psi$. We therefore write

$$
z
=
ai\left[
\sin\psi
+e\left(
-\alpha_a\cos M\sin\psi
+\beta_a\sin M\cos\psi
\right)
\right].
\tag{A13}
$$

The amplitude $i\geq0$ is the epicyclic inclination, measured in radians. Homogeneous first-order corrections can be absorbed into $i$ and the origin of $\psi$, so Equation (A13) introduces only the two required vertical integration constants.

To determine $\alpha_a$ and $\beta_a$, define $\mathcal L=d^2/dt^2+\nu_a^2$. Direct differentiation, using $\dot M=\kappa_a$ and $\dot\psi=\nu_a$, gives

$$
\begin{aligned}
\mathcal L(\cos M\sin\psi)
&=
-\kappa_a^2\cos M\sin\psi
-2\kappa_a\nu_a\sin M\cos\psi,\\
\mathcal L(\sin M\cos\psi)
&=
-\kappa_a^2\sin M\cos\psi
-2\kappa_a\nu_a\cos M\sin\psi.
\end{aligned}
\tag{A14}
$$

At order $ei$, the right side of Equation (A11) uses only the zeroth-order displacement $ai\sin\psi$. Substitution of Equation (A13), followed by matching the two trigonometric coefficients, therefore gives

$$
\begin{aligned}
\alpha_a\kappa_a^2-2\beta_a\kappa_a\nu_a&=-C_a,\\
2\alpha_a\kappa_a\nu_a-\beta_a\kappa_a^2&=0.
\end{aligned}
\tag{A15}
$$

The second equation gives $\beta_a=2\alpha_a\nu_a/\kappa_a$. Substituting this into the first equation determines

$$
\boxed{
\alpha_a
=
\frac{C_a}{4\nu_a^2-\kappa_a^2}
=
\frac{1+\frac{15}{2}q_a}{1+\frac{13}{2}q_a},
\qquad
\beta_a=\frac{2\alpha_a\nu_a}{\kappa_a}.
}
\tag{A16}
$$

The factor $\alpha_a$ measures the radial modulation of the vertical amplitude; $\beta_a$ measures the associated modulation of its phase.

Differentiating Equation (A13), without making another approximation, yields

$$
\begin{aligned}
v_z
=ai\bigg\{&
\nu_a\cos\psi\\
&+e(\alpha_a\kappa_a-\beta_a\nu_a)
       \sin M\sin\psi\\
&+e(\beta_a\kappa_a-\alpha_a\nu_a)
       \cos M\cos\psi
\bigg\}.
\end{aligned}
\tag{A17}
$$

Using the exact identity $\beta_a\kappa_a=2\alpha_a\nu_a$, define

$$
\begin{aligned}
\gamma_a&\equiv
\beta_a-\alpha_a\frac{\kappa_a}{\nu_a}
=
\alpha_a\left(
\frac{2\nu_a}{\kappa_a}-\frac{\kappa_a}{\nu_a}
\right),\\
F(M)&\equiv1-\alpha_a e\cos M,
&
G(M)&\equiv\beta_a e\sin M,\\
H(M)&\equiv-\gamma_a e\sin M,
&
K(M)&\equiv1+\alpha_a e\cos M.
\end{aligned}
\tag{A18}
$$

These four dimensionless functions describe the eccentric modulation of the vertical position and velocity. The forward coordinate map is then

$$
\boxed{
\begin{aligned}
\psi&=\lambda-\Omega,\\
z&=ai\,[F(M)\sin\psi+G(M)\cos\psi],\\
v_z&=a\nu_ai\,[H(M)\sin\psi+K(M)\cos\psi].
\end{aligned}}
\tag{A19}
$$

Equations (A18) and (A19) will be treated as exact definitions of the approximate coordinate map. Their physical accuracy is assessed separately in Section A.6.

For a circular orbit, an ascending crossing has $\psi=0$ and $\theta=\Omega$. For an eccentric orbit, $\Omega$ is a mean node: the periodic factors in Equation (A19) shift the actual crossing. Setting $z=0$ near $\psi=0$ gives $\psi_{\rm asc}=-\beta_a e\sin M+O(e^2,i^2)$, and hence

$$
\theta_{\rm asc}
=
\Omega+
\frac{2e}{\kappa_a}
(\Omega_a-\alpha_a\nu_a)\sin M
+O(e^2,i^2).
\tag{A20}
$$

Thus the longitude of an actual ascending crossing need not equal the mean nodal element at finite eccentricity. In an oblate potential the orbit is not confined to one fixed inclined plane.

A.4. Exact inversion of the vertical map

Suppose that the planar inversion has already supplied $a$, $e$, $M$, and $\lambda$. Define the dimensionless coordinates and inclination quadratures

$$
Z\equiv\frac{z}{a},
\qquad
W\equiv\frac{v_z}{a\nu_a},
\qquad
X\equiv i\sin\psi,
\qquad
Y\equiv i\cos\psi.
\tag{A21}
$$

Equation (A19) is the linear system

$$
\begin{pmatrix}Z\\W\end{pmatrix}
=
\begin{pmatrix}F&G\\H&K\end{pmatrix}
\begin{pmatrix}X\\Y\end{pmatrix}.
\tag{A22}
$$

Its determinant is

$$
D(M)
\equiv FK-GH
=
1-\alpha_a^2e^2\cos^2M
+\beta_a\gamma_ae^2\sin^2M.
\tag{A23}
$$

For the stable oblate branch, $\beta_a\gamma_a>0$. Consequently $D\geq1-\alpha_a^2e^2$, so the inverse exists for the nearly circular orbits considered here.

The exact inverse of Equation (A22) is

$$
\boxed{
X=\frac{KZ-GW}{D},
\qquad
Y=\frac{FW-HZ}{D}.
}
\tag{A24}
$$

The elements follow as

$$
\boxed{
i=\sqrt{X^2+Y^2},
\qquad
\psi=\operatorname{atan2}(X,Y),
\qquad
\Omega=\lambda-\psi\pmod{2\pi}.
}
\tag{A25}
$$

The argument order in $\operatorname{atan2}(X,Y)$ is intentional: $X$ is the sine quadrature and $Y$ the cosine quadrature.

Although the forward map retains eccentricity only through first order, its determinant contains $e^2$ terms. These must not be discarded when inverting the map. Keeping the full determinant makes the conversion algebraically exact; it does not assert second-order physical accuracy.

When $i=0$, both quadratures vanish and the node is undefined. The implementation should then retain a conventional node or work directly with $X$ and $Y$.

A.5. Drift algorithm and reversibility

The drift can be implemented in the following order.

1. Recover the planar elements.

Calculate the conserved axial angular momentum $h=rv_\theta$. The guiding radius follows from the same equation used in the planar code,

$$
h^2
=
\mu\left(a+\frac{3J_2R_p^2}{2a}\right),
\qquad
a
=
\frac{1}{2}
\left[
\frac{h^2}{\mu}
+
\sqrt{
\left(\frac{h^2}{\mu}\right)^2-6J_2R_p^2
}
\right].
\tag{A26}
$$

The plus sign selects the radially stable branch. Calculate the three frequencies from Equation (A8), and then calculate

$$
\begin{aligned}
E_s&\equiv e\sin M
=\frac{v_r}{a\kappa_a}\left(\frac{a}{r}\right)^3,\\
E_c&\equiv e\cos M
=\frac12\left[\left(\frac{a}{r}\right)^2-1\right],\\
e&=\sqrt{E_s^2+E_c^2},
\qquad
M=\operatorname{atan2}(E_s,E_c),\\
\lambda&=\theta-2\frac{\Omega_a}{\kappa_a}E_s,
\qquad
\tilde\omega=\lambda-\frac{\Omega_a}{\kappa_a}M.
\end{aligned}
\tag{A27}
$$

These are exact inverses of the supplied planar coordinate map. When $e=0$, choose any convenient $M$; $\lambda=\theta$ remains well defined.

2. Recover the vertical elements.

Calculate $\alpha_a$, $\beta_a$, and $\gamma_a$ from Equations (A16) and (A18), then evaluate $F,G,H,K,D$ at the recovered $M$. Use Equations (A24) and (A25) to recover $X,Y,i,\Omega$. This conversion is exact for Equation (A19).

3. Advance the elements.

During a drift of duration $\Delta t$, use

$$
\boxed{
\begin{aligned}
a'&=a,
&
e'&=e,
&
i'&=i,
&
\tilde\omega'&=\tilde\omega,\\
M'&=M+\kappa_a\Delta t,
&
\lambda'&=\lambda+\Omega_a\Delta t,\\
\Omega'&=\Omega+(\Omega_a-\nu_a)\Delta t,
&
\psi'&=\psi+\nu_a\Delta t.
\end{aligned}}
\tag{A28}
$$

The expression for $\psi'$ follows directly from $\psi'=\lambda'-\Omega'$.

For robust angle bookkeeping, either retain an unwrapped $M$ when using $\tilde\omega$, or advance $\lambda$ independently as shown above. Wrapping $M$ while leaving $\tilde\omega$ unchanged can otherwise introduce a spurious longitude jump.

It is also convenient to advance the vertical quadratures without evaluating a node:

$$
\begin{aligned}
X'&=X\cos(\nu_a\Delta t)+Y\sin(\nu_a\Delta t),\\
Y'&=Y\cos(\nu_a\Delta t)-X\sin(\nu_a\Delta t).
\end{aligned}
\tag{A29}
$$

This is Equation (A28) written in Cartesian inclination variables, and remains regular at $i=0$.

4. Reconstruct the coordinates.

Use the original planar expressions at $M'$. Equivalently, reconstruct the azimuth as $\theta'=\lambda'+2(\Omega_a/\kappa_a)e\sin M'$.

For the vertical coordinates, first recompute the modulation functions at the new anomaly $M'$, and then use

$$
\boxed{
\begin{aligned}
z'&=a\,[F(M')X'+G(M')Y'],\\
v_z'&=a\nu_a\,[H(M')X'+K(M')Y'].
\end{aligned}}
\tag{A30}
$$

It is important to remodulate at $M'$. Directly rotating $(z,v_z/\nu_a)$ as a constant-frequency oscillator would omit the retained $O(ei)$ terms.

To see why the drift is reversible, let $\mathcal T$ denote the complete coordinate-to-element transformation and $\mathcal S_{\Delta t}$ the element advance. The implemented coordinate drift is

$$
\mathcal D_{\Delta t}
=
\mathcal T^{-1}\mathcal S_{\Delta t}\mathcal T,
\qquad
\mathcal D_{-\Delta t}\mathcal D_{\Delta t}=I.
\tag{A31}
$$

The identity holds in exact arithmetic because the coordinate conversion is inverted exactly and the element advances change sign under $\Delta t\rightarrow-\Delta t$.

Thus the drift introduces no additional time-discretization truncation within this epicyclic model. It still has the physical approximation errors derived below, and an actual implementation has floating-point error. Exact reversibility should not be confused with exact motion in the full gravitational potential.

A.6. Accuracy and limits

A.6.1. Neglected terms

The approximation has three sources.

First, expansion of the vertical force omits terms cubic in $z$, giving vertical acceleration errors $O(an_a^2i^3)$. Second, expansion of its eccentric modulation omits $O(e^2)$ coefficients, giving errors $O(an_a^2e^2i)$. Third, retaining the supplied planar motion omits both its $O(e^2)$ residual and the $O(i^2)$ horizontal force from the vertical displacement.

More explicitly, substituting the reconstructed trajectory into the full $J_2$ equations gives

$$
\begin{aligned}
\ddot r-r\dot\theta^2+\frac{\partial\Phi}{\partial r}
&=
O\!\left(an_a^2[e^2+i^2]\right),\\
\frac{d}{dt}(r^2\dot\theta)&=0,\\
\ddot z+\frac{\partial\Phi}{\partial z}
&=
O\!\left(an_a^2[e^2i+i^3]\right).
\end{aligned}
\tag{A32}
$$

The coefficients implicit in these estimates are bounded away from marginal radial stability. In particular, the perturbation requires $e\alpha_a$, $e\beta_a$, and $e\gamma_a$ to be small. Since $\beta_a,\gamma_a$ grow as $\kappa_a$ becomes small, merely requiring $e\ll1$ is insufficient near $\kappa_a=0$.

The frequencies and mean precession rate also acquire corrections of relative order $e^2$ and $i^2$. Those omitted corrections can accumulate into secular phase errors over long integrations. Equation (A32) is a local residual estimate, not a uniform error bound for arbitrarily long times.

No powers of $q_a$ were discarded in deriving the displayed coefficients. Neglecting higher planetary harmonics is a separate assumption, not an error of order $e$ or $i$.

A.6.2. Circular limit

When $e=0$, the modulation functions become $F=K=1$ and $G=H=0$. Equations (A19) and (A25) reduce to

$$
\begin{aligned}
z&=ai\sin(\lambda-\Omega),\\
v_z&=a\nu_ai\cos(\lambda-\Omega),\\
i&=\sqrt{
\left(\frac{z}{a}\right)^2+
\left(\frac{v_z}{a\nu_a}\right)^2
},\\
\Omega&=
\lambda-
\operatorname{atan2}
\left(\frac{z}{a},\frac{v_z}{a\nu_a}\right).
\end{aligned}
\tag{A33}
$$

This is the expected vertical harmonic oscillator. Its forward and inverse transformations are exact, and its remaining physical errors arise from the omitted nonlinear inclination terms.

A.6.3. Nodal regression

Using the exact difference-of-squares identity, Equation (A12) becomes

$$
\boxed{
\dot\Omega
=
\Omega_a-\nu_a
=
-\frac{3n_a^2q_a}{\Omega_a+\nu_a}
=
-\frac32 n_aq_a+O(n_aq_a^2).
}
\tag{A34}
$$

Thus an oblate planet produces retrograde nodal precession. The last expression is only a small-$q_a$ check; the algorithm uses the unexpanded frequency difference.

A.6.4. Spherical-planet limit

When $J_2=0$,

$$
\Omega_a=\kappa_a=\nu_a=n_a,
\qquad
\alpha_a=1,
\qquad
\beta_a=2,
\qquad
\gamma_a=1,
\qquad
\dot\Omega=0.
\tag{A35}
$$

The vertical displacement then becomes

$$
z
=
ai\left[
(1-e\cos M)\sin(\lambda-\Omega)
+2e\sin M\cos(\lambda-\Omega)
\right].
\tag{A36}
$$

But the planar solution has $r=a(1-e\cos M)+O(ae^2)$ and $\theta=\lambda+2e\sin M$. Expanding the usual inclined-plane relation therefore gives

$$
\begin{aligned}
z
&=
ir\sin(\theta-\Omega)
+O\!\left(a[e^2i+i^3]\right),\\
v_z
&=
i\left[
v_r\sin(\theta-\Omega)
+v_\theta\cos(\theta-\Omega)
\right]
+O\!\left(an_a[e^2i+i^3]\right).
\end{aligned}
\tag{A37}
$$

The $i^3$ term accounts for replacing the exact geometrical factor $\tan i$ by $i$ when $r$ is cylindrical radius. Thus the solution reproduces a fixed tilt of the eccentric planar orbit through order $ei$, and the node is constant, as required in a spherical potential.

A.6.5. Planar limit

When $i=0$, Equation (A19) gives $z=v_z=0$ identically. The entire drift reduces to the supplied planar algorithm. The undefined node has no effect on the trajectory.

A.7. Validation status

The force expansion, coefficient matching, inverse, and limiting cases above provide analytic checks. Before code acceptance, the style guide additionally requires a SymPy verification and equation ledger, a blind re-derivation, and comparison with direct integration of Equation (A1). Those external checks have not been performed here.

The direct-integration comparison should verify the residual orders in Equation (A32), while a separate forward–backward test should recover the initial coordinates to floating-point precision. These test different properties: physical accuracy and algebraic reversibility.

Final results

- The circular vertical frequency is $\nu_a=n_a[1+(9/2)q_a]^{1/2}$.
- The vertical coordinate map is Equation (A19), with eccentric modulation factors from Equations (A16) and (A18).
- Its exact inverse is Equations (A24) and (A25); the full determinant must be retained.
- During a drift, $i$ is constant and the mean node advances by $(\Omega_a-\nu_a)\Delta t$. Advance the vertical quadratures and remodulate at the new $M$.
- The map is exactly reversible in exact arithmetic. Its physical acceleration errors are $O(e^2,i^2)$ horizontally and $O(e^2i,i^3)$ vertically.
