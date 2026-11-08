# derive: anthropic/claude-opus-5.5

Problem: `reviews/problems/vertical_motion.md` · 2026-11-08 · tokens in 5581, out 27678 · cost $0.5759

---

# Appendix A. Vertical Motion in the Epicyclic Drift Step

This appendix extends the epicyclic orbit elements of hahn25 (Appendix A) out of the planet's equatorial plane. The aim is to let the drift step of epi_int_lite advance particles that have a small vertical displacement $z$ and vertical velocity $v_z$, while keeping the property that makes that step reliable: the conversion from coordinates to elements and back is exact.

## A.1. Setting up the problem

During a drift step a ring particle moves only under the planet's gravity. The ring's self-gravity, the satellites and viscosity are applied in the code's separate "kick" steps, so we ignore them here. The planet is axisymmetric and is described by its mass $M_p$, radius $R_p$ and second zonal harmonic $J_2$.

We neglect the higher zonal harmonics $J_4, J_6, \ldots$ for the same reason hahn25 does: in the ring region their contributions are smaller than the $J_2$ term by further factors of $(R_p/a)^2$, and the in-plane solution we are extending does not contain them either. The particle is a test particle, and it is always outside the planet, where the multipole expansion of the planet's potential converges.

We use cylindrical coordinates $(r, \theta, z)$ centred on the planet, with $z = 0$ its equatorial plane. Here $r$ is the cylindrical radius, $\theta$ the longitude measured in the equatorial plane, and $v_r = \dot r$, $v_\theta = r\dot\theta$, $v_z = \dot z$.

We seek orbits that are nearly circular and nearly equatorial, so $e \ll 1$ and $i \ll 1$. The solution must reduce to the hahn25 solution when $z = v_z = 0$. This splits the problem into two parts, treated in turn:

- the vertical equation of motion and its oscillation frequency (Section A.2);
- a vertical solution built on top of the existing in-plane one (Sections A.3–A.4).

Its inverse and the drift step follow in Sections A.5–A.6, and its accuracy and its limits in Sections A.7–A.8.

Two named dimensionless ratios measure the small quantities other than $e$ and $i$. They are the Keplerian mean motion at $a$ and the oblateness parameter at $a$:

$$n_a = \sqrt{\frac{GM_p}{a^3}}, \qquad \epsilon_a = J_2\left(\frac{R_p}{a}\right)^2, \qquad (A1)$$

where $\epsilon_a$ is the fractional size of the planet's quadrupole potential at the ring. For Saturn's rings $\epsilon_a \sim 10^{-3}$. In this notation the in-plane frequencies of hahn25 are $\Omega_a = n_a(1 + \frac{3}{2}\epsilon_a)^{1/2}$ and $\kappa_a = n_a(1 - \frac{3}{2}\epsilon_a)^{1/2}$.

## A.2. The potential and the equations of motion

We start from the exterior potential of an axisymmetric planet truncated at $J_2$ (e.g. Murray & Dermott 1999, Chapter 6; equation number to be added):

$$\Phi(r, z) = -\frac{GM_p}{r_s}\left[1 - J_2\left(\frac{R_p}{r_s}\right)^2 P_2\!\left(\frac{z}{r_s}\right)\right], \qquad (A2)$$

where $r_s = (r^2 + z^2)^{1/2}$ is the distance from the planet's centre and $P_2(x) = (3x^2 - 1)/2$ is the Legendre polynomial. In the equatorial plane $P_2 = -\frac{1}{2}$, so the oblate planet's potential well there is deeper than a point mass's, as it should be.

Because the orbits are nearly equatorial, $|z| \ll r$, and we can expand Equation (A2) in the small quantity $z/r$. Using
$$r_s^{-1} = r^{-1}\left[1 - \tfrac{1}{2}(z/r)^2\right], \qquad r_s^{-3} = r^{-3}\left[1 - \tfrac{3}{2}(z/r)^2\right], \qquad z^2 r_s^{-5} = z^2 r^{-5},$$
each with errors of order $O(z^4/r^4)$, we obtain

$$\Phi(r, z) = \Phi_0(r) + \frac{1}{2}\nu^2(r)\,z^2 + O\!\left(\frac{GM_p z^4}{r^5}\right), \qquad (A3)$$

where

$$\Phi_0(r) = -\frac{GM_p}{r}\left[1 + \frac{1}{2}J_2\left(\frac{R_p}{r}\right)^2\right], \qquad (A4)$$

is the equatorial potential used by hahn25, and

$$\nu^2(r) = \left.\frac{\partial^2\Phi}{\partial z^2}\right|_{z=0} = \frac{GM_p}{r^3}\left[1 + \frac{9}{2}J_2\left(\frac{R_p}{r}\right)^2\right] \qquad (A5)$$

is the square of the vertical oscillation frequency. In Equation (A5) the "1" comes from the point mass. The factor $\frac{9}{2}$ collects $\frac{3}{2}$ from expanding $r_s^{-3}$ in the $J_2$ term and $3$ from the explicit $3z^2/r_s^2$ in $P_2$. Both contributions are positive, because the planet's equatorial bulge pulls a particle that is displaced out of the plane back toward it.

The equations of motion are $\ddot{\mathbf{r}} = -\nabla\Phi$. Since $\Phi$ is independent of $\theta$, the azimuthal equation conserves the specific angular momentum about the planet's spin axis exactly, even when $z \neq 0$:

$$h = r v_\theta = r^2\dot\theta = \text{constant}. \qquad (A6)$$

The radial and vertical equations, to lowest order in $z$, are

$$\ddot r - r\dot\theta^2 = -\frac{d\Phi_0}{dr} - \frac{1}{2}z^2\frac{d\nu^2}{dr} + O\!\left(\frac{GM_p z^4}{r^6}\right), \qquad (A7)$$

$$\ddot z = -\nu^2(r)\,z\left[1 + O(z^2/r^2)\right]. \qquad (A8)$$

Equation (A8) is the vertical equation of motion. It describes a harmonic oscillator whose frequency depends on the particle's instantaneous radius, and whose neglected terms are smaller by $O(z^2/r^2)$.

Equation (A7) shows the one way in which the vertical motion feeds back on the in-plane motion. The $z^2$ term changes the radial force by a fraction $O(z^2/r^2)$, so when $z = O(ia)$ the in-plane solution of hahn25 acquires errors of order $O(i^2)$ in addition to its $O(e^2)$ errors.

The vertical frequency at semimajor axis $a$ is $\nu_a \equiv \nu(a)$:

$$\nu_a = \sqrt{\frac{GM_p}{a^3}}\left[1 + \frac{9}{2}J_2\left(\frac{R_p}{a}\right)^2\right]^{1/2}. \qquad (A9)$$

Like $\Omega_a$ and $\kappa_a$, this is exact for the $J_2$-truncated potential; no expansion in $\epsilon_a$ has been made.

The three frequencies are tied together by an exact identity. Outside the planet $\Phi$ obeys Laplace's equation, which in the equatorial plane reads $r^{-1}\partial_r(r\,\partial_r\Phi_0) + \nu^2 = 0$. Writing $\Omega^2_a = a^{-1}\Phi_0'(a)$ and $\kappa_a^2 = \Phi_0''(a) + 3a^{-1}\Phi_0'(a)$, the first term equals $\kappa_a^2 - 2\Omega_a^2$, so

$$\nu_a^2 + \kappa_a^2 = 2\Omega_a^2. \qquad (A10)$$

This is the vacuum form of the familiar relation $\nu^2 + \kappa^2 = 2\Omega^2 + 4\pi G\rho$. It is satisfied by Equations (A1) and (A9) exactly, and it is used below to simplify the accuracy estimates.

To first order in $\epsilon_a$, the differences between the three frequencies are

$$\nu_a - \Omega_a = \Omega_a - \kappa_a = \frac{3}{2}\epsilon_a n_a\left[1 + O(\epsilon_a)\right]. \qquad (A11)$$

## A.3. The vertical solution for a circular orbit

We first find the solution in the limit $e = 0$, because it fixes the phase of the vertical oscillation. On a circular orbit $r = a$ and $\dot\theta = \Omega_a$ exactly. Equation (A8) is then a harmonic oscillator with constant frequency $\nu_a$, with the solution
$$z = a\sin i\,\sin(\nu_a t + \text{const}).$$

We want to write this phase in terms of the particle's longitude $\theta$ and a longitude of ascending node $\Omega$, as $\nu_a t + \text{const} = \theta - \Omega$. Since $\theta$ advances at $\Omega_a$, this is possible only if the node moves at the rate

$$\dot\Omega = \Omega_a - \nu_a = -\frac{3}{2}\epsilon_a n_a\left[1 + O(\epsilon_a)\right]. \qquad (A12)$$

The node therefore regresses, because the bulge makes the vertical oscillation faster than the orbital motion. Equation (A12) is the familiar nodal regression rate due to $J_2$; this is a first check of the solution.

## A.4. The vertical solution for an eccentric orbit

For $e \neq 0$ we have to decide how $z$ depends on the in-plane motion. The physical guide is the Keplerian limit $J_2 = 0$. In that limit an inclined orbit has
$$z = r\sin i\,\sin u \quad \text{exactly},$$
where $u$ is the argument of latitude. To first order in $i$, the argument of latitude equals the equatorial longitude measured from the node, $u = \theta - \Omega + O(i^2)$.

This suggests writing the vertical solution with the particle's actual radius $r$ and longitude $\theta$, taken from the in-plane solution, together with the regressing node of Equation (A12). The vertical elements are the inclination $i$ and the longitude of ascending node $\Omega$:

$$z = r\sin i\,\sin(\theta - \Omega), \qquad (A13)$$

$$v_z = \sin i\left[v_r\sin(\theta - \Omega) + w\cos(\theta - \Omega)\right], \qquad (A14)$$

where $r$, $\theta$, $v_r$ and $v_\theta$ are the in-plane coordinates given by hahn25's Equation (A9) in terms of $a$, $e$, $M$ and $\tilde\omega$, and

$$w = v_\theta - (\Omega_a - \nu_a)\,r = \frac{a^2\Omega_a}{r} + (\nu_a - \Omega_a)\,r. \qquad (A15)$$

The quantity $w = r(\dot\theta - \dot\Omega)$ is the particle's azimuthal velocity measured in a frame that rotates with the line of nodes. Equations (A13) and (A14) are thus the Keplerian expressions for $z$ and $\dot z$ written in that rotating frame.

Equation (A14) is the exact time derivative of Equation (A13) along the drift trajectory. To see this, we use three facts about that trajectory:
- the hahn25 solution satisfies $\dot r = v_r$ exactly, since differentiating $r = a(1 + 2e\cos M)^{-1/2}$ with $\dot M = \kappa_a$ gives $v_r$;
- it satisfies $\dot\theta = v_\theta/r = a^2\Omega_a/r^2$ exactly;
- during a drift $\dot\Omega = \Omega_a - \nu_a$ and $\dot i = 0$.

Equations (A13)–(A14) are therefore kinematically consistent. All of their error lies in how well they satisfy the dynamics, Equation (A8), which is examined in Section A.7.

The physical requirements fix what the algebra leaves open. The inclination must be non-negative, so $\sin i \geq 0$ with $0 \leq i \leq \pi/2$. $\Omega$ must be the ascending node, and indeed at $\theta = \Omega$ Equations (A13)–(A14) give $z = 0$ and $v_z = w\sin i > 0$, since $w > 0$. This holds because $\nu_a > \Omega_a$ for an oblate planet, which makes the second term of Equation (A15) positive.

We write the amplitude as $\sin i$ rather than $i$ so that Equation (A13) is the geometric Keplerian expression. The two differ at $O(i^3)$, which is beyond the accuracy of the solution (Section A.7).

We also keep $\Omega$ as the physical node, which regresses during a drift. The alternative is the device hahn25 uses for $\tilde\omega$, a constant element with the precession absorbed into the phase (here $\theta - \Omega$ would be replaced by $\theta + (\nu_a - \Omega_a)M/\kappa_a - \tilde\Omega$ with $\tilde\Omega$ fixed). That alternative is equally exact. We chose the physical node because it can be compared directly with nodal regression and bending-wave theory without first converting the element.

## A.5. The inverse

The inverse uses the fact that the in-plane conversion of hahn25 involves only $r$, $\theta$, $v_r$ and $v_\theta$. It is therefore unchanged when $z \neq 0$. In particular, $a$ is still obtained from $h = rv_\theta$, which Equation (A6) shows is exactly conserved in three dimensions.

Once $a$ is known, so are $\Omega_a$, $\kappa_a$, $\nu_a$ (Equation (A9)) and $w$ (Equation (A15)). Equations (A13) and (A14) are then two linear equations for $\sin i\,\sin(\theta - \Omega)$ and $\sin i\,\cos(\theta - \Omega)$. Multiplying Equation (A14) by $r$ and subtracting $v_r$ times Equation (A13) eliminates the sine term, and gives

$$\sin i\,\sin(\theta - \Omega) = \frac{z}{r}, \qquad \sin i\,\cos(\theta - \Omega) = \frac{r v_z - z v_r}{r w}. \qquad (A16)$$

Choosing the positive root of the sum of squares, as required by $\sin i \geq 0$, gives

$$\sin i = \left[\left(\frac{z}{r}\right)^2 + \left(\frac{r v_z - z v_r}{r w}\right)^2\right]^{1/2}, \qquad (A17)$$

$$\Omega = \theta - \operatorname{atan2}\!\left(\frac{z}{r},\ \frac{r v_z - z v_r}{r w}\right). \qquad (A18)$$

Equations (A17)–(A18) are the exact inverse of Equations (A13)–(A14). They need only $w \neq 0$, which always holds since $w > 0$ (Section A.4).

$\Omega$ is undefined when $\sin i = 0$, just as $\tilde\omega$ is undefined when $e = 0$. In that case any value reproduces $z = v_z = 0$ on the return conversion, so reversibility is not affected.

The numerator $r v_z - z v_r$ in Equation (A16) is, up to sign, the horizontal component of the particle's angular momentum per unit $r$ along the direction $\hat\theta$. In the Keplerian limit Equation (A16) therefore reduces to the familiar way of obtaining the node and inclination from the angular momentum vector.

## A.6. The drift step

During a drift of duration $\Delta t$, the in-plane elements evolve as in hahn25: $M \rightarrow M + \kappa_a\Delta t$, with $a$, $e$ and $\tilde\omega$ constant. The vertical elements evolve as

$$i \rightarrow i, \qquad \Omega \rightarrow \Omega + (\Omega_a - \nu_a)\,\Delta t. \qquad (A19)$$

The inclination is constant because the particle moves in a potential with no $\theta$-dependence and no time-dependence, apart from the terms neglected in Section A.7. The node regresses at the rate of Equation (A12).

The step does not advance $z$ and $v_z$ by integrating them. Instead, it recomputes them from Equations (A13)–(A14) using the advanced elements together with the advanced in-plane coordinates $r$, $\theta$ and $v_r$. The vertical phase therefore advances by
$$\Delta(\theta - \Omega) = \Delta\theta + (\nu_a - \Omega_a)\Delta t,$$
which reduces to $\nu_a\Delta t$ on a circular orbit.

As an algorithm, in the order the code runs it, the drift is:

1. Convert $(r, \theta, v_r, v_\theta)$ to $(a, e, M, \tilde\omega)$ with steps 1–4 of hahn25 Appendix A, computing $\Omega_a$ and $\kappa_a$ along the way. This step is exact.
2. Calculate $\nu_a$ from Equation (A9) and $w$ from Equation (A15). Then calculate $\sin i$ and $\Omega$ from Equations (A17)–(A18). This step is exact.
3. Advance $M \rightarrow M + \kappa_a\Delta t$ and $\Omega \rightarrow \Omega + (\Omega_a - \nu_a)\Delta t$, holding $a$, $e$, $\tilde\omega$ and $i$ fixed. This step is exact for the solution (A13)–(A14).
4. Convert back: calculate $r$, $\theta$, $v_r$ and $v_\theta$ from hahn25's Equation (A9), then $w$ from Equation (A15), and then $z$ and $v_z$ from Equations (A13)–(A14). This step is exact.

Steps 1, 2 and 4 are exact inverses of each other, so the round trip from coordinates to elements and back introduces no error of its own. The only departure of the drifted particle from its true trajectory is the extent to which Equations (A13)–(A14) and hahn25's Equation (A9) fail to satisfy the equations of motion (A7)–(A8). That is evaluated next.

## A.7. Accuracy

To find which terms the solution neglects, we substitute Equation (A13) into the vertical equation of motion (A8). We write the vertical phase as $\phi = \theta - \Omega$, so that $\dot\phi = h/r^2 + \delta$ with $\delta \equiv \nu_a - \Omega_a$. Differentiating Equation (A13) twice and using $d(r^2\dot\theta)/dt = 0$ (Equation (A6)), the residual $\mathcal{R} \equiv \ddot z + \nu^2(r)z$ is

$$\mathcal{R} = \sin i\left\{\left[E_r + D(r)\right]\sin\phi + 2\delta\,v_r\cos\phi\right\} + O(i^3 a n_a^2), \qquad (A20)$$

where

$$E_r = \ddot r - \frac{h^2}{r^3} + \frac{d\Phi_0}{dr} \qquad (A21)$$

is the residual of the planar radial equation of motion left by the hahn25 solution, and

$$D(r) = r\nu^2(r) - \frac{d\Phi_0}{dr} - \frac{2h\delta}{r} - r\delta^2. \qquad (A22)$$

The $O(i^3)$ remainder in Equation (A20) collects two contributions: the $z^3$ term of the vertical force (Equation (A8)), and the $O(i^2)$ in-plane error of Equation (A7) acting on $z$.

We take each term of Equation (A20) in turn.

**The $E_r$ term.** This is the in-plane error carried into the vertical motion through $r$. hahn25 shows it is $O(e^2 a\Omega_a^2)$. With $z \neq 0$ it also acquires the $O(i^2)$ term of Equation (A7).

**The $D(r)$ term.** This term measures how well a single constant node rate $\dot\Omega$ suits a particle whose radius oscillates. At $r = a$ it vanishes identically, since there $D(a) = a[\nu_a^2 - (\Omega_a + \delta)^2] = 0$.

Since $r - a = -ae\cos M + O(e^2 a)$, we can Taylor expand $D$ in the small quantity $(r - a)/a$. Using $\Phi_0'' = \kappa_a^2 - 3\Omega_a^2$ and Equation (A10), the expansion is
$$D(r) = \left[(1 + \gamma_a)\nu_a^2 - 2\Omega_a^2 + 4\Omega_a\nu_a\right](r - a) + O(e^2 a n_a^2),$$
where $\gamma_a = d\ln\nu^2/d\ln r|_a = -3 - 9\epsilon_a + O(\epsilon_a^2)$. The bracket is exactly zero in the Keplerian limit and equals $-9\epsilon_a n_a^2$ to first order in $\epsilon_a$. Therefore
$$D(r) = 9\epsilon_a n_a^2\, ae\cos M\left[1 + O(\epsilon_a)\right] + O(e^2 a n_a^2).$$

**The $2\delta v_r$ term.** This term is the Coriolis-like cross term $2\dot r\dot\phi$, whose part $2\dot r\dot\theta$ is cancelled exactly by angular momentum conservation. The constant node rate $\delta$ has no such partner. By Equation (A11) the term equals $3\epsilon_a n_a v_r = 3\epsilon_a n_a^2\,ae\sin M$ to first order in $\epsilon_a$.

Dividing by the size $a n_a^2\sin i$ of the vertical force, Equations (A13)–(A14) satisfy the vertical equation of motion to accuracy $O(i)$, with fractional errors of order

$$O(e^2), \qquad O(e\,\epsilon_a) = O\!\left(eJ_2(R_p/a)^2\right), \qquad O(i^2). \qquad (A23)$$

The in-plane solution, for its part, now has errors $O(e^2)$ and $O(i^2)$.

The cross term $O(e\epsilon_a)$ is the only error that is first order in $e$. It is smaller than the $O(e)$ error of the simpler choice $z = a\sin i\sin(\nu_a t + \text{const})$, which ignores the radial and azimuthal epicycles, by the factor $\epsilon_a \sim 10^{-3}$. Avoiding that larger error is why Equation (A13) uses the actual $r$ and $\theta$.

The solution also assumes that $\epsilon_a$ is small enough for the series of zonal harmonics to be truncated at $J_2$, with neglected terms $O(J_4(R_p/a)^4)$, as in hahn25.

The orders in Equation (A23), including the coefficients 9 and 3 of the $O(e\epsilon_a)$ terms, are to be checked against an orbit integration that makes none of these approximations. That integration should solve the full equations of motion in the potential (A2), varying $e$, $i$ and $\epsilon_a$ separately. Any disagreement should then be traced to one of the terms in Equation (A20).

## A.8. Limits

The solution passes four limit checks.

**Circular orbit ($e = 0$).** Here $r = a$, $v_r = 0$ and $w = a\nu_a$. Equations (A13)–(A14) become $z = a\sin i\,\sin\phi$ and $v_z = a\nu_a\sin i\,\cos\phi$, with $\dot\phi = \nu_a$. This is the exact solution of the linearized Equation (A8). Correspondingly, $E_r$, $D$ and $v_r$ all vanish in Equation (A20).

**Spherical planet ($J_2 = 0$).** Here $\nu_a = \Omega_a = \kappa_a$, so the node is fixed and $w = v_\theta$. Equations (A13)–(A14) reduce to the Keplerian $z = r\sin i\sin u$ and $\dot z = \sin i(\dot r\sin u + r\dot u\cos u)$ to first order in $i$. The residual (A20) reduces to $\sin i\,E_r\sin\phi$, so the vertical solution is exactly as accurate as the in-plane one.

**Equatorial orbit ($i = 0$).** Equations (A13)–(A14) give $z = v_z = 0$. Equation (A7) loses its $z^2$ term, and the scheme reduces exactly to that of hahn25.

**Precession rates.** By Equation (A11), the apsidal precession rate $\Omega_a - \kappa_a$ and the nodal regression rate $\Omega_a - \nu_a$ are equal and opposite to first order in $\epsilon_a$, which is the classical result for $J_2$. Using Equation (A10), their sum is exactly $2\Omega_a - \kappa_a - \nu_a = [(\Omega_a - \kappa_a)^2 + (\nu_a - \Omega_a)^2]/2\Omega_a = \frac{9}{4}\epsilon_a^2 n_a[1 + O(\epsilon_a)]$, which is positive and second order in $J_2$.

## Final results

1. **Vertical equation of motion.** $\ddot z = -\nu^2(r)z[1 + O(z^2/r^2)]$, with $\nu^2(r) = (GM_p/r^3)[1 + \frac{9}{2}J_2(R_p/r)^2]$ (Equations (A5), (A8)).

2. **Vertical frequency.** $\nu_a = (GM_p/a^3)^{1/2}[1 + \frac{9}{2}J_2(R_p/a)^2]^{1/2}$, which obeys $\nu_a^2 + \kappa_a^2 = 2\Omega_a^2$ exactly (Equations (A9), (A10)).

3. **Coordinates from elements.**
   - $z = r\sin i\,\sin(\theta - \Omega)$;
   - $v_z = \sin i[v_r\sin(\theta - \Omega) + w\cos(\theta - \Omega)]$, with $w = a^2\Omega_a/r + (\nu_a - \Omega_a)r$;
   - $r$, $\theta$ and $v_r$ are taken from the in-plane solution of hahn25 (Equations (A13)–(A15)).

4. **Exact inverse.**
   - $\sin i = [(z/r)^2 + ((rv_z - zv_r)/rw)^2]^{1/2}$;
   - $\Omega = \theta - \operatorname{atan2}(z/r,\ (rv_z - zv_r)/rw)$ (Equations (A17)–(A18)).

5. **Drift.**
   - $i$ is constant, and $\Omega \rightarrow \Omega + (\Omega_a - \nu_a)\Delta t$, which is nodal regression at $\approx -\frac{3}{2}J_2(R_p/a)^2 n_a$;
   - $z$ and $v_z$ are recomputed from Equations (A13)–(A14) with the advanced elements and in-plane coordinates;
   - the coordinate–element conversion is exactly reversible (Equation (A19), Section A.6).

6. **Accuracy.**
   - The vertical motion is accurate to $O(i)$, with fractional errors $O(e^2)$, $O(eJ_2(R_p/a)^2)$ and $O(i^2)$.
   - The in-plane motion has errors $O(e^2)$ and $O(i^2)$.
   - $J_4$ and higher harmonics are neglected (Equation (A23)).
