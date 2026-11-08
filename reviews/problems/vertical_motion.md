## Background

Epi_int_lite is an N-body code that simulates narrow planetary rings as streamlines of particles orbiting
an oblate planet of mass $M_p$, radius $R_p$ and second zonal harmonic $J_2$. All particles move in the
planet's equatorial plane, $z = 0$. Each timestep has a "drift" step, where a particle moves only under
the planet's gravity, and that drift is done by converting the particle's coordinates
$r, \theta, v_r, v_\theta$ to epicyclic orbit elements, advancing the elements, and converting back.

The in-plane solution it uses (Hahn et al. 2025, ApJ 986, 220, Appendix A) is the following. The mean
angular velocity and epicyclic frequency at semimajor axis $a$ are

$$\Omega_a = \sqrt{\frac{GM_p}{a^3}}\left[1 + \frac{3}{2}J_2\left(\frac{R_p}{a}\right)^2\right]^{1/2}, \qquad
\kappa_a = \sqrt{\frac{GM_p}{a^3}}\left[1 - \frac{3}{2}J_2\left(\frac{R_p}{a}\right)^2\right]^{1/2},$$

and in terms of the epicyclic elements $a$, $e$, $\tilde{\omega}$ (longitude of periapse) and $M$ (mean
anomaly) the particle's coordinates are

$$r = \frac{a}{\sqrt{1 + 2e\cos M}}, \qquad \theta = \frac{\Omega_a}{\kappa_a}(M + 2e\sin M) + \tilde{\omega},$$
$$v_r = ea\kappa_a\sin M\left(\frac{r}{a}\right)^3, \qquad v_\theta = \frac{a^2\Omega_a}{r}.$$

These conserve the specific angular momentum $h = rv_\theta = a^2\Omega_a$ exactly and solve the planar
equation of motion with errors of order $O(e^2)$. During a drift of duration $\Delta t$ the mean anomaly
advances as $M \rightarrow M + \kappa_a\Delta t$, with $a$, $e$ and $\tilde{\omega}$ constant. The equations
are inverted exactly: $a$ from $h = rv_\theta$ via $h^2 = a^4\Omega_a^2$; then $\Omega_a$, $\kappa_a$;
then $e\sin M = (v_r/a\kappa_a)(a/r)^3$ and $e\cos M = \frac{1}{2}[(a/r)^2 - 1]$; then
$\tilde{\omega} = \theta - (\Omega_a/\kappa_a)(M + 2e\sin M)$. Because the conversion in both directions
is exact, the drift step introduces no numerical error of its own.

Notation follows Goldreich & Tremaine: $\Omega_a$, $\kappa_a$ are frequencies at semimajor axis $a$; the
symbol $\Omega$ (no subscript) is reserved for the longitude of ascending node.

## The problem

Extend this solution into the vertical direction, so that particles can have a vertical displacement
$z \neq 0$ and vertical velocity $v_z \neq 0$, for orbits that are nearly circular and nearly equatorial.

1. Starting from the planet's gravitational potential, including $J_2$, derive the particle's vertical
   equation of motion to lowest order in $z$, and the frequency of its vertical oscillation at
   semimajor axis $a$.
2. Derive expressions for $z$ and $v_z$ analogous to the in-plane ones above, in terms of new orbit
   elements, the inclination $i$ and the longitude of ascending node $\Omega$, plus the in-plane elements.
3. Give the exact inverse: $i$ and $\Omega$ from the coordinates.
4. Say how $i$ and $\Omega$ change during a drift of duration $\Delta t$, and how the drift step should
   advance $z$ and $v_z$, so that the conversion stays exactly reversible.
5. State the accuracy of the result: which terms are neglected, and of what order in $e$ and $i$.
