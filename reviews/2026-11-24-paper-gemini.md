# paper: google/gemini-4-argon

Input: `paper/epi_int_3d.tex` · 2026-11-24 · tokens in 43783, out 10073 · cost $0.1883

---

### Finding 1 (Author's specific questions)
- **Type:** other
- **Severity:** high
- **Where:** Whole paper (`paper/epi_int_3d.tex`, Sections 1–4 and Appendices A–E)
- **Finding:** Responses to the author's three broad questions on simplifying the text and math, summarizing the physics of each derivation, and reorganizing the paper.
- **Suggestion:**

#### (1) Where descriptions, derivations and explanations can be made simpler, shorter, or more correct

1. **Section 2.2, Eqns. (10a)–(10b) and Section 3.4 (Vertical epicyclic solution and drift):**
   - *Problem:* Writing the unperturbed vertical harmonic oscillator in terms of the true longitude $\theta$ rather than the mean longitude $\lambda = (\Omega_a/\kappa_a)M + \tilde{\omega}$ means that differentiating Eqn. (10a) with respect to time gives $\dot{z} = a i (\dot{\theta} - \dot{\Omega})\cos(\theta - \Omega)$. Because $\dot{\theta} = \Omega_a[1 + 2e\cos M + \mathcal{O}(e^2)]$, the time derivative $\dot{z}$ differs from Eqn. (10b), $v_z = a\nu_a i\cos(\theta - \Omega)$, at order $\mathcal{O}(ei)$. That is consistent with the stated $\mathcal{O}(ei)$ accuracy, and makes exact spatial-to-element inversion trivial, but a careful reader who differentiates Eqn. (10a) will think Eqn. (10b) is missing a factor of $1 + 2e\cos M$.
   - *Suggested fix:* Add one clause after Eqn. (11): "Equations (10) replace the mean longitude $(\Omega_a/\kappa_a)M + \tilde{\omega}$ with $\theta$, which differs from it by $\mathcal{O}(e)$ and so alters $z$ and $v_z$ only at the already-neglected order $\mathcal{O}(ei)$."

2. **Section 3.1.1, Eqns. (13)–(14), and Appendix A, Eqns. (A1)–(A8) (Streamline wire geometry):**
   - *Problem:* Appendix A takes eight equations and 45 lines to find the component of a vector $\mathbf{d}_0$ perpendicular to a unit vector $\hat{\mathbf{t}}$. It introduces an unknown longitude step $\delta\theta'$ to the near point, writes the chord $\mathbf{d}_1(\delta\theta')$, and solves $\boldsymbol{\Delta}\cdot\hat{\mathbf{t}} = 0$ for $\delta\theta'$, even though Eqn. (A6), $\boldsymbol{\Delta} = \mathbf{d}_0 - (\mathbf{d}_0\cdot\hat{\mathbf{t}})\hat{\mathbf{t}}$, is an elementary vector projection that requires neither $\delta\theta'$ nor solving for the near point's coordinate first.
   - *Suggested fix:* In Appendix A, define $\mathbf{d}_0 = d_r\hat{\mathbf{r}} + d_z\hat{\mathbf{z}}$ and the unit tangent $\hat{\mathbf{t}} = (s_r\hat{\mathbf{r}} + \hat{\boldsymbol{\theta}} + s_z\hat{\mathbf{z}})/N$ from the streamline path $\mathbf{r}'(\theta')$. Then go straight to Eqn. (A6): because $\boldsymbol{\Delta}$ meets the straight wire at a right angle, $\boldsymbol{\Delta} = \mathbf{d}_0 - (\mathbf{d}_0\cdot\hat{\mathbf{t}})\hat{\mathbf{t}}$, and $\mathbf{d}_1 = -(\mathbf{d}_0\cdot\hat{\mathbf{t}})\hat{\mathbf{t}} = -Nc\,\hat{\mathbf{t}}$ with $c = (\mathbf{d}_0\cdot\hat{\mathbf{t}})/N = (d_r s_r + d_z s_z)/N^2$. This cuts Eqns. (A1), (A4), and the intermediate solve for $\delta\theta'$ while arriving at Eqns. (A7)–(A8) in half the space.

3. **Section 4.1.2, Eqns. (17)–(21) (Satellite vertical forcing):**
   - *Problem:* Equation (17) writes the vertical acceleration with $(z - z_s)$, the prose drops $z$ to get $A_z \simeq GM_s z_s / \Delta_s^3$, and then Eqn. (18) expands it in a double-sided sum over $m = -\infty$ to $\infty$ without showing the step $z_s = a_s i_s \sin(\theta_s - \Omega_s) = a_s i_s \sin(\nu_{a_s} t)$ or how the product of $\sin(\nu_{a_s} t)$ and the cosine series for $\Delta_s^{-3}$ becomes a single sine series.
   - *Suggested fix:* Since $z_s = a_s i_s \sin(\nu_{a_s}t)$ is independent of $\theta$, state that $\Delta_s^{-3}$ expands in a cosine series in $\phi = \theta - \theta_s$ with coefficients $b_{3/2}^{(m)}(\beta)/a_s^3$, and that combining $\cos m(\theta - \Omega_{a_s}t)$ with $\sin(\nu_{a_s}t)$ via $2\cos A\sin B = \sin(A+B) - \sin(A-B)$ gives Eqn. (18). Also drop $z$ in Eqn. (17) from the outset ("Since $|z| \ll |z_s|$, the satellite accelerates a ring particle vertically by $A_z \simeq GM_s z_s/\Delta_s^3$..."), which avoids defining a term that is immediately discarded.

4. **Section 4.1.4 and Section 4.2.4 (Wavenumber discrepancy in the finite lattice):**
   - *Problem:* Section 4.1.4 explains the simulated wavenumber being larger than Eqn. (27) by saying "This is because the simulated ring is a finite lattice of streamlines rather than the continuous sheet that Eqn. (27) assumes." Actually, two opposite effects are at work, and the equation tracker (`notes/equations.md`, ID B13) has the exact formula: (i) the finite lattice spacing $\Delta$ *weakens* the restoring gravity from $2\pi G\sigma|k|$ to $X(k\Delta) = (4G\sigma/\Delta)[\pi k\Delta/2 - (k\Delta)^2/4]$, which by itself *increases* $|k|$ for a given $D|x|$; and (ii) the finite radial width of the ring omits the gravity of streamlines beyond the inner and outer edges, $-2G\sigma(1/d_{\rm in} + 1/d_{\rm out})$, which *also* weakens the restoring force and raises $|k|$ for bending waves, whereas for density waves (where self-gravity is destabilizing) missing edge gravity and lattice discreteness act with opposite sign on the restoring force.
   - *Suggested fix:* Cite Eqn. (C1) ($X(\phi)$) and Appendix D in Section 4.1.4 (or write out ID B13 in Appendix D): in a discrete lattice of finite width, $2\pi G\sigma|k|$ in Eqn. (15) is replaced by $X(k\Delta) - 2G\sigma(1/d_{\rm in} + 1/d_{\rm out})$, which is smaller than $2\pi G\sigma|k|$ and therefore requires a slightly larger $|k|$ to match $D(x)$.

5. **Appendix C, Eqns. (C4)–(C5) and pages 7–8 (Resolution limit, viscous overstability, and the $m \simeq N_t/2$ sawtooth):**
   - *Problem:* Appendix C packs three distinct topics into one dense section: (i) the axisymmetric dispersion relation of the discrete streamline lattice with pressure and gravity (Eqns. C1–C4); (ii) the linear viscous overstability on the lattice (Eqn. C5); and (iii) an empirical investigation of an eccentric twisted ringlet that suffers a numerical $m \simeq N_t/2$ azimuthal sawtooth instability. Moreover, Eqn. (C5) is introduced without showing the viscous force equations it comes from, and the scaling $m_r \propto J_2 \Delta a^3 / \Delta\tilde{\omega}$ has a `%TODO: ... derive it properly` comment in the LaTeX source.
   - *Suggested fix:* If viscosity and pressure are deferred to a future paper (as Sections 3.2 and 3.3 state), keep only the cold gravitational lattice dispersion relation $X(\phi)$ and the Toomre/zigzag resolution criterion $\mathcal{R} < 2$ (Eqns. C1 and C4) here, and merge Appendix D into it (since vertical lattice gravity has the exact same sum $X(\phi)$, ID G5). Move the viscous overstability (Eqn. C5) and the viscous ringlet tests to the future viscosity/pressure paper, or give the viscous equations of motion (IDs O2, O4) explicitly if kept.

---

#### (2) Plain-language physics summaries of the derivations

1. **Vertical epicyclic motion (Section 2, Eqns. 1–12):**
   - *Plain-language summary:* An oblate planet's gravity pulls a slightly inclined ring particle back toward the equatorial plane with a linear restoring force, like a simple harmonic oscillator with vertical frequency $\nu_a$ slightly faster than the orbital frequency $\Omega_a$. Because the vertical restoring pull is slightly faster than the orbit, the particle crosses the ring plane a little early on each revolution, causing its orbital node to regress slowly while its horizontal and vertical motions remain decoupled to first order.
   - *Comparison to the paper:* The paper's presentation in Section 2 is already clear, brisk, and just as easy to follow as this summary.

2. **3D streamline wire gravity (Section 3.1 and Appendix A, Eqns. 13–14, A1–A8):**
   - *Plain-language summary:* Because neighboring streamlines in a narrow ring are much closer to a particle than the ring's radius of curvature, each perturbing streamline pulls on the particle like a straight, infinite wire of mass. The wire's pull points directly from the particle to the closest point on the wire—found by taking the radial and vertical offset to the wire at the particle's longitude and subtracting the component parallel to the wire's local 3D tilt—and that slight tilt gives the wire's gravity a small azimuthal component that transfers angular momentum between streamlines.
   - *Comparison to the paper:* The paper's derivation in Appendix A is harder to follow than this summary because it frames the geometry as a two-legged trip parameterized by an unknown longitude offset $\delta\theta'$ before using orthogonality, rather than projecting the crossing-point vector perpendicular to the wire's tangent vector directly, and because the crucial physical insight—that the wire's tilt is what allows streamlines to exchange angular momentum—is buried in Appendix B rather than stated in Section 3.1.

3. **2D planar limit and corrections to `epi_int_lite` (Appendix B, Eqns. B1–B2):**
   - *Plain-language summary:* When a streamline is perturbed by a multi-armed wave ($m \ge 2$), the wave pattern rotates at a different speed from the particles orbiting inside it, so a particle's velocity vector does not point along the instantaneous geometric shape of the streamline. Earlier versions of the code used the particle's velocity vector as the wire's tangent, which underestimated the wire's tilt—and therefore the gravitational torque between streamlines—by a factor of roughly $m$.
   - *Comparison to the paper:* The "Direction" paragraph in Appendix B explains this well, though burying it in Appendix B hides a major physical bug fix that is essential for simulating $m \ge 2$ spiral density and bending waves.

4. **Lattice dispersion relation and resolution limit (Appendix C, Eqns. C1–C5):**
   - *Plain-language summary:* Replacing a continuous disk with $N_r$ discrete streamlines spaced $\Delta$ apart softens self-gravity at short wavelengths and makes centered-difference pressure and viscosity completely blind to a grid-scale "zigzag" mode where alternating streamlines move in opposite directions. A cold or viscous ring simulation will stay free of grid-scale instabilities only if the streamline spacing $\Delta$ is kept coarser than the physical Toomre or viscous-overstability wavelength, which ironically means adding too many streamlines to a self-gravitating ringlet makes it go unstable.
   - *Comparison to the paper:* The paper's explanation of the zigzag mode (Eqns. C1–C4 and Fig. 8a) is very clear, but the second half of Appendix C (viscous overstability and the $m \simeq N_t/2$ sawtooth in twisted eccentric ringlets) jumps abruptly between linear axisymmetric theory and complex nonlinear numerical artifacts without enough setup to follow easily.

5. **Continuum bending-wave restoring force (Appendix D, Eqn. D1):**
   - *Plain-language summary:* When a flat sheet of streamlines is corrugated vertically into a wave of wavenumber $k$, each wire pulls neighboring wires back toward the mean plane with a vertical force proportional to their height difference divided by their radial separation squared. Integrating that $1/x^2$ pull across the whole sheet yields the familiar $-2\pi G\sigma|k|Z$ vertical restoring force of a continuous self-gravitating bending wave.
   - *Comparison to the paper:* Appendix D is only one paragraph and one equation, and is completely transparent.

6. **Forced and free spiral bending and density waves (Section 4, Eqns. 15–38):**
   - *Plain-language summary:* An inclined (or eccentric) moon pulls on the ring with a rotating $m$-armed spiral pattern that matches the ring particles' natural vertical (or radial) oscillation frequency at a resonant radius, launching a spiral wave whose wavelength shrinks as it propagates away from the resonance at the group velocity $c_g \propto G\sigma$. Far from resonance, the ring's motion is the sum of a non-propagating forced oscillation (the response massless particles would have) and the propagating self-gravitating wave that carries away the torque exerted by the moon.
   - *Comparison to the paper:* Sections 4.1.1–4.1.2 and 4.2.1–4.2.2 are concise and well-structured, following Shu (1984) closely while translating the equations into the paper's sign convention.

7. **Massless limit of the bending wave (Appendix E, Eqns. E1–E3):**
   - *Plain-language summary:* In steady state, lowering the ring's surface density $\sigma$ toward zero paradoxically makes the bending wave's amplitude grow as $\sigma^{-1/2}$, because a lighter ring propagates wave energy away more slowly ($c_g \propto \sigma$) and has to pile up a taller wave to carry away the satellite's torque. However, at any finite time $t$, the distance the wave has traveled from resonance ($|c_g|t \propto \sigma t$) shrinks to zero as $\sigma \rightarrow 0$, leaving behind only the massless forced response $L_V/x$ everywhere outside the resonance.
   - *Comparison to the paper:* Appendix E explains both the math (integration by parts of the Fresnel integral) and the physical resolution of the steady-state paradox admirably well.

---

#### (3) Organization of the paper

The overall progression of the main text (Section 1: Intro $\rightarrow$ Section 2: Vertical drift $\rightarrow$ Section 3: Ring forces $\rightarrow$ Section 4: Spiral wave tests) is logical and matches `hahn13` and `hahn25`. However, a few targeted structural changes would make the paper read much more cohesively:

1. **Swap Section 4.1 (Bending waves) and Section 4.2 (Density waves), or explain why Density Waves are in the paper at all:**
   - *Reason:* Right now, Appendix B reveals a crucial fix to the 2D in-plane gravity of `epi_int_lite` (using the geometric slope $s_r$ instead of the velocity ratio $v_r'/v_\theta'$, without which $m \ge 2$ Lindblad waves have the wrong torque by a factor of $m$), and Section 4.2 tests `epi_int_3d` on in-plane spiral density waves. Yet the title, abstract, and introduction frame the paper solely around vertical motion and spiral bending waves, without mentioning spiral density waves or the 2D slope fix. Either update the Title/Abstract/Intro to say the paper extends `epi_int_lite` to spiral density and bending waves (and present the 2D density-wave validation or the 2D slope fix explicitly), or keep Bending Waves first as the flagship result and introduce Section 4.2 as verifying that the 3D wire-gravity formulation also fixes `epi_int_lite`'s in-plane gravity for $m \ge 2$ Lindblad density waves (referencing Appendix B).
2. **Move Appendix D (Restoring force of a bending wave) immediately after Appendix B (or into Appendix C):**
   - *Reason:* Appendices A, B, and D are all short derivations of the wire gravity's properties (A: 3D wire geometry; B: 2D limit and `epi_int_lite` comparison; D: continuum integral of the vertical wire pull for a bending wave). Right now, the long and tangential Appendix C (Resolution limit and viscous overstability of eccentric ringlets) sits between Appendix B and Appendix D, breaking the flow from wire gravity to bending waves. Placing D right after B (so the order is: A. 3D Wire Geometry, B. Relation to 2D `epi_int_lite`, C. Continuum Restoring Force of a Bending Wave, D. Massless Limit of the Bending Wave, E. Resolution Limit of the Streamline Method) groups the appendices in the exact order they are cited in Sections 3 and 4.
3. **Trim or defer the unfinished parts (Sections 3.2, 3.3, and the second half of Appendix C):**
   - *Reason:* Sections 3.2 (Pressure) and 3.3 (Viscosity) currently say "will be completed later," while Appendix C analyzes pressure and viscous overstability in 2D eccentric ringlets using `epi_int_lite`'s equations. If this paper submits with only self-gravitating collisionless bending and density waves, fold Sections 3.2 and 3.3 into a single sentence at the start of Section 3 ("Pressure and viscosity, which are not used in the wave simulations of Section 4, are unchanged from `epi_int_lite` [Hahn et al. 2025]"), and save the viscous overstability and eccentric-ringlet sawtooth study in Appendix C for the follow-up paper where 3D viscosity and pressure are actually implemented and tested.

---

### Finding 2
- **Type:** error
- **Severity:** high
- **Where:** `paper/epi_int_3d.tex`, Section 4.2.4, lines 449–454: `"Over x < 0.01, outside the resonance zone, the wave's amplitude is within \dwAmpPercent\% of |Z_w|..."` and `"it is \dwKcontPercent\% larger than Eqn. (29) gives over x < 0.01"`
- **Finding:** At an inner Lindblad resonance, the spiral density wave propagates *outward* into the region $x > 0$ (as Section 4.2.1 states: "which lies outside the resonance, so $c_g > 0$ and $k > 0$... it travels outward", and Section 4.2.4 states: "The wave reaches $x = 0.01$ after the time $\tau_{0.01}$... when its front is at $x = \dwXfront$, short of the ring's outer edge at $x = \dwXout$"). Writing "Over $x < 0.01$, outside the resonance zone" twice in the second paragraph of Section 4.2.4 reverses the inequality sign (copy-pasted from `Over -0.01 < x` in Section 4.1.4): $x < 0.01$ includes the entire inner non-wave side $x < 0$ and the resonance zone $x \approx 0$, whereas the steady wave zone between the resonance zone and the reference point $x = 0.01$ is $0 < x < 0.01$ (or $x_{\rm res} < x < 0.01$).
- **Suggestion:** Change both occurrences of `Over $x < 0.01$, outside the resonance zone` in Section 4.2.4 to `Over $x < 0.01$ outside the resonance zone` $\rightarrow$ `Between the resonance zone and $x = 0.01$` (or `Over $0 < x < 0.01$, outside the resonance zone`).

---

### Finding 3
- **Type:** error
- **Severity:** medium
- **Where:** `paper/epi_int_3d.tex`, Section 4.2.2, line 400: `"As in Section 4.1.2, \code{epi\_int\_3d} keeps only this term, with $b_{1/2}^{(m)}$ replaced by its tangent line at $\beta = r_L/a_s$, and pulls each particle with $-\nabla\Phi_m$."` and Eqn. (34)
- **Finding:** If the code pulls each particle with $-\nabla\Phi_m$ in cylindrical coordinates where $\Phi_m(r,\theta,t) = -(GM_s/a_s)b_{1/2}^{(m)}(\beta)\cos m(\theta - \Omega_p t)$, then $A_r = -\partial\Phi_m/\partial r$ and $A_\theta = -(1/r)\partial\Phi_m/\partial\theta$. However, Section 4.1.2 did *not* replace $b_{3/2}^{(m)}(\beta)$ by its tangent line at $\beta = r_V/a_s$ (it kept $A_V$ constant at $\beta = r_V/a_s$), and more importantly, if $b_{1/2}^{(m)}(\beta)$ is replaced by its tangent line at $\beta_L = r_L/a_s$, then in $A_\theta = -(1/r)\partial\Phi_m/\partial\theta$, $b_{1/2}^{(m)}(\beta)$ varies linearly with $r$ across the ring rather than being evaluated at $\beta_L$. That is fine and accurate, but saying "As in Section 4.1.2" is slightly misleading because Section 4.1.2 evaluated the force amplitude at $\beta = r_V/a_s$ directly rather than linearizing the potential in $\beta$.
- **Suggestion:** Clarify to: `"As in Section \ref{sec:waves_theory}, \code{epi\_int\_3d} keeps only the resonant $m$th term, with $b_{1/2}^{(m)}(\beta)$ expanded to first order in $\beta - r_L/a_s$ so that its radial derivative in $-\nabla\Phi_m$ is constant across the narrow ring."`

---

### Finding 4
- **Type:** overstated
- **Severity:** medium
- **Where:** `paper/epi_int_3d.tex`, Section 4.1.3, lines 307–310: `"are in excellent agreement with Eqn. (28), to within \vresAgreePercent\%, and the differences are of order $x$, which Eqn. (28) neglects."` and `notes/equations.md` IDs B12 and B20
- **Finding:** According to `notes/equations.md` (IDs B12 and B20), the discrepancy of up to ~3.6% (or 1–4%, 13% beside $r_V$) between the simulated inclinations and Eqn. (28) comes from *three* distinct sources, not just $\mathcal{O}(x)$ terms in $D(x)$: (1) right beside $r_V$, the un-damped free oscillation remains at ~13% because the turn-on timescale $\tau$ is not adiabatic when $|\nu_a - \omega_m|\tau \lesssim 1$ (ID B11, B12); (2) the code's diagnostic `vertical_elements` computes $i = [(z/a)^2 + (v_z/a\nu_a)^2]^{1/2}$, whereas the forced motion oscillates at frequency $\omega_m \ne \nu_a$, introducing a systematic elliptic-integral factor $(2/\pi)E(1 - \omega_m^2/\nu_a^2) \simeq 1 - D/(4\nu_a^2)$ of order $\mathcal{O}(x)$ in the time-averaged inclination even when $D(x)$ is exact (ID B20); and (3) replacing $D(x)$ by $\mathcal{D}x$ and dropping $J_2$ in Eqn. (28).
- **Suggestion:** State briefly that the two streamlines straddling $r_V$ retain a small free oscillation because the satellite's turn-on time $\tau$ is shorter than the beat period $2\pi/|\nu_a - \omega_m|$ very close to resonance, while away from $r_V$ the $\mathcal{O}(x)$ difference arises because Eqn. (28) linearizes $D(x)$ and because Eqn. (12a) divides $v_z$ by $\nu_a$ rather than the forcing frequency $\omega_m$.

---

### Finding 5
- **Type:** unclear
- **Severity:** medium
- **Where:** `paper/epi_int_3d.tex`, Section 4.1.4, lines 354–356: `"On the nonwave side of the resonance, the forced response departs from linear theory's by \bwForcedMedianPercent\% typically (the median) and by up to \bwForcedPercent\%. That theory is for an infinite sheet, which has ring matter beyond this ring's outer edge, and that is a likely cause, though it has not been quantified."`
- **Finding:** On the nonwave side ($x > 0$) of the bending wave, the ring's equilibrium has no unperturbed vertical force from the missing exterior ring matter (since an uninclined ring lies at $z = 0$), and in the WKB regime the forced response $h \simeq L_V/x [1 + \mathcal{O}(\epsilon/x^2)]$ (Eqn. E2) is independent of the ring's self-gravity to leading order. A more direct cause of discrepancy on the nonwave side in a run where the satellite is switched on instantaneously at full mass at $t = 0$ is the undamped free vertical oscillation that was excited at $t = 0$ (which propagates inward at speed $|c_g|$, so on the nonwave side $x > 0$ the free oscillation excited at $x > 0$ travels toward $r_V$ and takes time to clear the outer edge), together with the linear approximation $D \simeq \mathcal{D}x$ in Eqn. (22).
- **Suggestion:** Check whether the discrepancy on the nonwave side is simply the $\mathcal{O}(x)$ linearization of $D(x)$ plus the initial free oscillation or edge truncation of the integral in Eqn. (D1), and state both possibilities rather than attributing it solely to missing exterior matter.

---

### Finding 6
- **Type:** missing
- **Severity:** medium
- **Where:** `paper/epi_int_3d.tex`, Section 3.1.2, line 183 (`%TODO: lambda = m_1 Omega_a/(2 pi v) assumes a stationary streamline`), and Appendix C, line 577 (`%TODO: this scaling is an order-of-magnitude estimate that fits the runs below; derive it properly`)
- **Finding:** Two unresolved `%TODO` comments remain in the LaTeX source, and the first one points to a genuine physical approximation that is never stated in the text: the code does not actually interpolate the longitudinal spacing of the particles to measure the instantaneous linear mass density $\lambda'(\theta')$ of a compressed or sheared streamline; rather, it computes $\lambda'$ from the continuity relation $\lambda = m_{\rm ring}\Omega_a / (2\pi N_r v)$ for a *stationary* streamline (which, as Appendix B "Direction" notes for velocities, is not exact when an $m \ge 2$ wave pattern rotates at $\Omega_p \ne \Omega_a$, though for linear tightly-wound waves $k r \gg m$ radial compression of streamline spacing dominates over azimuthal compression along a streamline).
- **Suggestion:** State explicitly in Section 3.1.2 how $\lambda'$ is computed for each particle before interpolation, and note why variations of $\lambda'$ along a streamline are negligible for tightly wound spiral waves ($|k|r \gg m$). In Appendix C, either derive the secular self-gravity precession scaling from Eqn. (C1) of Hahn et al. (2025) or remove the `%TODO`.

---

### Finding 7
- **Type:** unclear
- **Severity:** medium
- **Where:** `paper/epi_int_3d.tex`, Section 4.1.4, lines 379–384: `"The satellite pulls only vertically, so the ring does not feel this torque, and its total angular momentum $L_z$ can change only through the ring's gravity and the integrator's round-off error."`
- **Finding:** This sentence will puzzle many readers: how can the satellite exert a nonzero torque $T \simeq T_V$ on the bending wave (Fig. 4) while "the satellite pulls only vertically, so the ring does not feel this torque" and $L_z$ is conserved (Fig. 5)? The physical resolution is an artifact of the simplified forcing in Eqn. (21): a real inclined satellite has a 3D disturbing potential $\Phi_s(r,\theta,z,t) \simeq -A_V z \sin m(\theta - \Omega_p t)$, which exerts both a vertical force $A_z = -\partial\Phi_s/\partial z = -A_V \sin m(\theta - \Omega_p t)$ *and* a second-order azimuthal force $A_\theta = -(1/r)\partial\Phi_s/\partial\theta = -(m/r)A_V z \cos m(\theta - \Omega_p t)$ proportional to $i_s z \sim \mathcal{O}(i_s i)$. Because `epi_int_3d`'s vertical-resonance mode applies only $A_z$ and omits $A_\theta = -(1/r)\partial\Phi_s/\partial\theta$, the satellite injects vertical energy without changing the particles' $h_z = r v_\theta$, whereas in Section 4.1.4 $T = \sum \mu r A_\theta = -m A_V \sum \mu z \cos m(\theta - \Omega_p t)$ is evaluated post-hoc as the torque that $-\partial\Phi_s/\partial\theta$ *would* exert. Meanwhile, the wave still carries negative $z$-angular momentum flux inward because the ring's *own* 3D wire gravity includes the $A_{g,\theta} \propto -d_z s_z$ term (Appendix B), which transfers $L_z$ from the inner wave region back to the resonance zone, conserving total $L_z$.
- **Suggestion:** Explain this explicitly in two sentences: "Because the satellite's resonant potential $\Phi_s = A_V z \sin m(\theta - \Omega_p t)$ has an azimuthal gradient $-\partial\Phi_s/\partial\theta$ of second order, $\mathcal{O}(i_s i)$, \code{epi\_int\_3d} applies only its first-order vertical pull $A_z$, Eqn. (21), and evaluates the torque $T = -\sum \mu (\partial\Phi_s/\partial\theta)$ as a diagnostic. Consequently, the satellite applies no azimuthal kick in the simulation, and the ring's total vertical angular momentum $L_z$ is conserved up to the ring's internal forces and numerical error (Fig. 5), while the ring's own tangential gravity $A_{g,\theta} \propto -d_z s_z$ transports $L_z$ between streamlines."

---

### Finding 8
- **Type:** style
- **Severity:** low
- **Where:** Multiple locations in `paper/epi_int_3d.tex`
- **Finding:** Several minor breaches of the author's style guides (`joes_papers/writing_style.md` and `style/dataviz.md`):
  1. **Rule P62 / P63 (Every symbol in a figure is used in the text; check signs):** In Appendix B, line 500, the text writes `and $\cos\varphi = v_\theta'/v'$ (Fig. \ref{fig:wire_velocity})`, whereas the caption of Fig. 7 does not mention $\cos\varphi$, and in Fig. 6 (`wire_planar.pdf`) the caption defines $\mathbf{r}$ and $\mathbf{r}'$ as Positions from the planet ("the planet is far off to the left").
  2. **Rule P24 (Equation references):** In the captions of Figs. 6 and 8 (`paper/epi_int_3d.tex`, lines 469, 472, 608, 610), equation numbers are written without parentheses after `Eqn.` (`(Eqn. \ref{eqn:wire_d})`, `(Eqns. \ref{eqn:wire_2d})`, `(Eqn. \ref{eqn:zigzag_X}, solid red)`), whereas everywhere else in the paper and captions (Figs. 1–5, 7) they are written with parentheses: `Eqn. (\ref{eqn:...})`.
  3. **Rule P54 ("leading" and "trailing" in quotes at first use):** Section 4.1.1 (line 264) writes `and the wave is trailing \citep{S84}` without quotes around `'trailing'` at first use.
  4. **Rule V2 and V3 (One panel per figure; `\columnwidth` width):** Figures 1, 7 (`density_torque.pdf`), and 8 (`zigzag.pdf`) have two panels `(a)` and `(b)`, and Figs. 6, 7 (`wire_velocity.pdf`), and 8 use `\textwidth` (`0.5\textwidth`, `0.3\textwidth`, `0.9\textwidth` in `figure*`) rather than `\columnwidth` (though `vertical_resonance` and `zigzag` are noted in `dataviz.md` as deferred).
  5. **Rule P5 (Code's name starts a sentence capitalized):** Line 189 (`\code{Epi\_int\_lite}`), line 196 (`\code{Epi\_int\_3d}`), line 202 (`\code{Epi\_int\_3d}`), line 284 (`\code{Epi\_int\_3d}`), line 490 (`\code{Epi\_int\_lite}`), line 507 (`\code{Epi\_int\_3d}`), line 509 (`\code{Epi\_int\_lite}`), line 521 (`\code{Epi\_int\_3d}`) capitalize `\code{Epi\_int\_...}` correctly at the start of a sentence, except line 284 (`\code{Epi\_int\_3d} keeps only the resonant term`) is in the middle of a sentence after `so` and should be lowercase `\code{epi\_int\_3d}`.
  6. **Comment on line 647:** `\bibliographystyle{aasjournalv7} % v7.1 (June 2026) mislabels lone citations as "2006a"` and comments with dates in late 2026 (`2026-10-18`, `2026-11-12`): fine for internal notes, though clean up `%TODO` and draft comments before journal submission.
- **Suggestion:** Change `\code{Epi\_int\_3d}` to `\code{epi\_int\_3d}` on line 284, add parentheses around equation numbers in the captions of Figs. 6 and 8 (`Eqn.\ (\ref{...})`), and put `'trailing'` in quotes at first use in Section 4.1.1.

---

### Finding 9
- **Type:** unclear
- **Severity:** low
- **Where:** `paper/epi_int_3d.tex`, Appendix D, line 619: `"just as the same sum gives the $-2\pi G\sigma|k|$ term of a density wave in the plane."`
- **Finding:** In the dispersion relation as written in this paper—Eqn. (15) ($D = \nu_a^2 - \omega_m^2 = -2\pi G\sigma|k|$) and Eqn. (29) ($D = \kappa_a^2 - \omega_m^2 = +2\pi G\sigma|k|$)—the radial acceleration of a planar density wave $\xi_r = \xi_0 \cos kr$ is $A_{g,r} = +2\pi G\sigma|k|\xi_r$ (destabilizing, reducing $\omega_m^2 = \kappa_a^2 - 2\pi G\sigma|k|$), whereas for a bending wave $A_{g,z} = -2\pi G\sigma|k|Z$ (restoring, increasing $\omega_m^2 = \nu_a^2 + 2\pi G\sigma|k|$). Saying "gives the $-2\pi G\sigma|k|$ term of a density wave" blurs the opposite sign in Eqn. (29) vs. Eqn. (15).
- **Suggestion:** Change to: `"just as the same sum in the plane gives the $+2\pi G\sigma|k|\xi_r$ acceleration (and the $2\pi G\sigma|k|$ term of Eqn.\ (\ref{eqn:dw_disp})) of a spiral density wave."`

---

### Strongest parts

1. **Exact, reversible 3D drift step (Section 2):** The extension of `epi_int_lite`'s reversible epicyclic transformation to 3D via $(z, v_z/\nu_a) \leftrightarrow (i, \theta - \Omega)$ is clean, physically transparent, and introduces zero secular drift in the unperturbed vertical elements while keeping the horizontal and vertical drift steps decoupled.
2. **Identification and correction of the 2D wire-tangent bug for $m \ge 2$ waves (Appendix B and Section 3.1):** Recognizing that the particle velocity vector $\mathbf{v}'$ is not tangent to the streamline when the wave pattern rotates at $\Omega_p \ne \Omega_a$—and that using $\mathbf{v}'$ underestimates inter-streamline gravitational torques by a factor of $\sim m$ near an $m$-armed resonance—is an important physical and numerical insight that unlocks spiral wave simulations in streamline N-body codes.
3. **Comprehensive, quantitative wave benchmarks (Section 4 and Appendix E):** Testing both massless forced oscillations and self-gravitating propagating waves (amplitude, spatial waveform, wavenumber convergence with $N_r$, satellite torque, and conservation of $L_z$) for both vertical bending waves and in-plane Lindblad density waves provides a thorough validation of the code against Shu (1984) and Goldreich & Tremaine (1979), and Appendix E's resolution of the $\sigma \rightarrow 0$ limit is particularly sharp.

### Weakest parts

1. **Unfinished sections and scope mismatch (Abstract, Intro, Sections 3.2–3.3, Section 4.2):** The draft contains placeholder sections for 3D pressure and vertical viscosity ("will be completed later"), while including full spiral density wave tests (Section 4.2) and a 2D viscous overstability / numerical instability appendix (Appendix C) that are not motivated in the Title or Introduction.
2. **Appendix C's transition from clean lattice theory to unresolved numerical instability:** The first half of Appendix C (the discrete gravitational kernel $X(\phi)$ and zigzag blindness of centered differences) is a gem, but the second half introduces an unsolved $m \simeq N_t/2$ azimuthal sawtooth instability in twisted ringlets whose dependence on $N_t$ is non-monotonic and left unexplained.
3. **Verbose vector algebra in Appendix A:** The two-leg derivation of the perpendicular distance to a straight wire via $\delta\theta'$ takes much longer than projecting $\mathbf{d}_0$ orthogonal to $\hat{\mathbf{t}}$ directly.

### What a skeptical referee would push back on

1. **Why total ring angular momentum $L_z$ is conserved when a satellite exerts a torque $T_V$ on the ring (Section 4.1.4, Figs. 4 and 5):** A referee will immediately ask why the satellite exerts a steady torque $T \approx T_V$ in Fig. 4 while Fig. 5 shows $\Delta L_z = 0$. As noted in Finding 7, you must explain that the simulation applies only the first-order vertical force $A_z = -\partial\Phi_s/\partial z$ and omits the second-order azimuthal force $A_\theta = -(1/r)\partial\Phi_s/\partial\theta \sim \mathcal{O}(i_s i)$, while measuring $T$ diagnostically from the wave's phase and amplitude.
2. **The unresolved $m \simeq N_t/2$ azimuthal sawtooth instability in Appendix C:** If you tell a referee that a circular, non-gravitating, inviscid ringlet goes numerically unstable and crosses streamlines within 850 orbits when $N_t = 40$ or $320$ (because of 3-point Lagrange interpolation in longitude), the referee will ask why the bending and density wave simulations in Section 4 (with $N_r = 200$–$800$) can be trusted not to develop the same interpolation instability over longer timescales, and what $\mathcal{R}$ values the runs in Section 4 have. You should report the value of $\mathcal{R}$ for the bending and density wave runs in Section 4 (since $|k|\Delta$ is small, $\mathcal{R}$ is related to the shortest resolved wavelength at the inner/outer edge) and note why those runs remain completely smooth.
3. **Streamline linear density $\lambda'$ in non-stationary $m \ge 2$ waves (Section 3.1.2):** Because Appendix B makes a strong point that streamline shapes in $m \ge 2$ waves are not stationary in the inertial frame (so $\mathbf{v}'$ is not parallel to the streamline), a sharp referee will ask how $\lambda'(\theta)$ is computed when $m \ge 2$. If $\lambda'$ is still computed from $\lambda \propto 1/v$ as if the streamline were a stationary tube of constant mass flux, that is an $m=1$ secular approximation; you should state clearly in Section 3.1.2 that for tightly wound waves ($|k|r \gg m$), surface density variations come from the radial convergence of adjacent streamlines ($\Delta r$) rather than longitudinal variations of $\lambda'$ along a single streamline.
4. **Planet motion / indirect potential (Section 3.4):** Section 3.4 mentions that "Before and after each drift, small shifts of the particles' coordinates account for the planet's motion about the system's center of mass (Chambers 1999)." A referee will ask whether the planet moves vertically in response to the inclined ring or the analytic $m$-th harmonic satellite forcing (Eqn. 21), or whether the center-of-mass shift is zero/disabled for the analytic $m \ge 2$ Fourier satellite modes (as `notes/equations.md` ID B2 and L1 note, the indirect term only exists for $m = 0$ vertically and $m = 1$ horizontally). State in Section 3.4 or 4.1.2 that $m \ge 2$ modes do not displace the planet's center of mass.
