# Derivation style

How the author works out the physics: how a problem is set up, how approximations are made and checked,
and how a result is interpreted. Use it for every derivation in `notes/` and the paper's appendices.
Rules are numbered (D1, D2, …) so that reviews and the ledger can cite them. How equations sit in the
prose (introduced by what they are, symbols defined in a "where" clause, "Equation (n)" versus "Eqn. (n)")
is in `joes_papers/writing_style.md` ("Equations"); this guide is about the math itself.

**Sources.** Hahn (2007, ApJ 665, 856; "hahn07") §2, the bending-wave amplitude and dispersion relation
derived in the main text; Hahn & Spitale (2013; "hahn13") Appendix A, the viscous angular momentum flux;
Hahn et al. (2025; "hahn25") Appendices A–C, the epicyclic orbit elements, the gravitational angular
momentum flux and the checks against theory; with quotes from hahn08–hahn10 already gathered in
`joes_papers/writing_style.md`.

**Status.** Approved by the author 2026-11-08.

**D0. Derivations in `notes/` are done in this style from the start** (the author, 2026-11-08), not
written as working notes and restyled later. A derivation in `notes/` should read like a draft of the
paper's appendix, so putting it in the paper is mostly trimming. Derivations written before this guide
(e.g. `notes/phase2/streamline_gravity.md`) stay as they are.

## Setting up the problem

- **D1. Start from the physical system and say what is kept and what is ignored, and why.** "The disk is
  flat and thin and has a vertical half-width h and constant volume density ρ… The disk is assumed to be
  viscous, and its gravity is ignored here since this appendix is only interested in the angular momentum
  flux that is transported solely by viscosity" (hahn13 App. A). "For the moment we will assume that all
  rings are circular, noting that we will deal with the system's eccentricity evolution in a follow-up
  study" (hahn07).
- **D2. Split the problem into its physical parts and treat each in turn.** hahn07 writes the equation of
  motion as a sum, İ = İ|disk + İ|sat and ϖ̇ = ϖ̇|disk + ϖ̇|sat + ϖ̇|obl, gives each part its own subsection
  (Ring–Disk Evolution, Ring–Satellite Evolution), then sums them. Say up front which equation will be
  used for what: "Section 2.1 will use the first equation to solve for the wave amplitude… Section 2.2 will
  then use the other equation to solve for… the dispersion relation."
- **D3. Use the condition that defines the state being sought.** A steady wave: "The amplitude of a spiral
  bending wave that is in a steady state does not vary with time, so the disk inclinations obey İ(a) = 0";
  "a persistent spiral pattern must also rotate with a constant angular velocity" (hahn07). The
  solution follows from balancing terms: "When the wave is in steady state, the two I-excitation rates…
  are balanced, which yields the amplitude."
- **D4. Build from a standard starting point, cited by equation number**: the disturbing function
  (Brouwer & Clemence 1961; Murray & Dermott 1999), Euler's equation and the momentum flux density tensor
  (Landau & Lifshitz 1987), the epicyclic solution (hahn25 Eq. A9). For this project, Goldreich & Tremaine
  notation throughout; results taken from Shu are translated first (`CLAUDE.md`, "Units and notation").
- **D4a. Equations of motion and fluxes first.** Work from the equations of motion, forces and fluxes
  (hahn13, hahn25) by default; use orbit elements and the disturbing function (hahn07–hahn10) when that
  is more convenient, e.g. for secular precession rates, and say why (the author's preference, 2026-11-08).
- **D5. Physical quantities stay physical.** G, M, Ω, κ appear explicitly in derivations, even though the
  code sets G = M = 1. Quantities are made dimensionless by named ratios with a physical meaning:
  μ_d = πσa²/M "the ring's so-called normalized disk mass", x = (a′ − a)/a "the fractional distance",
  z = |k|ax "the distance from the ring edge in units of 2π wavelengths."

## Approximating

- **D6. Every approximation comes with its physical reason, in the same sentence, before it is made.**
  "As one might expect, this integral is dominated by the contributions from nearby annuli… Because of the
  steep dependence of b̃ on x′, we can replace the inclination I′(x′) and disk mass μ′_d(x′) with their values
  evaluated at the perturbed ring… and also pull them out of the integral" (hahn07). The pattern is
  *physical reason → "we can" or "which allows us to" → the step*.
- **D7. Name the small quantity and the order.** "a Taylor expansion of equation (30) in the small quantity
  x + Δ" (hahn07); "to first order in the small quantities Δa, Δe, and Δϖ̃ and assuming e₀ ≪ 1" (hahn25);
  "These equations satisfy the equation of motion (A1) to accuracy O(e) and with errors of order O(e²)"
  (hahn25).
- **D8. State the assumptions that need checking, and say where they will be checked.** "considerable
  progress can be made if we assume that k is constant over a wavelength… Also keep in mind that these
  results assumed that the wavenumber k varies little over a single wavelength; §2.2.2 will note when this
  approximation breaks down" (hahn07). The check comes later, in its own section or against a simulation
  ("Checking the Assumptions", hahn08).
- **D9. Use limits to simplify, and say which limit.** "In the |x′| ≪ 1 limit, the softened Laplace
  coefficient is…"; "We will be interested in a disk whose vertical thickness is small compared to the
  wavelength, so H ≪ 1, and A_H(z) ≃ sin(z)/z − Ci(z)" (hahn07).

## Results

- **D10. Collect the hard part into a dimensionless function**, give it a name and a one-line meaning,
  express it with standard functions where possible, and plot it. "the function A_H(z) is a dimensionless
  measure of the warped disk's perturbation of itself"; A_H and B_H in terms of Si and Ci (Abramowitz &
  Stegun 1972), all plotted in hahn07 Fig. 2.
- **D11. Check every result in its limits.** Downstream and at the edge: "For z ≫ 1, Ci(z) ≃ sin(z)/z −
  cos(z)/z² + O(z⁻³)… so A_H(z) ≃ cos(z)/z² downstream"; "B_H takes values of ½ ≤ B_H(z) ≤ 1 for z ≥ 0,
  with B_H(0) = ½ at the disk's inner edge"; thin and thick disks: "B_H^∞ ≃ 1. But if the disk is thick,
  H ≫ 1 and B_H^∞ ≃ H⁻¹, which indicates that the disk's ability to sustain a bending wave is weakened."
  A massless or planar limit that reproduces a known result is a check worth stating.
- **D12. Use physical requirements to fix what the math leaves open.** "So if I(z) is to be a finite
  constant, then the longitude offset must be φ_o = ±π/2… Of course, these inclinations must also be
  positive, so sin φ_o = ±1 = −s_k" (hahn07).
- **D13. Interpret each term of the result.** "The first term is the familiar differential precession that
  would occur if the disk were massless. The second term, which is proportional to the disk mass μ_d, is
  the additional precession that is due to the torque that the satellite exerts on the disk's spiral
  pattern"; "All terms in equation (32) are negative, so the disk precesses in a retrograde sense"; "The
  dispersion relation… thus tells us how the wavenumber |k(x)| must adjust throughout the disk" (hahn07).
- **D14. Simplify using exact identities and say which.** "but note that ∇x₂·Π₁ − Π₂·∇x₁ = Π₂₁ − Π₁₂ =
  σ₁₂ − σ₂₁ = 0 since the viscous stress tensor is symmetric (Equation (A11))" (hahn13).
- **D15. A result that the code uses is also written as an algorithm**: numbered steps in the order the
  code runs them, with where each formula comes from and what each step's accuracy is (hahn25 App. A,
  "1. Calculate the particle's specific angular momentum h = rv_θ and then solve… 2. Calculate Ω and κ
  using Equations (A7) and (A8)…"), ending with why it is accurate: "because steps 1–4 convert coordinates
  to elements exactly, epi_int_lite's drift step is not a significant source of numerical error."

## Checking (this project's rules)

- **D16. Every derived equation gets a sympy check** (`notes/sympy/`) and a row in `notes/equations.md`
  before it is used in code; the author hand-verifies it later (`CLAUDE.md`, "Rules"). The sympy script
  follows `style/code.md` and checks the steps in the order the derivation makes them, so a failure
  points to a step.
- **D17. Outside review is a blind re-derivation**: another model gets the problem as set up under D1–D3,
  not the answer, and the two are compared (`PLAN.md`, Phase 3b).
- **D18. A derivation is finished when it has been compared with a simulation** that makes none of its
  approximations (`style/testing.md`, T12), and any disagreement is traced to a named assumption (D8).

## One result, one derivation

- **D19. Derive one result, the one that will be used.** Variants (exact versus first order, two choices of
  tangent) go in a short note on why the chosen one was chosen, not side by side as G9–G12 are in
  `notes/phase2/streamline_gravity.md`.

## Claude's habits to avoid

From the derivation notes Claude wrote for this project (e.g. `notes/phase2/streamline_gravity.md`).

| Claude's habit | Instead |
|---|---|
| The result first, in bold, then its justification | setup, approximations with reasons, then the result (D1–D9) |
| Telegraphic fragments: "(1) Raising particle and wire together by z changes nothing (\|z\| ≪ r)." | full sentences with the physical reason first (D6) |
| Several versions of one result side by side | the one that is used, and why (D19) |
| Approximations made silently inside an equation (≃ with no order or small quantity stated) | the small quantity and order named (D7) |
| Code-unit shortcuts in the math (G = 1) | G, M, Ω kept explicit; named dimensionless ratios (D5) |
| No limits checked, terms not interpreted | limits (D11) and each term's meaning (D13) |

## Open questions for the author

*(None open. Settled 2026-11-08: the method (D4a); this guide applies to derivations in `notes/` from
the start (D0).)*

## From the author's edits

*(Each correction the author makes to a derivation becomes a rule here, so it isn't needed twice.)*
