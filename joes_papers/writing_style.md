# Writing style for the paper

The author's style, taken from his papers in this folder, plus the author's edits to this paper. Use it
when drafting text for `paper/epi_int_3d.tex`. Drafted from Hahn & Spitale (2013, ApJ 772, 122; "hahn13")
and Hahn et al. (2025, ApJ 986, 220; "hahn25") and approved by the author 2026-10-25. Extended 2026-10-31
from Hahn (2007, ApJ 665, 856; "hahn07"), Hahn (2008, ApJ 680, 1569; "hahn08"), Hahn, Spitale & Porco
(2009, ApJ 699, 686; "hahn09") and Hahn (2010, ApJ 719, 1699; "hahn10"). Ward & Hahn (2003, AJ 125, 3389)
is mostly Ward's prose and is not used.

Where the author's habits have changed since 2007, follow the recent papers (hahn13, hahn25). A habit
seen in only one paper is marked as such.

Rules are numbered (P1, P2, …) so that reviews and the ledger can cite them, as in the other guides in
`style/`. The numbers are permanent: a new rule takes the next free number, or a lettered one (P9a) to sit beside
a related rule, and a dropped rule's number is not reused. The same holds in every guide.

## Voice and tense

- **P1.** Mostly impersonal: "The following describes…", "This section…", "the code identifies…", "Here, …".
  "We" appears sparingly, for the authors' own judgments and choices: "we conclude", "we speculate",
  "we multiply time t by…". Not "I", and not "one" except in "allows one to…".
- **P2.** Derivation steps may be imperatives: "Begin by considering a planetary ring…" (hahn07), "Next,
  Fourier-expand the longitude-dependent factor…" (hahn09), "Now derive the threshold…" (hahn10),
  "Begin with…" (hahn25).
- **P3.** Present tense for what the code does and what figures and simulations show ("Figure 3 shows that…",
  "simulations show that…"). Past tense only for earlier work or completed runs ("was used to
  simulate"). Future tense for pointers ahead: "Section 2.2.3 will quantify when that assumption breaks
  down" (hahn08); "it will be shown below that this contribution is unimportant" (hahn07, hahn08).
- **P4.** "is used to" for a code, model or equation doing work: "The rings model is used to simulate…"
  (hahn07, hahn08); "Equations (31) and (32) are used to calculate…" (hahn10).
- **P5.** The code's name starts a sentence capitalized: "Epi_int_lite is a child of the epi_int N-body
  integrator…".
- **P6.** Code names are set in typewriter, as in hahn13 (hahn25 used roman): write `\code{epi\_int\_3d}`,
  a macro in the paper's preamble that can be switched back to roman in one line (the author's
  choice, 2026-10-26).

## Sentences and connectives

- **P7.** Long, chained sentences that carry the reasoning: "Because a ringlet is narrow, all particles are in
  close proximity to the nearby portions of all streamlines, which allows us to approximate a
  streamline as an infinitely long wire…". Physical reason first, then the approximation it licenses.
- **P8.** Characteristic openers, the most frequent first: "However,", "Consequently,", "Note that", "Note
  though that", "Also note that", "Evidently,", "For instance,", "Keep in mind (though) that", "Also keep
  in mind that", "Of course,", "The main point is that", "So,", "So if…", "Nonetheless,", "Rather,",
  "Indeed,", "Interestingly,", "With this in mind,", "As one might expect,", "Lastly,".
- **P9.** "Consequently," rather than "Therefore," (Consequently 2–14 times per paper; Therefore only in
  hahn13).
- **P10.** Sentences may start with "And" or "But": "But keep in mind that…" (hahn09).
- **P11.** Sentences that start with "Which" to extend the previous one: "Which is to be expected since…",
  "Which suggests that…". A recent habit (hahn25, 11 times); use it occasionally, not often. Likewise
  "Ditto for…" and "Otherwise," (hahn25 only).
- **P12.** "tells us that" for what an equation or result shows: "The dispersion relation … tells us how the
  wavenumber |k(x)| must adjust…" (hahn07; also hahn08, hahn10, hahn13, hahn25). "allows us to replace…"
  for a step an approximation licenses (hahn07, hahn08, hahn25).
- **P13.** "namely," and "provided…" to state a condition: "when the satellite's mass is sufficiently small,
  namely, when…" (hahn07). "This is due to…" / "which is due to…" to give a cause.
- **P14.** "i.e.," and "e.g.," inside sentences; "∼", "≃" and "O(e²)" in running text; "to lowest order in",
  "to first order in" for the order of an approximation.
- **P15.** Plain and a little informal where it helps: "apples-to-apples", "concave-downish", "aka", "merely
  for labeling purposes", "perfectly adequate", "a bit of fiction" (hahn07). Honest about limits: "this
  admittedly speculative scenario", "Keep in mind though that this does not mean…".
- **P16.** Heavy use of possessives for the object being described: "the ringlet's eccentricity gradient",
  "the particle's distance from the streamline".

## Equations

- **P17.** Equations are part of sentences, introduced by what they are: "draws a particle toward that
  streamline with acceleration A_g = 2Gλ/Δ, (1) where Δ is the particle's distance from the
  streamline."
- **P18.** Every new symbol is defined in a "where" clause right after the equation.
- **P19.** An earlier equation is named by its role, then its number: "the expected wave amplitude, Equation
  (47)", "the equation of motion (81)", "theory, Equation (C3)" (every paper).
- **P20.** Vectors are bold, including the unit vectors r̂, θ̂, ẑ, t̂ and Δ̂ (the author's choice, 2026-10-26).
- **P21.** Derivations are short and verbal; longer ones go to an appendix ("detailed in Appendix A"). The
  2007–2010 papers derived everything in the main text; hahn13 and hahn25 moved derivations to
  appendices, which this paper follows.
- **P22.** Each assumption is stated where it is made, and checked later: "considerable progress can be made if
  we assume that k is constant over a wavelength" (hahn07); "Also keep in mind that these results
  assumed that…" (hahn07); a subsection "Checking the Assumptions" (hahn08) or "Testing the
  Assumptions" (hahn10).
- **P23.** Footnotes may carry a sign or convention note: "Comparison of eq. (55) to the phase convention adopted
  in Shu (1984) shows that the signs of m and ω are reversed" (hahn08).
- **P24.** Equation references in this paper: "Equation (13)" / "Equations (13)–(15)" when the reference
  starts a sentence, "Eqn. (13)" / "Eqns. (13)" otherwise (the author's choice, 2026-10-25). Its
  journal papers use "Equation (19)" throughout, and "Eqn (B11)" in figure legends. (hahn07 and hahn08
  use "equation (46)", "eq. [46]" and "§2.2", the journal's house style then.)

## Comparing simulations to theory

- **P25.** Say what the code solves and which approximations it does not make, so that it is an independent
  check: "the model does so without making any of the wave assumptions invoked in §2.1. The model thus
  provides an independent check of the analytic results obtained above" (hahn07 and hahn08, nearly
  verbatim).
- **P26.** Say why the parameters were chosen: "Those parameters do not correspond to any real ring-satellite
  system; rather, these parameters were chosen to illustrate the results of §2 in the limit in which
  those results were obtained" (hahn07, hahn08).
- **P27.** Name the theory curve with its equation and inputs: "Also plotted is the expected wavenumber,
  equation (46), which compares favorably" (hahn07).
- **P28.** Agreement, from weakest to strongest: "compares favorably", "compares well", "in good agreement", "in
  excellent agreement", "agrees very well". The recent papers mostly say "excellent agreement" and
  "agrees very well"; "is overall very good, though modest discrepancies do exist for…" (hahn25).
- **P29.** Quantify: a fractional error in a conserved quantity, "conserves L_i with a fractional error of
  |ΔL_i/L_i| < 2×10⁻⁵" (hahn07); residuals as "the solutions' fractional errors, which … are all smaller
  than 4×10⁻¹²" (hahn09); "agree with Equations (21)–(23) to within 0.04%" (hahn10).
- **P30.** A disagreement: name the violated key assumption, say what still agrees, and give the size of the
  error. "This disagreement is a consequence of the simulated ringlets' evolution violating a key
  assumption of the … derivation" (hahn25); "Nonetheless, … we find equation (46) to be in excellent
  agreement with the wavenumber k(x) exhibited by the simulated wave. But equation (47) does
  overestimate the amplitude of this simulated wave by a factor of 4" (hahn07).
- **P31.** Explain every visible artifact: "the simulated curve gets a bit ragged at the disk's outer edge,
  which is where the bending wave is just starting to reflect" (hahn07); "Those wiggles are the short
  waves that are described in §2.2.2" (hahn08).
- **P32.** Use the run to confirm the setup's assumptions: "the disk's inclinations are constant, İ(a) = 0,
  which justifies our steady state assumptions, equation (3)" (hahn07).

## Citations

- **P33.** Author–year, either as a sentence subject ("P. Goldreich et al. (1995) show that…") or in a closing
  parenthetical list ("(P. Goldreich & S. Tremaine 1979a, 1979b; E. I. Chiang & P. Goldreich 2000)").
  "see" for pointers: "(see Appendix A)", "see J. M. Hahn & J. N. Spitale (2013) for additional
  details".
- **P34.** Specific equation or section numbers of other papers are cited when a result is used:
  "Equation (2.13) of J. E. Pringle (1981)", "Equation (30) of N. Borderies et al. (1983b)",
  "(and see Section 2.3.3 of … for details about how that is evaluated numerically)".
- **P35.** But not as a stack of parentheticals at the start of a paragraph or sentence (the author cut one
  from Section 3's opening, 2026-10-25). Put the citation where the specific fact is used.
- **P36.** Citations use the ApJ reference style as AASTeX 7 prints it, with first initials in the text
  ("J. M. Hahn & J. N. Spitale (2013)"), as in hahn25 (the author's choice, 2026-10-25).

## Figures and captions

- **P37.** Text cites figures as sentence subjects: "Figure 1 shows…", "Figure 20 plots…", "see Figures 5–6",
  "This is confirmed in Figure 5, which shows…" (hahn07), "The main point of this figure is to show
  that…" (hahn09).
- **P38.** In this paper: "Figure 1" when the reference starts a sentence, "Fig. 1" otherwise (the author's
  choice, 2026-10-25), as for equations.
- **P39.** Captions are full and explanatory: what is plotted against what, in which units, what each color or
  line style is, the run's parameters, and often what the reader should conclude ("…and show that
  higher-e ringlets settle into self-confinement sooner").
- **P40.** Common caption openers: "X plotted vs. Y" (hahn07, hahn09); "Rings model is used to simulate…"
  (hahn07); "Equation (27) is solved for…" (hahn10). A plot of a single function needs only one line:
  "Function C(|k₀|a) from eq. (49) plotted vs. the dimensionless wavenumber |k₀|a" (hahn07).
- **P41.** The theory curve is named with its equation and inputs: "The dashed line is the expected wave
  amplitude, eq. (47), with ε = 1" (hahn07); "simulation (dashed curves with dots) vs. theory (solid
  curves)" (hahn25).
- **P42.** Captions may point back for shared parameters ("Model parameters are otherwise identical to those
  adopted in Figure 2", hahn10), point to the text ("see Section 3.2.1 for details", hahn09), and give a
  numerical-quality figure ("The fractional error in this system's total angular momentum is
  |ΔL_i/L_i| < 5×10⁻⁵", hahn07).

## Structure

- **P43.** Introduction: the physical problem, prior work and its shortcomings, then "This … motivates this
  study, which uses…", and a roadmap ("Section 2 describes…, Section 3 then shows…"). A purpose
  statement may precede the roadmap: "Our purpose here is to determine whether…" (hahn08).
- **P44.** Readers may be told they can skip the derivations: "a reader who is uninterested in the many details
  and derivations can skip ahead to Section 3" (hahn09, hahn10).
- **P45.** When a companion paper has done similar mathematics, say so and be brief: "the derivations presented
  here will be succinct" (hahn08).
- **P46.** Method section explains the code in physical terms first, numerics second; code tests and
  derivations go to appendices ("Appendix C compares the evolution of epi_int_lite simulations to
  theoretical predictions…").
- **P47.** Natural units are stated once, with how to convert to physical units.
- **P48.** "Summary of Findings" ("Our main findings are as follows:") as a numbered list, each item a
  full-sentence finding with a pointer to its section or figure, and real-system numbers where they help
  ("would excite a bending wave having an initial wavelength of about 500 km", hahn07). The 2007–2010
  papers used prose summaries; the recent ones use the numbered list.
- **P49.** "Additional Follow-up Studies" as a numbered list with a short title per item.
- **P50.** Code availability stated plainly, with an invitation to contact the authors. Acknowledgments may be
  light-hearted. A sole author writes "The author thanks…" (hahn07, hahn08).

## Vocabulary

- **P51.** streamline, perturbing streamline, perturbed particle, ringlet, narrow eccentric ringlet, wire (of
  matter), linear density λ, tracer particles, drift step, velocity kick, natural units, nominal
  ringlet/model, self-gravity, secular perturbations.
- **P52.** Gravity of a streamline: "a straight and infinitely long wire of matter whose linear density is λ";
  Δ is "the particle's distance from the streamline" / "the separation between the particle and the
  streamline".
- **P53.** Particles without mass are "massless particles" (hahn13) or "tracer particles" (hahn25), not "test
  particles", which none of the author's papers use. A massless particle's forced motion is "the single
  particle solution" or "the single particle result" (hahn09); "forced eccentricity" (hahn09).
- **P54.** Waves: waves are "launched" and "propagate"; "downstream" / "far downstream" (hahn07, hahn08), with
  distance measured in wavelengths; "wave-excitation zone" (hahn07, hahn08); "wavenumber k",
  "dimensionless wavenumber |k|a", "wavelength λ = 2π/|k|", "dispersion relation", "group velocity",
  "long waves" / "short waves", "standing wave" (on reflection), "reflect"; "the fractional amplitude of
  the bending wave" (hahn07); "normalized disk mass" μ_d (hahn07, hahn08, hahn09). "Leading" and
  "trailing" go in quotes and are defined at first use: "Spiral waves having k < 0 are called 'leading'
  waves" (hahn07).
- **P55.** Resonances: x = (r − r_r)/r_r is "the particle's fractional distance from resonance"; D is "the
  particle's distance from resonance in frequency-squared units"; "the nominal resonance position,
  since exact resonance can be shifted further by the ring's internal forces"; "the nonwave side of the
  resonance" (all hahn09).

## From the author's edits to this paper

- **P56.** Prefer plain words to coined terms ("near point", "distance") and cut asides that another section
  already covers.
- **P57.** One idea per sentence where the sentence would otherwise need a parenthetical.
- **P58.** Lead to a result; don't open with it. A derivation should arrive at its conclusion ("…which means the
  near point lies very near the crossing point"), not state the inequality first and then justify it.
- **P59.** Take the simplest route to a result: an argument that |φ| ≪ 1 beats a chain of inequalities.
- **P60.** State a geometric fact and its consequence in one short clause: "these form a right triangle, so
  |d₀|² = Δ² + |d₁|²". Number such a result if a later argument uses it.
- **P61.** Define a quantity by what it means physically, tied to a direction or motion the reader knows: s is
  "the streamline's slope relative to circular motion in the ring plane, the θ̂ direction".
- **P62.** Every symbol in a figure is used in the text; drop the ones that aren't.
- **P63.** Check that each equation agrees with the figure that illustrates it (signs, which side is the
  hypotenuse).
- **P64.** Precise notation: a particle orbits at Ω_a (its mean angular velocity), not Ω.
- **P65.** When a result comes from the literature, state the solution with its defining parameters and cite the
  source's equation numbers; don't reproduce the equation of motion it solves.
- **P66.** Name a solution's parts by their roles: h = h_forced + h_wave, with h_wave given by cases (on the wave
  side and zero otherwise).
- **P67.** One defined quantity per display: "A_z = … (25), where A_V = … (26)", not a stacked block.
- **P68.** Write time-dependent quantities explicitly in t (λ_s = Ω_{a_s}t + λ_{s0}), so that what follows from
  them, such as a pattern speed, is visibly derived.
- **P69.** A test of the code is reported briefly and to the point: the theory it checks, the setup in a few
  sentences, one figure that compares simulation to the theory curve the text uses (only that curve),
  and the agreement in a number; the remaining detail stays in the test's notebook, linked in a footnote.
- **P70.** Theory first, then the tests it makes possible: a theory subsection ends by listing the tests ("gives
  us a few opportunities to test…"), and each test then gets its own subsection.
- **P71.** An approximate equation uses ≃ and says what it neglects in words ("where terms of order z⁴ are
  neglected"), rather than carrying a "+ 𝒪(z⁴)" term. Evaluate at a point with a bar, ∂Φ/∂r|_{z=0}, as with
  ν² = ∂²Φ/∂z²|_{z=0}. (2026-11-11)
- **P72.** Give a solution first and its parameters after it: "the particle's coordinates are (7), where the
  mean angular velocity and epicyclic frequency are (8)". Give a symbol a subscript when it is one component
  of a vector, or when the symbol is used elsewhere for something else (h_z, not h). (2026-11-12)
- **P73.** Number the lines of a multi-line equation (6a, 6b) when the text cites one line, and cite it by
  number, not "the second of Eqns. (6)". Don't display an equation that only restates an earlier one with a
  parameter evaluated; say so in words ("so ν(r) ≃ ν(a) in it") and go to its solution. Put the source of a
  displayed result after it: "the coordinates are (7) when written in terms of the epicyclic elements
  (Eqn. A9 of H25)". (2026-11-12)
- **P74.** Don't restate what is understood: leave "at r = a, z = 0" off quantities already subscripted a,
  and don't add a sentence verifying that a stated solution solves its equation. Caveats the argument doesn't
  use (e.g. "this holds only in a massless ring") are cut. (2026-11-12)
- **P75.** No comma before a "where" that defines the symbols of the inline expression just written: "at
  δθ′ = −c/r′ where c = (d₀·t̂)/N", "the angle φ where tan φ = s". Keep the comma when "where" names a place
  ("Near r_V, where …") or follows a displayed equation. (2026-11-12)
- **P76.** Prefer a short order-of-magnitude argument to a displayed derivation when only an inequality is
  needed: "displaces the streamline by about ea′, so s_r ∼ me" rather than writing out the mode. Say "minus its
  component" rather than "with its component removed". (2026-11-12)
- **P77.** When one displayed equation follows from another by a standard step, just name the step: "The Fourier
  expansion of Eqn. (25) in θ is (26)", not how the step is done ("expanding Δ⁻³ in Laplace coefficients") or a
  citation for a routine operation; define the symbols the result introduces after it. (2026-11-12)
- **P78.** Justify an approximation by what the work actually does, plainly: "The simulations here follow particles
  near a single mth resonance, so Epi_int_3d keeps only the resonant term", not "near one resonance only the
  particle's response to that one term matters". (2026-11-12)
