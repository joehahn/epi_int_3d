# Writing style for the paper

The author's style, taken from his papers in this folder, plus the author's edits to this paper. Use it
when drafting text for `paper/epi_int_3d.tex`. Drafted from Hahn & Spitale (2013, ApJ 772, 122; "hahn13")
and Hahn et al. (2025, ApJ 986, 220; "hahn25") and approved by the author 2026-10-25. Extended 2026-10-31
from Hahn (2007, ApJ 665, 856; "hahn07"), Hahn (2008, ApJ 680, 1569; "hahn08"), Hahn, Spitale & Porco
(2009, ApJ 699, 686; "hahn09") and Hahn (2010, ApJ 719, 1699; "hahn10"). Ward & Hahn (2003, AJ 125, 3389)
is mostly Ward's prose and is not used.

Where the author's habits have changed since 2007, follow the recent papers (hahn13, hahn25). A habit
seen in only one paper is marked as such.

## Voice and tense

- Mostly impersonal: "The following describes…", "This section…", "the code identifies…", "Here, …".
  "We" appears sparingly, for the authors' own judgments and choices: "we conclude", "we speculate",
  "we multiply time t by…". Not "I", and not "one" except in "allows one to…".
- Derivation steps may be imperatives: "Begin by considering a planetary ring…" (hahn07), "Next,
  Fourier-expand the longitude-dependent factor…" (hahn09), "Now derive the threshold…" (hahn10),
  "Begin with…" (hahn25).
- Present tense for what the code does and what figures and simulations show ("Figure 3 shows that…",
  "simulations show that…"). Past tense only for earlier work or completed runs ("was used to
  simulate"). Future tense for pointers ahead: "Section 2.2.3 will quantify when that assumption breaks
  down" (hahn08); "it will be shown below that this contribution is unimportant" (hahn07, hahn08).
- "is used to" for a code, model or equation doing work: "The rings model is used to simulate…"
  (hahn07, hahn08); "Equations (31) and (32) are used to calculate…" (hahn10).
- The code's name starts a sentence capitalized: "Epi_int_lite is a child of the epi_int N-body
  integrator…".
- Code names are set in typewriter, as in hahn13 (hahn25 used roman): write `\code{epi\_int\_3d}`,
  a macro in the paper's preamble that can be switched back to roman in one line (the author's
  choice, 2026-10-26).

## Sentences and connectives

- Long, chained sentences that carry the reasoning: "Because a ringlet is narrow, all particles are in
  close proximity to the nearby portions of all streamlines, which allows us to approximate a
  streamline as an infinitely long wire…". Physical reason first, then the approximation it licenses.
- Characteristic openers, the most frequent first: "However,", "Consequently,", "Note that", "Note
  though that", "Also note that", "Evidently,", "For instance,", "Keep in mind (though) that", "Also keep
  in mind that", "Of course,", "The main point is that", "So,", "So if…", "Nonetheless,", "Rather,",
  "Indeed,", "Interestingly,", "With this in mind,", "As one might expect,", "Lastly,".
- "Consequently," rather than "Therefore," (Consequently 2–14 times per paper; Therefore only in
  hahn13).
- Sentences may start with "And" or "But": "But keep in mind that…" (hahn09).
- Sentences that start with "Which" to extend the previous one: "Which is to be expected since…",
  "Which suggests that…". A recent habit (hahn25, 11 times); use it occasionally, not often. Likewise
  "Ditto for…" and "Otherwise," (hahn25 only).
- "tells us that" for what an equation or result shows: "The dispersion relation … tells us how the
  wavenumber |k(x)| must adjust…" (hahn07; also hahn08, hahn10, hahn13, hahn25). "allows us to replace…"
  for a step an approximation licenses (hahn07, hahn08, hahn25).
- "namely," and "provided…" to state a condition: "when the satellite's mass is sufficiently small,
  namely, when…" (hahn07). "This is due to…" / "which is due to…" to give a cause.
- "i.e.," and "e.g.," inside sentences; "∼", "≃" and "O(e²)" in running text; "to lowest order in",
  "to first order in" for the order of an approximation.
- Plain and a little informal where it helps: "apples-to-apples", "concave-downish", "aka", "merely
  for labeling purposes", "perfectly adequate", "a bit of fiction" (hahn07). Honest about limits: "this
  admittedly speculative scenario", "Keep in mind though that this does not mean…".
- Heavy use of possessives for the object being described: "the ringlet's eccentricity gradient",
  "the particle's distance from the streamline".

## Equations

- Equations are part of sentences, introduced by what they are: "draws a particle toward that
  streamline with acceleration A_g = 2Gλ/Δ, (1) where Δ is the particle's distance from the
  streamline."
- Every new symbol is defined in a "where" clause right after the equation.
- An earlier equation is named by its role, then its number: "the expected wave amplitude, Equation
  (47)", "the equation of motion (81)", "theory, Equation (C3)" (every paper).
- Vectors are bold, including the unit vectors r̂, θ̂, ẑ, t̂ and Δ̂ (the author's choice, 2026-10-26).
- Derivations are short and verbal; longer ones go to an appendix ("detailed in Appendix A"). The
  2007–2010 papers derived everything in the main text; hahn13 and hahn25 moved derivations to
  appendices, which this paper follows.
- Each assumption is stated where it is made, and checked later: "considerable progress can be made if
  we assume that k is constant over a wavelength" (hahn07); "Also keep in mind that these results
  assumed that…" (hahn07); a subsection "Checking the Assumptions" (hahn08) or "Testing the
  Assumptions" (hahn10).
- Footnotes may carry a sign or convention note: "Comparison of eq. (55) to the phase convention adopted
  in Shu (1984) shows that the signs of m and ω are reversed" (hahn08).
- Equation references in this paper: "Equation (13)" / "Equations (13)–(15)" when the reference
  starts a sentence, "Eqn. (13)" / "Eqns. (13)" otherwise (the author's choice, 2026-10-25). Its
  journal papers use "Equation (19)" throughout, and "Eqn (B11)" in figure legends. (hahn07 and hahn08
  use "equation (46)", "eq. [46]" and "§2.2", the journal's house style then.)

## Comparing simulations to theory

- Say what the code solves and which approximations it does not make, so that it is an independent
  check: "the model does so without making any of the wave assumptions invoked in §2.1. The model thus
  provides an independent check of the analytic results obtained above" (hahn07 and hahn08, nearly
  verbatim).
- Say why the parameters were chosen: "Those parameters do not correspond to any real ring-satellite
  system; rather, these parameters were chosen to illustrate the results of §2 in the limit in which
  those results were obtained" (hahn07, hahn08).
- Name the theory curve with its equation and inputs: "Also plotted is the expected wavenumber,
  equation (46), which compares favorably" (hahn07).
- Agreement, from weakest to strongest: "compares favorably", "compares well", "in good agreement", "in
  excellent agreement", "agrees very well". The recent papers mostly say "excellent agreement" and
  "agrees very well"; "is overall very good, though modest discrepancies do exist for…" (hahn25).
- Quantify: a fractional error in a conserved quantity, "conserves L_i with a fractional error of
  |ΔL_i/L_i| < 2×10⁻⁵" (hahn07); residuals as "the solutions' fractional errors, which … are all smaller
  than 4×10⁻¹²" (hahn09); "agree with Equations (21)–(23) to within 0.04%" (hahn10).
- A disagreement: name the violated key assumption, say what still agrees, and give the size of the
  error. "This disagreement is a consequence of the simulated ringlets' evolution violating a key
  assumption of the … derivation" (hahn25); "Nonetheless, … we find equation (46) to be in excellent
  agreement with the wavenumber k(x) exhibited by the simulated wave. But equation (47) does
  overestimate the amplitude of this simulated wave by a factor of 4" (hahn07).
- Explain every visible artifact: "the simulated curve gets a bit ragged at the disk's outer edge,
  which is where the bending wave is just starting to reflect" (hahn07); "Those wiggles are the short
  waves that are described in §2.2.2" (hahn08).
- Use the run to confirm the setup's assumptions: "the disk's inclinations are constant, İ(a) = 0,
  which justifies our steady state assumptions, equation (3)" (hahn07).

## Citations

- Author–year, either as a sentence subject ("P. Goldreich et al. (1995) show that…") or in a closing
  parenthetical list ("(P. Goldreich & S. Tremaine 1979a, 1979b; E. I. Chiang & P. Goldreich 2000)").
  "see" for pointers: "(see Appendix A)", "see J. M. Hahn & J. N. Spitale (2013) for additional
  details".
- Specific equation or section numbers of other papers are cited when a result is used:
  "Equation (2.13) of J. E. Pringle (1981)", "Equation (30) of N. Borderies et al. (1983b)",
  "(and see Section 2.3.3 of … for details about how that is evaluated numerically)".
- But not as a stack of parentheticals at the start of a paragraph or sentence (the author cut one
  from Section 3's opening, 2026-10-25). Put the citation where the specific fact is used.
- Citations use the ApJ reference style as AASTeX 7 prints it, with first initials in the text
  ("J. M. Hahn & J. N. Spitale (2013)"), as in hahn25 (the author's choice, 2026-10-25).

## Figures and captions

- Text cites figures as sentence subjects: "Figure 1 shows…", "Figure 20 plots…", "see Figures 5–6",
  "This is confirmed in Figure 5, which shows…" (hahn07), "The main point of this figure is to show
  that…" (hahn09).
- In this paper: "Figure 1" when the reference starts a sentence, "Fig. 1" otherwise (the author's
  choice, 2026-10-25), as for equations.
- Captions are full and explanatory: what is plotted against what, in which units, what each color or
  line style is, the run's parameters, and often what the reader should conclude ("…and show that
  higher-e ringlets settle into self-confinement sooner").
- Common caption openers: "X plotted vs. Y" (hahn07, hahn09); "Rings model is used to simulate…"
  (hahn07); "Equation (27) is solved for…" (hahn10). A plot of a single function needs only one line:
  "Function C(|k₀|a) from eq. (49) plotted vs. the dimensionless wavenumber |k₀|a" (hahn07).
- The theory curve is named with its equation and inputs: "The dashed line is the expected wave
  amplitude, eq. (47), with ε = 1" (hahn07); "simulation (dashed curves with dots) vs. theory (solid
  curves)" (hahn25).
- Captions may point back for shared parameters ("Model parameters are otherwise identical to those
  adopted in Figure 2", hahn10), point to the text ("see Section 3.2.1 for details", hahn09), and give a
  numerical-quality figure ("The fractional error in this system's total angular momentum is
  |ΔL_i/L_i| < 5×10⁻⁵", hahn07).

## Structure

- Introduction: the physical problem, prior work and its shortcomings, then "This … motivates this
  study, which uses…", and a roadmap ("Section 2 describes…, Section 3 then shows…"). A purpose
  statement may precede the roadmap: "Our purpose here is to determine whether…" (hahn08).
- Readers may be told they can skip the derivations: "a reader who is uninterested in the many details
  and derivations can skip ahead to Section 3" (hahn09, hahn10).
- When a companion paper has done similar mathematics, say so and be brief: "the derivations presented
  here will be succinct" (hahn08).
- Method section explains the code in physical terms first, numerics second; code tests and
  derivations go to appendices ("Appendix C compares the evolution of epi_int_lite simulations to
  theoretical predictions…").
- Natural units are stated once, with how to convert to physical units.
- "Summary of Findings" ("Our main findings are as follows:") as a numbered list, each item a
  full-sentence finding with a pointer to its section or figure, and real-system numbers where they help
  ("would excite a bending wave having an initial wavelength of about 500 km", hahn07). The 2007–2010
  papers used prose summaries; the recent ones use the numbered list.
- "Additional Follow-up Studies" as a numbered list with a short title per item.
- Code availability stated plainly, with an invitation to contact the authors. Acknowledgments may be
  light-hearted. A sole author writes "The author thanks…" (hahn07, hahn08).

## Vocabulary

- streamline, perturbing streamline, perturbed particle, ringlet, narrow eccentric ringlet, wire (of
  matter), linear density λ, tracer particles, drift step, velocity kick, natural units, nominal
  ringlet/model, self-gravity, secular perturbations.
- Gravity of a streamline: "a straight and infinitely long wire of matter whose linear density is λ";
  Δ is "the particle's distance from the streamline" / "the separation between the particle and the
  streamline".
- Particles without mass are "massless particles" (hahn13) or "tracer particles" (hahn25), not "test
  particles", which none of the author's papers use. A massless particle's forced motion is "the single
  particle solution" or "the single particle result" (hahn09); "forced eccentricity" (hahn09).
- Waves: waves are "launched" and "propagate"; "downstream" / "far downstream" (hahn07, hahn08), with
  distance measured in wavelengths; "wave-excitation zone" (hahn07, hahn08); "wavenumber k",
  "dimensionless wavenumber |k|a", "wavelength λ = 2π/|k|", "dispersion relation", "group velocity",
  "long waves" / "short waves", "standing wave" (on reflection), "reflect"; "the fractional amplitude of
  the bending wave" (hahn07); "normalized disk mass" μ_d (hahn07, hahn08, hahn09). "Leading" and
  "trailing" go in quotes and are defined at first use: "Spiral waves having k < 0 are called 'leading'
  waves" (hahn07).
- Resonances: x = (r − r_r)/r_r is "the particle's fractional distance from resonance"; D is "the
  particle's distance from resonance in frequency-squared units"; "the nominal resonance position,
  since exact resonance can be shifted further by the ring's internal forces"; "the nonwave side of the
  resonance" (all hahn09).

## From the author's edits to this paper

- Prefer plain words to coined terms ("near point", "distance") and cut asides that another section
  already covers.
- One idea per sentence where the sentence would otherwise need a parenthetical.
