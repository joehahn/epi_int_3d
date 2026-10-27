# Writing style for the paper

The author's style, taken from Hahn & Spitale (2013, ApJ 772, 122; "hahn13") and Hahn et al. (2025,
ApJ 986, 220; "hahn25"), plus the author's edits to this paper. Use it when drafting text for
`paper/epi_int_3d.tex`. Drafted and approved by the author 2026-10-25.

## Voice and tense

- Mostly impersonal: "The following describes…", "This section…", "the code identifies…", "Here, …".
  "We" appears sparingly, for the authors' own judgments and choices: "we conclude", "we speculate",
  "we multiply time t by…". Not "I", and not "one" except in "allows one to…".
- Present tense for what the code does and what figures and simulations show ("Figure 3 shows that…",
  "simulations show that…"). Past tense only for earlier work or completed runs ("was used to
  simulate").
- The code's name starts a sentence capitalized: "Epi_int_lite is a child of the epi_int N-body
  integrator…".
- Code names are set in typewriter, as in hahn13 (hahn25 used roman): write `\code{epi\_int\_3d}`,
  a macro in the paper's preamble that can be switched back to roman in one line (the author's
  choice, 2026-10-26).

## Sentences and connectives

- Long, chained sentences that carry the reasoning: "Because a ringlet is narrow, all particles are in
  close proximity to the nearby portions of all streamlines, which allows us to approximate a
  streamline as an infinitely long wire…". Physical reason first, then the approximation it licenses.
- Characteristic openers: "Note that", "Note though that", "Also note that", "Keep in mind (though)
  that", "The main point is that", "So,", "Consequently,", "Therefore,", "Nonetheless,", "Otherwise,",
  "Rather,", "Ditto for…", "Indeed,".
- Sentences that start with "Which" to extend the previous one: "Which is to be expected since…",
  "Which suggests that…". (A habit of the author's; use it occasionally, not often.)
- "i.e.," and "e.g.," inside sentences; "∼", "≃" and "O(e²)" in running text; "to lowest order in",
  "to first order in" for the order of an approximation.
- Plain and a little informal where it helps: "apples-to-apples", "concave-downish", "aka", "merely
  for labeling purposes", "perfectly adequate". Honest about limits: "this admittedly speculative
  scenario", "Keep in mind though that this does not mean…".
- Heavy use of possessives for the object being described: "the ringlet's eccentricity gradient",
  "the particle's distance from the streamline".

## Equations

- Equations are part of sentences, introduced by what they are: "draws a particle toward that
  streamline with acceleration A_g = 2Gλ/Δ, (1) where Δ is the particle's distance from the
  streamline."
- Every new symbol is defined in a "where" clause right after the equation.
- Derivations are short and verbal; longer ones go to an appendix ("detailed in Appendix A").
- Equation references in this paper: "Equation (13)" / "Equations (13)–(15)" when the reference
  starts a sentence, "Eqn. (13)" / "Eqns. (13)" otherwise (the author's choice, 2026-10-25). Its
  journal papers use "Equation (19)" throughout, and "Eqn (B11)" in figure legends.

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

- Text cites figures as sentence subjects: "Figure 1 shows…", "Figure 20 plots…", "see Figures 5–6".
- In this paper: "Figure 1" when the reference starts a sentence, "Fig. 1" otherwise (the author's
  choice, 2026-10-25), as for equations.
- Captions are full and explanatory: what is plotted against what, in which units, what each color or
  line style is, the run's parameters, and often what the reader should conclude ("…and show that
  higher-e ringlets settle into self-confinement sooner").

## Structure

- Introduction: the physical problem, prior work and its shortcomings, then "This … motivates this
  study, which uses…", and a roadmap ("Section 2 describes…, Section 3 then shows…").
- Method section explains the code in physical terms first, numerics second; code tests and
  derivations go to appendices ("Appendix C compares the evolution of epi_int_lite simulations to
  theoretical predictions…").
- Natural units are stated once, with how to convert to physical units.
- "Summary of Findings" ("Our main findings are as follows:") as a numbered list, each item a
  full-sentence finding with a pointer to its section or figure.
- "Additional Follow-up Studies" as a numbered list with a short title per item.
- Code availability stated plainly, with an invitation to contact the authors. Acknowledgments may be
  light-hearted.

## Vocabulary

- streamline, perturbing streamline, perturbed particle, ringlet, narrow eccentric ringlet, wire (of
  matter), linear density λ, tracer particles, drift step, velocity kick, natural units, nominal
  ringlet/model, self-gravity, secular perturbations.
- Gravity of a streamline: "a straight and infinitely long wire of matter whose linear density is λ";
  Δ is "the particle's distance from the streamline" / "the separation between the particle and the
  streamline".

## From the author's edits to this paper

- Prefer plain words to coined terms ("near point", "distance") and cut asides that another section
  already covers.
- One idea per sentence where the sentence would otherwise need a parenthetical.
