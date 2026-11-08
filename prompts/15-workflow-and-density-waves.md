# 15 — spiral density waves, style guides, outside review and a ledger

i'd like to make several changes to this project. dont make any changes until we've discussed. first, i'd like this project to simulate spiral density waves as well, with rigourous comparisons of theory to observations. also this project, agent-banking-watch, created a ledger dashboard to manage the development of a document, and i'd like to do similar to manage the development of this project's paper. the agent-banking-watch (ABW) project also used multiple ai models to verify and confirm that project's findings, should we do similar here? with multiple ai to: double check the math, check the code, check the paper, check the code tests? so that we can more efficiently do the math in my particular style rather than in this ai's style or that ai's style, and write in my voice. do the code development in my style (rather than ai's style) so I can spot-check the code at a glance. the jupyter notebooks in scientific repos show that i have a particular style to validate and test numerical work, and i'd like to preserve that style in this project's testing. if agreed, should different ai be used to solve the physics equations, confirm the math, write the results, write the code in my style, manage the tests? with the ledger providing visibility into all of the above? if so, should that db be called a ledger or something else? also should physics problem-solving and derivation have its own style guide? likewise for paper development, code dev, and testing? and for data visualization? should db include links to those style guides?

[Read ABW's README, `scripts/review.py` and its reviews. Proposed: density waves reuse the bending-wave setup (the Fourier satellite at an inner Lindblad resonance; planar, so `source/` and `fast/` can be compared), checked against the WKB dispersion relation, the Goldreich & Tremaine (1979) torque, Shu (1984) and convergence. A ledger that links each paper claim to its equation, sympy check, code, test and figure; keep the name "ledger". One author (Claude) and outside models as reviewers only, since several authors would mean several AI styles: blind re-derivation for math, diff-versus-equation for code, "would this test fail if the code were wrong?" for tests, claims and a referee read for the paper. Style guides for derivations, code, testing, data visualization and the paper, drawn from the author's pre-AI work, with numbered rules the ledger links to. Asked about Cassini data, restyling `fast/`, a review budget and linear-only.]

i'd also like to use ai to reduce the knock on edits, the goal is to complete this scientific project as correctly as possible with paper and visuals composed in my voice and style with fewest interventions by me. can you confirm that Claude is best or nearly best for "derives, writes code, writes tests and drafts prose, each against its style guide."? likewise should reviewers be GPT, Gemini and Grok or something else? regarding tests: comparing sim to theory usually requires a judgement call, lets lean on our expertise & judgement to decide if something does or doesnt agree, dont overthink this part. 1 we wont use real cassini data in current version of this project. 3 have ledger track expenses be sufficiently granular ie code dev versus math versus paper writing vs judging. 4 this project will only do linear waves. what about managing knockons? what does 2 mean?

[Claude for code, tests, repo work and prose (one consistent author working in the repo); for hard derivations GPT and Gemini reasoning models may be as good or better, so a bake-off on an already-verified derivation decides. Reviewers: GPT and Gemini, Grok optional for the paper; Perplexity not needed. Knock-ons: one source for every number (a macro file written by scripts, figures by script), dependencies with file hashes in the ledger so changes mark downstream items stale, a knock-on pass listing each fix for veto, and corrections turned into style rules. Costs by activity, outside reviews exact and Claude's split estimated. Question 2 meant restyling `fast/`, written in Claude's style.]

lets put reviewing fast code on our todo list, we dont have to tackle that now, can do that later. all sounds right

[`PLAN.md` draft v9: Phase 3b (style guides, one author and outside reviewers, fewer knock-ons, the ledger), Phase 3c (linear spiral density waves), a density-wave row in Phase 4 and a note that agreement with theory is a judgment call, Phase 5 and the decisions narrowed to linear waves only with no Cassini data, and restyling `fast/` on the to-do list.]

commit, then start the code and testing style guides

[Committed the plan. Drafted `style/code.md` (C1–C25) from `source/` and `style/testing.md` (T1–T21) from seven
pre-AI notebooks: the four ported tests and, from the epi_int_lite repo, the self-confining nominal run, the
number-of-streamlines convergence study and the equilibrium study. Each guide lists Claude's habits that break
its rules. Open questions: `exec` and `import *`, input checking, markdown cells, seaborn, restyling Claude's notebooks.]

the purpose of most tests are to: confirm sim output agrees with theory, show that output agrees with expectations, and that constants of the motion are conserved reasonably well. usually documents via to-the-point plots make the point at a glance. some of which are publication quality and go into journal articles or conference reports

[Added as the testing guide's statement of purpose, and T22: one point per plot, the key plot drawn as a candidate paper figure.]

2 follow my style 3 my notebooks often have a comment explaining its purpose, do that in my comments' style in voice. 4 plots should use seaborn and preserve my appearance. should we make a data vis style guide? 5 unsure we'll deal with that later. 1 i dont understand the question

[No input checks (C4); a purpose comment and no markdown (T3); the author's seaborn setup (T7); restyling Claude's notebooks on the to-do list. Explained question 1: `exec` of inputs.py and `import *` hide where names come from; recommended keeping `exec` everywhere and `import *` in notebooks, with explicit imports in new modules. Recommended a data-visualization guide covering notebook plots, paper figures, animations and the Pages site.]

agreed on 1, commit and do the data vis guide

[Committed the code and testing guides. Drafted `style/dataviz.md` from all 29 hahn25 figures (the author's
notebook plots placed directly in the paper) and the notebooks: seaborn's gray grid, one 10×6 panel at column
width, theory dashed in the simulation's color, the reference run thicker, word-then-symbol labels, normalized
axes, no titles on paper figures, PDFs made by script. Open questions: diagrams, color maps, restyling the
existing paper figures.]

1 unsure do whats best. 2 those claude-generated color maps are likely going to get deleted. 3 do later

[Diagrams are plain line drawings on white in the same font and palette (V17); no 2D color maps without
asking (V18); restyling `vertical_resonance` and `zigzag` on the to-do list.]

push. go ahead with the derivations guide

[Pushed. Drafted `style/derivations.md` from hahn07 §2 and the appendices of hahn13 and hahn25 (quotes checked
against the PDF pages): set up the physical system and split it into parts; every approximation with its physical
reason, small quantity and order; assumptions stated and checked later; dimensionless functions checked in their
limits; each term interpreted; results the code uses written as numbered steps. Open questions: whether `notes/`
derivations follow it from the start, and the preferred method.]

2 i prefer equations of motion and fluxes. but I'll occasionally use orbit elements and the disturbing function when that is more convenient 1 i dont understand

[Added as D4a. Explained question 1: derivations are written twice, as Claude's terse working notes in `notes/`
and again in the author's style for the paper; should the style apply from the start?]

the derivation in the notes should do the math in my preferred style

[D0: derivations in `notes/` follow the guide from the start, as drafts of the paper's appendices.]

push, then number the paper guide's rules

[Pushed. `joes_papers/writing_style.md`'s 70 rules numbered P1–P70, with a note that rule numbers are permanent
(new rules take the next free number or a lettered one, e.g. P9a); `tools/lint_paper.py` cites P9 and P53.]
