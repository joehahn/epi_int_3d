# 14 — a 5-page deck on using Claude Code for research

review this email and come up with a plan for a 5-page deck that will provide my space science collegauges a view into how i use AI for research. here is the email from [an SSI colleague]. possible title is 'Using Athropic Claude Code for research. page 1 might be everything the audience should know about CC in 3 or so bullets, assume some in audience are CC users and some are not. remaining pages could be about CC for literature search, CC for astrophysics & deriving equations of motion, CC for code development & testing, CC for analyzing and visualizing data, CC for drafting scientific publications. Focus will be about how I'm using for this project, but I want my findings to be generalized to the audience which will be phds in physics, astro, chemistry, geology, and education. each page should include 1+ visuals or scrreengrabs plus at most 2-3 short to-the-point bullets in my voice, mostly just to remind me what to discuss while on the page. also think about where to store this content, maybe in folder public/ai_for_research folder that will be private until finalyzed but might go public soon. what you think?

[pasted email, summarized here because this repo is public: an SSI AI LLM Working Group meeting asks for volunteers to give ~10-minute talks, a few slides plus Q&A, on how they use AI in their research.]

[Proposed a plan: five pages for a ten-minute slot, merging code and data analysis into one page; a visual from this project and 2–3 bullets per page; a closing generalization on every page; storage in a git-ignored talks/ folder (this repo is public) or a private claude.ai Slides artifact. Awaiting go-ahead.]

we are eventually going to talk about this project via linkedin post. so plan for other content that we will eventually want to expose to the public, plan accordinginly. should plan for this deck be put into .md file we can iterate on? talk to tomorro

[Yes. New git-ignored outreach/ folder (the repo is public; each piece is un-ignored when final): outreach/README.md lists the public pieces (this talk, a LinkedIn post that reuses the deck as a PDF carousel, a README section on how the project was built, the paper), the shared assets they draw on, and a pre-publication checklist; outreach/2026-10-28_ssi_talk/plan.md is the deck plan to iterate on.]

do outreach folder. dont put time/cost into deck. pop the .md plan into vscode

[The outreach/ folder is in place (git-ignored). Time and cost are out of the deck plan: page 1's visual is the terminal screenshot alone. The plan is open in VS Code.]

this is how i use CC, "  - It's an AI that works inside my project folder: reads files, writes and runs code, and I approve or redirect each step." but I havent used either of these:   - A rules file (CLAUDE.md) and memory carry my standards from one session to the next.   - Every prompt is saved, so the work is reproducible and auditable. so maybe they dont merit being discussed in deck? page 1 is about what the audience should know about CC, keeping in mind that probably half of the audience has not touched. What 2 topmost bullets should this audience know about CC that is relevant to their research in observational and theoritical astro and geology and planetary atmospheres and cratering and ionosphere/magnetospheric physics , chemistry

[Agreed, dropped both. Page 1's two bullets are now (1) it works inside your project folder on your own files, papers and data, and writes and runs code, with you approving each step, not a chat window to copy-paste from; (2) it's fast and usually right but sometimes confidently wrong, so pair it with checks you trust, and you stay the scientist. Updated in plan.md.]

I stay the scientist. can be expanded, i direct, i nudege CC out of rabbit holes, I gaurd against ai slop (think or a better phrase plze). what else do i do?

[Suggested "referee" for guarding against AI slop, plus what the record shows the user doing: choosing the problem and the approximations, setting standards, sanity-checking physics, asking "why?", hand-checking equations, deciding what is in and out, owning the voice, signing off. plan.md page 1 now has a third bullet: "I stay the scientist: I pick the problem and the physics, steer it out of rabbit holes, ask "why?" until the argument is transparent, and referee everything before it goes in the paper."]

do the above

[Dropped the "confidently wrong" bullet; page 1 now has two: works inside the project folder, and "I stay the scientist…".]

what kind of AI in "It's an AI that works inside my project folder:" and agent, a skill, something else?

[An agent: Anthropic calls Claude Code an agentic coding tool, i.e. the Claude LLM running in a loop with tools (read and edit files, run commands) in the terminal; skills are add-on instruction packages it can load, not what it is. plan.md now says "It's an AI agent that works inside my project folder".]

page 5 should say the CC is managing the latex, aastex, figure develpoment . page 4 should speak to using CC for numerical simulation, and test sim output against known theoretical predictions. drop "[an SSI-internal item] - Format: claude.ai Slides artifact (private link, exports .pptx/PDF) or a local .pptx." should Speedup factor (to be computed, not quoted from memory). be noted in code dev page?

[Done. Page 4 is now "Numerical simulations, tested against theory": the port and speedup, simulations checked against viscous spreading, libration (Borderies, Goldreich & Tremaine) and J2 nodal regression, and regression tests to 1e-11. Yes on the speedup, now computed: ~36× at 1000 streamlines (447 ms → 12 ms per step, notes/phase0b/driver.md), identical results. Page 5 is "Writing the paper": LaTeX/AASTeX, bibliography, figure scripts, compiling. Both open questions dropped.]

page 1 should say something about getting context from files, and how important that is

[Added as page 1's middle bullet: "Context is everything: before it answers, it reads my code, notes, papers and data, so its work fits my project instead of being a generic guess. The better my files, the better its work." Inserted after the user's own edited first bullet, which is unchanged.]

should It found 10 bugs in our older code refer to bugs and shortcomings? is there a better way to expose my dirty laundry? and   - It ported my old Python 2 code, then made it ~36× faster (447 ms → 12 ms per step at 1000 streamlines), with identical results (`notes/phase0b/driver.md`). was achieved by parallelization with which python package?

[Correction: 9 issues, not 10. Of them 3 are bugs (1, 4, 6), 1 an outdated test input (2), 2 approximations improved (7, 9), 1 a deliberate choice (8), 2 open limits of the numerical scheme (3, 5). Page 3 now reads "It reviewed my published code like a careful collaborator: 9 issues logged, from small bugs to approximations worth improving, and one matters for the spiral waves I'm after." The speedup is Numba: @njit compiles the loops (×6 for the gravity kernel alone) and prange runs them on 10 threads (×50 for the kernel; ×36 for a whole step). Page 4 now says "with Numba, which compiles Python loops and runs them in parallel".]

time to make deck?

[Built the 5-slide deck as a private claude.ai Slides artifact, https://claude.ai/artifact/PvHynoKBwd94cawaCktta3, following plan.md: (1) what Claude Code is, with a condensed real exchange from this session; (2) literature, with a bib entry and its SciX-check line; (3) physics, with paper Fig. 1, the user's "d1 << Delta isn't obvious" question and the G18/G19 sympy PASS lines; (4) simulations, with the libration-vs-theory plot from tests/libration; (5) writing, with page 3 of the draft and a real before/after edit. Speaker notes carry the full bullets. Sources copied to outreach/2026-10-28_ssi_talk/deck/.]

display deck outside of browser

[Rendered the five slides to outreach/2026-10-28_ssi_talk/Using_Claude_Code_for_Research.pdf with headless Chrome (from the deck sources and local images) and opened it in Preview.]

i tried to copy/paste first bullet but I'm just grabbing "An AI agent in co de a n wind ow"

[Preview mangles selection in multi-column Chrome PDFs; gave the bullet text.]

so fix the above, dont use crome when making pdf?

[Now built as a real PowerPoint file, outreach/2026-10-28_ssi_talk/Using_Claude_Code_for_Research.pptx, by make_pptx.py (python-pptx, real text boxes and speaker notes; Georgia, Helvetica Neue and Menlo, since the web fonts are not installed), and PowerPoint exports the PDF. Text extracts in order and copies cleanly; the bib entry keeps its indentation.]

it looks like AI made that deck. can you find examples of previous decks to use as style guides, they are likely in ~/Dropbox/academic/text/conferences
[and, mid-turn: "dda folder has something from 2026"; "epi_int_lite repo might have recent decks"; "dont use latex to make this deck, that is overkill"]

[The local conferences folder holds LaTeX/foils talks from 2010–2012 (pale yellow, Helvetica, lowercase bullets); dda/2026 holds only expenses, and the DDA 2026 poster is a colleague's. The joehahn/epi_int_lite repo has the recent PowerPoint decks, DDA 2025 and Space Jam 2023: 4:3, white, Calibri, a header line ending in a colon, nested round bullets in lowercase, matplotlib plots beside the text. outreach/slide_style.md records that style, and make_pptx.py now rebuilds the deck in it (PowerPoint, no LaTeX); PDF exported by PowerPoint.]

make background light yellow or light blue. each page has a lot of unused vertical whitespace. in page 1 put CC snapshot in box, maybe with white background?

[Light yellow (#FFFEE5, the old LaTeX talks' color) on every slide; larger header (28 pt) and bullets (20–24 pt) and bigger figures to fill the page; the CC exchange on page 1, the bib entry on page 2 and the sympy output on page 3 sit in white boxes with a thin grey border.]

 i think we want CC reads my code, notes, papers & data before it answers, so its work fits my project to also say that is also manages new code dev and paper & deck (and proposal drafting, eventually)

[Added a sub-bullet under "context is everything:" on page 1: "it also manages new code development, the paper & this deck (proposal drafting, eventually)".]

the a CC exchange from this project, condensed: example is trivial, can you find another example, maybe one that resulted in changes to derivation and equations, and can we somehow also fit before/after view of relevant part of paper? maybe via a click or popup if space is too limited?

[Page 1 now shows the exchange that changed Eqn. 21: the question about writing Delta^2 as d_r^2 + d_z^2 - (cN)^2, the answer (Pythagoras), the paper and sympy edits, PASS G16 and the compile; below it, crops of Eqn. 21 from the paper before and after (commit 9283d09). It fits on the slide, so no popup.]

what does Edit notes/sympy/… (new check) mean?

[Claude edited the symbolic-check script notes/sympy/streamline_gravity.py, adding a check (in G16) that the new Delta^2 form equals the old one; the path was truncated to fit. The line now reads "● Edit sympy script: verify new Δ²".]

i think that means that CC uses sympy to confirm the math, how does that work?

[Explained: sympy (Python computer algebra) rebuilds each result from the setup symbolically, eg G16 builds d0, t-hat and d1 as symbolic vectors, solves Delta.t = 0 for the near point, forms Delta^2, subtracts the paper's formula and asserts simplify() gives exactly 0, printing PASS; the check is an independent re-derivation from the definitions, not a re-read of the algebra; it confirms the algebra, not the physics assumptions, which the hand check and the user cover.]

should this be spelled out briefly in page 3?

[Yes. Page 3's first bullet now has two sub-bullets: "sympy (computer algebra) re-derives each equation independently, & must match exactly" and "sympy checks the algebra, I check the physics" (replacing "nothing goes into code until it passes sympy", now in the speaker notes with a fuller explanation). Bullets set to 19 pt so the slide fits.]

sympy = python package that lets one do symbolic math? should deck explain briefly how eqn 7 was derived, and solved?

[Yes, sympy is a Python package for symbolic math. Eqn. 7 is the vertical oscillator zddot = -nu_a^2 z (planet potential expanded to second order in z, evaluated at r = a; solved as z = a i sin(lambda - Omega), with the node regressing at Omega_a - nu_a). Recommended against a derivation in the 10-minute talk: one line in the speaker notes or a backup slide for Q&A. Proposed wording "sympy (a Python package for symbolic math)". Awaiting go-ahead.]

what derived Eqn 7 from newton's laws, cc or sympy or both?

[Both, in different roles. Prompt 01 asked Claude to derive the vertical motion; Claude did the derivation (Newton's second law with the oblate planet's potential, expand Phi to second order in z, evaluate at r = a, choose what to neglect, identify i and Omega). Claude then wrote notes/sympy/vertical_motion.py, in which sympy redoes the calculus mechanically (series expansion of Phi in z, d2Phi/dz2, differentiating the solution) and V6 sizes the neglected O(e i) term. Sympy did not choose the approximations.]

add that three-way split to page 3

[Page 3's first bullet is now "who does what:" with "CC: physical reasoning, approximations & derivation", "sympy (Python symbolic math): independently redoes the algebra & calculus" and "me: are the setup & approximations physically right?". To fit, the epi_int_lite review is now one bullet.]

what does Bash uv run python notes/sympy/… PASS G16: = d_r^2 + d_z^2 - (N c)^2 mean?

["Bash" = Claude ran a shell command; uv run python notes/sympy/streamline_gravity.py runs the sympy check script (uv manages the Python environment); "PASS G16: …" is its output, the check that Delta^2 = d_r^2 + d_z^2 - (Nc)^2 matches. Page 1's action lines now read in plain English: edits Eqn. 21 in the paper; adds a sympy check of the new Delta^2; runs the sympy checks: PASS G16: …; recompiles the paper.]

runs the sympy checks: PASS G16: Δ² = d_r² + d_z² − (Nc)² = you create sympy code to confirm mathematically that the revised eqn is correct?

[Yes: Claude wrote the sympy check (rebuild Delta^2 from the geometry, subtract the revised formula, simplify to exactly 0) and ran it. The slide now says "● writes sympy code that re-derives Δ² from the geometry, runs it:" above the PASS line.]

change to "Claude Code (CC) is an AI agent <CR> that works inside my project folder" and to "it reads my files, writes & runs code, <CR> and I approve each step" and revise "so its work fits my project" to sound like me, maybe something like 'so its outputs align with the project's goals' and to "and will  draft my next NASA proposal eventually" change “why?” to “why this and not that?”

[All made on page 1, with real line breaks inside the two bullets (make_pptx.py now turns \n into a line break): "…so its outputs align with the project’s goals", "it also manages new code development, the paper & this deck, and will draft my next NASA proposal eventually", and ask “why this and not that?”.]

change to "CC reads my code, notes, papers & data <CR> before it answers, <CR>so its outputs align with the project’s goals" and "I stay the scientist: • I pick the physics, steer CC out of rabbit holes, ask “why this and not that?”, and referee everything" should be expand to via anaology, ie CC=graduate student , i am the graduate student' research advisor, the student is very good and math and code and data analysis, but still needs to be prompted to derive this way and not that way, and to test test test, what you think?

[Line breaks made. The analogy is a good one for this audience, with two caveats worth a phrase: the student forgets between sessions unless it is written in the files, and it rarely says "I don't know". Proposed wording; page 1 is full, so it replaces the current "I stay the scientist" bullet rather than adding to it. Awaiting go-ahead.]

yes

[Page 1's last bullet is now "I’m the advisor, CC is my grad student:" with "very good at math, code & data analysis", "but needs direction: derive it this way, not that way, and test, test, test" and "it never tires, but forgets between sessions & rarely says “I don’t know”"; speaker notes expanded. To fit, the bullets start higher and are 19 pt.]

should Joe Hahn, Space Science Institute be another color dark blue? add date to 1st page?

[Both done: the name line and a right-aligned "SSI AI Working Group, Oct 28, 2026" on the same line, in dark blue (#1F3864).]

change to "it also manages all new code development, CR the paper & this deck, CR and will eventually draft my next funding proposal" and to "CC is very good at math" and to "derive it this way not that way CR and test, test, test and i'm not sure what the intent of "it never tires, but forgets between sessions & rarely says “I don’t know”" is

[Wording and line breaks made. The last sub-bullet marks where the grad-student analogy breaks: unlike a student, CC does not remember last week unless it is written in the project files, and it does not flag its own uncertainty. Offered a clearer rewording or dropping it to the speaker notes.]

maybe change to 'like a student it can drift, and rarely says i dont know' . also shink the left & right margins slightly so the middle gap is wider. and slide all of page 1's contnet up somewhat

[Done: last sub-bullet is now "like a student it can drift, and rarely says “I don’t know”"; page 1's content starts higher, the left and right margins are 0.2 in, the right column moved right for a wider middle gap, and the exchange box is shorter with the Eqn. 21 crops moved up under it.]

make left & right margin slightly wider so middle gap is slight ly narrower. did page 1 font get smaller?

[Margins now 0.3 in, gap narrower. Yes: page 1 bullets went from 20 to 19 pt (sub-bullets 16 to 15) when the grad-student bullet was added, to fit; nothing changed since.]

like a student it can drift, CR and rarely says “I don’t know”.

[Line break added.]

in page 2 can "Literature search & citations: • CC reads the papers in my project folder & finds what I need • eg the wire-gravit" be put in a left-aligned box, and nother box just to right contains the figure?

[Done: a framed box at upper left holds the header (26 pt) and the first bullet; a framed box to its right holds Fig. 3 cropped from papers/hahn13.pdf (page 6), captioned "Fig. 3 of Hahn & Spitale (2013)"; the citation bullets and the bib entry follow below.]

make upper left box be invisible and yellow background, and put new paper's figure in upper right

[Done: no frame on the upper left (header and first bullet sit on the yellow page); the upper-right box now holds the new paper's Fig. 1, captioned "the new paper’s Fig. 1, redrawn by CC from HS13".]

drop caption from right box, make box lines invisible, background yellow, and make figure as big as possible to fille the available space

[Done: caption and box removed; Fig. 1 re-rendered from paper/figures/wire_planar.pdf with its white made transparent (assets/fig1_transp.png), so it sits on the yellow page, enlarged to fill the upper right above the next bullet.]

i think we want to say but AI can invent citations. and what does "so CC checks every one against NASA ADS before it goes in the paper:" mean does CC query ADS to get the article's latex reference?

[Changed to "but AI can invent citations" and "so CC confirms each new citation in NASA ADS & records its ADS ID:". The session logs show only the ADS search endpoint (api.adsabs.harvard.edu/v1/search/query) was used: CC searched ADS for each new reference, confirmed it exists and its authors, year, journal, volume and page match, and recorded the bibcode in the "SciX check" line; it did not fetch the BibTeX from ADS's export endpoint. Offered to switch to pulling BibTeX straight from ADS.]

delete "& records its ADS ID:" and expand "so CC confirms each new citation in NASA ADS" to say CC confirmed  authors, year, journal, volume and page

[Now reads "so CC confirms each new citation’s authors, year, journal, volume & page in NASA ADS"; spacing adjusted so the bib box no longer touches it.]

shorten to "CC confirms each citation’s authors, year, journal, volume & page in NASA ADS" so it spans 1 line . and I would like to replace "any field: make it verify against the authoritative database" [and the bib entry] with a discussion of what CC knows internally (whats the word for that, it was trained on?) versus understanding extract from context docs. Do you recall our discussion of whether toomre grav instability was source of epi_int's zigzag, or whether viscous instability or viscous overstability or numeric instability? would it be worth discussing how much of granular astrophysics CC knows versus that gleaned from context papers or websearch?

[Shortened; it fits on one line. Terms: training (or pretraining) knowledge vs context. The zigzag discussion predates this session's memory, but the record shows it: prompts 07 and 09, notes/phase0b/nominal.md and Appendix C; the four viscous-overstability papers cited there (Schmit & Tscharnuter 1995, Schmidt et al. 2001, Latter & Ogilvie 2008, Lehmann et al. 2017) were not in papers/, so they came from training knowledge, and all four checked out in ADS. Recommended the replacement, with proposed text; awaiting go-ahead.]

yes make the swap

[Page 2's "any field" bullet and the bib box are replaced by "what CC knows from training vs what it reads in my folder:" with the three sub-bullets proposed; the bib entry's SciX line and the four overstability citations moved to the speaker notes.]

from training: the published literature means that CC is trained on widely referenced papers for toomre & viscous instability. but not my lightly reference hahn13 paper, correct?

[Mostly: Claude cannot reliably say which papers were in its training data (hahn13 is on arXiv, so it may have been), but widely discussed results are absorbed many times over while lightly cited papers are remembered thinly, if at all; putting hahn13 in papers/ made the question moot. Proposed rewording.]

yes make that change

[Page 2 now reads "from training: well-known results & widely cited papers, eg …" and "from my folder: my own lightly cited papers, code & notation, …".]

drop like equation numbers & signs , put "so the specifics come from the papers themselves" on its own line, and revise specifics come from the papers themselves to mean that the user copies the papers to project folder because CC gets its context from files

[Now "training knowledge is broad but fuzzy on details," then on its own line "so I copy the key papers into my project folder, where CC reads them".]

use bold for page titles

[Page headers on pages 2–5 are now bold; page 1's title keeps only "Claude Code" in bold.]

change "from my folder: my own lightly cited papers, code & notation, the specifics it would otherwise guess at" to "what's in  my project folder: my papers, code & notation, the specifics it would otherwise guess at" . also I think we documented a plan at the start of this project, should that be mentioned here?

[Changed. PLAN.md (prompt 03) holds the goal, phases and ground rules (record the prompts; verify, don't trust; human gates; git history as evidence; track time). Recommended mentioning it, on page 1 or 5 rather than page 2, with proposed wording; awaiting go-ahead.]

yes

[Page 5's closing advice now reads "advice: start with a written plan, insist on verification, keep a record".]

change to "from training: well-known results & widely cited papers, CR eg the viscous overstability papers it cited when we chased a numerical instability, CR none of which were in my folder; all four checked out in ADS" and to what’s in my project folder: my papers, code & notation, code output, notes , all the specifics that CC will need to advance the resaerch

[Both made on page 2, with the two line breaks in the first.]

"what’s in my project folder: my papers, code & notation, code output, notes, CR all the specifics that CC will need to advance the research"

[Line break added.]

does CC: physical reasoning, approximations & derivation mean that CC can solve PDE's ? or is it somehow parroting the derivation from training knowledge which include texts on such

[Neither extreme: it applies techniques learned from training to the problem at hand. Eqn. 7's result is textbook (the paper cites Borderies-Rappaport & Longaretti 1994 Eqn. 17), so close to recall adapted to the project's notation; the 3D near-point geometry (Eqns. 14–22) and the |d1| <= s Delta proof are new combinations of standard steps, not copied. It can solve standard ODEs/PDEs analytically and write numerical solvers, but can also produce plausible wrong answers on hard or novel ones, hence sympy and tests; Claude cannot introspect which steps were recall.]

can we replace the fig on page 5 with a screengrab of "Mostly recall, adapted to the project: …" and say that was obtained vie recall, and another screengrab of "New combinations of standard steps: …" that was AI applying 3D near-point geometry ? and delete "any field: let an independent check referee the algebra" [and the sympy PASS box]?
[after a clarifying question:] replace page 3's current figure with two screengrabs that focus on the relevant text

[Page 3: Fig. 1 replaced by two framed boxes holding Claude's two paragraphs, labeled "CC on how it derived Eqn. 7, recall:" and "…and the 3D near-point geometry, new work:"; the "any field" bullet and the PASS box deleted; bullets back to 20 pt.]

can the boxes on page 3 show screengrab of actualy paper content versus a paraphrase, with caption saying the ddot(z) eqn is solved via recall while the sequence of equations in 2nd box was obtained via stepwise resonaning about the geometry with reminder saying relevant figure is on previous page

[Page 3's right column now holds two framed crops of the paper PDF: Eqns. 7–8 with their text (page 2), captioned "Eqn. 7: solved via recall (a textbook result)", and Eqns. 19–21 with their text (page 3), captioned "Eqns. 19–21: stepwise reasoning about the geometry (figure on previous page)".]

i'm note sure what "me: are the setup physically right' means. i do think that bullet says that I confirm the equations by hand, suggest alternative derivations

[Now "me: I check the equations by hand & suggest alternative derivations".]

"me: I check the equations by hand & CR suggest alternative derivations" and i think we want to revise "sympy (Python symbolic math): independently redoes the algebra & calculus" to it says that CC calls sympy to confirm the algebra & calculus". does sympy also redo the math?

[Yes: Claude writes the sympy script (definitions, what to solve and compare); sympy, an ordinary program, carries out the algebra and calculus itself, exactly. Now "sympy (Python symbolic math): CC calls sympy, which redoes the algebra & calculus to confirm them" and "me: I check the equations by hand & / suggest alternative derivations".]

what does  redoes the algebra & calculus mean?

[Explained with Eqn. 7 (series expansion, d2Phi/dz2, differentiate the solution) and Eqn. 21 (build vectors, solve, dot product, simplify to zero); proposed clearer wording.]

yes

[Now "sympy (Python symbolic math): CC calls sympy, which re-derives each equation from its definitions & confirms the result".]

lets revise "eg I asked: “d₁ ≪ Δ isn’t obvious to me, how do we prove that transparently?” • CC found a short proof, put it in the paper, & added sympy checks" to illustrate me giving CC a small but non trivial task, with CC reasoning through it. Maybe replace "d₁ ≪ Δ isn’t obvious to me, how do we prove that transparently?" with something like "I suspect d₁ ≪ Δ, try to derive a simple proof ?

[Now "eg I asked CC: “I suspect d₁ ≪ Δ, try to derive a simple proof”", with sub-bullets "CC reasoned it through: |d₁| ≤ sΔ, where s ≪ 1 is the streamline’s slope" and "then put the proof in the paper & added sympy checks". Noted that the quote is now a paraphrase of the actual prompt.]

CC reasoned it through: |d₁| ≤ sΔ, <CR> where s ≪ 1 is the streamline’s slope and to "CC then put the proof in the paper & <CR> added sympy checks" and to "it found a few small bugs and two approximations worth improving"

[All three made on page 3; the epi_int_lite bullet now ends ": it found a few small bugs and two approximations worth improving" (issues 1, 4, 6 and 7, 9 of notes/epi_int_lite_issues.md).]

change to "worth improving upon" . and change page 4 title to something like "CC to manage numerical simuations that are tested against theory'

[Page 3 now ends "…two approximations worth improving upon"; page 4's title is "CC manages numerical simulations, tested against theory:" at 26 pt so it fits on one line.]

should "CC ported my old Python 2 code, then made it ~36× faster" state the python version that it got ported to? and what does "(compiled & parallel Python)" mean, is that more than a parallelization library?

[Yes, Python 3.14. Numba is a just-in-time compiler (Python loops to machine code; about 6x alone for the gravity kernel) that can also parallelize prange loops across cores (about 50x for the kernel on 10 threads, 36x per step).]

yes make both changes

[Page 4 now reads "CC ported my old Python 2 code to Python 3.14, then made it ~36× faster" and "using Numba, which compiles Python to machine code & runs it on all cores; identical results".]

we dont want to use the word parallelization? what does machine code mean?

[Parallelization is fine for this audience. Machine code: the CPU's own instructions; Fortran and C are compiled to it, plain Python is interpreted line by line; Numba compiles selected Python functions to it. Proposed "using Numba, which compiles Python like Fortran or C & parallelizes it across all CPU cores; identical results".]

yes make that change

[Done on page 4.]

any field: trust comes from known answers is vague. i think we want to say that trust in the derivations and code is earned by comparing simulations to well understood theory

[Page 4's last bullet is now "trust in the derivations & code is earned by comparing simulations to well-understood theory".]

last bullet can span the entire width of that page

[Page 4's "trust in the derivations & code…" bullet is now its own full-width text box below the plot.]

 comparing simula0ons to <CR> well- understood theory

[Line break before "well-understood theory".]

change "Wri-ng the paper:" to "Wri-ng the science:" . and i think "figure scripts" means that CC manages the code that generates the figures? regarding "CC learned my wriNng style from my last two papers" i think we also created a style guide which then resulted in writing that was even more in my voice, correct?

[Yes to both: CC writes and maintains the Python that makes each figure (eg paper/figures/wire_planar.py); notes/writing_style.md was drafted from hahn13 and hahn25 and grew a section of rules from the author's own edits, which CC reads before drafting. Page 5: title "Writing the science:", "CC manages the LaTeX/AASTeX, bibliography & the code that makes the figures", "CC drafted a style guide from my last two papers, refined it with my edits, & writes in my voice"; the paper image is a bit smaller so the bullets fit.]

make "advice: start with a written plan, insist on verification, keep a record" 1 line, and i think "written plan" should be about a written plan for the project. and what does keep a record mean, how is that useful?

[Now a one-line, full-width bullet at the bottom of page 5: "advice: start with a written project plan, insist on verification, keep a record". Explained the record (saved prompts, notes, equation tracker, git history) and why it helps; suggested a more specific wording; awaiting a choice.]

yes

[Advice now "advice: start with a written project plan, insist on verification, save every prompt & change".]
[Shortened to fit one line: "advice: start with a written project plan, verify everything, save every prompt & change".]
[Still wrapped; now "advice: write a project plan, verify everything, save every prompt & change" at 19 pt, which fits on one line.]

change to "CC also reviewed my predecessor code epi_int_lite:"

[Page 3 now reads "CC also reviewed my predecessor code / epi_int_lite: it found a few small bugs and two approximations worth improving upon".]

change to something like advice: start the research project by writing a project plan, verify everything along the way, and save every prompt & change

[Page 5's closing line is now that sentence, two lines at 20 pt across the full width.]

CC drafted a style guide from my last two papers, refined it with my edits, <CR> & writes in my voice

[Line break added.]

i think "who does what: • CC: physical reasoning, approximations & derivation" should be revised to state that I set up the problem to be solved, then CC does the physical and mathematical reasoning

[Page 3's "who does what:" now starts "me: I set up the problem to be solved", then "CC: then does the physical & mathematical reasoning"; bullets start higher and are 19 pt so the page fits.]

make "Joe Hahn, Space Science Institute " dark green

[Name line now dark green (#1E5631); the date line stays dark blue.]

## Prompt

use round numbers, change ~36× to ~35× or ~40×

[Changed ~36× to ~35× on page 4 (447/12.3 = 36.3, so 35 is the nearer round number); rebuilt the pptx and PDF.]

## Prompt

adjust page 5 margins so "after: …so it and its neighbor on either side bracket
θ." is one line

[Page 5: the bullet box now starts at x=0.3 in and is 5.85 in wide, and the paper image is shrunk from 4.8 to 4.5 in tall and moved to x=6.2 in. The "after:" line now fits on one line. Rebuilt the pptx and PDF.]

## Prompt

on page 5, slide image to left side and bullets to right

[Page 5: the paper image is now on the left (x=0.3 in, 4.5 in tall) and the bullets are on the right (x=3.9 in, 5.85 in wide). The "after:" line is still one line, and the advice bullet stays full width at the bottom. Rebuilt the pptx and PDF.]

## Prompt

slide title to right side. zoom in on image and slide it left so that ir runs off the left and upper edges of slide, i want that image's right column to be more visible

[Page 5: the title now starts at x=3.9 in, aligned with the bullets (header() gained an x argument). The paper image is zoomed to 6.6 in tall at (-1.8, -0.6) in, so it runs off the left and top edges. The right column (Eqns 19-22 and Sec. 3.2 "Computing the acceleration") fills the left third of the slide. Rebuilt the pptx and PDF.]

## Prompt

make image 33% bigger and slide right so that its right edge is almost adjacent to text box on its right

[Page 5: the paper image is 33% bigger (6.6 -> 8.8 in tall) at (-3.0, -0.6) in. Its right edge is at about 3.8 in, next to the text at 3.9 in, and it runs off the left, top and bottom edges. It is drawn first, so it sits behind the text. Because the image covered the full-width advice bullet, that bullet moved into the right column below the others, on 3 lines. Rebuilt the pptx and PDF.]

## Prompt

make image 33% bigger again, that will overflow into text on right which is ok since the image will go behind th text, with the text landing on the image's whitespace along its right margin

[Page 5: the paper image is 33% bigger again (8.8 -> 11.7 in tall) at (-4.55, -3.2) in. The right column of the paper (Eqn 21 through the end of Sec. 3.2) ends just left of the text at x=3.9 in. The paper's white right margin (about 0.7 in) lies behind the start of the bullets, and the yellow background resumes to its right. Rebuilt the pptx and PDF.]

## Prompt

slide image down so it top edge of image align with slide's top

[Page 5: the paper image now starts at y=0, so its top edge lines up with the top of the slide. It shows the running head and Eqns 19-22 down to the start of Sec. 3.2, and runs off the left and bottom edges. Rebuilt the pptx and PDF.]

## Prompt

move page 5 text to right a smidge and image to left a smidg so only bullets overlap in the image whitespace

[Page 5: the image moved 0.15 in left (x=-4.7 in) and the title and bullet boxes moved 0.15 in right (x=4.05 in, bullet boxes 5.75 in wide). The top-level bullet dots now sit on the white margin of the paper, and the text after them starts where the yellow resumes. The "after:" line is still one line. Rebuilt the pptx and PDF.]

## Prompt

display deck

[Opened Using_Claude_Code_for_Research.pptx in PowerPoint.]

## Prompt

display pdf

[Opened Using_Claude_Code_for_Research.pdf in Preview.]

## Prompt

pop log of prompts to vscode

[Opened prompts/14-ai-for-research-deck.md in VS Code.]

## Prompt

the audience liked the deck and the presentation. would it be worthwhile writing a LI post about AI for reasearch? For the purpose of demonstrating competency, and letting my linkedin network know that i understand the topic. and with an eye towards drawing traffic and future prosects to my LI profile and my website /Users/joehahn/Dropbox/prog/claude/jmh.datasciences . I've done similar for my other CC projects, they are stored in /Users/joehahn/Dropbox/prog/claude/jmh.datasciences/publications , what you think? should i make this deck pblically available via git? or upload to LI?

[Recommended doing all three: a LinkedIn post with the deck PDF attached as a document carousel; the deck committed to this public repo for a permanent link; and the post draft kept in jmh.datasciences/publications/epi-int-3d/linkedin/, following the convention there. Suggested angle: verification is how you trust AI output. Listed checks to run before going public. Nothing was drafted or published yet.]

## Prompt

draft publications/epi-int-3d/linkedin/00-ai-for-research.md in your website folder , add closing 6th page that notes my scientific affiliation but also that I am an AI consultant/builder at jmh-datasciences.com

[Added page 6, "About me:", to make_pptx.py. It covers the Ph.D. from Notre Dame; planetary scientist at SSI with 30+ papers; AI consultant & builder at JMH DataSciences (12+ years, Oracle AI CoE, using the approved wording from POSITIONING.md); the repo link; and jmh-datasciences.com plus the email in green. Rebuilt the deck as a 6-page PDF. A different presentation was open in PowerPoint, so the export now names the deck explicitly instead of using "active presentation". Drafted jmh.datasciences/publications/epi-int-3d/linkedin/00-ai-for-research.md, about 1,590 characters with no em dashes. It has a hook from the talk, the grad-student framing, four numbered checks, the 35x speedup, and a bridge to business AI ending on a question. The links go in the first comment. Also added an index entry to publications/README.md.]

## Prompt

change "Joe Hahn, Ph.D. in physics, University of Notre Dame / • planetary scientist at the Space Science Institute / • 30+ peer-reviewed astronomy & planetary-dynamics papers / • planetary rings: the code & paper in this talk are my current research" to Joe Hahn, planetary scientist at the Space Science Institute, link to my ssi page https://www.spacescience.org/bio.php?emp=JHAHN, child bullets refer to my research interets: planetary dynamics, planet formation, circustellar disks, zodiacal cloud

[Page 6: the first bullet is now "Joe Hahn, planetary scientist at the Space Science Institute". "Space Science Institute" is a hyperlink to the SSI bio page; runs() now takes an optional url, and the exported PDF keeps the link. Its child bullet is "research interests: planetary dynamics, planet formation, circumstellar disks & the zodiacal cloud". The Ph.D., 30+ papers and planetary-rings lines were removed. Rebuilt the pptx and PDF.]

## Prompt

add planetary rings and computational dynamics, change "12+ years building production AI & machine learning, eight of them at Oracle, the last three in its AI Center of Excellence" to "12+ years building production AI & machine learning, eight at Oracle, with  last three in its AI Center of Excellence" revise nect bullet to indicat that i now run my own independat AI consulting shop and change "& add senior hands to AI teams that have more work than people" to something like "& provide senior hands to AI teams that have more work than people and change "paper" to "draft paper"

[Page 6 changes:
- Research interests now begin with planetary rings and end with computational dynamics.
- The Oracle line reads "eight at Oracle, with the last three in its AI Center of Excellence" ("the" added).
- The next bullet reads "I now run my own independent AI consulting shop: I build AI tools that do work a machine should do, & provide senior hands to AI teams that have more work than people".
- "this project's code, draft paper & every prompt I gave CC".
Rebuilt the pptx and PDF.]

## Prompt

change to "also an AI consultant & builder at JMH DataSciences" and make "I build AI tools that do work a machine should do,
& provide senior hands to AI teams that have more work than people" sub bullets

[Page 6: the bullet now reads "also an AI consultant & builder at JMH DataSciences". The two lines "I build AI tools that do work a machine should do" and "& provide senior hands to AI teams that have more work than people" are now level-2 sub-bullets under "I now run my own independent AI consulting shop:", at 14 pt. Rebuilt the pptx and PDF.]

## Prompt

is 'research interests: planetary rings, planetary dynamics, planet formation,
circumstellar disks, the zodiacal cloud & computational dynamics' ordered correctly?

[Mostly. Rings first is right because the talk is about them. Suggested grouping the topics by setting (Solar System, then other stars) and putting computational dynamics last as the method: "planetary rings, planetary dynamics & the zodiacal cloud; planet formation & circumstellar disks; computational dynamics". Not edited yet.]

## Prompt

yes

[Page 6: the research interests now read "planetary rings, planetary dynamics & the zodiacal cloud; planet formation & circumstellar disks; computational dynamics". Rebuilt the pptx and PDF.]

## Prompt

what's next push deck to git?

[Proposed a commit plan and asked for the go-ahead:
- Un-ignore only the talk's PDF, make_pptx.py and the 7 assets it uses. Keep the pptx (speaker notes), plan.md, deck/ and the unused assets private, including hs13_fig3.png, an AAS figure.
- Add a README pointer.
- Flagged that prompt 14's first entry names the SSI organizer and mentions an SSI-internal arrangement; suggested redacting both before pushing.
- Then push, together with the 2 commits not yet pushed.]
