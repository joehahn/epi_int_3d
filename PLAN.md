# epi_int_3d — project plan (draft v9)

**Goal:** show how an AI coding agent (Claude Code, "CC") can speed up a real astrophysics
research project. The project extends the 2D epi_int_lite streamline N-body integrator
(Hahn & Spitale 2013; Hahn et al. 2025) to 3D, uses it to simulate spiral bending waves in
planetary rings, and drafts the resulting journal paper.

**The demo is the record of the work.** The README tells the story through the prompts that
drove each phase, what CC produced, and how the results were verified.

---

## Ground rules (apply to every phase)

1. **Record the final prompts.** Each prompt that drives the project is saved in its final,
   cleaned-up form as `prompts/NN-short-name.md`. No drafts, failed attempts or revision
   history. (The raw CC session logs, which hold the full history, stay local and are not
   part of the demo.)
2. **Verify, don't trust.** Every derived equation gets a symbolic check by CC (sympy), and a
   check by hand by you. Your checks can lag behind, so `notes/equations.md` tracks each
   equation's status: *derived → sympy-verified → hand-verified*. An equation may go into the
   code once it is sympy-verified; the paper is not submitted until all are hand-verified.
   Every citation CC introduces is checked against SciX (successor to NASA ADS, scixplorer.org).
3. **Human gates.** Each phase ends with a checkpoint where you review and approve before the
   next phase starts.
4. **Git history as evidence.** One branch or set of commits per phase, with readable messages,
   so a reader can follow the project's evolution.
5. **Track time and cost.** For each phase, record:
   - **CC time:** from each prompt to CC's final reply.
   - **Human time in session:** from CC's reply to your next prompt (reading, reviewing,
     thinking). Gaps longer than ~30 min count as "away" and are excluded.
   - **Human time off session:** e.g. hand-verifying equations, logged by you in
     `notes/human_time.md` (date, phase, task, hours).
   - **Token / $ cost** per phase; the $ cost is the API cost of those tokens (the logs'
     `totalCostUSD`), whatever plan actually paid for them.
   - Your estimate of how long the same work would have taken without AI.

   CC session logs (`~/.claude/projects/<project>/*.jsonl`) timestamp every message, so the
   first two, plus tokens and cost, are extracted by `tools/session_time.py` (phases are
   assigned by date in `notes/phases.csv`; a session's cost is known once it ends).

---

## Phase 0 — Setup and baseline

- GitHub repo `joehahn/epi_int_3d` (public since 2026-10-20), **starting fresh**, with no git
  connection to epi_int_lite. Start from the source on epi_int_lite's **`master` branch**
  (the hahn25 code); the Fourier-satellite forcing on the `lindblad_resonance` branch is
  ported separately in Phase 3 (see 2c). Cite
  https://github.com/joehahn/epi_int_lite and https://github.com/joehahn/epi_int (IDL) as
  origins in the README.
- **Port to modern Python.** Upstream is Python 2.7 (`print` statements, `execfile`). Target
  the newest CPython that Numba supports when the phase starts (likely 3.14), with NumPy 2.x.
  The port itself is a straight translation, with no speed changes yet, so it can be checked
  directly against the original. (When the `lindblad_resonance` forcing is ported later, its
  `scipy.misc.derivative` call, removed from modern SciPy, is replaced with the analytic
  Laplace-coefficient derivative.)
- Rerun a few existing `tests/` (e.g. `unperturbed`, `libration`, `spreading`) and confirm they
  reproduce the epi_int_lite results. **Decided:** the port is checked against the figures in
  epi_int_lite's test notebooks and against analytic expectations (e.g. unperturbed orbits,
  libration frequencies, viscous spreading rates), not against reruns of the Python 2 code,
  which can't run on this Mac without an x86 Python 2.7 install. These become the **2D
  regression suite**, and the straight port becomes the **reference implementation** that the
  faster code is checked against (to ~1e-12).
- Python environment managed by uv (`pyproject.toml`, `uv.lock`, project-local `.venv/`),
  Python 3.14 (set up early, with the vertical-motion write-up); Phase 0 adds Numba etc.
- Add `CLAUDE.md` (conventions, units, how to run tests), `prompts/`, `notes/`, and
  `tools/session_time.py` for time tracking.
- `papers/` holds the reference PDFs by others, `joes_papers/` the author's own (and the style guide
  drawn from them). Only openly licensed ones are committed (hahn25 is CC BY 4.0); the rest are
  git-ignored and listed in each folder's `README.md` with DOI and arXiv links.

**Deliverable:** a modern-Python epi_int_lite that reproduces its published results.

**Status (2026-10-17, prompt 05):** done. `source/` is the Python 3.14 port (changes in
`notes/port_changes.md`); it reproduces the notebooks' Python 2 results to 6–10 significant
figures (`notes/phase0/port_comparison.md`). The 2D regression suite is `uv run pytest`
(unperturbed, spreading, libration; golden files in `tests/golden/`). **Gate: approved 2026-10-17**
(the user compared each test notebook with its epi_int_lite original).

## Phase 0b — Performance

Done once, before the 3D work, so the 3D kernels are written in the fast form from the start.

**Target problem size** (for an m = 2–3 spiral bending wave):
- About 20 streamlines per wavelength and about 20 wave crests. The wavelength shrinks as 1/x
  with distance x from resonance, so the number of crests grows as x², and the shortest
  wavelength is at the far end of the simulated region.
- Nr ≈ 400 streamlines if their radial spacing adapts to the local wavelength, or ≈ 800 if the
  spacing is uniform (set by the shortest wavelength). Nt ≈ 20 particles per streamline, so
  about 8,000–16,000 particles in total.
- Gravity is ~Nr² × Nt ≈ 3–13 million streamline–particle interactions per kick. Rough
  estimate: a few hundredths of a second per kick on one core with Numba, a few thousandths
  on a 16–32-core VM. At ~10⁵–10⁶ steps per run, that is roughly minutes to an hour per run
  on a cloud VM, versus days for the current NumPy code. *(To be replaced by measured
  numbers from the benchmark.)*
- Adaptive (non-uniform) streamline spacing halves Nr, so it's worth supporting.

**Where the time goes now** (from reading `source/helper_fns.py`; to be confirmed by profiling):
- `ring_gravity` loops in Python over every other streamline, so the cost is O(Nr² × Nt)
  (Nr = streamlines, Nt = particles per streamline). Each pass builds about ten temporary
  arrays (`np.concatenate` shifts, Lagrange interpolation, a `sqrt`). This is the main cost,
  and bending waves will need many more streamlines than the narrow ringlets in hahn25.
- `velocity_kick` is called every step and re-sorts every streamline by longitude, recomputes
  orbit elements to get λ, and converts coordinates back and forth.
- Every output is kept in memory as a list and written only at the end, so long runs grow in
  memory and can't be restarted.

**Plan, in order of expected payoff:**
1. **Benchmark first.** A fixed set of problem sizes (e.g. Nr = 10, 100, 1000), profiled with
   cProfile/line_profiler, and timing recorded for the reference implementation. Every change
   below is measured against this and must match the reference results to ~1e-12.
2. **Compiled kernels with Numba.** Rewrite the gravity, pressure and viscosity sums as single
   `@njit` loops: no temporary arrays, the neighbor interpolation done inline. Expected 10–100×.
   No `fastmath`, so results stay reproducible.
3. **Multicore.** `prange` over target streamlines in those kernels, using all cores of one
   machine. Expected another ~N_cores×.
4. **Cheaper bookkeeping.** Keep particles sorted between steps (they barely move in longitude
   per step, so an insertion sort is nearly free) and avoid redundant coordinate conversions.
5. **Algorithmic speedup (research option).** Gravity from distant streamlines is smooth, so it
   might be computed with a cutoff plus a smooth-tail correction, or as an FFT convolution in
   streamline index, turning O(Nr²) into ~O(Nr log Nr). Only adopted if it passes the accuracy
   tests; otherwise left as future work.
6. **Output and restarts.** Write outputs as they happen (HDF5 or chunked `.npz`) and support
   checkpoint/restart, which matters for runs of 10⁵–10⁶ orbits.
7. **Parameter-study throughput.** Phase 5 is many independent runs, which is trivially
   parallel: run them at the same time across cores or cloud machines.

**Considered but not planned:** GPU (JAX/CuPy); Apple-silicon GPU support is still immature
and the arrays are probably too small to benefit. Free-threaded Python is not needed because
Numba releases the GIL. Higher-order symplectic integrators (fewer, larger steps) only if
profiling shows the step count, not the per-step cost, is the limit.

**Deliverable:** the same results as the reference implementation, with a speedup table in
the README (a good demo result in its own right).

**Status (2026-10-18, prompt 06):** step 1 (benchmark) done: `notes/phase0b/benchmark.md`.
Ring gravity is 99.6% of the time at Nr = 1000; s/step ≈ (Nr − 1)(33 μs + 0.018 μs Nr(Nt + 2)) on
one M5 core, so a bending-wave run of 10⁵–10⁶ steps takes 2 h–3 days with the current code. The
driver's streamline-crossing check misfires near θ = ±π for closely spaced eccentric
streamlines; the faster code must fix it. Closely spaced eccentric streamlines with Nt ≈ 20 go
unstable after 50–200 orbits (issue 3 in `notes/epi_int_lite_issues.md`), so the benchmark ringlet
is circular; this must be understood before Phase 5.

**Status (2026-10-18, prompt 07):** epi_int_3d's nominal ringlet, `tests/nominal/`
(`notes/phase0b/nominal.md`): 100 streamlines × 40 particles, width 0.0125 with tapered edges, mass
3×10⁻⁸, e = 10⁻⁴ with a π periapse twist, ν = 7×10⁻¹¹; smooth for 10⁴ orbits: self-gravity straightens
the twist and the body stays nearly aligned, as two sub-ringlets precessing slowly relative to each
other, with eccentric, librating edges. At 300 and 1000 streamlines it goes unstable, in ~190 and ~30 orbits: gravitational, not cured by
more particles per streamline or by ring pressure (Q ≤ 2), probably the cold ringlet's Toomre
instability once the streamlines resolve λ_crit. Ring-pressure equations checked against hahn13
(P1–P7); issues 4–5 logged, minor for the nominal ringlet; a pressure deep dive comes first if the
z ≠ 0 runs use pressure. Next: step 2, the Numba gravity kernel.

**Status (2026-10-19, prompt 08):** steps 2–3 done for gravity (`notes/phase0b/numba.md`):
`fast/kernels.py`, a Numba ring-gravity kernel with `prange` over target streamlines, and
`fast/epi_int_fast.py`, the reference driver with that kernel swapped in. Results are bit-for-bit the
reference's (benchmark states, regression suite, the full nominal run). Kernel ×50 on 10 threads;
whole driver ×7 at Nr = 100, ×32 at Nr = 1000 (14 ms/step at Nt = 20); the nominal run takes 3.3 min
instead of 24. The remaining time is NumPy bookkeeping (viscosity interpolation, coordinate
conversions, sorting), so step 4 is next. The nominal ringlet is now the zero-twist one, run for
5×10⁴ orbits (2.2 viscous timescales; 16 min with the fast driver): smooth, its body circularizes
and its outer edge stays eccentric. At 333 streamlines it crosses after ~145 orbits. Benchmarked at
Nr = 33, 100, 333 (`--ringlet nominal`): ×2.8, ×7.4, ×16.5, final states identical.

**Status (2026-10-19, evening):** the 333-streamline crossings are the m ≈ Nt/2 tangential sawtooth
along each streamline (issue 3) plus, with shear viscosity alone, a code bug that drops bulk viscosity
(issue 6, fixed opt-in in `tools/patches.py`); the viscous overstability theory (O1–O5) and the code's own
linear modes (`notes/numerics/viscous_stability.py`) are recorded. The nominal ringlet is now 1000 × 20
streamlines, width 0.00125, mass 3×10⁻¹¹ (λ_crit/2Δ ≈ 0.09), e = 10⁻⁵, zero twist, viscosity and
pressure off: smooth for 5×10⁴ orbits (2.8 h with the fast driver), self-gravity holding its twist
within 0.38 rad against J2's ~5 rad. The former nominal ringlet is frozen in
`benchmark/nominal_inputs.py` for the benchmark.

**Status (2026-10-20, prompt 10):** step 4 done (`notes/phase0b/driver.md`): `fast/epi_int_3d.py`, epi_int_3d's
own driver, with the bookkeeping (coordinate conversions, drift, sorting, linear density, crossing
check) as Numba kernels in `fast/bookkeeping.py`, bit-for-bit the reference. Only ×1.15–1.5 over the
gravity-kernel driver: at 1000 × 20 gravity is 8.7 of 12.3 ms per step, and the rest is the
transcendentals that bit-for-bit results require. Threads don't help the bookkeeping below ~50 000
particles with Numba's workqueue layer. Next: the vertical drift (Phase 2a) in this driver.

## Phase 1 — Literature review

- CC reads hahn13 and hahn25, plus key papers on epicyclic theory and bending waves:
  Borderies-Rappaport & Longaretti (1994); Shu, Cuzzi & Lissauer (1983); Shu (1984);
  Borderies, Goldreich & Tremaine on inclined/eccentric ringlets; and Cassini-era observations
  of the Mimas 5:3 bending wave. *(You are gathering the PDFs into `papers/`.)*
- **Notation.** Your treatment follows Goldreich & Tremaine (GT), while Shu's papers use a
  different convention: GT's m = −(Shu's m), and there may be other differences (the sign of
  the pattern speed, phase conventions, the definition of the resonance offset x, etc.).
  `notes/conventions.md` will map Shu ↔ GT symbol by symbol. Every equation taken from Shu is
  translated into GT notation before use, and the translation itself is sympy-checked and
  added to `notes/equations.md`.
- Output: `notes/literature.md`, covering what each paper provides, which equations we need,
  and the linear-theory predictions we can test against (bending-wave dispersion relation,
  wavelength vs. distance from resonance, damping length).

**Gate:** you confirm the summaries are accurate and nothing important is missing.

## Phase 2 — Theory

- **2a. Vertical epicyclic motion (largely done).** The first-order, exactly reversible
  extension of hahn25 Eq. A9 to (z, v_z) with new elements (i, Ω). Notation: Ω is the node, and Ω_a, κ_a, ν_a
  are the frequencies at semimajor axis a (as in Borderies-Rappaport & Longaretti 1994). A numerical check has
  already been done (round trip accurate to ~1e-16; errors scale as O(e², i², e·i)), and is
  kept as `notes/numerics/vertical_motion_check.py`. V1–V4 are sympy-verified
  (`notes/sympy/vertical_motion.py`) and written up condensed in the mock paper
  `paper/epi_int_3d.tex` (prompt 04). Remaining: your hand verification.
  **epi_int is treated as a 2D code** (Hahn & Spitale 2013 describe a planar model).
- **2b. Vertical self-gravity of streamlines.** Neither epi_int nor epi_int_lite has this;
  both are planar. The new physics: the vertical (and modified
  radial) acceleration that one inclined/warped streamline exerts on a particle in another.
  In the wire approximation this is roughly A ≈ −2Gλ(Δr, Δz)/(Δr² + Δz²). That needs to be
  derived carefully, including the regime λ_b ≲ A where in-plane and vertical forcing couple
  (hahn25 Sec. 7, item 3).
  **Status (2026-10-20, prompt 12):** derived and sympy-checked (G1–G6,
  `notes/phase2/streamline_gravity.md`, with a 3D version of hahn13 Fig. 3). Exact straight-wire pull
  2GλD/|D|²; to first order in the streamline's slopes, A_r = 2GλΔ_r/ρ², A_z = 2GλΔ_z/ρ²,
  ρ² = Δ_r² + Δ_z², valid at any Δ_z/Δ_r (so it covers the λ_b ≲ A regime). Summed over a sheet it gives
  the bending-wave term −2πGσ|k|Z, and on the lattice −X(φ)Z with P7's X. Taking the wire's direction from
  the velocity is part of the wire model and stays. epi_int_lite issue 7 (λ of the perturbed particle,
  not the perturbing streamline) is fixed in `fast/epi_int_3d.py`, opt-in via `fix_lambda = True`. Next, after
  your review: the kernel (Phase 3 step 2).
  **2026-10-24:** first the 2D kernel the way we want it (epi_int_lite issue 9): inputs `nearpoint`
  (the exact planar wire, (2Gλ/d_r)(1, −tan φ)) and `shape_slope` (tan φ from the shape, needed for
  m ≥ 2), G14. All three gravity fixes (with `fix_lambda`) are now epi_int_3d's default;
  `epi_int_lite_gravity = True` restores the reference's for the bit-for-bit regression, which stays as a
  software check. The physics is validated against theory instead: linear density and bending waves
  (needs 2c's satellite forcing), and meanwhile the m = 1 ringlet tests. Then the z ≠ 0 case.
- **2c. Satellite forcing.** Bending waves are launched by an inclined satellite at a vertical
  resonance. Existing code:
  - epi_int_lite `lindblad_resonance` branch: a "Fourier satellite", i.e. only the resonant
    m-th term of a circular, uninclined satellite's potential (horizontal forces only), with
    Laplace coefficients evaluated once at the ring's mean radius.

  Plan:
  - *Horizontal forcing:* port the Fourier satellite from the `lindblad_resonance` branch
    onto the master-based code (no rederivation).
  - *Vertical forcing:* rederive Shu's (1984) resonant vertical forcing term in GT notation
    and add it as a vertical Fourier satellite (inclined satellite, m-th vertical resonance).
  - *Cross-check:* also implement direct 3D gravity of a point-mass satellite (written new,
    written new). A direct inclined satellite
    contains all terms of the potential, so its response near the vertical resonance should
    match the Fourier version. That gives an independent test of the rederived forcing.
    Production runs use the Fourier satellite (cleaner, and it doesn't also excite nearby
    Lindblad resonances).
  - *Linear theory for testing:* rederive the bending-wave dispersion relation and the
    wavelength vs. x relation (λ_b ∝ 1/x) in GT notation, for comparison in Phase 4.
  - Satellite parameters are generic (not Saturn-specific), chosen to place the resonance
    conveniently within the simulated ring.
- **2d. Vertical treatment of pressure/viscosity.** Decide whether to model ring thickness,
  or keep a zero-thickness sheet with a vertical-velocity-dispersion parameter.
- **2e. Ring boundaries.** Bending waves travel toward the satellite (opposite to density
  waves) and would reflect off the ring's edge. Default approach: **stop each simulation
  before the wave front reaches the edge.** The simulated ring is made wide enough, with a
  margin beyond the ~20 crests of interest, that the wave is established over the measured
  region before the front (moving at the group velocity) arrives; the run length follows from
  the group velocity derived in 2c. *Nice-to-have, only if needed:* adapt epi_int's
  `radiative_bc` so waves can leave the region.

**Gate:** all Phase 2 equations are sympy-verified and listed in `notes/equations.md`; you
review the derivations (hand verification may continue in parallel with Phase 3).

## Phase 3 — Code: epi_int_3d

- Drift step: add (i, Ω) and nodal regression (Phase 2a).
- Kick step: vertical streamline gravity (2b), satellite horizontal + vertical forcing (2c),
  and the 2d terms, all written as Numba kernels (Phase 0b).
- A stop condition: end the run when the wave front nears the ring edge (2e).
- Optional non-uniform streamline spacing, to resolve the shrinking wavelength efficiently.
- Initial conditions, I/O and animation scripts extended for z and v_z.
- Keep a `z ≡ 0` path that reproduces the 2D results exactly.
- Remove the `fast_gravity` approximation; epi_int_3d uses exact streamline gravity only.

**Status (2026-10-20, prompt 11):** the drift step is done (`notes/phase3/vertical_drift.md`): z, v_z in
`fast/epi_int_3d.py`, the vertical drift as the rotation V8, inputs `inclination` and `node`, outputs
`z.npy`, `vz.npy`. The in-plane results stay the reference's bit for bit, inclined or not, and z ≡ 0 stays
exactly 0. Phase 4 tests done: z = 0 regression; single inclined test particles (`tests/inclined/`): i
constant to 2e-13, nodes at Ω_a − ν_a to 3e-11 rad; versus exact 3D J2 integration the nodal rate is
off by −4i² (= cos⁸i − 1, since a comes from h_z) plus O(e²). Next: vertical streamline gravity (2b).

**Status (2026-10-29):** step 2, vertical streamline gravity, is done (`notes/phase3/vertical_gravity.md`):
the 3D wire of the paper's Eqn. 25 as the kernel `ring_gravity_3d`, input `vertical_gravity` (default on);
z = 0 runs unchanged bit for bit. A rigidly inclined narrow ringlet's nodes stay locked (±0.012 rad) where
without it they spread 0.39 rad. Next: the quantitative inclined-ringlet test (nodal libration versus
theory), then satellite forcing (2c).

## Phase 3b — Workflow: style guides, outside review, the ledger (decided 2026-11-08)

**Goal:** finish the project as correctly as possible, with the math, code, tests, figures and paper
in your style and voice, and with as few interventions by you as possible. Modeled on
agent-banking-watch (ABW), adapted to a physics paper. Done in this order, before 3c, so that the
density waves are the first feature built through the whole pipeline.

1. **Style guides, one per activity**, each drawn from your pre-AI work the way
   `joes_papers/writing_style.md` was drawn from your papers, each with numbered rules (D3, C7, …)
   so that reviews and ledger cards can cite them, and each with a "From the author's edits" section
   that turns every correction you make into a rule, so the same correction isn't needed twice:
   - *Derivations* (D): from the appendices of the papers in `joes_papers/`.
   - *Code* (C): from epi_int_lite, i.e. `source/`.
   - *Testing and validation* (T): from the original epi_int_lite notebooks (commit `99a7753`).
   - *Data visualization* (V): from those notebooks' figures and the papers' figures.
   - *Paper* (P): the existing `joes_papers/writing_style.md`, rules numbered.
2. **One author, outside reviewers.** CC (Claude) derives, writes code, writes tests and drafts
   prose, each against its style guide: one author keeps the style consistent, and the author must
   be the one working in the repo (sympy, runs, tests, rendering). Models from other companies only
   review:
   - **Bake-off first.** GPT and Gemini (and Claude) each derive one already-verified result, e.g.
     the vertical drift (`notes/phase3/`), graded on correctness and closeness to your style. If
     another model clearly derives better, it drafts derivations and CC puts them in your style.
     *Done 2026-11-08* (`tools/review.py`, problem `reviews/problems/vertical_motion.md`, replies in
     `reviews/`, $2.55 in all). All three were correct. GPT (`openai/gpt-6-astra`) went furthest,
     solving for the O(e·i) term our V2–V8 drop; Claude (`anthropic/claude-opus-5.5`) used the
     Keplerian form z = r sin i sin(θ − Ω), leaving only an O(e·i·J₂(R_p/a)²) error; Gemini
     (`google/gemini-4-argon`) matched V2–V8. Ranked blind on style, Claude's was chosen (the author:
     "lets use b for now"), so CC stays the author of derivations, and GPT, the strongest on the math,
     is the main math reviewer. Both GPT's and Claude's vertical solutions remove most of the drift's
     O(e·i) error (periodic, at ν_a ± κ_a, not secular); the drift stays as is (the author, 2026-11-08),
     since that error is a fraction of a percent for the near-circular rings simulated here. Revisit for
     inclined eccentric ringlets.
   - **Reviewers:** GPT and Gemini as the main pair (math and code), Grok optionally as a third
     opinion on the paper. Models are picked from OpenRouter's list when the tooling is built, then
     pinned and logged. A script modeled on ABW's `scripts/review.py`; raw reviews kept unedited.
   - **Math:** blind re-derivation (the reviewer gets the problem, not CC's answer), then the two
     are compared. sympy stays the ground truth, and your hand check the final word.
   - **Code:** the diff plus the equation it implements (match, units, signs).
   - **Tests:** would this test fail if the code were wrong? Each reviewer names a plausible bug the
     test would miss.
   - **Paper:** every claim against the result behind it, plus a skeptical-referee read.
     Citations still go through SciX.
   - **Disputes:** where CC rejects or partly accepts a finding, a model not involved in it rules
     (ABW's `tiebreak.py`). Editorial decisions stay yours.
   - *Built 2026-11-08:* `tools/review.py` (`--derive`, `--code`, `--test`, `--paper`, each with
     its style guide; `--context`, `--diff`, `--focus`, `--dry-run`) and `tools/tiebreak.py`, sharing
     `tools/openrouter.py`. CC writes `reviews/<date>-<activity>-verdicts.json` after each review.
3. **Fewer knock-on edits.**
   - *One source for every number:* scripts write the paper's numbers into a LaTeX macro file the
     paper includes, and every figure is made by script, so nothing is hand-copied.
   - *Dependencies:* each ledger item records what it depends on (paragraph → equation → sympy
     check → code → test → figure) and those files' hashes; a change marks everything downstream stale.
   - *Knock-on pass:* after every change CC fixes the stale items and lists each fix with a
     one-line reason, so you can veto it.
   - *Built 2026-11-09:* `tools/chain.py` and the chain `notes/chain.json` (its state in
     `notes/chain_state.json`). Each node is a piece of a file: a table row of `notes/equations.md`, the
     sections of a sympy script, a Python function, a LaTeX equation or section, or a paragraph found by a
     phrase; whitespace and the Hand column don't count as edits. A node is edited, stale (something it
     needs changed), new, missing or ok; `--run` reruns the checks (sympy, pytest) of nodes that aren't ok.
     The chain covers every equation the paper rests on, 198 nodes, all confirmed 2026-11-10: the 57
     tracked equations (V, H, P, O, G, B; each needs the earlier ones it builds on) and their sympy checks;
     the drift, ring-gravity and satellite code and the test helpers; `test_vertical.py`, `test_kernels.py`;
     the inclined, vertical_resonance and bending_waves notebooks; the three figure scripts; every labeled
     equation and figure of the paper and 15 of its claims (paragraphs). Not in it yet: the paper's numbers
     (item 3 will make them nodes), the 2D regression suite, the explore/ runs behind Appendix B.
     `notes/sympy/viscous_overstability.py` prints True/False instead of asserting; its check command fails
     on a False.
   - *Built 2026-11-10:* `tools/paper_numbers.py` computes each quoted number from its run and writes
     `paper/numbers/<test>.tex` (`\newcommand`s, each with its source), which the paper `\input`s; each file is a
     chain node, so a rerun that changes a number marks the claims quoting it stale. First: Section 4.2 and
     Fig. vres (16 numbers, from `tests/vertical_resonance`); the PDF's text is unchanged. Still typed in:
     Appendices A and B, whose numbers come from git-ignored `benchmark/runs/` and `explore/` runs; they move
     to macros when those runs are scripted. The averaging start, 10^4 orbits, is set in both the figure
     script and `paper_numbers.py`.
4. **The ledger.** A claude.ai artifact dashboard, as in ABW, with its own database that the page
   and CC both read and write. Cards for each paper paragraph, equation, code unit, test and figure,
   linked along the dependency chain, so clicking a claim shows everything behind it and its status
   (drafted → sympy-checked → reviewed → hand-verified, or stale). Each card links to the style
   guide rules it is checked against and shows the reviewers' findings with CC's verdict, the
   tiebreak ruling, and a box for your note. `notes/equations.md` folds into it. **Costs** are
   broken down by activity (derive, code, test, write, review per model, tiebreak) plus your time
   from `notes/human_time.md`; outside reviews are exact, and CC's split within a session is an
   estimate (an extension of `tools/session_time.py`).
   - *Built 2026-11-10:* the ledger, https://claude.ai/artifact/L54vp4Q4cLf6Cm99Mu41JV (a claude.ai
     dashboard, private until shared). `tools/ledger.py` exports `notes/ledger/nodes.csv` (the chain, statuses,
     hand checks), `findings.csv` (verdicts, tiebreaks) and `costs.csv` (outside reviews per call, Claude
     Code per phase, Joe's hours); Claude uploads them. The page shows the chain by kind and status, a
     searchable list of items, and for the picked item everything it rests on, what rests on it, and the
     findings on its file; then the reviews and the costs. Notes go in as comments on the page. Not yet:
     Claude Code's cost split by activity within a phase (only per phase), style-rule links per card.
     `notes/phases.csv` gained phase 3b (from 2026-11-08 17:52 UTC).
   - *Rebuilt 2026-11-10 as the paper ledger,* https://claude.ai/artifact/KpgQvmghBVn7c4nTmFWxLZ, after the
     author's use of the ABW ledger: there he wrote line edits and questions in the note boxes and chose Use or
     Skip on suggested fixes, card by card, but accepted every keep/trim/cut call in one click and never asked
     about the source and reviewer panels. So the page is now the paper itself, one card per section or
     subsection (16), with the math typeset and the number macros filled in; each card lists its open issues
     (chain items not current, lint problems, findings of the latest paper review) with Use/Skip, and a note
     box; "Apply changes" sends the notes to Claude Code, which revises the paper and replies on each card.
     The chain table, the chart of items by kind, the hours chart and the findings table are gone; the chain
     feeds the cards instead, and the costs are one view. The first dashboard,
     https://claude.ai/artifact/L54vp4Q4cLf6Cm99Mu41JV, was deleted.

## Phase 3c — Spiral density waves (decided 2026-11-08)

*Deferred (the author, 2026-11-10): skipped for now; still planned.*

Linear only. The density-wave counterpart of `tests/bending_waves/`, in `tests/density_waves/`:
the Fourier satellite (`fast/satellite.py`) forcing an m ≥ 2 inner Lindblad resonance in a
self-gravitating ring. The motion is planar, so `source/` and `fast/epi_int_3d.py` can be compared.
Theory (derived in GT notation, sympy-checked, in `notes/equations.md`) and the comparisons:
- wavenumber k(r) against the WKB dispersion relation, including the onset at resonance;
- the angular momentum flux the wave carries against the Goldreich & Tremaine (1979) torque;
- the amplitude profile against Shu (1984), with the finite lattice's restoring force treated as B13
  does for bending waves;
- the wave travels outward from the ILR (bending waves travel inward), a quick sign check;
- convergence in the number of streamlines and the timestep.

Also a density wave in a non-gravitating, pressure-supported gas disk (the author, 2026-11-12), compared with
linear theory, as the test of how the code handles ring pressure. This needs the ring-pressure code
(`Q_ring`), which until then stays off (`Q_ring = -1`, CLAUDE.md) because it is untested; that rule is revised
when this test starts. The paper's Section 2.1 already lists pressure among the ring's internal forces.

Energy budget (the author, 2026-11-12): the bending-wave tests report only angular momentum (L_z is conserved
because the Fourier satellite pulls only along z). An energy check would need the driver to accumulate the
satellite's work and the ring's wire potential energy (E - W versus time); do that with the density waves, where
the satellite exerts in-plane forces.

## Phase 4 — Testing

| Test | Expectation |
|---|---|
| Round trip coords ↔ elements | machine precision |
| Single particle vs. exact 3D J₂ integration | errors O(e², i², e·i) |
| All 2D regression tests with z = 0 | identical to Phase 0 |
| Single inclined test particle | nodal regression at Ω_a − ν_a |
| Narrow inclined self-gravitating ringlet | rigid nodal precession / libration (vertical analog of the eccentric-ringlet test) |
| Satellite horizontal forcing | reproduces the `lindblad_resonance` branch / epi_int results for the same setup |
| Vertical Fourier satellite vs. direct inclined satellite | same forced vertical response near the resonance |
| Linear bending wave from a vertical resonance | wavelength vs. x and damping match the GT-notation linear theory (Phase 2c) |
| Wave-front timing | front position vs. time matches the group velocity; runs stop before it reaches the edge |
| Energy / angular momentum | conserved to the expected tolerance |

| Linear density wave from an inner Lindblad resonance (3c) | wavelength vs. x, torque and amplitude match linear theory; `source/` and `fast/epi_int_3d.py` agree |

Whether a simulation agrees with theory is a judgment call made by you and CC together, not a
pass/fail threshold, and not something the outside reviewers rule on (they check only that the
test code computes what it says it does).

**Gate:** all tests pass and the plots look right to you.

## Phase 5 — Science runs

**Scope (revised 2026-11-08): linear waves only**, bending and density. Saturn-specific
comparison with Cassini data is out of scope for this version of the project.

- Linear bending and density waves: validation and convergence (number of streamlines, timestep).
- *Dropped 2026-11-08:* the nonlinear bending-wave parameter study (λ_b ≲ A).
- **Compute: multicore VMs on Google Cloud.**
  - Software environment: the same `uv.lock` used on the Mac, installed on the VM with
    `uv sync --frozen` (optionally baked into a container image), so both run identical code.
  - Each run uses all cores of one VM through Numba's parallel loops. For the parameter study,
    either several runs share one large VM, or several smaller VMs run at once (the choice
    depends on the measured scaling with core count from Phase 0b).
  - Use compute-optimized VMs; consider spot (preemptible) VMs, which are much cheaper and
    are safe to use once checkpoint/restart exists (Phase 0b, item 6).
  - Outputs go to a Cloud Storage bucket; analysis and plots happen on the Mac.
  - CC drives this through the `gcloud` CLI. You authenticate once; **CC asks before creating
    or deleting any cloud resource**, and every VM is stopped when its runs finish. Cost is
    recorded per run for the time/cost table.
  - **A dedicated Google Cloud project** (e.g. `epi-int-3d`), linked to your existing billing
    account. This keeps its costs, permissions and resources separate from your other
    project, makes the per-phase cost easy to report, and lets everything be removed by
    deleting the project at the end. You create the project and link billing; CC does the rest.
  - **Budget: $100 for the whole project**, with Google Cloud budget alerts at 50%, 90% and
    100%. Note that a Google Cloud budget only sends alerts; it does not stop spending by
    itself. As a hard stop, VMs are created with a maximum run duration so they shut down
    automatically, and CC checks for running VMs at the end of every session.
  - Rough cost estimate (to be firmed up after the Phase 0b benchmark): a 32-vCPU
    compute-optimized VM costs roughly $1–2/hr on demand and a few tenths of a $/hr as a spot
    VM. A parameter study of a few dozen ~1-hour runs, plus validation runs, comes to roughly
    $10–30 on spot VMs, and somewhat more on demand. Disks and storage add a few $.

## Phase 6 — Paper draft

- Target journal: **ApJ** (same as hahn25).
- AASTeX draft generated by CC from the notes, derivations, test results and figures, in the
  style of hahn25. You own the science claims and the final text.
- CC checks citations against SciX and checks that every number in the text traces back to a
  run in the repo.

## Phase 7 — Demo write-up (README)

- Narrative organized by phase: the final prompts → what CC produced → how it was verified.
- Lessons learned: where AI accelerated the work, where it was wrong or overconfident, and what
  needed domain expertise.
- A time and cost table per phase, with cost as API cost: what the tokens would cost at API
  prices, labelled as such (decided 2026-10-20). Cloud compute, if any, is a separate line at
  its billed cost.

---

## Decisions so far

- Repo `joehahn/epi_int_3d`, private until ready, public since 2026-10-20 (MIT license); fresh git history; epi_int_lite referenced by URL.
- Target journal: ApJ.
- Science scope: linear bending and density waves only (revised 2026-11-08); no Saturn/Cassini data.
- Workflow (2026-11-08, Phase 3b): CC authors everything against per-activity style guides; GPT and
  Gemini review (Grok optional, paper only); a ledger dashboard tracks every item, its dependencies,
  reviews and costs by activity.
- Port to the newest Python that Numba supports, then optimize for speed (Phase 0b).
- Bending wave: m = 2–3, ~20 crests, ~20 streamlines per wavelength, ~20 particles per streamline.
- Math in Goldreich & Tremaine notation; Shu's results are rederived/translated (GT m = −Shu m).
- Code starts from epi_int_lite's `master` branch; the Fourier-satellite forcing is ported
  from the `lindblad_resonance` branch in Phase 3.
- Satellite horizontal forcing ported; vertical forcing rederived from Shu in GT notation and
  cross-checked against a newly written direct 3D point-mass satellite.
- Simulations run multicore on Google Cloud VMs in a dedicated GCP project (shared billing
  account); budget $100, scoped to that project (alerts + auto-shutdown).
- Ring edges: stop runs before the wave front reaches the edge; `radiative_bc` is a nice-to-have.
- Prompts: final versions only, no history.
- Equations: sympy-verified by CC, hand-verified by you (may lag); no separate numerical
  verification of derivations.
- Time tracking includes human time, from session timestamps plus your off-session log.
- Papers in `papers/`; only openly licensed PDFs are committed, the rest are linked by DOI/arXiv.
- `fast_gravity = False` and `Q_ring = -1` (pressure off) in all runs: the fast-gravity
  approximation was never useful, and the pressure code was never used in hahn25 nor rigorously
  tested. `fast_gravity` is removed from epi_int_3d (Phase 3).
- The README reports Claude Code's cost as API cost (tokens at API prices), not what the plan
  that paid for it cost (2026-10-20).

## Open questions to settle

- **To do, next: rerun everything after the vertical-phase change** (decided 2026-11-12). The vertical phase is now
  measured from the particle's longitude θ, z = a i sin(θ − Ω), rather than from the mean longitude λ
  (`fast/bookkeeping.vertical_drift`, `fast/epi_int_3d.py` `incline` and `vertical_elements`; notes/equations.md V2, V3, V8;
  paper Sections 2.2, 2.3, 3.4). On `tests/inclined` it measured no loss of accuracy against the exact 3D motion: the z
  error is set by the in-plane O(e²) longitude error either way. Nothing has been rerun since: regenerate the vertical
  golden files (`UPDATE_GOLDEN=1 uv run pytest tests/test_vertical.py`, after checking the change is the expected O(ei)),
  run `uv run pytest`, rerun every test with vertical motion (`tests/inclined`, `tests/vertical_resonance`,
  `tests/bending_waves`, and any other that sets an inclination or satellite) with `tools/run_simulations.py`,
  re-execute their notebooks, remake the paper's figures and `tools/paper_numbers.py` macros, and check the paper's quoted
  numbers. The in-plane results and the 2D regression suite are unchanged (the vertical drift doesn't touch r, θ).
- **To do: an outside review of Section 3.1.1's derivation** (asked 2026-11-13). No outside model has reviewed the
  perturbing streamline's geometry (equations G1-G19, `notes/sympy/streamline_gravity.py`); the 2026-11-08 reviews covered
  only the vertical motion and `fast/bookkeeping.py`. Write `reviews/problems/streamline_geometry.md` with the setup and the
  results but not our derivation, then `tools/review.py --derive` on it, asking for the shortest route to those results.
  It costs money: show the `--dry-run` estimate and ask first. Record verdicts, tiebreak the disputed ones.
- **To do, later: review `fast/` against the code style guide** (decided 2026-11-08). `fast/` was
  written in CC's own style before the guide existed. Restyle the parts you read
  (`fast/epi_int_3d.py`, `fast/satellite.py`) so they can be spot-checked at a glance; the Numba
  kernels and bookkeeping probably stay as they are, since the regression suite proves them bit for
  bit. The regression suite makes any restyling safe.
- **Undecided, later: restyle the notebooks Claude wrote** (e.g. `tests/bending_waves/bending_waves.ipynb`)
  against `style/testing.md`, or leave them as they are (2026-11-08).
- **To do, later: restyle the paper's data figures** `paper/figures/vertical_resonance.py` and
  `zigzag.py` against `style/dataviz.md` (2026-11-08).

- **To do, later: fix the tangential zigzag** (decided 2026-10-19, not analyzed further for now). At
  333 streamlines the nominal ringlet's streamlines cross through a particle-scale sawtooth along each
  streamline, azimuthal wavenumber m ≈ Nt/2, probably numerical (the Lagrange interpolation between
  neighbouring streamlines' particles); it needs neither eccentricity, self-gravity nor viscosity.
  Bending-wave runs (closely spaced streamlines, Nt ≈ 20) must not suffer it. See issue 3 in
  `notes/epi_int_lite_issues.md`; candidate fixes: filter the m = Nt/2 component, or a different
  interpolation stencil, after linearizing the interpolation to find the mechanism.
- **To do, later: conserve angular momentum in the ring-gravity kicks** (2026-10-23). The wire forces
  are not reciprocal: a few parts in 10⁶ of the streamlines' torques fail to cancel in a twisted
  ringlet, a secular |ΔL|/L of about 2×10⁻¹⁶ per orbit for the nominal ringlet at 100 streamlines
  (`notes/numerics/angular_momentum.md`). Neither λ′ nor the shape slope fixes it. A fix would compute
  each streamline pair's interaction once and apply equal and opposite torques.
- **To do, later: prune `explore/` outputs** once we are sure they are not needed (decided 2026-10-19).
  `explore/` (~5 GB) and `.venv/` are excluded from Dropbox sync (`com.dropbox.ignored`). Candidates:
  `ladder/`, `nominal/`, `taper/`, `pressure/` (~2.5 GB, the finished parameter search in
  `notes/phase0b/nominal.md`) and `twist0_n333_*` (~1.9 GB, the 333-streamline crash study); their
  launch scripts and notes stay, so any run can be regenerated.
- **Idea (user, 2026-10-19): a free-wave variant of the nominal ringlet**, possibly more interesting
  than the aligned eccentric one. Start the nominal ringlet circular
  with only its innermost streamline eccentric (e ~ 10⁻⁵) and let self-gravity carry the disturbance
  across it as a free (unforced) m = 1 spiral density wave; once z ≠ 0 works, do the same with an
  inclined innermost streamline for a free m = 1 spiral bending wave. A quick WKB reading, to be
  derived and sympy-checked before use: for m = 1, (Ω_p − Ω)² = κ² − 2πGσ|k| gives
  Ω_p − ϖ̇(a) ≈ πGσ|k|/κ > 0, so density waves live where ϖ̇(a) < Ω_p, outside the radius where the
  pattern speed equals the local precession rate; the bending-wave analogue,
  (Ω_p − Ω)² = ν_z² + 2πGσ|k|, gives Ω_p − Ω̇_node(a) ≈ −πGσ|k|/ν_z < 0, also outside its resonance
  since |Ω̇_node| falls outward. With the pattern speed set near the inner edge's rate, both would
  propagate outward. Related: hahn13 Sec. 3.3, the B-ring edge's free m = 1 mode, whose ILR lies
  inside the edge. The ringlet may need to be wider, in wavelengths, to watch the wave travel.
