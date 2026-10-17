# epi_int_3d — project plan (draft v8)

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
     `notes/human_time.md` (date, task, hours).
   - **Token / $ cost** per phase.
   - Your estimate of how long the same work would have taken without AI.

   CC session logs (`~/.claude/projects/<project>/*.jsonl`) timestamp every message, so the
   first two can be extracted by a small script, `tools/session_time.py`, written in Phase 0.

---

## Phase 0 — Setup and baseline

- GitHub repo `joehahn/epi_int_3d` (private now, public later), **starting fresh**, with no git
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
  reproduce the Python 2 / hahn25 results. These become the **2D regression suite**, and the
  straight port becomes the **reference implementation** that the faster code is checked against.
- Add `CLAUDE.md` (conventions, units, how to run tests), `prompts/`, `notes/`, and
  `tools/session_time.py` for time tracking.
- `papers/` holds the reference PDFs. Only openly licensed ones are committed (hahn25 is
  CC BY 4.0); the rest are git-ignored and listed in `papers/README.md` with DOI and arXiv links.

**Deliverable:** a modern-Python epi_int_lite that reproduces its published results.

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
  extension of hahn25 Eq. A9 to (z, v_z) with new elements (i, Ω̂). A numerical check has
  already been done (round trip accurate to ~1e-16; errors scale as O(e², i², e·i)).
  Remaining: sympy verification, entry in `notes/equations.md`, LaTeX write-up.
  Note: epi_int (IDL) already tracks z, v_z, I, N with the same first-order formulas, using
  vertical phase ω̃ + M − N. Ours uses ω̃ + (Ω/κ)M − Ω̂, which makes the vertical phase advance
  at exactly ν; the two agree to first order.
- **2b. Vertical self-gravity of streamlines.** Neither code has this: epi_int's
  `ring_gravity` explicitly assumes z = 0 and sets A_z = 0. The new physics: the vertical (and modified
  radial) acceleration that one inclined/warped streamline exerts on a particle in another.
  In the wire approximation this is roughly A ≈ −2Gλ(Δr, Δz)/(Δr² + Δz²). That needs to be
  derived carefully, including the regime λ_b ≲ A where in-plane and vertical forcing couple
  (hahn25 Sec. 7, item 3).
- **2c. Satellite forcing.** Bending waves are launched by an inclined satellite at a vertical
  resonance. Existing code:
  - epi_int_lite `lindblad_resonance` branch: a "Fourier satellite", i.e. only the resonant
    m-th term of a circular, uninclined satellite's potential (horizontal forces only), with
    Laplace coefficients evaluated once at the ring's mean radius.
  - epi_int (IDL) `sat_grav`: the same Fourier option, plus fully 3D direct satellite gravity
    (point-mass satellites with z, giving A_z).

  Plan:
  - *Horizontal forcing:* port the Fourier satellite from the `lindblad_resonance` branch
    onto the master-based code (no rederivation).
  - *Vertical forcing:* rederive Shu's (1984) resonant vertical forcing term in GT notation
    and add it as a vertical Fourier satellite (inclined satellite, m-th vertical resonance).
  - *Cross-check:* also port epi_int's direct 3D satellite gravity. A direct inclined satellite
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

- Drift step: add (i, Ω̂) and nodal regression (Phase 2a).
- Kick step: vertical streamline gravity (2b), satellite horizontal + vertical forcing (2c),
  and the 2d terms, all written as Numba kernels (Phase 0b).
- A stop condition: end the run when the wave front nears the ring edge (2e).
- Optional non-uniform streamline spacing, to resolve the shrinking wavelength efficiently.
- Initial conditions, I/O and animation scripts extended for z and v_z.
- Keep a `z ≡ 0` path that reproduces the 2D results exactly.

## Phase 4 — Testing

| Test | Expectation |
|---|---|
| Round trip coords ↔ elements | machine precision |
| Single particle vs. exact 3D J₂ integration | errors O(e², i², e·i) |
| All 2D regression tests with z = 0 | identical to Phase 0 |
| Single inclined test particle | nodal regression at Ω − ν |
| Narrow inclined self-gravitating ringlet | rigid nodal precession / libration (vertical analog of the eccentric-ringlet test) |
| Satellite horizontal forcing | reproduces the `lindblad_resonance` branch / epi_int results for the same setup |
| Vertical Fourier satellite vs. direct inclined satellite | same forced vertical response near the resonance |
| Linear bending wave from a vertical resonance | wavelength vs. x and damping match the GT-notation linear theory (Phase 2c) |
| Wave-front timing | front position vs. time matches the group velocity; runs stop before it reaches the edge |
| Energy / angular momentum | conserved to the expected tolerance |

**Gate:** all tests pass and the plots look right to you.

## Phase 5 — Science runs

**Scope (decided): linear validation plus a nonlinear parameter study.** Saturn-specific
comparison with Cassini data is out of scope for this paper.

- Linear bending waves: validation and convergence (number of streamlines, timestep).
- **Nonlinear spiral bending waves** (λ_b ≲ A). The open problem flagged in hahn25: does a
  shrinking bending wave go nonlinear and drive in-plane density waves?
- Parameter study: vary the forcing amplitude (hence wave amplitude A) and the ring surface
  density σ₀ (hence wavelength λ_b), and map where the transition to nonlinearity occurs and
  whether density waves are launched.
- **Compute: multicore VMs on Google Cloud.**
  - Software environment pinned in a container image (or a pinned conda/uv environment), so
    the VM and the Mac run identical code.
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
- A time and cost table per phase.

---

## Decisions so far

- Repo `joehahn/epi_int_3d`, private until ready, then public; fresh git history; epi_int_lite referenced by URL.
- Target journal: ApJ.
- Science scope: linear validation + nonlinear parameter study (no Saturn/Cassini comparison).
- Port to the newest Python that Numba supports, then optimize for speed (Phase 0b).
- Bending wave: m = 2–3, ~20 crests, ~20 streamlines per wavelength, ~20 particles per streamline.
- Math in Goldreich & Tremaine notation; Shu's results are rederived/translated (GT m = −Shu m).
- Code starts from epi_int_lite's `master` branch; the Fourier-satellite forcing is ported
  from the `lindblad_resonance` branch in Phase 3.
- Satellite horizontal forcing ported; vertical forcing rederived from Shu in GT notation and
  cross-checked against epi_int's direct 3D satellite gravity.
- Simulations run multicore on Google Cloud VMs in a dedicated GCP project (shared billing
  account); budget $100, scoped to that project (alerts + auto-shutdown).
- Ring edges: stop runs before the wave front reaches the edge; `radiative_bc` is a nice-to-have.
- Prompts: final versions only, no history.
- Equations: sympy-verified by CC, hand-verified by you (may lag); no separate numerical
  verification of derivations.
- Time tracking includes human time, from session timestamps plus your off-session log.
- Papers in `papers/`; only openly licensed PDFs are committed, the rest are linked by DOI/arXiv.

## Open questions to settle

*None at the moment.*
