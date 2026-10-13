# 03 — Draft the project plan

I want to turn this project into a demo of using AI to accelerate research. The README will be
about the prompts that instruct Claude Code to:

1. review the relevant papers,
2. derive the equations for vertical motion,
3. start from the epi_int_lite code at https://github.com/joehahn/epi_int_lite and adapt it to
   track the particles' |z| > 0 motion (the new code is called **epi_int_3d**),
4. test the code,
5. use the code to simulate spiral bending waves, and
6. draft the journal paper.

Draft a high-level plan for this project in `PLAN.md`; we will iterate on it. Plan only — do
not start any of the work. Requirements:

**Repository**
- A new public GitHub repo `joehahn/epi_int_3d`, starting fresh with no git connection to
  epi_int_lite; refer to the epi_int_lite URL as needed. Start private and make it public later.
- Start from epi_int_lite's `lindblad_resonance` branch, which has the satellite forcing.
  The older epi_int code (IDL) is at https://github.com/joehahn/epi_int.
- Reference papers go in `papers/`. Only openly licensed papers (e.g. hahn25, CC BY 4.0) are
  committed; link the others via DOI and arXiv.

**Working rules**
- Save the final version of each prompt in `prompts/`, with no revision history.
- Claude Code verifies every derived equation with sympy; I verify them by hand. My checks
  will lag, so track which equations still need my verification. No separate numerical
  verification of derivations.
- Track time and cost per phase, including my (human) development time, using the
  timestamps in the Claude Code session logs where possible.

**Code**
- Port the Python 2.7 code to the most suitable modern Python, aiming to shrink execution
  time as much as possible (parallelization or other numerical techniques).
- Simulations will run on multiple cores in a VM in a dedicated Google Cloud project, with a
  budget of $100.

**Science**
- Simulate a low-m (m = 2 or 3) spiral bending wave: ~20 streamlines per wavelength, ~20 wave
  crests, ~20 particles per streamline. The wavelength shrinks as 1/x with distance x from
  resonance.
- The wave is launched by a satellite at a vertical resonance. The satellite's horizontal
  forcing is derived in hahn13 and implemented in epi_int; its vertical forcing and the
  wavelength formula are derived in Shu's papers and will be rederived here.
- My mathematical treatment follows Goldreich & Tremaine, which differs from Shu's: GT's m is
  −1 × Shu's m, and there may be other differences.
- Scope: linear validation plus a nonlinear parameter study (no Saturn/Cassini comparison).
- Stop each simulation before the wave reaches the ring edge; a radiative boundary condition
  is a nice-to-have.
- Target journal: ApJ, the same as hahn25.
