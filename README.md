# epi_int_3d

**This project:** An AI coding agent (Claude Code) drives a real astrophysics research project end to
end: reviewing the literature, deriving and symbolically checking the equations, porting and extending
a published N-body code, testing it, running the simulations and drafting the journal paper, with
every equation, citation and number checked along the way. The science: spiral bending waves in
planetary rings.

**Author:** Joseph M. Hahn, Ph.D.  
Research Scientist, [Space Science Institute](https://www.spacescience.org/bio.php?emp=JHAHN) (planetary ring dynamics) · jhahn@SpaceScience.org  
Independent AI and machine learning consultant, [JMH DataSciences](https://jmh-datasciences.com) · [LinkedIn](https://www.linkedin.com/in/hahnjoe/) · joe.hahn@jmh-datasciences.com

**License:** [MIT](LICENSE), [details](#license)  
**Project:** Started 2026-Oct-11 · branch `main`  
Built end-to-end with [Claude Code](https://claude.com/claude-code).

**For researchers:** epi_int_3d extends the streamline N-body code epi_int_lite (Hahn et al. 2025,
ApJ, 986, 220) into the vertical direction, to simulate spiral bending waves in self-gravitating
planetary rings. Its faster drivers reproduce the published code bit for bit, every derived equation
has a symbolic (sympy) check, and the code, tests, notes and the draft paper are all open.

**For R&D and industry:** this repo shows that an AI agent can carry a long, open-ended technical
investigation when a domain expert builds the checks around it: symbolic algebra for every derivation,
every citation checked against the literature database, regression against a reference
implementation, and every prompt and hour of work logged ([prompts/](prompts/), [PLAN.md](PLAN.md)).
It is the same pattern as my other repos, applied to research rather than to a business decision.

---

*Work in progress — Phase 0b: epi_int_lite has been ported to Python 3 (the reference
implementation, in `source/`), and faster drivers (`fast/`) give its results bit for bit.
The vertical motion has begun: inclined orbits drift correctly in the planet's potential.*

epi_int_3d will extend the 2D [epi_int_lite](https://github.com/joehahn/epi_int_lite) streamline
N-body integrator for gravitating planetary rings (Hahn et al. 2025, ApJ, 986, 220; `source/`
started as a copy of its master branch, commit 5e0942a) into the vertical direction, and use
it to simulate spiral bending waves.
The project also serves as a demo of using an AI coding agent (Claude Code) to accelerate
research: reviewing papers, deriving equations, porting and extending code, testing,
running simulations, and drafting the journal paper.

- Plan: [PLAN.md](PLAN.md)
- Prompts that drive the project: [prompts/](prompts/)
- Equation verification status: [notes/equations.md](notes/equations.md)
- Reference papers: [papers/](papers/)
- Python 3 port of epi_int_lite, and how it compares with the Python 2 results:
  [notes/port_changes.md](notes/port_changes.md),
  [notes/phase0/port_comparison.md](notes/phase0/port_comparison.md)
- Issues found in epi_int_lite: [notes/epi_int_lite_issues.md](notes/epi_int_lite_issues.md)

Setup (macOS or Linux):
- Python: [uv](https://docs.astral.sh/uv/) manages the interpreter and packages, pinned in
  `pyproject.toml`, `.python-version` and `uv.lock`, and installed in this folder's `.venv/`
  only. Run scripts with `uv run python <script>`; `uv sync` rebuilds the environment.
- LaTeX: the paper compiles with [tectonic](https://tectonic-typesetting.github.io/)
  (`brew install tectonic`; `cd paper && tectonic epi_int_3d.tex`). The AASTeX v7.0.2 class
  and bibliography style are kept in `paper/`.

### Tests:

Each tests folder holds a test's inputs.py, its output/, and a Jupyter notebook that compares that
output to theory:

- unperturbed: 51 eccentric streamlines without ring forces should stay on their initial orbits, [check_orbits.ipynb](tests/unperturbed/check_orbits.ipynb) (compare to [the same test run with epi_int_lite](https://github.com/joehahn/epi_int_lite/blob/master/tests/unperturbed/check_orbits.ipynb))
- spreading: a viscous ringlet's surface density should spread like the analytic solution, [spreading.ipynb](tests/spreading/spreading.ipynb) (compare to [the same test run with epi_int_lite](https://github.com/joehahn/epi_int_lite/blob/master/tests/spreading/spreading.ipynb))
- libration: 20 narrow and wide ringlets of differing mass should librate at the rate predicted by
  Borderies, Goldreich & Tremaine, [libration_frequency.ipynb](tests/libration/libration_frequency.ipynb) (compare to [the same test run with epi_int_lite](https://github.com/joehahn/epi_int_lite/blob/master/tests/libration/libration_frequency.ipynb))
- nominal: epi_int_3d's nominal ringlet, 1000 streamlines, eccentric with aligned periapses, inviscid,
  should stay smooth for 5×10⁴ orbits while self-gravity keeps J2 from winding it up, [nominal.ipynb](tests/nominal/nominal.ipynb)
  (about 3 h with `--fast` below, days with the reference driver; see [notes/phase0b/nominal.md](notes/phase0b/nominal.md))
- inclined: inclined test particles should keep their inclinations while their nodes regress at the
  rate set by J2, [inclined.ipynb](tests/inclined/inclined.ipynb) (epi_int_3d only; see [notes/phase3/vertical_drift.md](notes/phase3/vertical_drift.md))

To run all 24 epi_int_3d test simulations in parallel (about 3 h with `--fast`, set by the nominal ringlet;
run the nominal ringlet with `--fast`, since the reference driver would take days), or just one, or one
with the faster driver (same results, bit for bit):

    uv run python tools/run_simulations.py
    uv run python tools/run_simulations.py tests/spreading
    uv run python tools/run_simulations.py tests/nominal --fast

then use Jupyter to inspect a test's output, e.g.

    cd tests/unperturbed
    uv run jupyter lab check_orbits.ipynb

and click Run > Run All Cells.

The regression suite reruns those tests (some shortened) and checks their results in about 40 s:

    uv run pytest

### Performance:

Time per timestep on an Apple M5 laptop, for a ringlet of Nr circular streamlines with 20 particles each
(`uv run python benchmark/benchmark.py`, details in [notes/phase0b/benchmark.md](notes/phase0b/benchmark.md)):

| code | Nr = 10 | Nr = 100 | Nr = 1000 |
|---|---:|---:|---:|
| epi_int_lite, Python 3 port (1 core) | 0.49 ms | 6.3 ms | 450 ms |
| + Numba gravity kernel, 10 threads (`fast/epi_int_fast.py`) | 0.38 ms | 0.88 ms | 14 ms |
| + own driver, bookkeeping kernels (`fast/epi_int_3d.py`) | 0.30 ms | 0.73 ms | 12 ms |

and for the former nominal ringlet (eccentric, viscous, 40 particles per streamline; `--ringlet nominal`):

| code | Nr = 33 | Nr = 100 | Nr = 333 |
|---|---:|---:|---:|
| epi_int_lite, Python 3 port (1 core) | 1.8 ms | 12.4 ms | 105 ms |
| + Numba gravity kernel, 10 threads (`fast/epi_int_fast.py`) | 0.66 ms | 1.7 ms | 6.3 ms |
| + own driver, bookkeeping kernels (`fast/epi_int_3d.py`) | 0.49 ms | 1.1 ms | 4.3 ms |

The faster drivers give bit-for-bit the same results ([notes/phase0b/numba.md](notes/phase0b/numba.md),
[notes/phase0b/driver.md](notes/phase0b/driver.md)).

## About the author

I am Joseph M. Hahn, Ph.D., an independent AI and machine learning consultant. Through **JMH DataSciences** I build AI and machine learning systems that automate real decisions. Before going independent I spent eight years inside Oracle's AI Center of Excellence delivering AI systems for enterprise clients in manufacturing, oil and gas, public sector, and retail, and before that four years building machine learning systems on large data platforms.

I am also a planetary dynamicist: a Research Scientist at the [Space Science Institute](https://www.spacescience.org/bio.php?emp=JHAHN) (Ph.D., University of Notre Dame, 1997), working on the dynamics of planetary rings; my most recent paper is Hahn et al. (2025, ApJ, 986, 220), on self-confining eccentric ringlets, which describes epi_int_lite. For academic correspondence and citations, use jhahn@SpaceScience.org.

This repo applies an AI agent to research rather than to a business decision. The agent does the work, while the scaffolding keeps it honest: every derived equation gets a symbolic check, every citation is checked against the literature database, and every faster code path must reproduce the reference implementation bit for bit. If that shape matches a problem in your business, the work I do and what it costs are at [jmh-datasciences.com](https://jmh-datasciences.com).

**Related work:**

- [**epi_int_lite**](https://github.com/joehahn/epi_int_lite). The 2D streamline N-body code for planetary rings that this project extends (Hahn et al. 2025, ApJ, 986, 220).
- [**geo-herd-rider**](https://github.com/joehahn/geo-herd-rider). An AI agent reads a continuous feed of unstructured news and makes a routine judgment call a person would otherwise make by hand, while deterministic code handles everything that can be computed rather than judged.
- [**diplomacy-A2A**](https://github.com/joehahn/diplomacy-A2A). Seven AI agents play the seven nations of *Diplomacy*, the classic negotiation board game, forming alliances, bargaining, and betraying each other over the A2A protocol. Playtested with several LLMs powering the agents, including different LLMs head to head.

## License

[MIT](LICENSE): anyone may use, change and redistribute this code and writing, commercially or not,
as long as the copyright notice is kept.

**Citing this work.** If you use the code or results in published research, please cite Hahn et al.
(2025, ApJ, 986, 220) and this repository ([CITATION.cff](CITATION.cff); GitHub's "Cite this
repository" button gives the reference).

**Third-party content.** The AASTeX class and bibliography style in `paper/` belong to the American
Astronomical Society and are distributed under their own terms; `papers/hahn25.pdf` is published by
the AAS under CC BY 4.0. They are not covered by the MIT license above.
