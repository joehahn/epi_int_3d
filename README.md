# epi_int_3d

*Work in progress — Phase 0: epi_int_lite has been ported to Python 3 (the reference
implementation, in `source/`).*

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
- nominal: epi_int_3d's nominal ringlet, eccentric with a π periapse twist, should stay smooth and
  coherent for 10⁴ orbits while self-gravity keeps J2 from winding it up, [nominal.ipynb](tests/nominal/nominal.ipynb)
  (about 30 min; see [notes/phase0b/nominal.md](notes/phase0b/nominal.md))

To run all 23 epi_int_3d test simulations in parallel (about 30 min), or just one:

    uv run python tools/run_simulations.py
    uv run python tools/run_simulations.py tests/spreading

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
| epi_int_lite, Python 3 port (1 core) | 0.47 ms | 6.3 ms | 450 ms |

Joe Hahn, joe.hahn@jmh-datasciences.com
