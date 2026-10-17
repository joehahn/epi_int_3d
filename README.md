# epi_int_3d

*Work in progress — Phase 0: epi_int_lite has been ported to Python 3 (the reference
implementation, in `source/`).*

epi_int_3d will extend the 2D epi_int_lite streamline N-body integrator for gravitating
planetary rings into the vertical direction, and use it to simulate spiral bending waves.
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

Setup (macOS or Linux):
- Python: [uv](https://docs.astral.sh/uv/) manages the interpreter and packages, pinned in
  `pyproject.toml`, `.python-version` and `uv.lock`, and installed in this folder's `.venv/`
  only. Run scripts with `uv run python <script>`; `uv sync` rebuilds the environment.
- Tests: `uv run pytest` runs the 2D regression suite (about 40 s). A simulation runs from
  its test folder, e.g. `cd tests/spreading && uv run python ../../source/epi_int_lite.py`.
- LaTeX: the paper compiles with [tectonic](https://tectonic-typesetting.github.io/)
  (`brew install tectonic`; `cd paper && tectonic epi_int_3d.tex`). The AASTeX v7.0.2 class
  and bibliography style are kept in `paper/`.

Built on:
- epi_int_lite (Python): https://github.com/joehahn/epi_int_lite — Hahn et al. 2025, ApJ, 986, 220.
  `source/` started as a copy of its master branch (commit 5e0942a).
- epi_int (IDL): https://github.com/joehahn/epi_int — Hahn & Spitale 2013, ApJ, 772, 122

Joe Hahn, joe.hahn@jmh-datasciences.com
