# epi_int_3d

Extends the 2D epi_int_lite streamline N-body integrator for gravitating planetary rings
(Hahn et al. 2025) into the vertical direction, to simulate spiral bending waves. It is also a
demo of AI-assisted research, so the record of the work matters: see `PLAN.md` (phases,
decisions, ground rules) and `README.md`.

## Layout

- `source/`: epi_int_lite ported to Python 3 — the **reference implementation**. It is a straight
  translation of the Python 2 code (changes in `notes/port_changes.md`); don't optimize or
  restructure it. Faster code (Phase 0b) and the 3D code (Phase 3) are checked against it.
- `tests/`: each test folder holds an `inputs.py`, its epi_int_lite notebook ported to Python 3
  (executed on the port's output) and a git-ignored `output/`. The Python 2 notebooks, with
  their original outputs, are git commit `99a7753`. `tests/analysis.py` reproduces the notebooks'
  calculations for `tests/test_regression.py`, the 2D regression suite.
- `notes/`: `epi_int_lite_issues.md` (bugs found in epi_int_lite; log each new one there),
  `equations.md` (equation verification tracker), `sympy/` (symbolic checks),
  `numerics/`, `phases.csv` and `human_time.md` (time tracking), per-phase write-ups.
- `benchmark/`: the Phase 0b benchmark (`benchmark.py`, `line_profile.py`, the benchmark ringlet's
  `inputs.py`), its results (`results/`) and the reference final states that faster code must
  match (`reference/`); scratch runs go in the git-ignored `runs/`.
- `paper/`: the mock ApJ paper (AASTeX 7.0.2). `papers/`: reference PDFs (git-ignored except hahn25).
- `prompts/`: the final version of each prompt that drove the project, `NN-short-name.md`.
- `tools/`: `session_time.py` (time/token/cost per phase), `compare_notebooks.py` (Python 2
  versus Python 3 notebook outputs), `run_simulations.py` (runs the test simulations in parallel),
  `patches.py` (run-time changes to `source/`: the fixed crossing check, edge taper, edge-flux fix).
- `fast/`: the faster code (Phase 0b). `kernels.py` holds Numba kernels; `epi_int_fast.py` runs the
  reference driver with ring gravity replaced by the kernel, bit-for-bit the same results.
- `tests/nominal/`: epi_int_3d's own nominal ringlet (not from epi_int_lite); its `inputs.py`
  applies `tools/patches.py`. `explore/`: scratch parameter studies (`run.py`, `diagnose.py`),
  outputs git-ignored.

## Running

- Python is managed by uv (Python 3.14, `pyproject.toml`, `uv.lock`, project-local `.venv/`).
  Always run through it: `uv run python ...`, `uv run pytest`.
- A simulation runs from its test folder, which supplies `inputs.py`:
  `cd tests/spreading && uv run python ../../source/epi_int_lite.py`. Parameters can be
  overridden with `-m '{"total_number_of_outputs": 20}'`. Outputs go to `output/`.
- Faster driver: `cd tests/spreading && uv run python ../../fast/epi_int_fast.py` (all cores; set
  `NUMBA_NUM_THREADS` for fewer). `EPI_INT_DRIVER=fast/epi_int_fast.py uv run pytest` checks it.
- Regression suite: `uv run pytest` (about 40 s). It compares final states with
  `tests/golden/*.npz` to relative precision 1e-11. Regenerate those (`UPDATE_GOLDEN=1 uv run
  pytest`) only when the reference implementation itself is meant to change.
- Notebooks read their test's git-ignored `output/`; make it with `uv run python
  tools/run_simulations.py [tests/<test>]` (skips runs that have `output/`). Then
  `uv run jupyter lab <notebook>`, or `uv run jupyter nbconvert --to notebook --execute
  --inplace <notebook>`. Each notebook's first cell puts `source/` on the import path.
- Paper: `cd paper && tectonic epi_int_3d.tex`.
- Symbolic checks: `uv run python notes/sympy/vertical_motion.py`.

## Units and notation

- Code units: G = 1, planet mass = 1, the ringlet's inner streamline starts at r = 1, so the
  orbit period there is ≈ 2π. Typical tests use planet radius Rp = 0.5 and J2 = 0.01.
- Cylindrical coordinates r, θ (`t` in code), z; velocities v_r, v_θ (`vr`, `vt`), v_z.
- Paper notation follows Goldreich & Tremaine: Ω_a, κ_a, ν_a are the mean angular velocity,
  epicyclic and vertical frequencies at semimajor axis a; Ω is the longitude of ascending node
  (`node` in code); M_p, R_p are the planet's mass and radius; M is the mean anomaly.
  Shu's papers use other conventions (e.g. GT's m = −Shu's m, vertical frequency μ).

## Rules

- Every derived equation gets a sympy check before it is used in code, and is tracked in
  `notes/equations.md`; the user hand-verifies them later.
- Check every new citation against SciX (scixplorer.org). The user's bibliography,
  `paper/jmh_bibliography.bib`, stays local and is never committed; entries it lacks go in
  `paper/cc_added.bib`.
- epi_int (IDL) is a 2D code: don't cite or port its unpublished vertical (z) code.
- Ask before creating or deleting any cloud resource.
- Always run with `fast_gravity = False` and `Q_ring = -1`: the fast-gravity approximation and
  the ring-pressure code are untested and unused (`fast_gravity` is dropped in epi_int_3d).
