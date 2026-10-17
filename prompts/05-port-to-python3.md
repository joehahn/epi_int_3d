# 05 — Port epi_int_lite to Python 3

Start Phase 0 of `PLAN.md`: bring the epi_int_lite code into this repo and port it to modern
Python, as the reference implementation that later, faster versions are checked against.

- Copy `source/` and the `unperturbed`, `libration` and `spreading` tests from the master
  branch of https://github.com/joehahn/epi_int_lite into this repo, with no git history. Cite
  that repo as the origin in the README.
- Port the code from Python 2.7 to the Python 3.14 environment already set up with uv, using
  NumPy 2.x. Make it a straight translation: no speedups and no changes to the algorithms,
  so the port can be compared directly with the original. Fix anything that modern Python or
  NumPy no longer supports, and list each such fix.
- Make random initial conditions reproducible with a seed.
- Check the port by rerunning the three tests and comparing against the figures in
  epi_int_lite's test notebooks, and against analytic expectations: unperturbed orbits that
  stay on their initial orbits, the libration frequency of a narrow eccentric ringlet, and
  the viscous spreading rate. Port epi_int_lite's test notebooks to Python 3 too, and rerun
  them on the new code's output, so that the same analysis code produces both sets of results.
  Show me the comparisons.
- Turn those checks into a pytest regression suite, the **2D regression suite**, that runs
  quickly with `uv run pytest`.
- Add a `CLAUDE.md` with the project's conventions, units, notation and how to run the code
  and tests.
