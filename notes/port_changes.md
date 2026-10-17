# Python 2 → 3 port of epi_int_lite: changes

The reference implementation in `source/` is epi_int_lite's master branch
(https://github.com/joehahn/epi_int_lite, commit 5e0942a) translated to Python 3.14 and NumPy
2.x. The untouched Python 2.7 import is commit `24f214a`, so `git diff 24f214a -- source tests`
shows every change. The algorithms are unchanged, and there are no speedups.

## Changes needed by Python 3

| Change | Where |
|---|---|
| `print x` statements → `print(x)` | `source/epi_int_lite.py` (20), `source/helper_fns.py` (5), every `inputs.py` (1) |
| `execfile('inputs.py')` → `exec(open('inputs.py').read())` | `source/epi_int_lite.py` |
| `dict.iteritems()` → `dict.items()` | every `inputs.py` (command-line parameter overrides) |
| `/` → `//` where an integer division is used as an array index (`/` returns a float in Python 3) | `helper_fns.peri_corotate()`, `fluxes.py` (2) |

No NumPy 2 incompatibilities were found: the tests run with `python -W error`, which turns
any deprecation warning into an error.

## Other changes

- **`tests/unperturbed/inputs.py` predated the current driver.** It defined
  `confine_edges`, which the driver no longer reads, and lacked `confine_inner_edge` and
  `confine_outer_edge`, which it requires; as committed in epi_int_lite it stops with a
  `NameError`. It now sets both to `False`, the old setting.
- **Optional random seed.** A new input `random_seed` (default `None`: unseeded, as before) seeds
  NumPy's random number generator in `initialize_streamline()`, which only the `log-e` initial
  orbits use. `tests/unperturbed/inputs.py` sets `random_seed = 1`, so its random eccentricities
  are reproducible; they differ from those of the notebook's unseeded Python 2 run.

## Not ported

- The test notebooks. Their calculations are reproduced in `tests/analysis.py` and used by the
  regression suite (`tests/test_regression.py`) and by `tools/compare_port.py`; their figures
  and printed outputs are kept in each test's `reference/` folder.
- The other tests (`equilibrium`, `edge_torques`, `timestep`, `particles_per_streamline`,
  `rates`, `profile`), the hahn25 production runs (`runs/`) and the paper (`paper/`).

## Agreement with Python 2

The port reproduces the notebooks' Python 2 results to 6–10 significant figures, not bit for
bit: the runs differ at the level of floating-point roundoff (different compilers, math library
and CPU since 2018), which accumulates over 10⁵ timesteps. See `notes/phase0/port_comparison.md`.
