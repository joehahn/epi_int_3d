# Code style

The author's style for Python code, taken from his pre-AI code so that new code reads like his and
can be spot-checked at a glance. Use it for all new code in `fast/`, `tools/`, `tests/` and `explore/`.
Rules are numbered so that reviews and the ledger can cite them (C1, C2, …).

**Sources.** epi_int_lite (2017–2018), i.e. this repo's `source/` (`helper_fns.py`, `epi_int_lite.py`,
`libration.py`, `fluxes.py`) and the tests' `inputs.py` files, read before the Python 3 port
(`git show 24f214a:source/...`). Notebook code follows the same rules; how notebooks test things is in
`style/testing.md`.

**Status.** Approved by the author 2026-11-08. `fast/` predates this guide and is
written in Claude's own style; restyling it is on the to-do list in `PLAN.md`.

## Files and layout

- **C1. Header block.** Every file starts with a comment block: the file name, a "by" line with
  author, email and date, and one line saying what the file does. Notes on provenance follow.
  ```python
  #helper_fns.py
  #
  #by Joe Hahn, joe.hahn@jmh-datasciences.com, 28 September 2017.
  #
  #helper functions used by epi_int_lite.py
  ```
  Code that Claude writes for the author says so: `#by Joe Hahn and Claude Code, ..., 8 November 2026.`
- **C2. Few, flat modules of plain functions.** One helper module per purpose (`helper_fns.py`,
  `libration.py`, `fluxes.py`), each a list of top-level functions. No classes, no inheritance, no
  decorators beyond what Numba needs, no packages within packages. State that must persist is held in
  a plain dict (`monitor`, `initial_orbits`).
- **C3. A driver script that reads top to bottom.** The driver (`epi_int_lite.py`) is a script, not a
  `main()` function: read inputs, print them, initialize, loop, save. Each step of the loop has its
  own comment line, so the algorithm can be read off the comments alone.
- **C4. Inputs are module-level variables in `inputs.py`.** One commented group per physical idea,
  with the sign convention that turns a feature off stated in the comment:
  ```python
  #ring kinematic shear and bulk viscosity, set < 0 to turn off
  shear_viscosity = 1.0e-13
  bulk_viscosity = 1.0*shear_viscosity
  ```
  Options with several choices are a small dict with a `'shape'`-style key (`initial_orbits`).
  Code doesn't check its inputs or raise errors for settings it can't handle; the comment states what
  is supported (the author's choice, 2026-11-08).
- **C4a. Reading inputs and importing helpers** (the author's choice, 2026-11-08). Every script and
  notebook reads its parameters by running `inputs.py`, as `source/` does
  (`with open('inputs.py') as fp: exec(fp.read())`). Notebooks import helpers with
  `from helper_fns import *`. New code modules list what they import,
  `from helper_fns import Omega, Kappa, coords2elem`, so a reader can see where each name comes from,
  and two modules that define the same name (`Omega` is in both `helper_fns.py` and `fast/satellite.py`)
  can't silently shadow each other.
- **C5. The driver prints every input** it reads, one per line, as `print('name =', value)`, so the
  run's log records its parameters.

## Comments

- **C6. A one-line comment above every function** says what it computes, in lowercase, with no
  space after `#` and no final period: `#epicyclic frequency`, `#convert coordinates to orbit elements`.
  No docstrings.
- **C7. Comments name the physics and its source**: `#evolve system using Chamber's (1993) 2nd order
  kick-drift-kick scheme`, `#BGT libration frequency`. For this project, cite the equation's label in
  `notes/equations.md` or the paper: `#vertical drift (V8)`.
- **C8. Comments mark each block inside a function**: `#in ring interior`, `#along ring edges`,
  `#at outer streamline`. A trailing comment explains a single line or a value:
  `random_seed = None     #inputs.py may set this to make random initial orbits reproducible`.
- **C9. Candor is fine in comments**: "buyer beware", "this works for 2-streamline...", "i dont this this
  plot/calcuation is correct..." (verbatim, `nominal.ipynb`). A known limitation is stated where it applies, not hidden.

## Names

- **C10. Physics arrays get short symbols** from the equations: `r, t, vr, vt` (t is the longitude θ),
  `a, e, wt, M` (wt is ϖ), `z, vz`, `Omg, Kap` for Ω and κ, `lambda_` for λ (underscore because
  `lambda` is reserved).
- **C11. Parameters and counts get long, plain snake_case names** that read as English:
  `number_of_streamlines`, `particles_per_streamline`, `total_ring_mass`, `timesteps_per_output`,
  `viscous_timescale`, `fractional_angular_momentum_change`.
- **C12. Suffixes carry meaning**, used consistently:
  - `_inner`, `_outer`, `_mean` (or `_avg`), `_middle` for which streamline;
  - `_init` or `_0` for an initial value (`da_init`, `e0`), and `_0` for the central planet (`r_0`);
  - `_now` for the slice at one output time (`r_now = r[tidx]`);
  - `_bc` for barycentric, `w` for longitude-wrapped arrays (`rw`, `tw`);
  - `_list` for a list being accumulated, a trailing `z` or `s` for a time series (`rz`, `r_bcs`);
  - `_cgs`, `_km`, `_sec`, `_gm` for physical units.
- **C13. Index names say what they index**: `streamline_idx`, `pidx` (particle), `tidx` or `t_idx`
  (output time), `sidx` (streamline), `idx` for a boolean mask: `idx = (dwt > np.pi)`.
- **C14. Unit constants are spelled out locally** where used: `GM = 1.0`, `mass_0 = 1.0`,
  `twopi = 2.0*np.pi`, rather than imported from a constants module.

## Arithmetic

- **C15. Products and powers written out.** `a*a`, `ar2*ar`, `c*c`; `**` only for fractional or large
  powers (`**0.75`, `**(1.0/3.0)`). Floats carry a decimal point: `2.0`, `1.5*J2`, `1.0e-7`.
- **C16. Repeated subexpressions get their own name** before the formula that uses them, so the
  formula matches the paper's equation term by term:
  ```python
  a2 = a*a
  Ra2 = (Rp*Rp)/a2
  Kappa2 = (GM/a2/a)*(   1.0 - (1.5*J2)*Ra2 - Ar*(a2*(3.0/GM))   )
  ```
- **C17. Spacing.** No spaces around `*` and `/`; spaces around `+`, `-`, `=` and comparisons. Three
  spaces inside the outer parentheses set off the main bracketed factor of a formula, as above and in
  `a = Rp*(   c + np.sqrt(c*c - 1.5*J2)   )`.
- **C18. Parenthesized conditions**: `if (G_ring > 0.0):`, `while (number_of_outputs < total_number_of_outputs):`,
  `if ((c > 0.0) or (shear_viscosity > 0.0)):`. Boolean masks too: `idx = (dr < 0)`.

## Functions and data flow

- **C19. Arguments in, arrays out.** Functions take what they need as positional arguments, physical
  parameters first (`J2, Rp, ...`) then arrays (`r, t, vr, vt`), and return a tuple, unpacked by name at
  the call: `a, e, wt, M = coords2elem(J2, Rp, r, t, vr, vt)`. The same names on both sides of the call.
- **C20. One statement per line.** No semicolons, no chained one-liners, no lambdas, no nested
  comprehensions. A simple comprehension is fine (`[de.max() for de in delta_e]`). Long calls break
  with `\` and an indented continuation.
- **C21. Accumulate in a loop, then convert.** `x_list += [x]` inside the loop, `x = np.array(x_list)`
  after it, rather than preallocating or vectorizing at the expense of readability.
- **C22. Optional behavior via keyword flags with defaults**: `Kappa(J2, Rp, a, Ar=0.0, kappa_squared=False)`,
  `interpolate_fn(t, f, n, interpolate=True)`.
- **C23. Speed tricks are commented with their reason**: `#shift 2D array ... is significantly faster
  than .roll()`. Numba kernels follow every rule above that Numba allows (names, comments, spacing,
  written-out arithmetic); their loops are the formula written index by index.

## Output

- **C24. Plain `print` with the name, then the value**: `print('dt =', dt)`,
  `print('max fractional angular momentum change = ', x.max())`. Units go in the name:
  `print('viscous_timescale (orbits) = ', ...)`. Use `np.round(x, 6)` for rounding, not f-string format
  specs.
- **C25. Outputs as `.npy` files plus a pickled dict** (`monitor.pkl`) in `output/`, written by one
  `save_output()` and read back by its mirror `restore_output()`.

## Claude's habits to avoid

These appear in code Claude wrote for this project before this guide existed (e.g. `fast/satellite.py`,
`tests/analysis.py`, the `bending_waves` notebook). Each one breaks a rule above.

| Claude's habit | Instead |
|---|---|
| Docstrings (`"""Return a test's parameters..."""`) | a `#` comment line (C6) |
| Classes holding derived constants (`class Satellite`) | functions and a dict (C2) |
| One-letter containers: `p`, `o`, `d`, `k` for parameters, outputs, loop items | named variables (C10–C13) |
| `lambda` functions, `;` on one line, several statements per line | one statement per line (C20) |
| f-strings with format specs: `f'{x:.4g}'` | `print('x =', x)` (C24) |
| `**2`, `**3` | `x*x`, `x*x*x` (C15) |
| `pathlib`, `try/finally`, type checks and `raise` for every bad input | plain paths and a comment stating the convention (C4) |
| Dense formulas that inline every subexpression | named intermediates (C16) |

## Open questions for the author

1. **Python 3 features.** None of these are in `source/`, being Python 2: f-strings, `pathlib`, `dataclasses`.
   The guide rules them out. Is that right, or are some welcome?

## From the author's edits

*(Each correction the author makes to code becomes a rule here, so it isn't needed twice.)*
