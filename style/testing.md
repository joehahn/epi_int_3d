# Testing and validation style

How the author validates numerical work: a Jupyter notebook per test that runs the simulation's output
through the physics it should obey, and shows by plots and printed ratios whether it does. Use it for
every new test and validation notebook. Rules are numbered (T1, T2, …) so that reviews and the ledger can
cite them. Code in notebooks follows `style/code.md`.

**Sources.** epi_int_lite's notebooks (2017–2018). The four ported into this repo, as they were before
the port (`git show 99a7753:<path>`): `tests/unperturbed/check_orbits.ipynb`, `tests/spreading/spreading.ipynb`,
`tests/libration/libration_frequency.ipynb`, `tests/libration/narrow_ring/ring_mass_1.0e-9/libration.ipynb`.
And from https://github.com/joehahn/epi_int_lite: `runs/self_confining/nominal/nominal.ipynb`,
`runs/self_confining/number_of_streamlines/number_of_streamlines.ipynb` (a convergence study),
`tests/equilibrium/equilibrium.ipynb` (simulations versus a theory curve).

**Status.** Approved by the author 2026-11-08. Plotting rules here are the ones tied to
testing; the look of figures will be set by the data-visualization guide (`style/dataviz.md`, to come).

## What a test is for

Most tests do one of three things (the author, 2026-11-08):
1. confirm that the simulation's output agrees with theory;
2. show that the output agrees with expectations;
3. show that the constants of the motion are conserved reasonably well.

Each is documented by to-the-point plots that make the point at a glance. Some of those plots are
publication quality and go into journal articles or conference reports, so a test's key plot is drawn
as a candidate paper figure from the start (T2, T22).

**Judging agreement.** Whether a simulation agrees with theory is a judgment the author and Claude make
together from the plots and the printed ratios, not a pass/fail threshold (`PLAN.md`, Phase 4). The
pytest suites (`tests/test_regression.py` etc.) are a separate software check that results haven't
changed; they are not where the physics is validated.

## The test folder

- **T1. One folder per test**, holding `inputs.py`, the notebook named after the test
  (`spreading/spreading.ipynb`), and the run's `output/`. A family of runs is a parent folder with one
  subfolder per parameter value, named by that value (`narrow_ring/ring_mass_1.0e-9/`,
  `number_of_streamlines/14/`), plus a notebook in the parent that compares them.
- **T2. Paper figures are saved by the notebook** that makes them:
  `plt.savefig('figs/libration.pdf', bbox_inches='tight')`.

## The notebook's opening cells

- **T3. Header cell**: a comment block with the notebook's name, author and date, and one line saying
  what it shows or checks:
  ```python
  #check_orbits.ipynb
  #by Joe Hahn, jmh.datasciences@gmail.com, 18 October 2017.
  #evolve a handful of unperturbed particles having 1e-8 < e < 1e-2 particles for 10^5 orbits
  #and show that the fractional variations in a and e are tiny
  ```
  The expected result is stated up front, in that line. This comment, in the author's lowercase
  comment voice, is the notebook's statement of purpose; there are no markdown cells (the author's
  choice, 2026-11-08).
- **T4. Restore cell**: read `inputs.py`, `restore_output()`, convert to elements with `coords2elem()`,
  then print the array shape, the parameters that matter for this test, the range of the key quantity
  (`e.min()`, `e.max()`), the execution time and the `monitor` dict. The printout is the record of
  which run the notebook looked at.
- **T5. A timescale cell**: compute the physical timescale the test is about, print it (in code units
  and in orbits), and rescale time by it: `times_viscous = times/viscous_timescale`. Plots then use
  t/τ_ν, t/τ_e, t/τ_dyn or orbits (often in units of 10⁴ orbits), never raw code time.
- **T6. A physical-units cell where it helps intuition**: the ringlet's mass, width, orbit period and
  timescales in cgs, km and years for a Saturn-like ring, so that code units mean something.
- **T7. A plot-prep cell** with seaborn, which sets the look of every plot (the author's choice,
  2026-11-08):
  ```python
  #prep plots
  %matplotlib inline
  import seaborn as sns
  import matplotlib.pyplot as plt
  color_seq = plt.rcParams['axes.prop_cycle'].by_key()['color']
  from matplotlib import rcParams
  sns.set(font_scale=1.5, font='DejaVu Sans')
  import numpy as np
  import pandas as pd
  ```
  `sns.set(font_scale=1.1, ...)` for the larger square phase-space plots.

## One check per cell

- **T8. Each cell computes or plots one thing**, and its first line says what:
  `#plot delta_a versus time`, `#compute system's total_angular_momentum in barycentric reference frame`,
  `#equilibrium dwt and de`. No markdown cells; the comment lines carry the narrative.
- **T9. Quantities are normalized** to their initial value or natural scale before plotting:
  Δa/Δa₀, (a − a_mean)/Δa₀, e/e₀, ΔL/L₀, r/r_mid, σ/σ̄. Tiny quantities are scaled by a power of ten,
  and the label says so: `factor = 1.0e-7`, `y = delta_L/factor`, `ax.set_ylabel(r'$\Delta L/L_0$    $(10^{-7})$')`.
- **T10. Print the number behind each key plot**: `print('max delta_L = ', delta_L.max())`,
  `print('q.std()/q.mean() = ', ...)`. Equilibrium values are the median after the transient:
  `dwt_eq = np.median(dwt[idx])`.

## Kinds of check

- **T11. Conservation and precision.** Plot the fractional change of what should be conserved
  (|Δa|, |Δe|, |Δh|, ΔL/L) over the run, on log axes when it spans decades, and against the parameter it
  might depend on (`#plot each streamline's max|de| versus initial e`). Print the maximum.
- **T12. Simulation versus theory, overlaid.** The simulation as a solid line (or markers for a family
  of runs), the theory as a dashed line in the same color, labeled `'simulation'` and `'theory'` (or the
  theory's name, `'BGT'`). A theoretical constant is a horizontal line, `y = x*0 + np.sqrt(3.0)/2.0`, with
  its value as the label. Theory curves are computed in the notebook from the formula, with a helper
  function defined there if needed (`e_prime_eqn`, solved with `fsolve`).
- **T13. The ratio, printed.** For a single number the theory predicts, print simulation, theory and
  their ratio: `print('omega/omega_bgt = ', omega/omega_bgt)`.
- **T14. Fits for rates and frequencies.** A rate is a linear fit, `coeff = np.polyfit(times, phase_unwrapped, 1)`,
  with the fit overplotted on the data (`'linear fit'`) and the angle unwrapped first (`unwrap_angle`).
- **T15. Convergence by overlay.** Runs at different resolution (number of streamlines, particles per
  streamline, timestep) plotted on one axis, labeled by the parameter value, the lowest-resolution run
  drawn thicker as the reference. Then one summary quantity versus resolution, with a fit if it trends.
- **T16. Parameter scans as a table.** Loop over the run folders, collect each run's summary into a dict,
  build `df = pd.DataFrame(data)`, show `df.tail()` and sorted subsets of columns, then plot the summary
  versus the parameter (log axes for masses) with the theory curve across the whole range.
- **T17. Test a helper inline before using it.** A cell `#test the above` calls a new function on one
  set of values and shows the result, before the function is used in a plot.
- **T18. Keep what was learned.** Results from earlier trials can stay as comments in the cell that
  produced them (`#2 streamlines max_delta_L=1.6e-7`, `#3 streamlines max_delta_L=2.1e-7`).

## Plots, as used in tests

- **T19.** `f, ax = plt.subplots(figsize=(10, 6))`, one axis per figure (8×8 for phase-space plots).
  Then `x = ...`, `y = ...`, `p = ax.plot(x, y, linestyle='-', linewidth=2, label=...)`, every keyword
  spelled out, and each `ax.set_...` call assigned to `p` so the notebook prints nothing extra.
- **T20. Axis labels give the symbol and the unit or scale, separated by four spaces**:
  `'time t    (orbits)'`, `r'time t    ($10^4$ orbits)'`, `r'$\Delta\tilde{\omega}$    (radians)'`.
  Labels name the quantity in words where a symbol would be cryptic: `r'radial width $\Delta a/\Delta a_0$'`.
- **T21.** A title only when it adds something the axes don't (`'e variations vs time'`); legends
  placed explicitly when the default hides data (`loc=(0.79, 0.1)`).
- **T22. One point per plot, readable at a glance.** Each plot answers the question in its cell's first
  comment line (T8) and nothing else. Axis limits are set to frame that answer (`ax.set_ylim(1, 6)`),
  not left at matplotlib's defaults when those hide it. A plot that makes the test's main point is drawn
  so it could go into a paper unchanged: normalized axes (T9), labeled units (T20), theory overlaid (T12),
  and saved to `figs/` (T2).

## Claude's habits to avoid

These appear in notebooks Claude wrote for this project (e.g. `tests/bending_waves/bending_waves.ipynb`).

| Claude's habit | Instead |
|---|---|
| Long markdown cells explaining the physics and the result | a comment line per cell (T8); the explanation goes in `notes/` or the paper |
| A setup cell that computes many quantities at once and prints them with f-strings | one quantity per cell, printed by name (T8, T10) |
| Several plots' worth of setup on one line (`ax.set_xlabel(...); ax.set_ylabel(...)`) | one call per line (T19) |
| Pass/fail statements in prose ("agrees to within 1.3%") with no ratio printed | the ratio printed (T13), the judgment made together |
| Colormaps and color cycles chosen per plot (`plt.cm.viridis(...)`) | the default color sequence, `color_seq` (T7) |

## Open questions for the author

*(None open. Settled 2026-11-08: a purpose comment, no markdown (T3); seaborn (T7); restyling the
notebooks Claude wrote is deferred, on the to-do list in `PLAN.md`.)*

## From the author's edits

*(Each correction the author makes to a test notebook becomes a rule here, so it isn't needed twice.)*
