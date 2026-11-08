# Data-visualization style

How the author's figures look: notebook plots, paper figures and animations. Rules are numbered (V1, V2, …)
so that reviews and the ledger can cite them. How a plot is used in a test (one point per plot, theory
overlaid, ratios printed) is in `style/testing.md`; captions and how the text cites figures are in
`joes_papers/writing_style.md` ("Figures and captions"); plotting code follows `style/code.md`.

**Sources.** Hahn et al. (2025, ApJ 986, 220; "hahn25"), whose 29 figures are the author's notebook plots
placed directly in the paper, and the epi_int_lite notebooks listed in `style/testing.md`. Hahn & Spitale
(2013; "hahn13") and the earlier papers were made with IDL (black and white, boxed axes, gray bands for
uncertainty); where they differ, follow hahn25.

**Status.** Approved by the author 2026-11-08.

## The look

- **V1. seaborn's gray grid**, set once in the plot-prep cell (T7):
  `sns.set(font_scale=1.5, font='DejaVu Sans')`. That gives the gray background with white grid lines,
  no box, and seaborn's default ("deep") palette. Colors come from it in order, as
  `color_seq = plt.rcParams['axes.prop_cycle'].by_key()['color']`; no colormaps for line plots.
- **V2. One panel per figure.** `f, ax = plt.subplots(figsize=(10, 6))`; `figsize=(8, 8)` with
  `font_scale=1.1` for a square phase-space plot (I_y versus I_x). Every hahn25 figure has one panel;
  comparisons go on one axis, not in a grid of subplots.
- **V3. In the paper, each figure is one column wide**, `\includegraphics[width=\columnwidth]`. The 10×6
  figure at font scale 1.5 shrinks to text-sized labels at that width, which is why those values are used.

## Lines and markers

- **V4. Data lines** `linestyle='-', linewidth=2`. A family of many curves (one per output time or per
  run) uses `linewidth=1` or `alpha=0.7–0.8` so they stay distinct.
- **V5. The reference run stands out.** In a family of runs, the nominal or lowest-resolution run is
  drawn thicker (`linewidth=3`–`6`, `markersize=15`) in its own color (hahn25 Figs. 17, 18, 20, 21).
- **V6. Theory is dashed** (`linestyle='dashed'`) in the color of the simulation it goes with; a
  constant is a dotted horizontal line, `y = x*0 + np.sqrt(3.0)/2.0`, labeled with its value; a fit is
  labeled `'linear fit'`.
- **V7. Discrete runs are markers joined by lines**: `marker='o', markersize=7–9, linestyle='-'`.
  Points from a parameter scan are markers alone, with the theory as a line across the whole range.

## Axes, labels and legends

- **V8. Axis labels are words, then the symbol**, lowercase: `'radial width $\Delta a/\Delta a_0$'`,
  `r'time $t/\tau_\nu$'`, `'eccentricity e'`, `r'surface density $\sigma/\sigma_0$'`,
  `r'true anomaly $\varphi/\pi$'`. A unit or scale follows after four spaces, in parentheses:
  `r'$\Delta L/L_0$    $(10^{-7})$'`, `r'time t    ($10^4$ orbits)'`.
- **V9. Axes are dimensionless and normalized**: time in the timescale of the physics (t/τ_ν, orbits),
  angles in units of π, distances and densities relative to a natural scale (Δa/Δa₀, σ/σ₀). Two quantities
  of different size can share an axis if one is rescaled and the label says so: `'Δϖ/10    and    Δe'`
  (hahn25 Fig. 4).
- **V10. Log axes** for anything spanning decades (|Δe| over 10⁵ orbits, ring mass in a parameter scan).
- **V11. Limits frame the point** of the plot (`ax.set_ylim(1, 6)`), not matplotlib's defaults.
- **V12. Legends inside the axes**, framed (seaborn's default), with short entries: `'outer'`, `'mean'`,
  `'inner'`; the parameter's value alone for a family of runs (`'2'`, `'3'`, … streamlines; `'0.002'`);
  `r'$t/\tau_\nu$ = 0.3'` for output times, with a second quantity in the same entry if it helps
  (`"e' = 0.011"`, hahn25 Fig. 9). Placed explicitly (`loc=(0.79, 0.1)`) when the default covers data. Where a legend would be
  crowded, curves can be labeled with text on the plot instead (hahn25 Fig. 13).
- **V13. No titles on paper figures**: the caption says what is plotted. Notebook diagnostics may have a
  short title (`'e variations vs time'`).

## Paper figures

- **V14. A paper figure is the notebook plot**, made by the notebook or a script beside it, saved as PDF
  with `plt.savefig('figs/name.pdf', bbox_inches='tight')`, and never edited by hand. In this repo the
  paper's figures are in `paper/figures/`, one script per figure that writes `name.pdf` (and a `name.png`
  preview), so a rerun regenerates every figure (`PLAN.md`, Phase 3b, "one source for every number").
- **V15. Every symbol and line in the figure appears in the caption or text**, and the caption follows
  `joes_papers/writing_style.md`: what is plotted versus what, in which units, what each color is.

## Animations

- **V16. Only when the result needs one** (an evolving wave a static plot can't show); static plots
  otherwise, and ask before making one. Each frame lists the run's parameters, with `dt = <value>` just
  below Nr and Nt. Written as a GIF and as a `docs/` page with play and pause (GitHub Pages), and opened
  in Firefox.

## Claude's habits to avoid

From figures Claude made for this project before this guide (`paper/figures/vertical_resonance.py`,
`zigzag.py`, the `bending_waves` notebook and animations).

| Claude's habit | Instead |
|---|---|
| matplotlib's default white style, boxed axes | seaborn's gray grid (V1) |
| Small figures sized for the column (`figsize=(4.6, 3.6)`) | 10×6 at font scale 1.5, shrunk to the column (V2, V3) |
| Simulation as black dots, theory as a solid line | simulation solid or markers, theory dashed in the same color (V4–V7) |
| Colormaps for line families (`plt.cm.viridis(...)`) | `color_seq` in order (V1) |
| Several panels with `gridspec` | one panel per figure (V2) |
| Symbol-only or bracketed labels (`inclination ⟨i⟩`) | words, then the symbol (V8) |
| `plt.tight_layout()` and PNG output for the paper | `savefig('….pdf', bbox_inches='tight')` (V14) |

## Diagrams and color maps

- **V17. Diagrams are plain line drawings on white**, without the gray grid, since a grid behind a
  geometric sketch is clutter and there are no data values to read off it. They match the plots
  otherwise: DejaVu Sans, the same palette (`color_seq`), lowercase word-then-symbol labels (V8), one
  column wide. hahn25 has no diagrams; this rule is Claude's call, which the author left open (2026-11-08).
- **V18. No 2D color maps by default.** hahn25 has none, and the color maps Claude made for the
  bending-wave runs will probably be deleted (the author, 2026-11-08). Show a 2D field as a family of
  line cuts (V4, V5) where that makes the point; ask the author before making a color map.

## Open questions for the author

*(None open. Settled 2026-11-08: diagrams (V17), color maps (V18); restyling the existing paper figures
`vertical_resonance` and `zigzag` is deferred, on the to-do list in `PLAN.md`.)*

## From the author's edits

*(Each correction the author makes to a figure becomes a rule here, so it isn't needed twice.)*
