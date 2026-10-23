# Angular momentum conservation of the ring-gravity kicks

Measured 2026-10-23. Scripts: `dL.py` and `torque.py` (session scratch, not kept). Runs:
`benchmark/runs/lamfix/{off,on}` and `benchmark/runs/Lgrav/{on,off}` (git-ignored).

## What was measured

The nominal ringlet, at 100 streamlines × 20 particles, inviscid, dt = 0.5, `fast/epi_int_3d.py`.

| run | total \|ΔL\|/L | rate |
|---|---|---|
| 5000 orbits, reference λ (`fix_lambda` off) | 1.08×10⁻¹² | secular, −2×10⁻¹⁶ per orbit |
| 5000 orbits, `fix_lambda` on | 1.12×10⁻¹² | the same |
| 500 orbits, ring gravity on | 1.0×10⁻¹³ | −2.1×10⁻¹⁶ per orbit |
| 500 orbits, ring gravity off (`G_ring = 0`) | 9×10⁻¹⁵ | +1.8×10⁻¹⁷ per orbit |

The drift step conserves each particle's h = r v_θ exactly, so the drift comes from the kicks, and
nearly all of it from ring gravity. The torques don't cancel because the wire forces are not
reciprocal. It is not roundoff.

**Fraction of the torques that fails to cancel**, Σ m r A_θ / Σ |m r A_θ|, with the forces evaluated on
one snapshot:

| λ | wire direction | initial state | final state (5000 orbits) |
|---|---|---|---|
| perturbed particle's (epi_int_lite) | velocity (epi_int_lite) | −1×10⁻¹⁸ | +3.1×10⁻⁶ |
| perturbed particle's | shape (3-point Lagrange slope) | +2×10⁻¹⁴ | +6.2×10⁻⁶ |
| perturber's own (`fix_lambda`) | velocity | −1×10⁻¹⁸ | −6.9×10⁻⁷ |
| perturber's own | shape | +2×10⁻¹⁴ | +2.5×10⁻⁶ |

The initial state is aligned and symmetric, so its torques cancel almost exactly. The final state is
twisted.

## What it means

- A few parts in 10⁶ of the streamlines' torques fail to cancel. That matches the ~10⁻⁵–10⁻⁶ the
  author recalled.
- The cause is the wire model, not the details of the local force. A particle of streamline a feels a′
  as a straight wire through a′'s crossing point at that particle's longitude. The particles of a′
  feel a through other points, with other straight-line fits. So the two pulls are never exactly
  equal and opposite.
- Neither λ′ (issue 7) nor the shape slope fixes it. λ′ helps this snapshot a few times but leaves
  the 5000-orbit drift unchanged. The shape slope, taken from the derivative of the 3-point
  polynomial, is noisier than the interpolated values and makes the snapshot slightly worse.
- Making conservation exact would mean building it in: compute each pair interaction once and apply
  equal and opposite torques to the two streamlines. That is a design change, logged as a to-do in
  PLAN.md.

Caveats: one ringlet, one snapshot, and the crude 3-point slope.
