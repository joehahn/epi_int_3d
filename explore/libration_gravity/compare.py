#compare.py: libration diagnostics (tests/analysis.py libration_metrics, as in libration_frequency.ipynb) of the
#runs made by run_all.py, epi_int_3d's ring gravity ('new') versus epi_int_lite's ('lite'), and 'lite' versus
#the reference's own outputs in tests/libration/*/output (source/epi_int_lite.py), where those exist.
#Run:  uv run python explore/libration_gravity/compare.py
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO/'tests'))
from analysis import load, libration_metrics


def max_dr(new, lite):
    """Largest |r_new - r_lite| over outputs and particles, after aligning each streamline's particles: the driver
    re-indexes a particle that crosses theta = +-pi, and the two runs can do so one output apart."""
    worst = 0.0
    for k in range(lite['r'].shape[0]):
        for s in range(lite['r'].shape[1]):
            worst = max(worst, min(np.abs(np.roll(new['r'][k, s], sh) - lite['r'][k, s]).max() for sh in (-1, 0, 1)))
    return worst

print(f"{'run':28s} {'w/w_bgt lite':>12s} {'w/w_bgt new':>12s} {'dw/w':>9s} {'dI/I':>9s} {'d(de)/de':>9s} "
      f"{'max|dr|/width':>13s} {'lite=ref':>8s}")
for f in sorted((REPO/'tests/libration').glob('*_ring/ring_mass_*')):
    rel = f.relative_to(REPO/'tests/libration')
    p, lite = load(f, HERE/rel/'lite'/'output')
    _, new = load(f, HERE/rel/'new'/'output')
    ml, mn = libration_metrics(p, lite), libration_metrics(p, new)
    dr = max_dr(new, lite)/p['radial_width']
    ref = f/'output'
    same = '-'
    if ref.is_dir():
        _, o = load(f)
        same = 'yes' if (o['r'].shape == lite['r'].shape and np.array_equal(o['r'], lite['r'])
                         and np.array_equal(o['vt'], lite['vt'])) else 'NO'
    print(f"{str(rel):28s} {ml['omega_over_omega_bgt']:12.6f} {mn['omega_over_omega_bgt']:12.6f} "
          f"{mn['omega']/ml['omega'] - 1:9.1e} {mn['I_avg']/ml['I_avg'] - 1:9.1e} {mn['de_eq']/ml['de_eq'] - 1:9.1e} "
          f"{dr:13.1e} {same:>8s}")
