#run_all.py: the libration tests (tests/libration) with fast/epi_int_3d.py, twice each: 'new' (epi_int_3d's default
#ring gravity, epi_int_lite issues 7 and 9 fixed) and 'lite' (epi_int_lite_gravity = True). Outputs (git-ignored) in
#explore/libration_gravity/<ring>/<mass>/<new|lite>/output.
#Run:  uv run python explore/libration_gravity/run_all.py
import json, os, subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def run(job):
    folder, g = job
    out = HERE/folder.relative_to(REPO/'tests/libration')/g/'output'
    if out.is_dir():
        return f'skip {out}'
    m = {'output_folder': str(out)}
    if g == 'lite':
        m['epi_int_lite_gravity'] = True
    subprocess.run(['uv', 'run', 'python', str(REPO/'fast/epi_int_3d.py'), '-m', json.dumps(m)], cwd=folder,
        env=dict(os.environ, NUMBA_NUM_THREADS='1'), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    return f'done {out}'


jobs = [(f, g) for f in sorted((REPO/'tests/libration').glob('*_ring/ring_mass_*')) for g in ('new', 'lite')]
with ThreadPoolExecutor(8) as ex:
    for msg in ex.map(run, jobs):
        print(msg, flush=True)
