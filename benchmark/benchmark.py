#!/usr/bin/env python3
"""Phase 0b benchmark: time an epi_int_lite-style driver on the benchmark ringlet.

The ringlet (benchmark/inputs.py) is the wide 1e-9 libration ringlet with Nt = 20 particles per
streamline and viscosity on. As the number of streamlines Nr grows, radial_width and
total_ring_mass grow with it, so the streamline spacing (1e-5) and surface density stay fixed.

Each case runs the driver twice, for n1 and n2 timesteps, as a separate process from a scratch
folder under benchmark/runs/; the time per step is the difference in wall-clock time divided by
n2 - n1, so start-up, initialization and output cancel. The n2 run's final state is the
reference that faster code must reproduce.

With --ringlet nominal it runs epi_int_3d's former nominal ringlet (benchmark/nominal_inputs.py: eccentric,
tapered edges, viscous, Nt = 40) at Nr = 33, 100 and 333 streamlines, with its width and mass fixed, so the streamline
spacing shrinks as Nr grows. (At Nr = 333 it goes unstable after ~1800 steps, far beyond the 250 timed.)

Run:  uv run python benchmark/benchmark.py                  #time the reference implementation
      uv run python benchmark/benchmark.py --profile        #also cProfile each case
      uv run python benchmark/benchmark.py --driver <path> --label <name> --check
                                                            #time other code, compare with reference
      uv run python benchmark/benchmark.py --ringlet nominal [--driver ... --label ... --check]
"""
import argparse
import json
import os
import platform
import pstats
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

BENCH = Path(__file__).resolve().parent
REPO = BENCH.parent
sys.path.insert(0, str(REPO/'source'))
from helper_fns import restore_output, wrap_ring, interpolate_fn

SPACING = 1.0e-5            #streamline spacing in a
SIGMA_MASS = 1.0e-9/1.0e-3  #total_ring_mass/radial_width, the wide 1e-9 libration ringlet's

#(Nr, Nt, n1, n2): n2 - n1 timesteps are timed
CASES = [(10, 20, 200, 1200), (30, 20, 200, 1200), (100, 20, 100, 600), (300, 20, 50, 250),
         (1000, 20, 50, 250), (100, 10, 100, 600), (100, 40, 100, 600)]
REFERENCE_SIZES = {(10, 20), (100, 20), (1000, 20)}
NOMINAL_CASES = [(33, 40, 200, 1200), (100, 40, 100, 600), (333, 40, 50, 250)]
NOMINAL_INPUTS = BENCH/'nominal_inputs.py'          #the nominal ringlet before 2026-10-19
NOMINAL_WIDTH = 0.0125
GOLDEN_RTOL = 1e-12


def case_name(Nr, Nt, ringlet='benchmark'):
    return ('nominal_' if ringlet == 'nominal' else '') + f'nr{Nr:04d}_nt{Nt:03d}'


def params(Nr, Nt, steps, ringlet='benchmark'):
    if ringlet == 'nominal':
        return {'number_of_streamlines': Nr, 'particles_per_streamline': Nt,
                'timesteps_per_output': steps, 'total_number_of_outputs': 1, 'output_folder': 'output'}
    width = SPACING*(Nr - 1)
    return {'number_of_streamlines': Nr, 'particles_per_streamline': Nt, 'radial_width': width,
            'total_ring_mass': SIGMA_MASS*width, 'timesteps_per_output': steps,
            'total_number_of_outputs': 1, 'output_folder': 'output'}


def run(driver, Nr, Nt, steps, profile=None, ringlet='benchmark'):
    """Run the driver for `steps` timesteps; return wall time and the run folder."""
    folder = BENCH/'runs'/f'{case_name(Nr, Nt, ringlet)}_{steps}'
    shutil.rmtree(folder, ignore_errors=True)
    folder.mkdir(parents=True)
    shutil.copy(NOMINAL_INPUTS if ringlet == 'nominal' else BENCH/'inputs.py', folder/'inputs.py')
    cmd = [sys.executable]
    if profile:
        cmd += ['-m', 'cProfile', '-o', str(profile)]
    cmd += [str(driver), '-m', json.dumps(params(Nr, Nt, steps, ringlet))]
    env = dict(os.environ, PYTHONPATH=str(Path(driver).resolve().parent))
    t0 = time.perf_counter()
    with open(folder/'output.log', 'w') as log:
        subprocess.run(cmd, cwd=folder, stdout=log, stderr=subprocess.STDOUT, env=env, check=True)
    return time.perf_counter() - t0, folder


def final_state(folder):
    r, t, vr, vt, times, monitor = restore_output(str(folder/'output'))
    return r[-1], t[-1], vr[-1], vt[-1], monitor


def min_separation(r, t, spacing=SPACING):
    """Smallest radial gap between adjacent streamlines over the initial spacing, wrapping in
    longitude (the driver's own crossing check does not wrap, and misfires near t = +-pi)."""
    rw = wrap_ring(r, longitude=False)
    tw = wrap_ring(t, longitude=True)
    dr = (rw - interpolate_fn(tw, rw, -1, interpolate=True))[1:, 1:-1]
    return dr.min()/spacing


def compare(state, ref):
    """Largest relative difference from the reference final state (angles modulo 2 pi)."""
    out = {}
    for k in ('r', 'vr', 'vt'):
        out[k] = float(np.abs(state[k] - ref[k]).max()/np.abs(ref[k]).max())
    dt = np.angle(np.exp(1j*(state['t'] - ref['t'])))
    out['t'] = float(np.abs(dt).max()/np.pi)
    return out


def machine():
    def sysctl(key):
        try:
            return subprocess.run(['sysctl', '-n', key], capture_output=True, text=True).stdout.strip()
        except OSError:
            return ''
    return {'cpu': sysctl('machdep.cpu.brand_string') or platform.processor(),
            'cores': os.cpu_count(), 'performance_cores': sysctl('hw.perflevel0.physicalcpu'),
            'memory_GB': round(int(sysctl('hw.memsize') or 0)/2**30),
            'os': platform.platform(), 'python': platform.python_version(),
            'numpy': np.__version__}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--driver', default=str(REPO/'source'/'epi_int_lite.py'))
    parser.add_argument('--label', default='reference')
    parser.add_argument('--profile', action='store_true', help='cProfile the n2 run of each case')
    parser.add_argument('--check', action='store_true', help='compare final states with reference')
    parser.add_argument('--cases', help='comma-separated Nr:Nt, e.g. 100:20,1000:20')
    parser.add_argument('--ringlet', choices=['benchmark', 'nominal'], default='benchmark')
    args = parser.parse_args()
    nominal = args.ringlet == 'nominal'
    cases = NOMINAL_CASES if nominal else CASES
    if args.cases:
        wanted = {tuple(map(int, c.split(':'))) for c in args.cases.split(',')}
        cases = [c for c in CASES if c[:2] in wanted]
    results = {'label': args.label, 'ringlet': args.ringlet, 'driver': os.path.relpath(Path(args.driver).resolve(), REPO),
               'date': time.strftime('%Y-%m-%d'), 'machine': machine(), 'cases': []}
    (BENCH/'results').mkdir(exist_ok=True)
    (BENCH/'reference').mkdir(exist_ok=True)
    print(f"{'Nr':>5} {'Nt':>4} {'s/step':>10} {'us/(Nr^2 Nt)':>13} {'min gap':>8}")
    for Nr, Nt, n1, n2 in cases:
        t1, _ = run(args.driver, Nr, Nt, n1, ringlet=args.ringlet)
        prof = BENCH/'runs'/f'{case_name(Nr, Nt, args.ringlet)}.prof' if args.profile else None
        t2, folder = run(args.driver, Nr, Nt, n2, profile=prof, ringlet=args.ringlet)
        if prof:     #profiling slows the run, so time it again without
            t2, folder = run(args.driver, Nr, Nt, n2, ringlet=args.ringlet)
        per_step = (t2 - t1)/(n2 - n1)
        r, t, vr, vt, monitor = final_state(folder)
        state = dict(r=r, t=t, vr=vr, vt=vt)
        case = {'Nr': Nr, 'Nt': Nt, 'n1': n1, 'n2': n2, 'wall_n1_s': t1, 'wall_n2_s': t2,
                's_per_step': per_step,
                'min_gap': float(min_separation(r, t, NOMINAL_WIDTH/(Nr - 1) if nominal else SPACING)),
                'driver_crossing_flag': monitor.get('streamline_crossing_timestep'),
                'nan': monitor.get('nan_timestep')}
        ref_file = BENCH/'reference'/f'{case_name(Nr, Nt, args.ringlet)}_{n2}.npz'
        if args.label == 'reference' and (nominal or (Nr, Nt) in REFERENCE_SIZES):
            np.savez(ref_file, **state)
        if args.check and ref_file.exists():
            case['max_rel_diff'] = compare(state, np.load(ref_file))
            case['matches_reference'] = all(v < GOLDEN_RTOL for v in case['max_rel_diff'].values())
        results['cases'].append(case)
        print(f'{Nr:5d} {Nt:4d} {per_step:10.4g} {1e6*per_step/(Nr*Nr*Nt):13.4g} '
              f"{case['min_gap']:8.3f}" + (f"  matches reference: {case['matches_reference']}"
                                            if 'matches_reference' in case else ''), flush=True)
        if prof:
            print('    top functions by own time:')
            st = pstats.Stats(str(prof))
            rows = sorted(st.stats.items(), key=lambda kv: -kv[1][2])[:6]
            for (fname, line, func), (cc, nc, tt, ct, callers) in rows:
                print(f'      {tt:8.3f} s  {nc:8d} calls  {Path(fname).name}:{func}')
    out = BENCH/'results'/f"{args.label}{'_nominal' if nominal else ''}_{machine()['cpu'].replace(' ', '-')}.json"
    out.write_text(json.dumps(results, indent=1) + '\n')
    print(f'results written to {out.relative_to(REPO)}')


if __name__ == '__main__':
    main()
