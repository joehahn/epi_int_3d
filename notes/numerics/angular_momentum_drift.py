#angular_momentum_drift.py
#
#Fractional drift of the total angular momentum, and the L exchanged among streamlines, in the
#benchmark/runs/lamfix/{off,on} runs (notes/numerics/angular_momentum.md).
#
#Run from the repo root:  uv run python notes/numerics/angular_momentum_drift.py

import sys, pickle, runpy, numpy as np
sys.path.insert(0, 'source'); sys.path.insert(0, 'tests')
from helper_fns import planeto2bary
for run in ('off', 'on'):
    d = f'benchmark/runs/lamfix/{run}/'
    p = runpy.run_path(d + 'inputs.py')
    mon = pickle.load(open(d + 'output/monitor.pkl', 'rb'))
    p.update(mon['modified_params'])
    r, t, vr, vt, times = (np.load(d + f'output/{x}.npy') for x in ('r', 't', 'vr', 'vt', 'times'))
    Nr, Nt = r.shape[1:]
    mp = p['total_ring_mass']/(Nr*Nt)
    L = []
    for k in range(len(times)):
        r_bc, t_bc, vr_bc, vt_bc, r_0, t_0, vr_0, vt_0 = planeto2bary(p['total_ring_mass'], r[k], t[k], vr[k], vt[k])
        L.append(np.sum(mp*r_bc*vt_bc) + 1.0*r_0*vt_0)
        Lring = np.sum(mp*r[k]*vt[k])
    L = np.array(L)
    Lr = np.array([np.sum(mp*r[k]*vt[k]) for k in range(len(times))])
    print(run, f'Nr={Nr} Nt={Nt} orbits={times[-1]/(2*np.pi):.0f}',
          f'max|dL/L| total={np.abs((L-L[0])/L[0]).max():.2e}',
          f'max|dL_ring/L_ring|={np.abs((Lr-Lr[0])/Lr[0]).max():.2e}',
          f'final dL_ring/L_ring={(Lr[-1]-Lr[0])/Lr[0]:+.2e}', f'mass ratio={p["total_ring_mass"]:.1e}')
print('internal exchange vs net drift (ring L, units of L_ring):')
for run in ('off', 'on'):
    d = f'benchmark/runs/lamfix/{run}/'
    p = runpy.run_path(d + 'inputs.py'); p.update(pickle.load(open(d + 'output/monitor.pkl', 'rb'))['modified_params'])
    r, vt, times = (np.load(d + f'output/{x}.npy') for x in ('r', 'vt', 'times'))
    Nr, Nt = r.shape[1:]; mp = p['total_ring_mass']/(Nr*Nt)
    Lsl = (mp*r*vt).sum(axis=2)                       # per-streamline L(t)
    L0 = Lsl.sum(axis=1)[0]
    exch = np.abs(Lsl - Lsl[0]).sum(axis=1).max()/L0  # total |dL| moved among streamlines
    drift = np.abs(Lsl.sum(axis=1) - L0).max()/L0
    print(f'  {run}: max sum|dL_streamline| = {exch:.2e}, max |net dL| = {drift:.2e}, ratio = {drift/exch:.1e}')
