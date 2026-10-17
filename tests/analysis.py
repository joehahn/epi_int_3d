#analysis.py
#
#Diagnostics of the unperturbed, spreading and libration tests, reproducing the calculations
#in epi_int_lite's test notebooks (check_orbits.ipynb, spreading.ipynb, libration.ipynb).
#Used by the regression suite (test_regression.py) and by tools/compare_port.py.

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO/'source'))
from helper_fns import (coords2elem, restore_output, get_lambda, delta_f, surface_density,
    planeto2bary, unwrap_angle)
from libration import orbit_deltas, calculate_Hq


def read_inputs(test_folder, modified_params=None):
    """Execute a test's inputs.py and return its parameters as a dict."""
    params = {}
    code = (Path(test_folder)/'inputs.py').read_text()
    saved_argv = sys.argv
    sys.argv = ['inputs.py']        #keep inputs.py's argparse from reading pytest's arguments
    try:
        exec(code, params)
    finally:
        sys.argv = saved_argv
    params.update(modified_params or {})
    return params


def load(test_folder, output_folder=None, modified_params=None):
    """Return a test's parameters and its output r, t, vr, vt, times, monitor, a, e, wt, M."""
    p = read_inputs(test_folder, modified_params)
    folder = Path(output_folder) if output_folder else Path(test_folder)/p['output_folder']
    r, t, vr, vt, times, monitor = restore_output(str(folder))
    a, e, wt, M = coords2elem(p['J2'], p['Rp'], r, t, vr, vt)
    return p, dict(r=r, t=t, vr=vr, vt=vt, times=times, monitor=monitor, a=a, e=e, wt=wt, M=M)


def unperturbed_metrics(p, o):
    """Largest change in each particle's a, e and h = r*vt over the run (check_orbits.ipynb)."""
    h = o['r']*o['vt']
    return {'max_delta_a': np.abs(o['a'] - o['a'][0]).max(),
            'max_delta_e': np.abs(o['e'] - o['e'][0]).max(),
            'max_delta_h': np.abs(h - h[0]).max(),
            'e_min': o['e'][0].min(), 'e_max': o['e'][0].max()}


def spreading_angular_momentum(p, o):
    """Fractional change in the system's total barycentric angular momentum (spreading.ipynb)."""
    particle_mass = p['total_ring_mass']/(p['number_of_streamlines']*p['particles_per_streamline'])
    L = []
    for k in range(len(o['times'])):
        r_bc, t_bc, vr_bc, vt_bc, r_0, t_0, vr_0, vt_0 = planeto2bary(p['total_ring_mass'],
            o['r'][k], o['t'][k], o['vr'][k], o['vt'][k])
        L.append(np.sum(particle_mass*r_bc*vt_bc) + 1.0*r_0*vt_0)
    L = np.array(L)
    return (L - L[0])/L[0]


def spreading_profile(p, o, t_visc):
    """Simulated and analytic surface density at time t/tau_nu (spreading.ipynb, cell 11).

    Returns r, sigma_sim/sigma_bar and sigma_analytic/sigma_bar at the streamlines, where the
    analytic profile is the notebook's exp(-(x-1)^2/tau)/x^(3/4)/sqrt(4 pi tau), x = r/R."""
    times = o['times']
    tau_nu = p['radial_width']**2/(12*p['shear_viscosity'])
    k = np.where(times/tau_nu >= t_visc)[0][0]
    R = o['r'][0].mean()
    sd_init = p['total_ring_mass']/np.pi/R**2
    r2, t2, vr2, vt2 = o['r'][k], o['t'][k], o['vr'][k], o['vt'][k]
    lambda_ = get_lambda(p['total_ring_mass'], p['number_of_streamlines'], p['J2'], p['Rp'],
        r2, t2, vr2, vt2)
    sd = surface_density(lambda_, delta_f(r2, t2)).mean(axis=1)/sd_init
    r1 = r2.mean(axis=1)
    x = r1/R
    tau = 12.0*p['shear_viscosity']*times[k]/R**2
    sd_exp = np.exp(-(x - 1)**2/tau)/x**0.75/np.sqrt(4*np.pi*tau)
    return times[k]/tau_nu, r1, sd, sd_exp


def spreading_metrics(p, o, t_visc=(1.0, 2.5, 5.0, 10.0, 20.0, 50.0)):
    dL = spreading_angular_momentum(p, o)
    out = {'max_delta_L': np.abs(dL).max(), 'e_max': o['e'].max()}
    tau_nu = p['radial_width']**2/(12*p['shear_viscosity'])
    for tv in t_visc:
        if o['times'][-1]/tau_nu < tv:
            continue
        tv_actual, r1, sd, sd_exp = spreading_profile(p, o, tv)
        #compare interior streamlines; the edge streamlines' one-sided densities are cruder
        out[f'sigma_rel_err_t{tv:g}'] = np.abs(sd[1:-1]/sd_exp[1:-1] - 1).max()
    return out


def libration_metrics(p, o):
    """Libration of a narrow two-streamline ringlet (libration.ipynb)."""
    times = o['times']
    a_inner, a_outer, a_mean, da, e_inner, e_outer, e_mean, de, wt_inner, wt_outer, dwt = \
        orbit_deltas(times, o['r'], o['t'], o['vr'], o['vt'], o['a'], o['e'], o['wt'],
            p['J2'], p['Rp'])
    idx = times > 0
    dwt_eq = np.median(dwt[idx])
    de_eq = np.median(de[idx])
    H, q, e_prime, wt_prime = calculate_Hq(a_mean, e_mean, da, de, dwt)
    H_avg = H[idx].mean()
    Ix = (de - de_eq)*a_mean/da
    Iy = e_mean*a_mean*(dwt - dwt_eq)/da
    I = np.sqrt(Ix**2 + Iy**2)
    #libration frequency of Borderies, Goldreich & Tremaine, and from the simulation's phase
    n = np.sqrt(1.0/(a_mean**3).mean())
    omega_bgt = p['total_ring_mass']*((a_mean/da).mean()**2)*n*H_avg/np.pi
    phase = unwrap_angle(np.arctan2(Iy, Ix))
    omega = np.polyfit(times, phase, 1)[0]
    return {'omega': omega, 'omega_bgt': omega_bgt, 'omega_over_omega_bgt': omega/omega_bgt,
            'H_mean': H.mean(), 'I_avg': I[idx].mean(), 'q_std_over_mean': q.std()/q.mean(),
            'dwt_eq': dwt_eq, 'de_eq': de_eq, 'e_prime': e_prime[idx].mean(),
            'wt_prime': wt_prime[idx].mean(), 'e_min': o['e'].min(), 'e_max': o['e'].max()}
