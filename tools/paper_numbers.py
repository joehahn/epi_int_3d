#paper_numbers.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 10 November 2026.
#
#compute the numbers the paper quotes from the simulations and inputs that produce them, and write them
#as LaTeX macros into paper/numbers/<test>.tex, which paper/epi_int_3d.tex \inputs, so that no number is
#copied into the paper by hand. Each file is a node of the dependency chain (notes/chain.json), so a rerun
#that changes a number marks the claims that quote it as stale. Phase 3b of PLAN.md.
#
#usage, from the repo's top folder, after the test's simulation has run
#(uv run python tools/run_simulations.py tests/vertical_resonance):
#   uv run python tools/paper_numbers.py                       all tests
#   uv run python tools/paper_numbers.py vertical_resonance    one test
#(tests/bending_waves also reads its half-width resolution runs, if made: bending_waves.ipynb Section 7)
#
#macro names are letters only (TeX's rule), prefixed by the test: \vresNstreamlines. Numbers with a power
#of ten, e.g. $4\times10^{-11}$, must be used in math mode; plain integers and decimals work anywhere.

import argparse
import datetime
import math
import os
import sys

#paths
repo_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
numbers_folder = os.path.join(repo_folder, 'paper', 'numbers')
sys.path.insert(0, os.path.join(repo_folder, 'source'))
sys.path.insert(0, os.path.join(repo_folder, 'fast'))
sys.path.insert(0, os.path.join(repo_folder, 'tests'))
import numpy as np
from vertical_analysis import (load3d, simulated_h, bending_wave, ring_sigma, wave_wavenumber, nu_a,
    resonance_D, satellite_torque, ring_angular_momentum, late_outputs, resolution_runs, X_STEADY, tau_001,
    resonance_zone, wave_front, front_speed)
from satellite import Satellite
import density_analysis as dw

today = datetime.date.today().isoformat()


#x as LaTeX, to the given significant digits, rounded up if it is an upper bound;
#powers of ten below 1e-3 or from 1e4 up, and a mantissa of 1 dropped: 10^{4}, 4\times10^{-11}
def latex_number(x, digits=2, up=False):
    if (x == 0.0):
        return '0'
    exponent = int(math.floor(math.log10(abs(x))))
    scale = 10.0**(exponent - digits + 1)
    if (up):
        mantissa_digits = math.ceil(abs(x)/scale - 1.0e-9)
    else:
        mantissa_digits = round(abs(x)/scale)
    value = math.copysign(mantissa_digits*scale, x)
    #rounding up can carry into the next power of ten
    exponent = int(math.floor(math.log10(abs(value))))
    if ((exponent >= -3) and (exponent < 4)):
        decimals = max(0, digits - 1 - exponent)
        text = f'{value:.{decimals}f}'                #keeps trailing zeros: 1.00, not 1.0
        if (decimals == 0):
            text = str(int(round(value)))
        return text
    mantissa = value/10.0**exponent
    mantissa_text = str(round(mantissa, digits - 1)).rstrip('0').rstrip('.')
    if (mantissa_text == '1'):
        return '10^{' + str(exponent) + '}'
    return mantissa_text + '\\times10^{' + str(exponent) + '}'


#the wavenumber of a free wave on a uniform lattice of streamlines spaced Delta apart, from |D| = X(k Delta), X of
#the paper's resolution appendix (P7; B13 without its missing-edge term): (4 G sigma/Delta)(pi phi/2 - phi^2/4) = |D|,
#phi = k Delta, whose root below pi is phi = pi - (pi^2 - |D| Delta/(G sigma))^(1/2)
def lattice_wavenumber(D_abs, sigma, spacing, G=1.0):
    return (np.pi - np.sqrt(np.pi**2 - D_abs*spacing/(G*sigma)))/spacing


#R = 2 pi^2 G sigma/(kappa^2 Delta), the paper's resolution parameter (P7), with kappa at the ring's middle
def resolution_R(p, G=1.0):
    from helper_fns import Kappa
    spacing = p['radial_width']/(p['number_of_streamlines'] - 1)
    kappa = Kappa(p['J2'], p['Rp'], 1.0 + 0.5*p['radial_width'])
    return 2.0*np.pi**2*G*ring_sigma(p)/(kappa**2*spacing)


#the largest sawtooth along any streamline at output k, the alternating part (-1)^j of its particles' radii and
#heights, sorted by longitude, after a fit of harmonics m <= 3 is removed (the m ~ N_t/2 disturbance of the
#resolution appendix)
def sawtooth(o, k=-1):
    largest = 0.0
    for key in ('r', 'z'):
        if key not in o:
            continue
        for j in range(o['t'].shape[1]):
            i = np.argsort(o['t'][k, j])
            t, y = o['t'][k, j, i], o[key][k, j, i]
            A = np.column_stack([np.ones_like(t)] + [f(m*t) for m in range(1, 4) for f in (np.cos, np.sin)])
            c, *_ = np.linalg.lstsq(A, y, rcond=None)
            res = y - A@c
            largest = max(largest, abs(np.mean(res*(-1.0)**np.arange(len(res)))))
    return largest


#the massless ring at a vertical resonance, Section 4.2 and Fig. vres
def vertical_resonance_numbers():
    test_folder = os.path.join(repo_folder, 'tests', 'vertical_resonance')
    p, o = load3d(test_folder)
    sat = Satellite(p['satellite'], p['J2'], p['Rp'])
    orbits = o['times']/(2.0*np.pi)
    numbers = []
    #the run's inputs
    numbers += [('vresJtwo', latex_number(p['J2'], 1), 'inputs.py J2')]
    numbers += [('vresRp', latex_number(p['Rp'], 1), 'inputs.py Rp')]
    numbers += [('vresNstreamlines', str(p['number_of_streamlines']), 'inputs.py number_of_streamlines')]
    numbers += [('vresNparticles', str(p['particles_per_streamline']), 'inputs.py particles_per_streamline')]
    numbers += [('vresRmax', latex_number(1.0 + p['radial_width'], 3), 'inputs.py 1 + radial_width')]
    numbers += [('vresM', str(sat.m), 'inputs.py satellite m')]
    numbers += [('vresIs', latex_number(sat.inclination, 1), 'inputs.py satellite inclination')]
    numbers += [('vresAs', latex_number(p['satellite']['r'], 4), 'inputs.py satellite r')]
    numbers += [('vresRV', latex_number(sat.r_V, 3), 'fast/satellite.py r_V (B4)')]
    numbers += [('vresMs', latex_number(sat.mass_final, 1), 'inputs.py satellite mass_final')]
    numbers += [('vresTau', latex_number(sat.time_grow, 1), 'inputs.py satellite time_grow')]
    numbers += [('vresTauOrbits', latex_number(sat.time_grow/(2.0*np.pi), 2), 'time_grow/2 pi')]
    numbers += [('vresOrbits', latex_number(orbits[-1], 2), 'the run\'s length in orbits')]
    #the average over the late run, as in paper/figures/vertical_resonance.py
    late = late_outputs(p, o)
    numbers += [('vresAvgOrbits', latex_number(orbits[-1] - orbits[late][0], 1),
        'orbits averaged over, after 6 growth times (late_outputs) to the end')]
    i_avg = o['inc'][late].mean(axis=(0, 2))
    x_avg = ((o['a'][late] - sat.r_V)/sat.r_V).mean(axis=(0, 2))
    i_forced = sat.mass_final*sat.inclination*sat.beta_V*sat.beta_V*sat.b/(   6.0*(sat.m - 1)*np.abs(x_avg)   )
    ratio = i_avg/i_forced
    largest_departure = np.abs(ratio - 1.0).max()
    print('simulation/B12 =', np.round(ratio.min(), 4), 'to', np.round(ratio.max(), 4))
    numbers += [('vresAgreePercent', latex_number(100.0*largest_departure, 1, up=True),
        'largest |<i>/B12 - 1| in percent, rounded up')]
    #D = nu_a^2 - omega_m^2 (Eqn. 25) at each streamline's semimajor axis: it vanishes between the two streamlines
    #nearest r_V, which have the largest inclinations (Fig. 2a)
    sat_unused, D = resonance_D(p, sat.r_V*(1 + x_avg))
    j = np.where(np.diff(np.sign(D)) != 0)[0][0]
    peak = set(np.argsort(i_avg)[-2:])
    assert peak == {j, j + 1}, 'the largest inclinations are not at the streamlines either side of D = 0'
    numbers += [('vresXnear', latex_number(np.abs(x_avg[[j, j + 1]]).max(), 2),
        'largest |x| of the two streamlines either side of D = 0, which have the largest <i>')]
    return numbers


#the comparisons of tests/bending_waves/bending_waves.ipynb Sections 3, 4 and 7 for one run, behind the wave front
#at the run's end: the median fractional excess of the simulated wavenumber over the continuum's, the lattice's and
#the lattice-with-edges' (B7, B13), and the median wave amplitude over |Z_w| (B8)
def bending_wave_comparisons(p, o):
    sat, x, h = simulated_h(p, o, -1)
    sigma = ring_sigma(p)
    A_V = sat.mass_final*sat.A_V_per_mass
    Z_w = A_V*np.sqrt(sat.r_V/(sigma*sat.D_V))
    #where the wave is steady at the run's end, X_STEADY < x outside the resonance zone, as in the notebook
    x_lo = X_STEADY
    x_hi = -resonance_zone(p, sat)
    h_forced = 1j*A_V/(sat.D_V*x)
    h_wave = h - h_forced
    x_mid = 0.5*(x[1:] + x[:-1])
    k_sim = np.angle(h_wave[1:]/h_wave[:-1])/np.diff(x)/sat.r_V
    idx = (x_mid < x_hi) & (x_mid > x_lo)
    k_cont = wave_wavenumber(p, sat, x_mid[idx])
    behind = (x < x_hi) & (x > x_lo)
    comparisons = {}
    comparisons['k_cont'] = np.median(k_sim[idx]/k_cont - 1.0)
    spacing = p['radial_width']/(p['number_of_streamlines'] - 1)
    k_lat = lattice_wavenumber(sat.D_V*np.abs(x_mid[idx]), sigma, spacing)
    comparisons['k_lattice'] = np.median(k_lat/k_cont - 1.0)
    comparisons['amplitude'] = np.median(np.abs(h_wave[behind]))/Z_w
    return sat, x, h, comparisons


#the self-gravitating ring's bending wave, Section 4.3 and Fig. bwave
def bending_waves_numbers():
    test_folder = os.path.join(repo_folder, 'tests', 'bending_waves')
    p, o = load3d(test_folder)
    sat, x, h, comparisons = bending_wave_comparisons(p, o)
    sigma = ring_sigma(p)
    A_V = sat.mass_final*sat.A_V_per_mass
    Z_w = A_V*np.sqrt(sat.r_V/(sigma*sat.D_V))
    orbits = o['times']/(2.0*np.pi)
    numbers = []
    #the run's inputs
    numbers += [('bwNstreamlines', str(p['number_of_streamlines']), 'inputs.py number_of_streamlines')]
    numbers += [('bwNparticles', str(p['particles_per_streamline']), 'inputs.py particles_per_streamline')]
    numbers += [('bwRmax', latex_number(1.0 + p['radial_width'], 3), 'inputs.py 1 + radial_width')]
    numbers += [('bwRingMass', latex_number(p['total_ring_mass'], 1), 'inputs.py total_ring_mass')]
    numbers += [('bwSigma', latex_number(sigma, 2), 'ring_sigma, the uniform surface density')]
    numbers += [('bwM', str(sat.m), 'inputs.py satellite m')]
    numbers += [('bwIs', latex_number(sat.inclination, 1), 'inputs.py satellite inclination')]
    numbers += [('bwMs', latex_number(sat.mass_final, 2), 'inputs.py satellite mass_final')]
    numbers += [('bwAs', latex_number(p['satellite']['r'], 4), 'inputs.py satellite r')]
    numbers += [('bwRV', latex_number(sat.r_V, 4), 'fast/satellite.py r_V (B4)')]
    across = 100.0*(sat.r_V - 1.0)/p['radial_width']
    numbers += [('bwAcrossPercent', latex_number(across, 2), 'r_V\'s place across the ring, percent from the inner edge')]
    numbers += [('bwOrbits', latex_number(orbits[-1], 2), 'the run\'s length in orbits')]
    #tau_0.01, the wave's travel time from r_V to x = -0.01 at the group velocity (B10)
    tau = tau_001(p, sat)
    numbers += [('bwTauOrbits', latex_number(tau/(2.0*np.pi), 2), 'tau_0.01 in orbits (B10)')]
    numbers += [('bwRunFraction', latex_number(o['times'][-1]/tau, 2), 'run length over tau_0.01')]
    x_in = (1.0 - sat.r_V)/sat.r_V
    numbers += [('bwXin', latex_number(x_in, 2), 'the ring inner edge\'s x')]
    x_front = wave_front(p, o, -1)
    numbers += [('bwXfront', latex_number(x_front, 2), 'the measured wave front\'s x at the run\'s end')]
    numbers += [('bwFrontSpeed', latex_number(front_speed(p, o), 2), 'front_speed: the front\'s mean speed over 0.25 < t/tau_0.01 < 1.25, in units of |c_g|')]
    #resolution and linearity at the inner edge, from the continuum wavenumber (B7)
    k_in = sat.D_V*abs(x_in)/(2.0*np.pi*sigma)
    spacing = p['radial_width']/(p['number_of_streamlines'] - 1)
    numbers += [('bwPerWavelength', latex_number(2.0*np.pi/k_in/spacing, 2), 'streamlines per wavelength at the inner edge (B7)')]
    numbers += [('bwkZ', latex_number(k_in*Z_w, 1), 'largest |k Z_w|, at the inner edge')]
    #crests between the measured front at the run's end and r_V: the continuum phase D_V r_V x^2/(4 pi G sigma), over 2 pi
    crests = sat.D_V*sat.r_V*x_front*x_front/(4.0*np.pi*sigma)/(2.0*np.pi)
    numbers += [('bwCrests', latex_number(crests, 2), 'wavelengths between the front at the end and r_V')]
    #agreement with linear theory
    numbers += [('bwAmpPercent', latex_number(100.0*abs(comparisons['amplitude'] - 1.0), 1), 'median |h_wave|/|Z_w| - 1, percent')]
    numbers += [('bwKcontPercent', latex_number(100.0*comparisons['k_cont'], 2), 'median k/k_continuum - 1, percent (B7)')]
    numbers += [('bwKlatticePercent', latex_number(100.0*comparisons['k_lattice'], 2),
        'median k_lattice/k_continuum - 1 over the same stretch, percent (lattice_wavenumber)')]
    numbers += [('bwR', latex_number(resolution_R(p), 2), 'resolution_R, 2 pi^2 G sigma/(kappa^2 Delta)')]
    numbers += [('bwSawtooth', latex_number(sawtooth(o)/Z_w, 1, up=True), 'largest sawtooth at the run\'s end over |Z_w|, rounded up')]
    #the same ring at 400 and 1600 streamlines (resolution_runs, Fig. 4): k converges toward the continuum's
    for (p_res, o_res), name in zip(resolution_runs(test_folder)[::2], ('Coarse', 'Fine')):
        c_res = bending_wave_comparisons(p_res, o_res)[3]
        k_res = c_res['k_cont']
        numbers += [('bwKlatticePercent' + name, latex_number(100.0*c_res['k_lattice'], 2),
            f"{p_res['number_of_streamlines']} streamlines: median k_lattice/k_continuum - 1, percent")]
        numbers += [('bwKcontPercent' + name, latex_number(100.0*k_res, 2),
            f"{p_res['number_of_streamlines']} streamlines: median k/k_continuum - 1, percent (B7)")]
        numbers += [('bwNstreamlines' + name, str(p_res['number_of_streamlines']), 'resolution_runs')]
        numbers += [('bwFrontSpeed' + name, latex_number(front_speed(p_res, o_res), 2),
            f"{p_res['number_of_streamlines']} streamlines: front_speed, in units of |c_g|")]
    epsilon_half = np.sqrt(2.0*np.pi*sigma/(sat.r_V*sat.D_V))
    outside = (x > 2.0*epsilon_half)
    h_theory = bending_wave(sat, sigma, x)
    forced_excess = (np.abs(h - h_theory)/np.abs(h_theory))[outside].max()
    numbers += [('bwForcedPercent', latex_number(100.0*forced_excess, 1), 'largest |h - B8|/|B8| outside the resonance zone, percent')]
    forced_median = np.median((np.abs(h - h_theory)/np.abs(h_theory))[outside])
    numbers += [('bwForcedMedianPercent', latex_number(100.0*forced_median, 1), 'median |h - B8|/|B8| outside the resonance zone, percent')]
    #the satellite's torque against linear theory, B9, over the run's last half (paper/figures/torque.py)
    torque_theory = -np.pi**2*sat.m*sigma*sat.r_V**2*A_V**2/sat.D_V
    torque = np.array([satellite_torque(p, o, k).sum() for k in range(len(o['times']))])/torque_theory
    late = torque[len(torque)//2:]
    numbers += [('bwTorqueRatio', latex_number(late.mean(), 3), 'mean T/T_V over the last half of the run (B9)')]
    numbers += [('bwTorqueScatter', latex_number(100.0*np.abs(late - 1.0).max(), 2, up=True), 'largest |T/T_V - 1| over the last half, percent, rounded up')]
    #the ring's angular momentum about z (paper/figures/angular_momentum.py)
    L = ring_angular_momentum(p, o)
    numbers += [('bwDeltaL', latex_number(np.abs((L - L[0])/L[0]).max(), 2, up=True), 'largest |Delta L_z/L_z| over the run, rounded up')]
    print('amplitude =', np.round(comparisons['amplitude'], 4), '   crests =', np.round(crests, 2))
    return numbers


#the massless ring at a Lindblad resonance, Section 4.2.3 (no figure: tests/lindblad_resonance/lindblad_resonance.ipynb test 9)
def lindblad_resonance_numbers():
    test_folder = os.path.join(repo_folder, 'tests', 'lindblad_resonance')
    p, o = dw.load3d(test_folder)
    sat = dw.satellite(p)
    orbits = o['times']/(2.0*np.pi)
    numbers = []
    numbers += [('lresM', str(sat.m), 'inputs.py satellite m')]
    numbers += [('lresMs', latex_number(sat.mass_final, 1), 'inputs.py satellite mass_final')]
    numbers += [('lresAs', latex_number(p['satellite']['r'], 4), 'inputs.py satellite r')]
    numbers += [('lresRL', latex_number(sat.r_L, 3), 'fast/satellite.py r_L (L2)')]
    numbers += [('lresOrbits', latex_number(orbits[-1], 2), 'the run\'s length in orbits')]
    #each streamline's eccentricity, averaged over its particles and the late run (late_outputs), against L9's e_f
    late = dw.late_outputs(p, o)
    numbers += [('lresAvgOrbits', latex_number(orbits[-1] - orbits[late][0], 1),
        'orbits averaged over, after 6 growth times (late_outputs) to the end')]
    e_avg = o['e'][late].mean(axis=(0, 2))
    x_avg = ((o['a'][late] - sat.r_L)/sat.r_L).mean(axis=(0, 2))
    e_f = (sat.mass_final/1.0)*sat.beta_L*abs(sat.beta_L*sat.db + 2*sat.m*sat.b)/(3.0*(sat.m - 1)*np.abs(x_avg))
    ratio = e_avg/e_f
    print('simulation/L9 =', np.round(ratio.min(), 4), 'to', np.round(ratio.max(), 4))
    numbers += [('lresAgreePercent', latex_number(100.0*np.abs(ratio - 1.0).max(), 1, up=True),
        'largest |<e>/e_f - 1| in percent (L9), rounded up')]
    #D = kappa_a^2 - omega_m^2 at each streamline's semimajor axis changes sign between the two streamlines that have
    #the largest eccentricities
    D = dw.resonance_D(p, sat.r_L*(1 + x_avg))[1]
    j = np.where(np.diff(np.sign(D)) != 0)[0][0]
    assert set(np.argsort(e_avg)[-2:]) == {j, j + 1}, 'the largest eccentricities are not either side of D = 0'
    numbers += [('lresXnear', latex_number(np.abs(x_avg[[j, j + 1]]).max(), 2),
        'largest |x| of the two streamlines either side of D = 0, which have the largest <e>')]
    return numbers


#the self-gravitating ring's density wave, Section 4.2.4 and Figs. dwave, dwtorque
def density_waves_numbers():
    test_folder = os.path.join(repo_folder, 'tests', 'density_waves')
    p, o = dw.load3d(test_folder)
    sat, x, xi = dw.simulated_xi(p, o, -1)
    sigma = dw.ring_sigma(p)
    Z_w = dw.wave_amplitude(p, sat)
    orbits = o['times']/(2.0*np.pi)
    numbers = []
    numbers += [('dwM', str(sat.m), 'inputs.py satellite m')]
    numbers += [('dwMs', latex_number(sat.mass_final, 1), 'inputs.py satellite mass_final')]
    numbers += [('dwAs', latex_number(p['satellite']['r'], 4), 'inputs.py satellite r')]
    numbers += [('dwRL', latex_number(sat.r_L, 4), 'fast/satellite.py r_L (L2)')]
    across = 100.0*(sat.r_L - 1.0)/p['radial_width']
    numbers += [('dwAcrossPercent', latex_number(across, 2), 'r_L\'s place across the ring, percent from the inner edge')]
    tau = dw.tau_001(p, sat)
    numbers += [('dwTauOrbits', latex_number(tau/(2.0*np.pi), 2), 'tau_0.01 in orbits (L6)')]
    numbers += [('dwOrbits', latex_number(orbits[-1], 2), 'the run\'s length in orbits')]
    numbers += [('dwRunFraction', latex_number(o['times'][-1]/tau, 2), 'run length over tau_0.01')]
    x_out = (1.0 + p['radial_width'] - sat.r_L)/sat.r_L
    numbers += [('dwXout', latex_number(x_out, 2), 'the ring outer edge\'s x')]
    x_front = dw.wave_front(p, o, -1)
    numbers += [('dwXfront', latex_number(x_front, 2), 'the measured wave front\'s x at the run\'s end')]
    k_out = dw.wave_wavenumber(p, sat, x_out)
    numbers += [('dwkZ', latex_number(k_out*Z_w, 1), 'largest |k Z_w|, at the outer edge (L5, L7)')]
    #where the wave is steady at the end: outside the resonance zone, x < X_STEADY
    st = dw.steady(p, sat, x)
    xi_wave = xi - dw.L_L(sat)/x
    numbers += [('dwAmpPercent', latex_number(100.0*abs(np.median(np.abs(xi_wave[st]))/Z_w - 1.0), 1),
        'median |xi - L_L/x|/|Z_w| - 1 over the steady stretch, percent')]
    x_mid, k_sim = dw.wavenumber_profile(x, xi_wave, sat.r_L)
    st_mid = dw.steady(p, sat, x_mid)
    numbers += [('dwKcontPercent', latex_number(100.0*np.median(k_sim[st_mid]/dw.wave_wavenumber(p, sat, x_mid[st_mid]) - 1.0), 2),
        'median k/k_continuum - 1 over the steady stretch, percent (L5)')]
    k_lat = lattice_wavenumber(sat.D_L*np.abs(x_mid[st_mid]), sigma, p['radial_width']/(p['number_of_streamlines'] - 1))
    numbers += [('dwKlatticePercent', latex_number(100.0*np.median(k_lat/dw.wave_wavenumber(p, sat, x_mid[st_mid]) - 1.0), 2),
        'median k_lattice/k_continuum - 1 over the steady stretch, percent (lattice_wavenumber)')]
    numbers += [('dwSawtooth', latex_number(sawtooth(o)/Z_w, 1, up=True), 'largest sawtooth at the run\'s end over |Z_w|, rounded up')]
    numbers += [('dwNparticles', str(p['particles_per_streamline']), 'inputs.py particles_per_streamline')]
    assert (k_sim[st_mid] > 0).all(), 'the wave is not trailing (k > 0) everywhere in the steady stretch'
    inside = x < -2.0*np.sqrt(2.0*np.pi*sigma/(sat.r_L*sat.D_L))
    xi_theory = dw.density_wave(sat, sigma, x)
    dev = np.abs(xi - xi_theory)/np.abs(xi_theory)
    numbers += [('dwForcedMedianPercent', latex_number(100.0*np.median(dev[inside]), 1),
        'median |xi - L7|/|L7| inside the resonance zone (x < -2 eps^1/2), percent')]
    numbers += [('dwForcedPercent', latex_number(100.0*dev[inside].max(), 1),
        'largest |xi - L7|/|L7| inside the resonance zone, percent')]
    #the satellite's torque against linear theory, L8, over the run's last half
    torque = np.array([dw.satellite_torque(p, o, k).sum() for k in range(len(o['times']))])
    ratio = torque/dw.torque_theory(p, sat)
    late = ratio[len(ratio)//2:]
    numbers += [('dwTorqueRatio', latex_number(late.mean(), 3), 'mean T/T_L over the last half of the run (L8)')]
    numbers += [('dwTorqueScatter', latex_number(100.0*np.abs(late - 1.0).max(), 2, up=True),
        'largest |T/T_L - 1| over the last half, percent, rounded up')]
    return numbers


#write one test's macros
def write_numbers(test, numbers):
    path = os.path.join(numbers_folder, test + '.tex')
    lines = ['%' + test + '.tex', '%', '%written by tools/paper_numbers.py; do not edit by hand', '']
    for name, value, source in numbers:
        lines += ['\\newcommand{\\' + name + '}{' + value + '}' + '   %' + source]
        print(name, '=', value)
    with open(path, 'w') as fp:
        fp.write('\n'.join(lines) + '\n')
    print('wrote', os.path.relpath(path, repo_folder))


#the tests that supply numbers, and their functions
producers = {'vertical_resonance': vertical_resonance_numbers, 'bending_waves': bending_waves_numbers,
             'lindblad_resonance': lindblad_resonance_numbers, 'density_waves': density_waves_numbers}

#parse the command line
parser = argparse.ArgumentParser()
parser.add_argument('tests', type=str, nargs='*', default=list(producers.keys()), help='tests to compute numbers for')
args = parser.parse_args()

#compute and write each test's numbers
os.makedirs(numbers_folder, exist_ok=True)
for test in args.tests:
    print('test =', test)
    numbers = producers[test]()
    write_numbers(test, numbers)
