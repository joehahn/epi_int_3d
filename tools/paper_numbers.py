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
from vertical_analysis import (load3d, forced_vertical, simulated_h, bending_wave, ring_sigma, wave_wavenumber, nu_a,
    resonance_frequencies, satellite_torque, ring_angular_momentum)
from satellite import Satellite

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
        text = str(round(value, decimals))
        if (decimals == 0):
            text = str(int(round(value)))
        return text
    mantissa = value/10.0**exponent
    mantissa_text = str(round(mantissa, digits - 1)).rstrip('0').rstrip('.')
    if (mantissa_text == '1'):
        return '10^{' + str(exponent) + '}'
    return mantissa_text + '\\times10^{' + str(exponent) + '}'


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
    numbers += [('vresOrbits', latex_number(orbits[-1], 1), 'the run\'s length in orbits')]
    #the average over the late run, as in paper/figures/vertical_resonance.py
    late_orbits = 1.0e4
    late = (orbits > late_orbits)
    numbers += [('vresAvgOrbits', latex_number(orbits[-1] - late_orbits, 1), 'orbits averaged over, from 1e4 to the end')]
    i_avg = o['inc'][late].mean(axis=(0, 2))
    x_avg = ((o['a'][late] - sat.r_V)/sat.r_V).mean(axis=(0, 2))
    off_resonance = (np.abs(x_avg) > 1.0e-6)
    i_forced = sat.mass_final*sat.inclination*sat.beta_V*sat.beta_V*sat.b/(   6.0*(sat.m - 1)*np.abs(x_avg[off_resonance])   )
    ratio = i_avg[off_resonance]/i_forced
    largest_departure = np.abs(ratio - 1.0).max()
    print('simulation/B12 =', np.round(ratio.min(), 4), 'to', np.round(ratio.max(), 4))
    numbers += [('vresAgreePercent', latex_number(100.0*largest_departure, 1, up=True),
        'largest |<i>/B12 - 1| in percent, rounded up')]
    #each particle's z against the exact solution B11, as in the notebook
    sat_unused, exact = forced_vertical(p, o)
    z_error = np.abs(o['z'] - exact['z']).max()
    print('max |z - B11| =', z_error)
    numbers += [('vresZerror', latex_number(z_error, 1, up=True), 'max |z - B11|, rounded up')]
    #nu_a^2 - omega_m^2 from each streamline's forced amplitude, against the exact value (paper/figures/resonance_location.py)
    sat_unused, x, simulated, exact_value = resonance_frequencies(p, o)
    off = (np.abs(x) > 1.0e-6)
    departure = np.abs(simulated[off]/exact_value[off] - 1.0).max()
    print('largest |simulated/exact nu^2 - omega_m^2 - 1| =', departure)
    numbers += [('vresDpercent', latex_number(100.0*departure, 1, up=True),
        'largest |(nu^2 - omega_m^2) from forced amplitude/exact - 1| in percent, rounded up')]
    return numbers


#the comparisons of tests/bending_waves/bending_waves.ipynb Sections 3, 4 and 7 for one run, behind the wave front
#at the run's end: the median fractional excess of the simulated wavenumber over the continuum's, the lattice's and
#the lattice-with-edges' (B7, B13), and the median wave amplitude over |Z_w| (B8)
def bending_wave_comparisons(p, o):
    sat, x, h = simulated_h(p, o, -1)
    sigma = ring_sigma(p)
    A_V = sat.mass_final*sat.A_V_per_mass
    Z_w = A_V*np.sqrt(sat.r_V/(sigma*sat.D_V))
    x_in = (1.0 - sat.r_V)/sat.r_V
    #behind the front at the run's end, as in the notebook
    x_lo = 0.56*x_in
    x_hi = 0.125*x_in
    h_forced = 1j*A_V/(sat.D_V*x)
    h_wave = h - h_forced
    x_mid = 0.5*(x[1:] + x[:-1])
    k_sim = np.angle(h_wave[1:]/h_wave[:-1])/np.diff(x)/sat.r_V
    idx = (x_mid < x_hi) & (x_mid > x_lo)
    k_cont, k_lat, k_edge = wave_wavenumber(p, sat, x_mid[idx])
    behind = (x < x_hi) & (x > x_lo)
    comparisons = {}
    comparisons['k_cont'] = np.median(k_sim[idx]/k_cont - 1.0)
    comparisons['k_lat'] = np.median(k_sim[idx]/k_lat - 1.0)
    comparisons['k_edge'] = np.median(k_sim[idx]/k_edge - 1.0)
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
    nu_V = nu_a(p['J2'], p['Rp'], sat.r_V)
    orbits = o['times']/(2.0*np.pi)
    numbers = []
    #the run's inputs
    numbers += [('bwNstreamlines', str(p['number_of_streamlines']), 'inputs.py number_of_streamlines')]
    numbers += [('bwNparticles', str(p['particles_per_streamline']), 'inputs.py particles_per_streamline')]
    numbers += [('bwRmax', latex_number(1.0 + p['radial_width'], 4), 'inputs.py 1 + radial_width')]
    numbers += [('bwRingMass', latex_number(p['total_ring_mass'], 1), 'inputs.py total_ring_mass')]
    numbers += [('bwSigma', latex_number(sigma, 2), 'ring_sigma, the uniform surface density')]
    numbers += [('bwM', str(sat.m), 'inputs.py satellite m')]
    numbers += [('bwIs', latex_number(sat.inclination, 1), 'inputs.py satellite inclination')]
    numbers += [('bwMs', latex_number(sat.mass_final, 2), 'inputs.py satellite mass_final')]
    numbers += [('bwAs', latex_number(p['satellite']['r'], 4), 'inputs.py satellite r')]
    numbers += [('bwRV', latex_number(sat.r_V, 5), 'fast/satellite.py r_V (B4)')]
    across = 100.0*(sat.r_V - 1.0)/p['radial_width']
    numbers += [('bwAcrossPercent', latex_number(across, 2), 'r_V\'s place across the ring, percent from the inner edge')]
    numbers += [('bwOrbits', latex_number(orbits[-1], 2), 'the run\'s length in orbits')]
    #the wave's crossing time from r_V to the inner edge at the group velocity (B10)
    c_g = np.pi*sigma/nu_V
    crossing_time = (sat.r_V - 1.0)/c_g
    numbers += [('bwCrossOrbits', latex_number(crossing_time/(2.0*np.pi), 2), 'crossing time to the inner edge in orbits (B10)')]
    numbers += [('bwRunFraction', latex_number(o['times'][-1]/crossing_time, 2), 'run length over the crossing time')]
    #resolution and linearity at the inner edge, from the continuum wavenumber (B7)
    x_in = (1.0 - sat.r_V)/sat.r_V
    k_in = sat.D_V*abs(x_in)/(2.0*np.pi*sigma)
    spacing = p['radial_width']/(p['number_of_streamlines'] - 1)
    numbers += [('bwPerWavelength', latex_number(2.0*np.pi/k_in/spacing, 2), 'streamlines per wavelength at the inner edge (B7)')]
    numbers += [('bwkZ', latex_number(k_in*Z_w, 1), 'largest |k Z_w|, at the inner edge')]
    #crests between the front at the run's end and r_V: the continuum phase D_V r_V x^2/(4 pi G sigma), over 2 pi
    x_front = -c_g*o['times'][-1]/sat.r_V
    crests = sat.D_V*sat.r_V*x_front*x_front/(4.0*np.pi*sigma)/(2.0*np.pi)
    numbers += [('bwCrests', latex_number(crests, 2), 'wavelengths between the front at the end and r_V')]
    #agreement with linear theory
    numbers += [('bwAmpPercent', latex_number(100.0*abs(comparisons['amplitude'] - 1.0), 1), 'median |h_wave|/|Z_w| - 1, percent')]
    numbers += [('bwKcontPercent', latex_number(100.0*comparisons['k_cont'], 2), 'median k/k_continuum - 1, percent (B7)')]
    numbers += [('bwKlatPercent', latex_number(100.0*comparisons['k_lat'], 2), 'median k/k_lattice - 1, percent (B13)')]
    numbers += [('bwKedgePercent', latex_number(100.0*comparisons['k_edge'], 1), 'median k/k_lattice+edges - 1, percent (B13)')]
    epsilon_half = np.sqrt(2.0*np.pi*sigma/(sat.r_V*sat.D_V))
    outside = (x > 2.0*epsilon_half)
    h_theory = bending_wave(sat, sigma, x)
    forced_excess = (np.abs(h - h_theory)/np.abs(h_theory))[outside].max()
    numbers += [('bwForcedPercent', latex_number(100.0*forced_excess, 1), 'largest |h - B8|/|B8| outside the resonance zone, percent')]
    radial = np.abs(o['r'] - o['r'][0]).max()
    numbers += [('bwRadial', latex_number(radial, 1), 'largest radial excursion of any particle')]
    numbers += [('bwSpacing', latex_number(spacing, 2), 'streamline spacing')]
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
    #the resolution runs of a ring half as wide (the notebook's Section 7), if they have been made
    narrow = {'radial_width': 0.004, 'total_ring_mass': 1.0e-9, 'number_of_streamlines': 270, 'number_of_particles': 8640,
        'satellite': dict(p['satellite'], r=2.090314689118), 'output_folder': 'output_narrow'}
    doubled = dict(narrow, number_of_streamlines=540, number_of_particles=17280, output_folder='output_narrow_ns540')
    if (os.path.exists(os.path.join(test_folder, 'output_narrow')) and os.path.exists(os.path.join(test_folder, 'output_narrow_ns540'))):
        p_narrow, o_narrow = load3d(test_folder, modified_params=narrow)
        comparisons_narrow = bending_wave_comparisons(p_narrow, o_narrow)[3]
        p_doubled, o_doubled = load3d(test_folder, modified_params=doubled)
        comparisons_doubled = bending_wave_comparisons(p_doubled, o_doubled)[3]
        numbers += [('bwNarrowKcontPercent', latex_number(100.0*comparisons_narrow['k_cont'], 2), 'half-width ring, 270 streamlines: k/k_continuum - 1, percent')]
        numbers += [('bwDoubledKcontPercent', latex_number(100.0*comparisons_doubled['k_cont'], 2), 'half-width ring, 540 streamlines: k/k_continuum - 1, percent')]
        numbers += [('bwDoubledKedgePercent', latex_number(100.0*comparisons_doubled['k_edge'], 1), 'half-width ring, 540 streamlines: k/k_lattice+edges - 1, percent')]
    else:
        print('the half-width runs output_narrow and output_narrow_ns540 were not found: their numbers are left out')
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
producers = {'vertical_resonance': vertical_resonance_numbers, 'bending_waves': bending_waves_numbers}

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
