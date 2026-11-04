#!/usr/bin/env python

#inputs.py
#
#by Joe Hahn, joe.hahn@jmh-datasciences.com, 18 October 2017.
#
#tests/bending_waves: a self-gravitating ring of uniform surface density whose outer part straddles a
#satellite's m = 2 inner vertical resonance at r_V = 1.0032, 80% of the way across the ring (epi_int_3d,
#4 November 2026). The satellite is at full mass from t = 0, and the spiral bending wave it launches at r_V
#propagates inward at the group velocity pi G sigma/nu (notes/equations.md B10), reaching the inner edge
#after about 4100 orbits; the run stops before then. The wave's amplitude, wavelength and phase are compared
#with Shu (1984) (notes/equations.md B7-B8, paper Section 4.3). Sized in notes/phase2c/bending_wave_theory.md
#section 3: about 11 streamlines per wavelength at the inner edge, 32 particles per streamline (16 per arm),
#and max |k Z| ~ 0.04. Run with fast/epi_int_3d.py.
#
#define input input parameters

#set number of streamlins and particles per streamline (minumum is 3)
number_of_streamlines = 270
particles_per_streamline = 32
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.5
timesteps_per_output = 628
total_number_of_outputs = 70

#ring radial width assuming circular orbits
radial_width = 0.004

#total ring mass
total_ring_mass = 1.0e-9

#ring's gravitation constant is usually G_ring=1 but set G_ring < 0 to turn off ring gravity.
#Also set fast_gravity=False since there is very little speed benefit when the fast_gravity approximation is used
G_ring = 1.0
fast_gravity = False

#ring kinematic shear and bulk viscosity, set < 0 to turn off
shear_viscosity = -1.0
bulk_viscosity = shear_viscosity

#add fictitious torques at inner and/or outer streamlines, to oppose any radial spreading
confine_inner_edge = False
confine_outer_edge = False

#ring pressure scales with Toomre's Q_ring, set Q_ring < 0 to turn off
Q_ring = -1.0

#oblateness parameters, set J2=0 to turn off
Rp = 0.5
J2 = 0.01

#choose ringlet's initial orbits
initial_orbits = {
    'shape':'circular',
}

#epi_int_3d only: a Fourier satellite (fast/satellite.py), only the resonant m-th term of its potential, set
#mass_final < 0 to turn off. Its radius puts the inner vertical resonance at r_V = 1.0032.
satellite = {
    'resonance': 'vertical',          #'vertical' (inner vertical resonance) or 'horizontal' (inner Lindblad)
    'm': 2,                           #number of spiral arms, m >= 2 for 'vertical'
    'mass_final': 2.5e-7,             #satellite mass/planet mass once grown
    'time_grow': -1.0,                #mass = mass_final*(1 - exp(-time/time_grow)); <= 0: mass_final from t = 0
    'r': 2.090314689118,              #satellite's orbit radius a_s (circular)
    't': 0.0,                         #its initial mean longitude lambda_s
    'inclination': 0.01,              #i_s
    'node': 0.0,                      #its initial longitude of ascending node Omega_s
}

#output folder
output_folder = 'output'

#the wrapped streamline-crossing check (tools/patches.py)
import os, sys
repo = os.path.abspath('.')
while not os.path.exists(os.path.join(repo, 'pyproject.toml')):
    repo = os.path.dirname(repo)
sys.path.insert(0, os.path.join(repo, 'tools'))
import patches
patches.apply()

#parse any optional commandline modifications to the above, note these modifications will be ignored by Jupyter
import argparse
import json
parser = argparse.ArgumentParser() 
parser.add_argument('-m', '--modified_params', type=str, dest='modified_params', required=False)
args, unknown_args = parser.parse_known_args()
modified_params = None
if (args.modified_params):
    modified_params_str = args.modified_params
    modified_params = json.loads(modified_params_str)
    print('modified_params =', modified_params)
    for key, val in modified_params.items():
        exec(key + '=val')
