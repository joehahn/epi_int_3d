#!/usr/bin/env python

#inputs.py
#
#by Joe Hahn, joe.hahn@jmh-datasciences.com, 18 October 2017.
#
#tests/density_waves: the ring of tests/bending_waves (1 <= r <= 1.02, total mass 3e-8, so sigma = 2.36e-7, 800
#streamlines of 32 particles) with its inner part straddling a satellite's m = 2 inner Lindblad resonance at
#r_L = 1.005, so x = (r - r_L)/r_L runs from -0.005 (its inner edge) to +0.0149 (epi_int_3d, 23 November 2026). The
#satellite, in the ring plane, is at full mass from t = 0, and the spiral density wave it launches at r_L propagates
#outward at the group velocity pi G sigma/kappa (notes/equations.md L6), reaching x = 0.01 at
#tau_0.01 = 0.01 r_L kappa/(pi G sigma), about 2130 orbits. The run stops at 1.5 tau_0.01, as tests/bending_waves does.
#The satellite's mass keeps |k xi_w| = 0.01 at the outer edge (L5, L7). Run with fast/epi_int_3d.py.
#
#define input input parameters

#set number of streamlins and particles per streamline (minumum is 3)
number_of_streamlines = 800
particles_per_streamline = 32
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.5
timesteps_per_output = 672
total_number_of_outputs = 61

#ring radial width assuming circular orbits
radial_width = 0.02

#total ring mass
total_ring_mass = 3.0e-8

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
#mass_final < 0 to turn off. Its radius puts the inner Lindblad resonance at r_L = 1.005.
satellite = {
    'resonance': 'horizontal',        #'vertical' (inner vertical resonance) or 'horizontal' (inner Lindblad)
    'm': 2,                           #number of spiral arms, m >= 2
    'mass_final': 2.0e-10,            #satellite mass/planet mass once grown
    'time_grow': -1.0,                #mass = mass_final*(1 - exp(-time/time_grow)); <= 0: mass_final from t = 0
    'r': 1.5902293522815265,          #satellite's orbit radius a_s (circular, in the ring plane)
    't': 0.0,                         #its initial mean longitude lambda_s
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
