#!/usr/bin/env python

#inputs.py
#
#by Joe Hahn, joe.hahn@jmh-datasciences.com, 18 October 2017.
#
#define input input parameters
#
#Exploration (epi_int_3d, Phase 0b): candidate nominal ringlet for this project. Width 0.01,
#e = 1e-3 with periapses twisted by pi from inner to outer edge, viscosity chosen to double the
#width over 1e4 orbits, and a ring mass large enough that self-gravity, not J2, controls the
#twist (it librates instead of winding up). Runs override these values with -m (see launch.sh and
#notes/phase0b/nominal.md).

#set number of streamlins and particles per streamline
number_of_streamlines = 10
particles_per_streamline = 40
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.5
timesteps_per_output = 629
total_number_of_outputs = 200

#ring radial width assuming circular orbits
radial_width = 1.0e-2

#total ring mass
total_ring_mass = 3.0e-8

#ring's gravitation constant is usually G_ring=1 but set G_ring < 0 to turn off ring gravity.
#Also set fast_gravity=True for approximate gravity that is 2x faster and almost as accurate.
G_ring = 1.0
fast_gravity = False

#ring kinematic shear and bulk viscosity, set < 0 to turn off
shear_viscosity = 6.0e-11
bulk_viscosity = shear_viscosity

#add fictitious torques at inner and/or outer streamlines, to oppose any radial spreading
confine_inner_edge = False
confine_outer_edge = False

#ring pressure scales with Toomre's Q_ring, set Q_ring < 0 to turn off
Q_ring = -1.0

#oblateness parameters
Rp = 0.5
J2 = 0.01

#choose ringlet's initial orbits
initial_orbits = {
    'shape':'eccentric',
    'e':0.001,
    'e_prime':0.0,
    'w_prime':0.315730       #periapse twist of pi across the ringlet: w_prime = pi e a_mean/radial_width
}

#output folder
output_folder = 'output'

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
