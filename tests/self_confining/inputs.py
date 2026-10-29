#!/usr/bin/env python3

#inputs.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 18 October 2017.
#
#define input input parameters
#
#epi_int_lite's nominal self-confining ringlet (its runs/self_confining/nominal, Hahn et al. 2025), ported to
#Python 3. Run with fast/epi_int_3d.py and its default ring gravity (tools/run_simulations.py does); nominal.ipynb
#compares the result with the original epi_int_lite run's figures.

#set number of streamlins and particles per streamline
number_of_streamlines = 2
particles_per_streamline = 241
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.5
timesteps_per_output = 1700
total_number_of_outputs = 1000

#ring radial width assuming circular orbits
radial_width = 0.0001

#total ring mass
total_ring_mass = 1.0e-10

#ring's gravitation constant is usually G_ring=1 but set G_ring < 0 to turn off ring gravity.
#Also set fast_gravity=False since there is very little speed benefit when the fast_gravity approximation is used
G_ring = 1.0
fast_gravity = False

#ring kinematic shear and bulk viscosity, set < 0 to turn off
shear_viscosity = 1.0e-13
bulk_viscosity = 1.0*shear_viscosity

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
    'e':0.01,
    'e_prime':0.0,
    'w_prime':0.0
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
