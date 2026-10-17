#!/usr/bin/env python

#inputs.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 18 October 2017.
#
#define input input parameters

#set number of streamlins and particles per streamline (minumum is 3)
number_of_streamlines = 51
particles_per_streamline = 3
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.5
timesteps_per_output = 1260
total_number_of_outputs = 1000

#ring radial width assuming circular orbits
radial_width = 0.1

#total ring mass
total_ring_mass = 0.0

#ring's gravitation constant is usually G_ring=1 but set G_ring < 0 to turn off ring gravity.
#Also set fast_gravity=True for approximate gravity that is 2x faster and almost as accurate.
G_ring = -1.0
fast_gravity = False

#ring kinematic shear and bulk viscosity, set < 0 to turn off
shear_viscosity = -1.0
bulk_viscosity = shear_viscosity

#add fictitious torque at inner and outer streamlines, to oppose any radial spreading
confine_edges = False

#ring pressure scales with Toomre's Q_ring, set Q_ring < 0 to turn off
Q_ring = -1.0

#oblateness parameters, set J2=0 to turn off
Rp = 0.5
J2 = 0.01

#choose ringlet's initial orbits
initial_orbits = {
    'shape':'log-e',
    'e':(1.0e-8, 1.0e-2),
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
    print 'modified_params =', modified_params
    for key, val in modified_params.iteritems():
        exec(key + '=val')
