#!/usr/bin/env python

#inputs.py
#
#epi_int_3d's nominal ringlet, for testing and as the starting point for the bending-wave runs (see
#notes/phase0b/nominal.md): a narrow eccentric ringlet of 1000 streamlines with aligned periapses, light
#enough to be well resolved (lambda_crit/2Delta ~ 0.09) yet massive enough that its self-gravity, not J2,
#controls its periapse twist, with tapered edges. Viscosity and pressure are off, as in the first
#bending-wave runs. It stays smooth for 5e4 orbits, its twist oscillating within ~0.5 rad where J2 alone
#would wind it up by ~5 rad. About 2.8 h with fast/epi_int_fast.py. tools/run_simulations.py runs it with
#fast/epi_int_3d.py and its default ring gravity (since 2026-10-24; notes/phase0b/nominal.md used the reference's).

#set number of streamlins and particles per streamline
number_of_streamlines = 1000
particles_per_streamline = 20
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs (5e4 orbits, one output per 100)
dt = 0.5
timesteps_per_output = 1258
total_number_of_outputs = 500

#ring radial width assuming circular orbits
radial_width = 0.00125

#total ring mass
total_ring_mass = 3.0e-11

#ring's gravitation constant is usually G_ring=1 but set G_ring < 0 to turn off ring gravity.
G_ring = 1.0
fast_gravity = False

#ring kinematic shear and bulk viscosity, set < 0 to turn off
shear_viscosity = 0.0
bulk_viscosity = 0.0

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
    'e':1.0e-5,
    'e_prime':0.0,
    'w_prime':0.0            #no initial periapse twist
}

#output folder
output_folder = 'output'

#epi_int_lite has no edge taper: apply it, and the longitude-wrapped streamline-crossing check
#(issue 1 in notes/epi_int_lite_issues.md), from tools/patches.py; source/ is unchanged
import os, sys
repo = os.path.abspath('.')
while not os.path.exists(os.path.join(repo, 'pyproject.toml')):
    repo = os.path.dirname(repo)
sys.path.insert(0, os.path.join(repo, 'tools'))
import patches
patches.apply(taper=0.2)

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
