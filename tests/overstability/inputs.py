#!/usr/bin/env python

#inputs.py
#
#by Joe Hahn and Claude Code, joe.hahn@jmh-datasciences.com, 27 November 2026.
#
#tests/overstability: the viscous overstability of axisymmetric epicyclic oscillations. The ring of tests/density_waves
#(1 <= r <= 1.02, total mass 3e-8, 800 streamlines of 32 particles, J2 = 0.01) without the satellite, kicked at t = 0
#into one small axisymmetric epicyclic wave of lattice wavenumber phi (radians per streamline): streamline j displaced
#radially by xi_0 cos(phi j) and kicked by v_r = kappa xi_0 sin(phi j) (tools/patches.py seed_mode). The wave grows or
#decays at the rate of notes/equations.md O7 on the lattice (K = sin^2 phi/Delta^2, x = X(phi)/kappa^2, P7). This run
#has shear viscosity only, so linear theory says it grows; tools/run_simulations.py also makes runs at other phi, with
#nu_b = nu_s and nu_b = 2 nu_s, and without ring gravity. Kappa falls outward, so the wave winds up, its wavenumber
#rising by 1.5 kappa t/r: about 1.5% over this run, which is kept short (100 orbits) by a viscosity ten times
#tests/density_waves'. Analysis: tests/overstability_analysis.py.
#
#define input input parameters

#set number of streamlins and particles per streamline (minumum is 3)
number_of_streamlines = 800
particles_per_streamline = 32
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs: an output every 2 orbits, for 100 orbits
dt = 0.5
timesteps_per_output = 25
total_number_of_outputs = 50

#ring radial width assuming circular orbits
radial_width = 0.02

#total ring mass
total_ring_mass = 3.0e-8

#ring's gravitation constant is usually G_ring=1 but set G_ring < 0 to turn off ring gravity.
#Also set fast_gravity=False since there is very little speed benefit when the fast_gravity approximation is used
G_ring = 1.0
fast_gravity = False

#ring kinematic shear and bulk viscosity, set < 0 to turn off
shear_viscosity = 1.0e-11
bulk_viscosity = -1.0

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

#the seeded wave: its lattice wavenumber phi (radians per streamline) and initial radial amplitude xi_0 in units of
#the streamline spacing
seed_phi = 0.5*3.141592653589793
seed_amplitude = 3.0e-4

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

#the seed, after any modifications
patches.seed_mode(seed_phi, seed_amplitude)
