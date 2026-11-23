#!/usr/bin/env python

#inputs.py
#
#by Joe Hahn, joe.hahn@jmh-datasciences.com, 18 October 2017.
#
#tests/lindblad_resonance: a zero-mass ring of test particles straddling a satellite's m = 2 inner Lindblad
#resonance at r_L = 1.02, midway between the two middle streamlines, while the satellite's mass grows slowly: the
#in-plane counterpart of tests/vertical_resonance (epi_int_3d, 23 November 2026). Each particle starts on a circular
#orbit and settles into the forced epicycle xi_r = -Psi/(r D), D = kappa_a^2 - omega_m^2, an eccentricity
#e_f = |Psi|/(r_L^2 D_L |x|) near r_L (notes/equations.md L3, L9), plus a small free epicycle that never decays.
#Run with fast/epi_int_3d.py.
#
#define input input parameters

#set number of streamlins and particles per streamline (minumum is 3)
number_of_streamlines = 80                     #even, so that r_L = 1.02 lies midway between two streamlines
particles_per_streamline = 6
number_of_particles = number_of_streamlines*particles_per_streamline

#set timestamp, timesteps per output, and total number of outputs
dt = 0.5
timesteps_per_output = 5000
total_number_of_outputs = 400                 #t = 1e6, 10 growth times

#ring radial width assuming circular orbits
radial_width = 0.04

#total ring mass
total_ring_mass = 0.0

#ring's gravitation constant is usually G_ring=1 but set G_ring < 0 to turn off ring gravity.
#Also set fast_gravity=False since there is very little speed benefit when the fast_gravity approximation is used
G_ring = -1.0
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
#mass_final < 0 to turn off. Its radius puts the inner Lindblad resonance at r_L = 1.02, mid-ring. Its mass keeps the
#forced radial displacement of the two streamlines beside r_L (|x| = 2.5e-4) at 2e-5, well under their spacing.
satellite = {
    'resonance': 'horizontal',        #'vertical' (inner vertical resonance) or 'horizontal' (inner Lindblad)
    'm': 2,                           #number of spiral arms, m >= 2
    'mass_final': 1.0e-8,             #satellite mass/planet mass once grown
    'time_grow': 1.0e5,               #mass = mass_final*(1 - exp(-time/time_grow))
    'r': 1.6141146120822665,          #satellite's orbit radius a_s (circular, in the ring plane)
    't': 0.0,                         #its initial mean longitude lambda_s
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
