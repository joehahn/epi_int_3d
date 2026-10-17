#!/usr/bin/env python

#epi_int_lite.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 27 August 2017.
#
#this simulates the dynamical evolution of self-gravitating planetary rings.
#A portion of this code was written while drinking and singing at
#the Water Tank karaoke bar in northwest Austin TX, so buyer beware.

#read input parameters
import numpy as np
execfile('inputs.py')
print 'number_of_streamlines =', number_of_streamlines
print 'particles_per_streamline =', particles_per_streamline
print 'dt =', dt
print 'timesteps_per_output = ', timesteps_per_output
print 'total_number_of_outputs =', total_number_of_outputs
print 'radial_width =', radial_width
print 'total_ring_mass =', total_ring_mass
print 'ring gravitation constant =', G_ring
print 'fast_gravity =', fast_gravity
print 'shear_viscosity =', shear_viscosity
print 'bulk_viscosity =', bulk_viscosity
print 'confine_inner_edge =', confine_inner_edge
print 'confine_outer_edge =', confine_outer_edge
print "Toomre's Q_ring =", Q_ring
print 'Rp =', Rp
print 'J2 =', J2
print 'initial_orbits =', initial_orbits
print 'output_folder =', output_folder

#initialize orbits
from helper_fns import *
r, t, vr, vt, c, monitor = initialize_streamline(number_of_streamlines, particles_per_streamline,
    radial_width, total_ring_mass, G_ring, fast_gravity, shear_viscosity, bulk_viscosity, 
    confine_inner_edge, confine_outer_edge, Q_ring, J2, Rp, initial_orbits)

#prep for main loop
timestep = 0
number_of_outputs = 0
(rz, tz, vrz, vtz, timestepz) = ([r], [t], [vr], [vt], [timestep])

#evolve system using Chamber's (1993) 2nd order kick-drift-kick scheme 
print 'evolving system...'
while (number_of_outputs < total_number_of_outputs):
    #kick velocities forwards by timestep +dt/2
    r, t, vr, vt = velocity_kick(J2, Rp, G_ring, shear_viscosity, bulk_viscosity, c, total_ring_mass, number_of_streamlines, \
        r, t, vr, vt, dt/2.0, fast_gravity, confine_inner_edge, confine_outer_edge)
    timesteps_since_output = 0
    while (timesteps_since_output < timesteps_per_output):
        #kick coordinates to account for central body's motion about center of mass
        r, t, vr, vt = coordinate_kick(dt/2, total_ring_mass, r, t, vr, vt)
        #convert coordinates to elements
        a, e, wt, M = coords2elem(J2, Rp, r, t, vr, vt)
        #advance mean anomaly during drift step
        M = drift(a, M, J2, Rp, dt)
        #convert orbit elements to coordinates
        r, t, vr, vt = elem2coords(J2, Rp, a, e, wt, M)
        #kick coordinates to account for central body's motion about center of mass
        r, t, vr, vt = coordinate_kick(dt/2, total_ring_mass, r, t, vr, vt)    
        #kick velocities
        r, t, vr, vt = velocity_kick(J2, Rp, G_ring, shear_viscosity, bulk_viscosity, c, total_ring_mass, number_of_streamlines, \
            r, t, vr, vt, dt, fast_gravity, confine_inner_edge, confine_outer_edge)
        #updates
        timestep += 1
        timesteps_since_output += 1
        #update streamline monitor
        monitor = monitor_streamlines(monitor, r, t, timestep)
    #kick velocities backwards by timestep -dt/2
    r, t, vr, vt = velocity_kick(J2, Rp, G_ring, shear_viscosity, bulk_viscosity, c, total_ring_mass, number_of_streamlines, \
        r, t, vr, vt, -dt/2.0, fast_gravity, confine_inner_edge, confine_outer_edge)
    #save output
    rz, tz, vrz, vtz, timestepz = store_system(rz, tz, vrz, vtz, timestepz, r, t, vr, vt, total_ring_mass, timestep)
    number_of_outputs += 1
    #update display as needed
    if (20*number_of_outputs%total_number_of_outputs == 0):
         continue_sim = update_display(number_of_outputs, total_number_of_outputs, dt, timestep, monitor)
         if (continue_sim == False):
             break

#save results
timez = np.array(timestepz)*dt
save_output(rz, tz, vrz, vtz, timez, modified_params, monitor, output_folder)
print 'execution time (minutes) = ', (monitor['current_time'] - monitor['start_time'])/60.0
