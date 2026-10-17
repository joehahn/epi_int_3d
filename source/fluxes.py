#fluxes.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 6 June 2018
#
#compute viscous angular momentum flux and luminosity through ringlet's central streamline

import numpy as np
from helper_fns import get_lambda, accelerations

#compute angular momentum flux and luminosity
def calculate_angular_momentum_flux(total_ring_mass, number_of_streamlines, particles_per_streamline, J2, Rp, G_ring, 
        shear_viscosity, bulk_viscosity, c, r, t, vr, vt, wt, times, fast_gravity, confine_inner_edge, confine_outer_edge):
    angular_momentum_flux_list = []
    angular_momentum_luminosity_list = []
    for t_idx, _ in enumerate(times):
        r_now = r[t_idx]
        t_now = t[t_idx]
        vr_now = vr[t_idx]
        vt_now = vt[t_idx]
        wt_now = wt[t_idx]
        #streamlines' linear density 
        lambda_ = get_lambda(total_ring_mass, number_of_streamlines, J2, Rp, r_now, t_now, vr_now, vt_now)
        #preserve only two middlemost streamlines, to force acceleations() to compute one-sided accelerations
        middle_index = number_of_streamlines/2
        rc = r_now[middle_index-1:middle_index+1,:]
        tc = t_now[middle_index-1:middle_index+1,:]
        vrc = vr_now[middle_index-1:middle_index+1,:]
        vtc = vt_now[middle_index-1:middle_index+1,:]
        wtc = wt_now[middle_index-1:middle_index+1,:]
        #accelerations
        Ar, At = accelerations(lambda_, G_ring, shear_viscosity, bulk_viscosity, c, rc, tc, vrc, vtc, fast_gravity, confine_inner_edge, confine_outer_edge)
        #angular_momentum_flux
        angular_momentum_flux = lambda_*rc*At
        angular_momentum_flux = angular_momentum_flux[-1]
        angular_momentum_flux_list += [angular_momentum_flux]
        #torque=angular momentum luminosity
        m1 = total_ring_mass/(number_of_streamlines*particles_per_streamline)
        torque = m1*rc*At
        torque = torque[-1]
        torque = torque.sum()
        angular_momentum_luminosity_list += [torque]
    #convert lists to arrays
    angular_momentum_flux = np.array(angular_momentum_flux_list)
    angular_momentum_luminosity = np.array(angular_momentum_luminosity_list)
    #preserve central streamline's history
    rc = r[:, middle_index-1:middle_index+1,:]
    rc = rc[:, -1, :]
    tc = t[:, middle_index-1:middle_index+1,:]
    tc = tc[:, -1, :]
    wtc = wt[:, middle_index-1:middle_index+1,:]
    wtc = wtc[:, -1, :]
    return angular_momentum_flux, angular_momentum_luminosity, rc, tc, wtc

#compute energy flux and luminosity
def calculate_energy_flux(total_ring_mass, number_of_streamlines, particles_per_streamline, J2, Rp, G_ring, 
        shear_viscosity, bulk_viscosity, c, r, t, vr, vt, wt, times, fast_gravity, confine_inner_edge, confine_outer_edge):
    energy_flux_list = []
    energy_luminosity_list = []
    for t_idx, _ in enumerate(times):
        r_now = r[t_idx]
        t_now = t[t_idx]
        vr_now = vr[t_idx]
        vt_now = vt[t_idx]
        wt_now = wt[t_idx]
        #streamlines' linear density 
        lambda_ = get_lambda(total_ring_mass, number_of_streamlines, J2, Rp, r_now, t_now, vr_now, vt_now)
        #preserve only two middlemost streamlines, to force acceleations() to compute one-sided accelerations
        middle_index = number_of_streamlines/2
        rc = r_now[middle_index-1:middle_index+1,:]
        tc = t_now[middle_index-1:middle_index+1,:]
        vrc = vr_now[middle_index-1:middle_index+1,:]
        vtc = vt_now[middle_index-1:middle_index+1,:]
        wtc = wt_now[middle_index-1:middle_index+1,:]
        #accelerations
        Ar, At = accelerations(lambda_, G_ring, shear_viscosity, bulk_viscosity, c, rc, tc, vrc, vtc, fast_gravity, confine_inner_edge, confine_outer_edge)
        #angular_momentum_flux
        energy_flux = lambda_*(Ar*vrc + At*vtc)
        energy_flux = energy_flux[-1]
        energy_flux_list += [energy_flux]
        #work=energy luminosity
        m1 = total_ring_mass/(number_of_streamlines*particles_per_streamline)
        work = m1*(Ar*vrc + At*vtc)
        work = work[-1]
        work = work.sum()
        energy_luminosity_list += [work]
    #convert lists to arrays
    energy_flux = np.array(energy_flux_list)
    energy_luminosity = np.array(energy_luminosity_list)
    return energy_flux, energy_luminosity
