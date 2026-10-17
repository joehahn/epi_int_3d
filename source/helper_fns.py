#helper_fns.py
#
#by Joe Hahn, jmh.datasciences@gmail.com, 28 September 2017.
#
#helper functions used by epi_int_lite.py

#angular frequency
def Omega(J2, Rp, a, Ar=0.0):
    GM = 1.0
    a2 = a*a
    Ra2 = (Rp*Rp)/a2
    Omega2 = (GM/a2/a)*(   1.0 + (1.5*J2)*Ra2 - Ar*(a2/GM)   )
    return np.sqrt(Omega2)

#epicyclic frequency
def Kappa(J2, Rp, a, Ar=0.0, kappa_squared=False):
    GM = 1.0
    a2 = a*a
    Ra2 = (Rp*Rp)/a2
    Kappa2 = (GM/a2/a)*(   1.0 - (1.5*J2)*Ra2 - Ar*(a2*(3.0/GM))   )
    if (kappa_squared):
        return Kappa2
    else:
        return np.sqrt(Kappa2)

#adjust angles to live between -Pi and Pi
import numpy as np
def adjust_angle(angle):
    idx = angle > np.pi
    angle[idx] -= 2.0*np.pi
    idx = angle < -np.pi
    angle[idx] += 2.0*np.pi    
    return angle

#unwrap angle that lives between -Pi and Pi
def unwrap_angle(angle):
    angle_unw = angle.copy()
    twopi = 2.0*np.pi
    for idx in range(1, len(angle)):
        delta = angle_unw[idx] - angle_unw[idx - 1]
        if (delta < -np.pi):
            angle_unw[idx:] += twopi 
        if (delta > np.pi):
            angle_unw[idx:] -= twopi 
    return angle_unw

#convert cylindrical coordinates and velocities to cartesian
def rt2xy(r, t, vr, vt):
    sint = np.sin(t)
    cost = np.cos(t)
    x = r*cost
    y = r*sint
    vx = vr*cost - vt*sint
    vy = vr*sint + vt*cost
    return x, y, vx, vy

#convert cartesian coordinates and velocities to cylindrical
def xy2rt(x, y, vx, vy):
    r = np.sqrt(x*x + y*y)
    t = np.arctan2(y, x)
    sint = np.sin(t)
    cost = np.cos(t)
    vr =  vx*cost + vy*sint
    vt = -vx*sint + vy*cost
    return r, t, vr, vt

#order particles in each streamline by their longitudes
def sort_particles(r, t, vr, vt):
    for streamline_idx in range(len(t)):
        longitude_idx = t[streamline_idx].argsort()
        r[streamline_idx] = r[streamline_idx][longitude_idx]
        t[streamline_idx] = t[streamline_idx][longitude_idx]
        vr[streamline_idx] = vr[streamline_idx][longitude_idx]
        vt[streamline_idx] = vt[streamline_idx][longitude_idx]
    return r, t, vr, vt

#drift step advances M
def drift(a, M, J2, Rp, dt):
    return M + Kappa(J2, Rp, a)*dt

#shift 2D array 1 element to right when n=1 and left when -1, is significantly faster than .roll()
def sidestep(x, n):
    Ny, Nx = x.shape
    if (n > 0):
        left = x[:,-1].reshape(Ny, 1)
        right = x[:,0:-1]
    else:
        left = x[:, 1:]
        right = x[:, 0].reshape(Ny, 1)
    return np.concatenate((left, right), axis=1)

#shift 2D array vertically n rows, is significantly faster than .roll()
def advance(x, n):
    lower = x[-n:]
    upper = x[:-n]
    return np.concatenate((lower, upper), axis=0)

#use lagrange polynomial to evaluate function f that is evaluated n adjacent
#streamlines away and sampled at longitude t
def interpolate_fn(t, f, n, interpolate=True):
    if (interpolate):
        t1 = advance(t, -n)
        t0 = sidestep(t1,  1)
        t2 = sidestep(t1, -1)
        f1 = advance(f, -n)
        f0 = sidestep(f1,  1)
        f2 = sidestep(f1, -1)
        f_n = lagrange_poly_fit(t0, t1, t2, f0, f1, f2, t)
    else:
        #skip lagrange interpolation
        f_n = np.roll(f, (-n, 0), axis=(0,1))
    return f_n

#fit 2nd order lagrange polynomial to data (x0,y0),(x1,y1),(x2,y2) & interpolate y(x)
def lagrange_poly_fit(x0, x1, x2, y0, y1, y2, x):
    dx0 = x - x0
    dx1 = x - x1
    dx2 = x - x2
    dx10 = x1 - x0
    dx20 = x2 - x0
    dx21 = x2 - x1
    l0 =  (dx1*dx2)/(dx10*dx20)
    l1 = -(dx0*dx2)/(dx10*dx21)
    l2 =  (dx0*dx1)/(dx20*dx21)
    y = y0*l0 + y1*l1 + y2*l2
    return y

#wrap the ring's coordinate array about in longitude
def wrap_ring(c, longitude=False):
    Nr, Nt = c.shape
    left = c[:, -1].copy().reshape(Nr, 1)
    right = c[:, 0].copy().reshape(Nr, 1)
    if (longitude):
        twopi = 2*np.pi
        left -= twopi
        right += twopi
    cw = np.concatenate((left, c, right), axis=1)
    return cw

#compute ring surface density
def surface_density(lambda_, dr):
    sd = lambda_/dr
    return sd

#calculate (radial distance between exterior streamline and interior streamline)/2
def delta_f(f, t):
    #in ring interior
    f_plus  = interpolate_fn(t, f,  1)
    f_minus = interpolate_fn(t, f, -1)
    if (f.shape[0] > 2):
        df = (f_plus - f_minus)/2
    else:
        df = np.zeros_like(f)
    #along ring edges
    df[0] = f_plus[0] - f[0]
    df[-1] = f[-1] - f_minus[-1]
    return df

#calculation derivative df/dr
def df_dr(delta_f, delta_r):
    return delta_f/delta_r

#acceleration due to ring self-gravity
def ring_gravity(lambda_, G_ring, r, t, vr, vt, fast_gravity):
    if  (fast_gravity == True):
        v = np.sqrt(vr*vr + vt*vt)
        cos_phi = vt/v
        sin_phi = vr/v 
    Ar = 0
    At = 0
    two_G_lambda = 2.0*G_ring*lambda_
    Nr, Nt = r.shape
    for shft in range(1, Nr):
        dr = interpolate_fn(t, r, -shft) - r
        A = two_G_lambda/dr
        if (fast_gravity == False):
            vri = interpolate_fn(t, vr, -shft)
            vti = interpolate_fn(t, vt, -shft)
            vi = np.sqrt(vri*vri + vti*vti)
            cos_phi = vti/vi
            sin_phi = vri/vi
        Ar += A*cos_phi
        At -= A*sin_phi
    return Ar, At

#acceleration due to pressure P
def A_P(lambda_, sd, P, t, delta_P, delta_r):
    dPdr = df_dr(delta_P, delta_r)
    #acceleration in ring interior
    A = -dPdr/sd
    #at inner streamline
    A[0] = -P[0]/lambda_[0]
    #at outer streamline, interpolated from neighbor streamline delta_L = 4.56e-7
    P_outer = interpolate_fn(t, P, 1)
    P_outer = P_outer[-2]
    A[-1] = P_outer/lambda_[-1]
    return A

#acceleration due to ring pressure
def ring_pressure(c, lambda_, sd, r, t, vr, vt, v, delta_r):
    #pressure
    P = (c*c)*sd
    delta_P = delta_f(P, t)
    #acceleration
    A = A_P(lambda_, sd, P, t, delta_P, delta_r)
    #radial and tangential components
    Ar =  A*(vt/v)
    At = -A*(vr/v)
    return Ar, At

#angular acceleration due to ring viscosity
def ring_shear_viscosity(shear_viscosity, lambda_, sd, r, t, vr, vt, v, delta_r):
    w = vt/r
    delta_w = delta_f(w, t)
    dw_dr = df_dr(delta_w, delta_r)
    #viscous pseudo-pressure
    P = -(shear_viscosity*sd)*r*dw_dr
    delta_P = delta_f(P, t)
    #tangential viscous acceleration
    At = A_P(lambda_, sd, P, t, delta_P, delta_r)
    return At

#radial acceleration due to ring viscosity
def ring_bulk_viscosity(shear_viscosity, bulk_viscosity, lambda_, sd, r, t, vr, vt, v, delta_r):
    delta_vr = delta_f(vr, t)
    dvr_dr = df_dr(delta_vr, delta_r)
    #viscosity coefficients
    nu_1 = bulk_viscosity
    nu_2 = bulk_viscosity
    if (shear_viscosity > 0):
        nu_1 += 4.0*shear_viscosity/3.0
        nu_2 -= 2.0*shear_viscosity/3.0
    #viscous pseudo-pressure
    P = -(nu_1*sd)*dvr_dr - (nu_2*sd)*(vr/r)
    delta_P = delta_f(P, t)
    #radial viscous acceleration
    Ar = A_P(lambda_, sd, P, t, delta_P, delta_r)
    return Ar

#additional tangential acceleration due to edge-torques
def edge_torques(r, At, confine_inner_edge, confine_outer_edge):
    idxs = []
    if (confine_inner_edge): idxs += [0]
    if (confine_outer_edge): idxs += [-1]
    for idx in idxs:
        specific_torque = (r[idx]*At[idx]).mean()
        At[idx] -= specific_torque/r[idx]
    return At

#calculate radial and tangential accelerations due to ring gravity, pressure, viscosity
def accelerations(lambda_, G_ring, shear_viscosity, bulk_viscosity, c, r, t, vr, vt, fast_gravity, confine_inner_edge, confine_outer_edge):
    #wrap ring around in longitude
    rw = wrap_ring(r, longitude=False)
    tw = wrap_ring(t, longitude=True)
    vrw = wrap_ring(vr, longitude=False)
    vtw = wrap_ring(vt, longitude=False)
    vw = np.sqrt(vrw*vrw + vtw*vtw)
    lw = wrap_ring(lambda_, longitude=False)
    Ar = 0
    At = 0
    #acceleration due to streamline gravity
    if (G_ring > 0.0):
        A = ring_gravity(lw, G_ring, rw, tw, vrw, vtw, fast_gravity)
        Ar += A[0]
        At += A[1]
    #acceleration due to streamline pressure and viscosity
    if ((c > 0.0) or (shear_viscosity > 0.0)):
        delta_rw = delta_f(rw, tw)
        sdw = surface_density(lw, delta_rw)
        if (c > 0.0):
            A = ring_pressure(c, lw, sdw, rw, tw, vrw, vtw, vw, delta_rw)
            Ar += A[0]
            At += A[1]
        if (shear_viscosity > 0.0):
            At += ring_shear_viscosity(shear_viscosity, lw, sdw, rw, tw, vrw, vtw, vw, delta_rw)
        if (bulk_viscosity > 0.0):
            Ar += ring_bulk_viscosity(shear_viscosity, bulk_viscosity, lw, sdw, rw, tw, vrw, vtw, vw, delta_rw)
    #add confinement torque at ring edges, as needed
    if (confine_inner_edge or confine_outer_edge):
        At = edge_torques(rw, At, confine_inner_edge, confine_outer_edge)
    #drop left and right edges from Ar,At
    if (type(Ar) != int):
        Ar = Ar[:, 1:-1]
    if (type(At) != int):
        At = At[:, 1:-1]
    return Ar, At

#compute streamline's linear density
def get_lambda(total_ring_mass, number_of_streamlines, J2, Rp, r, t, vr, vt):
    a, e, wt, M = coords2elem(J2, Rp, r, t, vr, vt)
    Omg = Omega(J2, Rp, a)
    v = np.sqrt(vr*vr +vt*vt)
    factor = total_ring_mass/(number_of_streamlines*2*np.pi)
    lambda_ = factor*(Omg/v)
    return lambda_

#velocity kicks due to ring gravity, viscosity, pressure
def velocity_kick(J2, Rp, G_ring, shear_viscosity, bulk_viscosity, c, total_ring_mass, number_of_streamlines, \
        r, t, vr, vt, dt, fast_gravity, confine_inner_edge, confine_outer_edge):
    #convert mixed-center coordinates to planetocentric
    r, t, vr, vt = mixed2planeto(total_ring_mass, r, t, vr, vt)
    #order particles by longitude
    r, t, vr, vt = sort_particles(r, t, vr, vt)
    #compute streamlines' linear density 
    lambda_ = get_lambda(total_ring_mass, number_of_streamlines, J2, Rp, r, t, vr, vt)
    #radial acceleration due to ring gravity, viscosity, pressure
    Ar, At = accelerations(lambda_, G_ring, shear_viscosity, bulk_viscosity, c, r, t, vr, vt, fast_gravity, confine_inner_edge, confine_outer_edge)
    #kick velocity
    vr += Ar*dt
    vt += At*dt
    #convert planetocentric coordinates to mixed-center coordinates
    r, t, vr, vt = planeto2mixed(total_ring_mass, r, t, vr, vt)
    return r, t, vr, vt

#kick coordinates to account for central body's motion about center of mass
def coordinate_kick(dt, total_ring_mass, r, t, vr, vt):
    mass_0 = 1.0
    factor = dt*total_ring_mass/mass_0
    x, y, vx, vy = rt2xy(r, t, vr, vt)
    x += factor*vx.mean()
    y += factor*vy.mean()
    r, t, vr, vt = xy2rt(x, y, vx, vy)
    return r, t, vr, vt

#convert orbit elements to coordinates
def elem2coords(J2, Rp, a, e, wt, M, Ar=0.0, sort_particle_longitudes=True):
    e_sin_M = e*np.sin(M)
    e_cos_M = e*np.cos(M)
    r = a/np.sqrt(1.0 + 2.0*e_cos_M)
    Omg = Omega(J2, Rp, a, Ar=Ar)
    Kap = Kappa(J2, Rp, a, Ar=Ar)
    t = adjust_angle(   (Omg/Kap)*(M + 2.0*e_sin_M) + wt   )
    ra = r/a
    ra3 = ra*ra*ra
    vr = (a*Kap*ra3)*e_sin_M
    vt = a*a*Omg/r
    return r, t, vr, vt

#convert coordinates to orbit elements
def coords2elem(J2, Rp, r, t, vr, vt, Ar=0.0):
    GM = 1.0
    h = r*vt
    c = (h*h)/(2.0*GM*Rp)
    a = Rp*(   c + np.sqrt(c*c - 1.5*J2)   )
    Omg = Omega(J2, Rp, a, Ar=Ar)
    Kap = Kappa(J2, Rp, a, Ar=Ar)
    ar = a/r
    ar2 = ar*ar
    ar3 = ar2*ar
    e_cos_M = (ar2 - 1.0)/2.0
    aK = a*Kap
    e_sin_M = vr*(ar3/aK)
    e = np.sqrt(e_sin_M*e_sin_M + e_cos_M*e_cos_M)
    M = np.arctan2(e_sin_M, e_cos_M)
    wt = adjust_angle(   t - (Omg/Kap)*(M + 2.0*e_sin_M)   )
    return a, e, wt, M

#convert planetocentric coordinates and velocities to barycentric
def planeto2bary(total_ring_mass, r, t, vr, vt):
    mass_0 = 1.0
    mass_total = mass_0 + total_ring_mass
    N_particles = r.size
    factor = total_ring_mass/mass_total/N_particles
    x, y, vx, vy = rt2xy(r, t, vr, vt)
    x_0  = -factor*x.sum()
    y_0  = -factor*y.sum()
    vx_0 = -factor*vx.sum()
    vy_0 = -factor*vy.sum()
    x_bc  =  x +  x_0
    y_bc  =  y +  y_0
    vx_bc = vx + vx_0
    vy_bc = vy + vy_0
    r_bc, t_bc, vr_bc, vt_bc = xy2rt(x_bc, y_bc, vx_bc, vy_bc)
    r_0, t_0, vr_0, vt_0 = xy2rt(x_0, y_0, vx_0, vy_0)
    return r_bc, t_bc, vr_bc, vt_bc, r_0, t_0, vr_0, vt_0

#convert barycentric coordinates and velocities to planetocentric
def bary2planeto(total_ring_mass, r_bc, t_bc, vr_bc, vt_bc):
    mass_0 = 1.0
    N_particles = r_bc.size
    factor = total_ring_mass/mass_0/N_particles
    x_bc, y_bc, vx_bc, vy_bc = rt2xy(r_bc, t_bc, vr_bc, vt_bc)
    x_0_bc  = -factor*x_bc.sum()
    y_0_bc  = -factor*y_bc.sum()
    vx_0_bc = -factor*vx_bc.sum()
    vy_0_bc = -factor*vy_bc.sum()
    x = x_bc - x_0_bc
    y = y_bc - y_0_bc
    vx = vx_bc - vx_0_bc
    vy = vy_bc - vy_0_bc
    r, t, vr, vt = xy2rt(x, y, vx, vy)
    return r, t, vr, vt

#convert planetocentric r,v to planetocentic r and barycentric v
def planeto2mixed(total_ring_mass, r, t, vr, vt):
    r_bc, t_bc, vr_bc, vt_bc, r_0, t_0, vr_0, vt_0 = planeto2bary(total_ring_mass, r, t, vr, vt)
    return r, t, vr_bc, vt_bc

#convert planetocentric r,v to planetocentic r and barycentric v
def mixed2planeto(total_ring_mass, r, t, vr_bc, vt_bc):
    r_bc, t_bc, vr_ignore, vt_ignore, r_0, t_0, vr_0, vt_0 = planeto2bary(total_ring_mass, r, t, vr_bc, vt_bc)
    r_ignore, t_ignore, vr, vt = bary2planeto(total_ring_mass, r_bc, t_bc, vr_bc, vt_bc)
    return r, t, vr, vt

#append current r,t,vr,vt,a,timestep to lists rz,tz etc
def store_system(rz, tz, vrz, vtz, timestepz, r, t, vr, vt, total_ring_mass, timestep):
    rc = r.copy()
    tc = t.copy()
    vrc = vr.copy()
    vtc = vt.copy()
    rc, tc, vrc, vtc = mixed2planeto(total_ring_mass, rc, tc, vrc, vtc)
    rc, tc, vrc, vtc = sort_particles(rc, tc, vrc, vtc)
    rz.append(rc)
    tz.append(tc)
    vrz.append(vrc)
    vtz.append(vtc)
    timestepz.append(timestep)
    return rz, tz, vrz, vtz, timestepz

#check for streamine crossing and nans
def monitor_streamlines(monitor, r, t, timestep):
    monitor['current_timestep'] = timestep
    #check for nan in r
    nan_timestep = monitor['nan_timestep']
    if ((np.isnan(r).any() == True) and (nan_timestep == None)):
        print 'nan coordinate at timestep = ', timestep
        monitor['nan_timestep'] = timestep
    #check for streamline crossing where dr[1] = streamline 1's radial distance relative to streamline 0
    dr = r - interpolate_fn(t, r, -1, interpolate=True)
    dr[0] = dr[1]
    idx = (dr < 0)
    streamline_crossing_timestep = monitor['streamline_crossing_timestep']
    if ((idx.sum() > 0) and (streamline_crossing_timestep == None)):
        print 'steamlines cross at timestep = ', timestep
        monitor['streamline_crossing_timestep'] = timestep
    return monitor

#save orbit element arrays in files
import pickle
import os
def save_output(r, t, vr, vt, times, modified_params, monitor, output_folder):
    cmd = 'mkdir -p ' + output_folder
    q = os.system(cmd)
    np.save(output_folder + '/r.npy', r)
    np.save(output_folder + '/t.npy', t)
    np.save(output_folder + '/vr.npy', vr)
    np.save(output_folder + '/vt.npy', vt)
    np.save(output_folder + '/times.npy', times)
    monitor['modified_params'] = modified_params
    with open(output_folder + '/monitor.pkl', 'wb') as fp:
        pickle.dump(monitor, fp, protocol=pickle.HIGHEST_PROTOCOL)

#restore orbit elements from files
def restore_output(output_folder):
    r = np.load(output_folder + '/r.npy')
    t = np.load(output_folder + '/t.npy')
    vr = np.load(output_folder + '/vr.npy')
    vt = np.load(output_folder + '/vt.npy')
    times = np.load(output_folder + '/times.npy')
    with open(output_folder + '/monitor.pkl', 'rb') as fp:
        monitor = pickle.load(fp)
    return r, t, vr, vt, times, monitor

#print eta as needed, and end simulation if monitor says something bad happened to a self_interacting ringlet
import time as tm
def update_display(number_of_outputs, total_number_of_outputs, dt, timestep, monitor):
    monitor['current_time'] = int(tm.time())
    exec_time_min = (monitor['current_time'] - monitor['start_time'])/60.0
    eta_min = int((total_number_of_outputs - number_of_outputs)*exec_time_min/number_of_outputs)
    print 'time = ' + str(timestep*dt) + \
        '    number of outputs = ' + str(number_of_outputs) + \
        '    number of orbits = ' + str(int(timestep*dt/2.0/np.pi)) + \
        '    eta (minutes) = ', eta_min
    continue_sim = True
    if (monitor['self_interacting'] == True):
        for key in ['streamline_crossing_timestep', 'nan_timestep']:
            if (monitor[key]):
                print 'sim terminated at timestep = ' + str(timestep)
                continue_sim = False
    return continue_sim

#initialize streamlines
def initialize_streamline(number_of_streamlines, particles_per_streamline, radial_width,
    total_ring_mass, G_ring, fast_gravity, shear_viscosity, bulk_viscosity,
    confine_inner_edge, confine_outer_edge, Q_ring, J2, Rp, initial_orbits):
    
    #initialize particles in circular orbits
    a_streamlines = np.linspace(1.0, 1.0 + radial_width, num=number_of_streamlines)
    a_list = []
    for sma in a_streamlines:
        a_list.append(np.zeros(particles_per_streamline) + sma)
    a = np.array(a_list)
    e = np.zeros_like(a)
    wt = np.zeros_like(a)
    #particles anomalies are uniformly spaced
    M_streamline = np.linspace(-np.pi, np.pi, num=particles_per_streamline, endpoint=False)
    M_list = [M_streamline]*number_of_streamlines
    M = np.array(M_list)
    
    #modify initial orbits as needed
    if (initial_orbits['shape'] == 'circular'):
        pass
    if (initial_orbits['shape'] == 'eccentric'):
        e_init = initial_orbits['e']
        adeda = initial_orbits['e_prime']
        a_mean = a_streamlines.mean()
        e = e_init + adeda*(a - a_mean)/a_mean
        aedwtda = initial_orbits['w_prime']
        wt += aedwtda*(a - a_mean)/a_mean/e_init
    if (initial_orbits['shape'] == 'breathing mode'):
        e_init = initial_orbits['e']
        e[:] = e_init
        wt = M.copy()
        M[:] = 0.0
    if (initial_orbits['shape'] == 'log-e'):
        #streamlines' e is lograthmically distributed between initial_e[0] & initial_e[1] with random M,wt
        wt = np.random.uniform(low=-np.pi, high=np.pi, size=a.shape)
        M = np.random.uniform(low=-np.pi, high=np.pi, size=a.shape)
        e = np.zeros_like(M)
        e_init = initial_orbits['e']
        for idx in range(number_of_streamlines):
            e[idx] += np.exp( np.random.uniform(low=np.log(e_init[0]), high=np.log(e_init[1])) )
    
    #tweak longitude of periapse away from common value so that eccentric streamlines are closed loops
    Omg = Omega(J2, Rp, a)
    Kap = Kappa(J2, Rp, a)
    wt = wt - (Omg/Kap - 1)*M
    
    #ring coordinates
    r, t, vr, vt = elem2coords(J2, Rp, a, e, wt, M)
    r, t, vr, vt = sort_particles(r, t, vr, vt)
    
    #compute streamlines' linear density 
    lambda_ = get_lambda(total_ring_mass, number_of_streamlines, J2, Rp, r, t, vr, vt)
    
    #ring sound speed c
    c = 0.0
    if (Q_ring > 0.0):
        delta_r = delta_f(r, t)
        sd = surface_density(lambda_, delta_r)
        G = 1.0
        c = (Q_ring*np.pi*G*sd/Omg).mean()
    
    #adjust vt to compensate for ring's radial accelerations
    Ar, At = accelerations(lambda_, G_ring, shear_viscosity, bulk_viscosity, c, r, t, vr, vt, fast_gravity, confine_inner_edge, confine_outer_edge)
    rAr = r*Ar
    for idx in range(number_of_streamlines):
        rAr[idx] = rAr[idx].mean()
    vt = np.sqrt(vt*vt - rAr)
    
    #convert planetocentric coordinates to mixed-center coordinates
    r, t, vr, vt = planeto2mixed(total_ring_mass, r, t, vr, vt)
    
    #ringlet is self_interacting if gravity or viscosity or pressure is turned on
    self_interacting = False
    if ((G_ring > 0) or (shear_viscosity > 0) or (bulk_viscosity > 0) or (Q_ring > 0)):
        self_interacting = True
        
    #this dict is used to track execution time and when streamlines cross or nan is generated
    start_time = int(tm.time())
    monitor = {'start_time':start_time, 'current_time':start_time, 'current_timestep':None, 'streamline_crossing_timestep':None, 
        'nan_timestep':None, 'self_interacting':self_interacting}

    return r, t, vr, vt, c, monitor

#recompute coordinates in coordinate system that co-rotates with ringlet's middle streamline's peri
def peri_corotate(r, t, vr, vt, wt):
    number_of_streamlines = r.shape[0]
    s_idx = (number_of_streamlines - 1)/2
    r_middle_streamline = r[s_idx]
    t_idx = np.argmin(r_middle_streamline)
    wt_middle_streamline = wt[s_idx]
    wt_min = wt_middle_streamline[t_idx]
    tw = adjust_angle(t - wt_min)
    wts = adjust_angle(wt - wt_min)
    rs, ts, vrs, vts = sort_particles(r, tw, vr, vt)
    return rs, ts, vrs, vts, wts
