# Outputs of epi_int_lite tests/spreading/spreading.ipynb (master 5e0942a), extracted verbatim.
# Figures from cell N are saved as cellNN.png in this folder.

## cell 1
```python
#restore output
from helper_fns import *
execfile('inputs.py')
r, t, vr, vt, times, monitor = restore_output(output_folder)
a, e, wt, M = coords2elem(J2, Rp, r, t, vr, vt)
print 'r.shape = ', r.shape
print 'total_ring_mass = ', total_ring_mass
print 'shear_viscosity = ', shear_viscosity
print 'bulk_viscosity = ', bulk_viscosity
print 'initial_orbits = ', initial_orbits
print 'confine_inner_edge = ', confine_inner_edge
print 'confine_outer_edge = ', confine_outer_edge
print 'e.min() = ', e.min()
print 'e.max() = ', e.max()
print 'execution time (min) = ', (monitor['current_time'] - monitor['start_time'])/60.0
print 'monitor = ', monitor
```
```
r.shape =  (1001, 15, 3)
total_ring_mass =  1e-10
shear_viscosity =  1e-13
bulk_viscosity =  1e-13
initial_orbits =  {'shape': 'circular'}
confine_inner_edge =  False
confine_outer_edge =  False
e.min() =  0.0
e.max() =  5.435438028278856e-08
execution time (min) =  17.1166666667
monitor =  {'current_time': 1732624654, 'streamline_crossing_timestep': None, 'start_time': 1732623627, 'modified_params': None, 'nan_timestep': None, 'current_timestep': 834000, 'self_interacting': True}
```

## cell 2
```python
#time in units of viscous timescale
viscous_timescale = (radial_width**2)/(12*shear_viscosity)
print 'viscous_timescale = ', viscous_timescale
times_viscous = times/viscous_timescale
```
```
viscous_timescale =  8333.33333333
```

## cell 5
```python
#plot streamline's mutual a-evolution
twopi = 2*np.pi
sns.set(font_scale=1.5, font='DejaVu Sans')
f, ax = plt.subplots(figsize=(10, 6))
x = times_viscous
y = (a_outer - a_avg[0])/delta_a[0]
p = ax.plot(x, y, linestyle='-', linewidth=2, alpha=1, label='outer')
y = (a_avg - a_avg[0])/delta_a[0]
p = ax.plot(x, y, linestyle='-', linewidth=2, alpha=1, label='mean')
y = (a_inner - a_avg[0])/delta_a[0]
p = ax.plot(x, y, linestyle='-', linewidth=2, alpha=1, label='inner')
p = ax.set_xlabel(r'time $t/\tau_\nu$')
p = ax.set_ylabel('(a - a$_{mean})/\Delta a_0$')
p = ax.legend(loc=(0.79, 0.1))
```
```
[figure: cell05.png]
```

## cell 6
```python
#plot delta_a versus time
delta_a = a_outer - a_inner
print 'initial delta_a = ', delta_a[0]
x = times_viscous
y = delta_a/delta_a[0]
sns.set(font_scale=1.5, font='DejaVu Sans')
f, ax = plt.subplots(figsize=(10, 6))
p = ax.plot(x, y, linestyle='-', linewidth=2)
p = ax.set_ylabel(r'radial width $\Delta a/\Delta a_0$')
p = ax.set_xlabel(r'time $t/\tau_\nu$')
```
```
initial delta_a =  0.00010000000000010001
[figure: cell06.png]
```

## cell 7
```python
#plot streamline's e-evolution
twopi = 2*np.pi
sns.set(font_scale=1.5, font='DejaVu Sans')
f, ax = plt.subplots(figsize=(10, 6))
x = times_viscous
y = e_outer/1.0e-7
p = ax.plot(x, y, linestyle='-', linewidth=1, alpha=0.8, label='outer')
y = e_avg/1.0e-7
p = ax.plot(x, y, linestyle='-', linewidth=1, alpha=0.8, label='mean')
y = e_inner/1.0e-7
p = ax.plot(x, y, linestyle='-', linewidth=1, alpha=0.8, label='inner')
p = ax.set_ylabel(r'eccentricity e   ($10^{-7}$)')
p = ax.set_xlabel(r'time $t/\tau_\nu$')
p = ax.legend()
```
```
[figure: cell07.png]
```

## cell 8
```python
#plot r(theta)
twopi = 2*np.pi
f, ax = plt.subplots(figsize=(10, 6))
a_avg = a.mean(axis=2)                    #average a across longitude
a_0 = a_avg[0].mean()
da_0 = (a_avg[0][-1] - a_avg[0][0])/2
for pidx in range(number_of_streamlines):
    da = (a_avg[:, pidx] - a_0)/da_0
    p = ax.plot(times_viscous, da, linestyle='-', linewidth=2)
#ax.set_xlim(-0.05, 1.05)
yrng = 0.015
#ax.set_ylim(-yrng, yrng)
#ax.set_xlim(-0.2, 10.2)
p = ax.set_xlabel(r'time $t/\tau_\nu$')
p = ax.set_ylabel('a - a$_{mid}$')
```
```
[figure: cell08.png]
```

## cell 9
```python
#compute system's total_angular_momentum in barycentric reference frame
mass_0 = 1.0
particle_mass = total_ring_mass/(number_of_streamlines*particles_per_streamline)
ring_angular_momentum = []
planet_angular_momentum = []
for idx, tm in enumerate(times):
    r_now = r[idx]
    t_now = t[idx]
    vr_now = vr[idx]
    vt_now = vt[idx]
    r_bc, t_bc, vr_bc, vt_bc, r_0, t_0, vr_0, vt_0 = planeto2bary(total_ring_mass, r_now, t_now, vr_now, vt_now)
    l_0 = mass_0*r_0*vt_0
    l_ring = + np.sum(particle_mass*r_bc*vt_bc)
    ring_angular_momentum += [l_ring]
    planet_angular_momentum += [l_0]
ring_angular_momentum = np.array(ring_angular_momentum)
print 'mean ring_angular_momentum = ', ring_angular_momentum.mean()
planet_angular_momentum = np.array(planet_angular_momentum)
print 'mean planet_angular_momentum = ', planet_angular_momentum.mean()
total_angular_momentum = ring_angular_momentum + planet_angular_momentum
delta_L = (total_angular_momentum - total_angular_momentum[0])/total_angular_momentum[0]
print 'max delta_L = ', delta_L.max()
```
```
mean ring_angular_momentum =  1.0018981614430117e-10
mean planet_angular_momentum =  1.8001613665332788e-37
max delta_L =  1.1319086049183388e-07
```

## cell 10
```python
#plot delta_total_angular_momentum over time
#2 streamlines max_delta_L=1.6e-7
#3 streamlines max_delta_L=2.1e-7
#4 streamlines max_delta_L=2.4e-5
#5 streamlines max_delta_L=1.9e-5 particles_per_streamline=3,11,31
#6 streamlines max_delta_L=2.5e-5
#9 streamlines max_delta_L=2.1e-5
#15 streamlines max_delta_L=1.7e-5
factor = 1.0e-7
sns.set(font_scale=1.5, font='DejaVu Sans')
f, ax = plt.subplots(figsize=(10, 6))
x = times_viscous
y = delta_L/factor
p = ax.plot(x, y, linestyle='-', linewidth=2, label='total')
p = ax.set_ylabel(r'$\Delta L/L_0$    $(10^{-7})$')
p = ax.set_xlabel(r'time $t/\tau_\nu$')
```
```
[figure: cell10.png]
```

## cell 11
```python
#plot sigma(r) at selected times
plot_times = [1.0, 2.5, 5.0, 10.0, 20.0, 50.0]
R = r[0].mean()
sd_init = total_ring_mass/np.pi/(R**2)
f, ax = plt.subplots(figsize=(10, 6))
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
sns.set(font_scale=1.5, font='DejaVu Sans')
for idx, plot_time in enumerate(plot_times):
    t_idx = np.where(times_viscous >= plot_time)[0][0]
    r_mean = r[t_idx].mean(axis=1)
    dr_in = r_mean[1] - r_mean[0]
    dr_out = r_mean[-1] - r_mean[-2]
    r_padded = np.array([r_mean[0] - dr_in] + r_mean.tolist() + [r_mean[-1] + dr_out])
    x = r_padded/R
    tau = 12.0*shear_viscosity*times[t_idx]/(R**2)
    arg = ((x - 1)**2)/tau
    sd_exp = np.exp(-arg)/(x**0.75)/np.sqrt(2*twopi*tau)
    lbl = r'$t/\tau_\nu=$ ' + str(times_viscous[t_idx])[0:3]
    if (lbl[-1] == '.'): lbl=lbl.strip('.')
    ax.plot(r_padded, sd_exp, linestyle='-', linewidth=2, alpha=0.7, color=colors[idx], label=lbl)
    r_2D = r[t_idx]
    t_2D = t[t_idx]
    vr_2D = vr[t_idx]
    vt_2D = vt[t_idx]
    lambda_ = get_lambda(total_ring_mass, number_of_streamlines, J2, Rp, r_2D, t_2D, vr_2D, vt_2D)
    delta_r = delta_f(r_2D, t_2D)
    sd_2D = surface_density(lambda_, delta_r)
    r_1D = r_2D.mean(axis=1)
    sd_1D = sd_2D.mean(axis=1)
    dr_in = r_1D[1] - r_1D[0]
    dr_out = r_1D[-1] - r_1D[-2]
    r_padded = np.array([r_1D[0] - dr_in] + r_1D.tolist() + [r_1D[-1] + dr_out])
    sd_padded = np.array([0.0] + sd_1D.tolist() + [0.0])
    ax.plot(r_padded, sd_padded/sd_init, marker='o', markersize=3, linestyle='-.', linewidth=1, alpha=1.0, 
        color=colors[idx])
p = ax.set_xlabel(r'radial distance $r/r_{mid}$')
p = ax.set_ylabel(r'surface density $\sigma/\bar{\sigma}$')
p = ax.xaxis.set_major_locator(plt.MaxNLocator(5))
p = plt.legend(loc='center right')
plt.savefig('figs/sigma_spreading.pdf', bbox_inches='tight')
```
```
[figure: cell11.png]
```

