# Outputs of epi_int_lite tests/libration/libration_frequency.ipynb (master 5e0942a), extracted verbatim.
# Figures from cell N are saved as cellNN.png in this folder.

## cell 2
```python
#restore output
parent_folders = [
    'narrow_ring/ring_mass_1.5e-11/', 'narrow_ring/ring_mass_2.4e-11/', 'narrow_ring/ring_mass_3.7e-11/', 
    'narrow_ring/ring_mass_6.5e-11/', 'narrow_ring/ring_mass_1.3e-10/', 'narrow_ring/ring_mass_2.5e-10/', 
    'narrow_ring/ring_mass_5.0e-10/', 'narrow_ring/ring_mass_1.0e-9/', 'narrow_ring/ring_mass_2.0e-9/', 
    'wide_ring/ring_mass_1.0e-10/', 'wide_ring/ring_mass_2.0e-10/', 'wide_ring/ring_mass_3.3e-10/', 
    'wide_ring/ring_mass_5.0e-10/', 'wide_ring/ring_mass_7.0e-10/', 'wide_ring/ring_mass_1.0e-9/', 
    'wide_ring/ring_mass_1.5e-9/', 'wide_ring/ring_mass_2.0e-9/', 'wide_ring/ring_mass_3.3e-9/', 
    'wide_ring/ring_mass_5.0e-9/', 'wide_ring/ring_mass_8.5e-9/']
data = []
from helper_fns import *
from libration import *
sns.set(font_scale=1.5, font='DejaVu Sans')
sns.set(font_scale=1.1, font='DejaVu Sans')
f_narrow, ax_narrow = plt.subplots(figsize=(8, 8))
f_wide, ax_wide = plt.subplots(figsize=(8, 8))
for idx, parent_folder in enumerate(parent_folders):
    #read simulation output
    input_file = parent_folder + "inputs.py"
    execfile(input_file)
    d = {'parent_folder':parent_folder, 'radial_width':radial_width, 'total_ring_mass':total_ring_mass, 
        'Rp':Rp, 'J2':J2} 
    the_output_folder = parent_folder + output_folder
    print 'reading output in folder = ', the_output_folder
    r, t, vr, vt, times, lambda0 = restore_output(the_output_folder)
    a, e, wt, M = coords2elem(J2, Rp, r, t, vr, vt)
    a_inner, a_outer, a_mean, da, e_inner, e_outer, e_mean, de, wt_inner, wt_outer, dwt = \
        orbit_deltas(times, r, t, vr, vt, a, e, wt, J2, Rp)
    H, q, e_prime, wt_prime = calculate_Hq(a_mean, e_mean, da, de, dwt)
    d['ring_width'] = 'narrow'
    if ('wide' in parent_folder):
        d['ring_width'] = 'wide'
    d['a'] = a_mean.mean()
    d['e'] = e_mean.mean()
    d['da'] = da.mean()
    d['de_eq'] = de.mean()
    d['dwt_eq'] = dwt.mean()
    d['q'] = q.mean()
    d['H'] = H.mean()
    d['e_prime'] = e_prime.mean()
    d['wt_prime'] = wt_prime.mean()
    #BGT libration frequency
    GM = 1.0
    M_central = 1.0
    a3 = d['a']**3
    a_over_da = d['a']/d['da']
    n = np.sqrt(GM/a3)
    omega_bgt = (d['total_ring_mass']/M_central)*(a_over_da**2)*n*d['H']/np.pi
    d['omega_bgt']  = omega_bgt    
    d['T_lib_bgt'] = 2.0*np.pi/d['omega_bgt']
    #compute libration phase & frequency
    Ix = (de - d['de_eq'] )*d['a']/d['da']
    Iy = d['e']*d['a']*(dwt - d['dwt_eq'])/d['da']
    I = np.sqrt(Ix**2 + Iy**2)
    d['I'] = I.mean()
    d['Delta_I_over_I'] = (I.max() - I.min())/d['I']
    d['Delta_q_over_q'] = (q.max() - q.min())/d['q']
    d['Delta_H_over_H'] = (H.max() - H.min())/d['H']
    d['Iy_over_Ix'] = Iy.max()/Ix.max()
    phase = np.arctan2(Iy, Ix)
    phase_unwrapped = unwrap_angle(phase)
    coeff = np.polyfit(times, phase_unwrapped, 1)
    phase_fit = coeff[1] + coeff[0]*times
    omega = coeff[0]
    d['omega'] = omega
    d['T_lib'] = 2.0*np.pi/d['omega']
    data += [d]
    #plot
    x = Ix
    y = Iy
    label = d['ring_width'] + r'    $m_r$ = ' + str(d['total_ring_mass'])
    ax = ax_narrow
    if (d['ring_width'] == 'wide'):
        ax = ax_wide
    linewidth = 1
    alpha = 0.7
    if (d['Delta_H_over_H'] > 0.01):
        linewidth = 6
        alpha = 1
    p = ax.plot(x, y, linestyle='-', linewidth=linewidth, label=label, alpha=alpha)
    p = ax.set_title('ring_width = ' + d['ring_width'])
    p = ax.set_xlabel(r'$I_x = (\Delta e - \Delta e_{eq})a/\Delta a$')
    p = ax.set_ylabel(r'$I_y = ea\Delta\tilde{\omega}/\Delta a$')
    xylim = (-0.02, 0.02)
    p = ax.set_xlim(xylim)
    p = ax.set_ylim(xylim)
    p = ax.legend(loc=(1.1, 0))
df = pd.DataFrame(data)
```
```
reading output in folder =  narrow_ring/ring_mass_1.5e-11/output
reading output in folder =  narrow_ring/ring_mass_2.4e-11/output
reading output in folder =  narrow_ring/ring_mass_3.7e-11/output
reading output in folder =  narrow_ring/ring_mass_6.5e-11/output
reading output in folder =  narrow_ring/ring_mass_1.3e-10/output
reading output in folder =  narrow_ring/ring_mass_2.5e-10/output
reading output in folder =  narrow_ring/ring_mass_5.0e-10/output
reading output in folder =  narrow_ring/ring_mass_1.0e-9/output
reading output in folder =  narrow_ring/ring_mass_2.0e-9/output
reading output in folder =  wide_ring/ring_mass_1.0e-10/output
reading output in folder =  wide_ring/ring_mass_2.0e-10/output
reading output in folder =  wide_ring/ring_mass_3.3e-10/output
reading output in folder =  wide_ring/ring_mass_5.0e-10/output
reading output in folder =  wide_ring/ring_mass_7.0e-10/output
reading output in folder =  wide_ring/ring_mass_1.0e-9/output
reading output in folder =  wide_ring/ring_mass_1.5e-9/output
reading output in folder =  wide_ring/ring_mass_2.0e-9/output
reading output in folder =  wide_ring/ring_mass_3.3e-9/output
reading output in folder =  wide_ring/ring_mass_5.0e-9/output
reading output in folder =  wide_ring/ring_mass_8.5e-9/output
[figure: cell02.png]
[figure: cell02.png]
```

## cell 3
```python
df.tail()
```
```
    Delta_H_over_H  Delta_I_over_I  Delta_q_over_q         H         I  \
15        0.007492        0.172379        0.098978  0.518846  0.011060   
16        0.001468        0.384567        0.038808  0.509446  0.003077   
17        0.003736        0.148574        0.390285  0.502480  0.015690   
18        0.000404        0.255192        0.194651  0.500523  0.003627   
19        0.000034        0.371260        1.551052  0.500009  0.003372   

    Iy_over_Ix    J2   Rp         T_lib     T_lib_bgt    ...         e  \
15    1.010921  0.01  0.5  24207.059076  25387.294418    ...     0.005   
16    1.150970  0.01  0.5  18876.885983  19394.608985    ...     0.005   
17    1.028900  0.01  0.5  11797.083701  11929.966539    ...     0.005   
18    1.133126  0.01  0.5   7849.984665   7911.284574    ...     0.005   
19    1.103572  0.01  0.5   4630.268711   4668.869618    ...     0.005   

     e_prime     omega  omega_bgt                parent_folder         q  \
15  0.219315  0.000260   0.000247  wide_ring/ring_mass_1.5e-9/  0.219461   
16  0.157039  0.000333   0.000324  wide_ring/ring_mass_2.0e-9/  0.157055   
17  0.079553  0.000533   0.000527  wide_ring/ring_mass_3.3e-9/  0.080338   
18  0.037145  0.000800   0.000794  wide_ring/ring_mass_5.0e-9/  0.037233   
19  0.003429  0.001357   0.001346  wide_ring/ring_mass_8.5e-9/  0.004360   

    radial_width  ring_width total_ring_mass  wt_prime  
15         0.001        wide    1.500000e-09  0.000378  
16         0.001        wide    2.000000e-09  0.000124  
17         0.001        wide    3.300000e-09  0.000563  
18         0.001        wide    5.000000e-09  0.000022  
19         0.001        wide    8.500000e-09  0.000180  

[5 rows x 24 columns]
```

## cell 4
```python
#plot Delta_H_over_H vs total_ring_masses
sns.set(font_scale=1.5, font='DejaVu Sans')
f, ax = plt.subplots(figsize=(10, 6))
ring_widths = df.ring_width.unique()
for idx, ring_width in enumerate(ring_widths):
    jdx = (df.ring_width == ring_width)
    dfs = df[jdx]
    x = dfs.total_ring_mass
    y = dfs.Delta_H_over_H
    color = color_seq[idx]
    if (idx ==1):
        color = color_seq[idx + 1]
    p = ax.plot(x, y, marker='o', markersize=8, linestyle='-', linewidth=2, label=ring_width, color=color)
p = ax.set_xlabel(r'ringlet mass $m_r$')
p = ax.set_ylabel(r'$\Delta H/H$')
p = ax.set_xlim(1.0e-11, 1.0e-8)
p = ax.set_ylim(1.0e-5, 2)
p = ax.set_xscale('log')
p = ax.set_yscale('log')
p = plt.legend()
```
```
[figure: cell04.png]
```

## cell 5
```python
#plot Delta_I_over_I vs total_ring_masses
sns.set(font_scale=1.5, font='DejaVu Sans')
f, ax = plt.subplots(figsize=(10, 6))
ring_widths = df.ring_width.unique()
for idx, ring_width in enumerate(ring_widths):
    jdx = (df.ring_width == ring_width)
    dfs = df[jdx]
    x = dfs.total_ring_mass
    y = dfs.Delta_I_over_I
    color = color_seq[idx]
    if (idx ==1):
        color = color_seq[idx + 1]
    p = ax.plot(x, y, marker='o', markersize=8, linestyle='-', linewidth=2, label=ring_width, color=color)
p = ax.set_xlabel(r'ringlet mass $m_r$')
p = ax.set_ylabel(r'$\Delta I/I$')
p = ax.set_xlim(1.0e-11, 1.0e-8)
p = ax.set_ylim(-0.1, 1.2)
p = ax.set_xscale('log')
p = plt.legend()
```
```
[figure: cell05.png]
```

## cell 6
```python
#confirm lower-mass ringlet have high Delta_I_over_I
idx = df.ring_width == 'narrow'
dfs = df[idx]
cols = ['total_ring_mass', 'Delta_I_over_I', 'Delta_H_over_H', 'Delta_q_over_q', 'I', 'q', 'H', 'e_prime']
dfs.sort_values('Delta_I_over_I', ascending=False)[cols]
```
```
   total_ring_mass  Delta_I_over_I  Delta_H_over_H  Delta_q_over_q         I  \
0     1.500000e-11        0.980383        0.142761        0.019474  0.018148   
1     2.400000e-11        0.756769        0.094842        0.026952  0.018276   
2     3.700000e-11        0.628604        0.035678        0.019979  0.010065   
3     6.500000e-11        0.433985        0.013835        0.019719  0.006496   
4     1.300000e-10        0.306081        0.003896        0.020192  0.003586   
7     1.000000e-09        0.302031        0.000247        0.196464  0.002850   
5     2.500000e-10        0.251453        0.002882        0.058572  0.005288   
8     2.000000e-09        0.229301        0.000039        1.362296  0.004010   
6     5.000000e-10        0.217160        0.001910        0.198671  0.007970   

          q         H   e_prime  
0  0.922651  1.877579  0.922476  
1  0.849876  1.244576  0.849707  
2  0.748900  0.907790  0.748849  
3  0.571044  0.669004  0.571021  
4  0.339078  0.547714  0.339067  
7  0.028973  0.500317  0.028903  
5  0.178548  0.512287  0.178507  
8  0.004730  0.500010  0.003280  
6  0.079835  0.502415  0.079634
```

## cell 7
```python
#confirm lower-mass ringlet have high Delta_I_over_I
idx = df.ring_width == 'wide'
dfs = df[idx]
cols = ['total_ring_mass', 'Delta_I_over_I', 'Delta_H_over_H', 'Delta_q_over_q', 'I', 'q', 'H', 'e_prime']
dfs.sort_values('Delta_I_over_I', ascending=False)[cols]
```
```
    total_ring_mass  Delta_I_over_I  Delta_H_over_H  Delta_q_over_q         I  \
9      1.000000e-10        1.081308        0.202324        0.019515  0.020365   
10     2.000000e-10        0.702399        0.064805        0.019847  0.013086   
11     3.300000e-10        0.591216        0.030522        0.021008  0.009686   
12     5.000000e-10        0.434086        0.014336        0.020154  0.006680   
16     2.000000e-09        0.384567        0.001468        0.038808  0.003077   
19     8.500000e-09        0.371260        0.000034        1.551052  0.003372   
13     7.000000e-10        0.279583        0.007040        0.018493  0.004577   
14     1.000000e-09        0.276725        0.003424        0.018466  0.003222   
18     5.000000e-09        0.255192        0.000404        0.194651  0.003627   
15     1.500000e-09        0.172379        0.007492        0.098978  0.011060   
17     3.300000e-09        0.148574        0.003736        0.390285  0.015690   

           q         H   e_prime  
9   0.943540  2.282090  0.943316  
10  0.840987  1.200162  0.840900  
11  0.711838  0.836522  0.711789  
12  0.573650  0.671161  0.573625  
16  0.157055  0.509446  0.157039  
19  0.004360  0.500009  0.003429  
13  0.453388  0.593216  0.453374  
14  0.333162  0.545892  0.333153  
18  0.037233  0.500523  0.037145  
15  0.219461  0.518846  0.219315  
17  0.080338  0.502480  0.079553
```

## cell 8
```python
#plot T_lib vs total_ring_masses
sns.set(font_scale=1.5, font='DejaVu Sans')
f, ax = plt.subplots(figsize=(10, 6))
ring_widths = df.ring_width.unique()
for idx, ring_width in enumerate(ring_widths):
    jdx = (df.ring_width == ring_width)
    dfs = df[jdx]
    x = dfs.total_ring_mass
    y = dfs.T_lib/(2*np.pi)
    color = color_seq[idx]
    if (idx ==1):
        color = color_seq[idx + 1]
    p = ax.plot(x, y, marker='o', markersize=8, linestyle='-', linewidth=2, label=ring_width, color=color)
    y = dfs.T_lib_bgt/(2*np.pi)
    p = ax.plot(x, y, linestyle='dashed', linewidth=2, label='', color=color)
p = ax.set_xlabel(r'ringlet mass $m_r$')
p = ax.set_ylabel(r'libration period $T_{lib}$    (orbits)')
p = ax.set_xlim(8.0e-12, 2.0e-8)
p = ax.set_ylim(5.0e2, 5.0e4)
p = ax.set_xscale('log')
p = ax.set_yscale('log')
p = plt.legend()
plt.savefig('figs/libration.pdf', bbox_inches='tight')
```
```
[figure: cell08.png]
```

