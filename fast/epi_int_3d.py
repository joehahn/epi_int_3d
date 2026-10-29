#!/usr/bin/env python3
"""epi_int_3d's driver: epi_int_lite's streamline integrator with vertical motion.

It runs the same kick-drift-kick loop as the reference driver, source/epi_int_lite.py, from the
same inputs.py, but its own bookkeeping: the ring's state lives in arrays that are updated in place,
the coordinate conversions, drift, sorting and linear density are Numba kernels
(fast/bookkeeping.py), and ring gravity is the Numba kernel of fast/kernels.py.

Vertical motion (Phase 3): each particle also has z and vz. The drift moves them in the planet's
potential (notes/equations.md V1-V4, coded as the rotation V8), so inclined orbits keep their
inclination and their nodes regress at Omega_a - nu_a; z and vz also take part in the barycentric
bookkeeping (coordinate kick, mixed <-> planetocentric velocities). Ring gravity is the 3D wire
(fast/kernels.py ring_gravity_3d, the paper's Eqn. 23): each streamline pulls along the perpendicular
to its near point, in r, theta and z, and A_z kicks vz. A ring at z = vz = 0 stays there and gets the
planar kernel's in-plane result bit for bit. With epi_int_lite_gravity = True (below) there are no
vertical ring forces, ring gravity acts as if every particle were at z = 0, and the in-plane results
are the reference's bit for bit, with or without inclination.

Optional inputs (inputs.py or -m): inclination i (default 0) and node (default 0), given to every
particle at the start, with z = a i sin(lambda - node), vz = a nu_a i cos(lambda - node) (V2).
Ring gravity has three fixes of epi_int_lite's (notes/epi_int_lite_issues.md issues 7 and 9), each
an optional input, on by default:
- fix_lambda: each perturbing streamline pulls with its own linear density, interpolated to the
  perturbed particle's longitude, rather than the perturbed particle's;
- nearpoint: the pull is 2 G lambda/Delta with Delta the perpendicular distance to the wire, rather
  than the radial gap at the particle's longitude (fast/kernels.py, notes/equations.md G14);
- shape_slope: the wire's direction comes from the perturbing streamline's shape rather than its
  velocity, which m >= 2 patterns need.
Optional input vertical_gravity (default: nearpoint and shape_slope): the 3D wire, which is the
nearpoint and shape_slope wire out of the plane and so requires both; False keeps the planar kernel.
Optional input epi_int_lite_gravity (default False): True turns all four off, unless set
explicitly, giving epi_int_lite's ring gravity, as the regression tests against the reference use.
Outputs z.npy and vz.npy (planetocentric, like the others) join the reference's output files.
Pressure and viscosity, when on, are still the reference's NumPy code (helper_fns.accelerations),
and so is anything that inputs.py patches through tools/patches.py: initialization (the edge
taper), the accelerations (fix_viscosity, fix_edge) and the crossing check, which runs as a kernel
when it is the patched (wrapped) one and as the reference's code otherwise.

Run it like the reference driver, from a folder holding inputs.py:
    cd tests/nominal && uv run python ../../fast/epi_int_3d.py [-m '{"total_number_of_outputs": 20}']
The gravity kernel uses all cores; set NUMBA_NUM_THREADS to use fewer.
"""
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO/'source'))
sys.path.insert(0, str(REPO/'fast'))
import helper_fns
import kernels
import bookkeeping


def nu_a(J2, Rp, a):
    """Vertical frequency at semimajor axis a (notes/equations.md V1)."""
    return np.sqrt((1.0/a**3)*(1.0 + 4.5*J2*(Rp/a)**2))


def vertical_elements(J2, Rp, r, t, z, vr, vt, vz):
    """Planetocentric coordinates to elements a, e, wt, M, inclination i and node (V3)."""
    a, e, wt, M = helper_fns.coords2elem(J2, Rp, r, t, vr, vt)
    lam = wt + (helper_fns.Omega(J2, Rp, a)/helper_fns.Kappa(J2, Rp, a))*M
    nu = nu_a(J2, Rp, a)
    inc = np.sqrt((z/a)**2 + (vz/(a*nu))**2)
    node = helper_fns.adjust_angle(lam - np.arctan2(z, vz/nu))
    return a, e, wt, M, inc, node


def restore_output(output_folder):
    """helper_fns.restore_output plus z and vz (zeros for output written without them)."""
    r, t, vr, vt, times, monitor = helper_fns.restore_output(output_folder)
    folder = Path(output_folder)
    z = np.load(folder/'z.npy') if (folder/'z.npy').exists() else np.zeros_like(r)
    vz = np.load(folder/'vz.npy') if (folder/'vz.npy').exists() else np.zeros_like(r)
    return r, t, z, vr, vt, vz, times, monitor


gravity_kernel = kernels.ring_gravity_parallel    #main() sets it from fix_lambda, nearpoint, shape_slope
gravity_kernel_3d = None                          #and this one when vertical_gravity is on
vertical = {'zw': None, 'Az': None}               #Ring.velocity_kick hands over the wrapped z, gets A_z back


def ring_gravity(lambda_, G_ring, r, t, vr, vt, fast_gravity):
    """Drop-in for helper_fns.ring_gravity; epi_int_3d always has fast_gravity = False. With
    vertical_gravity the 3D kernel also returns A_z, which is left in vertical['Az']."""
    if fast_gravity:
        raise ValueError('fast_gravity is not supported: set fast_gravity = False')
    if gravity_kernel_3d is None or vertical['zw'] is None:    #None during initialization, at z = 0
        return gravity_kernel(lambda_, G_ring, r, t, vr, vt)
    Ar, At, vertical['Az'] = gravity_kernel_3d(lambda_, G_ring, r, t, vertical['zw'])
    return Ar, At


PARALLEL_MIN_PARTICLES = 50000     #below this the bookkeeping kernels run serially: with Numba's workqueue
                                   #threading layer a parallel loop costs ~60 us to launch, more than these
                                   #kernels take at 20 000 particles (rt2xy 0.15 ms, xy2rt 0.38, drift 0.34)


class Kernels:
    """bookkeeping's kernels, as attributes, serial or parallel."""
    def __init__(self, table):
        self.__dict__.update(table)


class Ring:
    """The ring's state r, t, vr, vt (mixed-center coordinates, as in the reference driver), its
    work arrays, and the reference's steps, each with the reference's arithmetic."""

    def __init__(self, r, t, vr, vt, p):
        self.bk = Kernels(bookkeeping.PARALLEL if r.size >= PARALLEL_MIN_PARTICLES else bookkeeping.SERIAL)
        self.p = p
        self.shape = r.shape
        self.N = r.size
        self.r, self.t, self.vr, self.vt = (np.array(x, dtype=float, order='C') for x in (r, t, vr, vt))
        self.z, self.vz = np.zeros(self.shape), np.zeros(self.shape)
        new = lambda: np.empty(self.shape)
        self.X, self.Y, self.VX, self.VY, self.R, self.T, self.S, self.C, self.lam = (new() for _ in range(9))
        self.a, self.e, self.wt, self.M, self.NU = (new() for _ in range(5))
        self.rw, self.tw, self.zw = (np.empty((self.shape[0], self.shape[1] + 2)) for _ in range(3))
        f = lambda x: x.reshape(-1)          #flat views for the elementwise kernels
        self.f = {k: f(getattr(self, k)) for k in ('r', 't', 'vr', 'vt', 'z', 'vz', 'X', 'Y', 'VX', 'VY', 'R',
                                                   'T', 'S', 'C', 'lam', 'a', 'e', 'wt', 'M', 'NU')}
        m = p['total_ring_mass']
        J2, Rp = p['J2'], p['Rp']
        self.consts = (1.5*J2, Rp*Rp, Rp, 2.0*1.0*Rp)
        self.J2_45 = 4.5*J2
        self.f_bary = m/(1.0 + m)/self.N                      #planeto2bary's factor
        self.f_planeto = m/1.0/self.N                         #bary2planeto's factor
        self.f_lambda = m/(p['number_of_streamlines']*2*np.pi)

    def rt2xy(self, r, t, vr, vt):
        F = self.f
        self.bk.rt2xy(r, t, vr, vt, F['X'], F['Y'], F['VX'], F['VY'], F['S'], F['C'])

    def coordinate_kick(self, dt):
        """helper_fns.coordinate_kick."""
        F = self.f
        self.rt2xy(F['r'], F['t'], F['vr'], F['vt'])
        factor = dt*self.p['total_ring_mass']/1.0
        self.bk.shift(F['X'], F['Y'], factor*self.VX.mean(), factor*self.VY.mean())
        self.bk.xy2rt(F['X'], F['Y'], F['VX'], F['VY'], F['r'], F['t'], F['vr'], F['vt'], F['S'], F['C'])
        self.bk.shift1(F['z'], factor*self.vz.mean())

    def drift(self, dt):
        """coords2elem, drift and elem2coords."""
        F = self.f
        J2_15, Rp2, Rp, two_GM_Rp = self.consts
        self.bk.drift(F['r'], F['t'], F['vr'], F['vt'], dt, J2_15, Rp2, Rp, two_GM_Rp,
                 F['a'], F['e'], F['wt'], F['M'], F['S'], F['C'])
        self.bk.vertical_drift(F['a'], F['z'], F['vz'], dt, self.J2_45, Rp2, F['NU'], F['S'], F['C'])

    def mixed2planeto(self):
        """helper_fns.mixed2planeto: barycentric velocities to planetocentric."""
        F = self.f
        self.rt2xy(F['r'], F['t'], F['vr'], F['vt'])
        self.bk.shift(F['X'], F['Y'], -self.f_bary*self.X.sum(), -self.f_bary*self.Y.sum())
        self.bk.polar(F['X'], F['Y'], F['R'], F['T'])                         #r_bc, t_bc
        self.rt2xy(F['R'], F['T'], F['vr'], F['vt'])
        x0, y0 = -self.f_planeto*self.X.sum(), -self.f_planeto*self.Y.sum()
        vx0, vy0 = -self.f_planeto*self.VX.sum(), -self.f_planeto*self.VY.sum()
        self.bk.shift4(F['X'], F['Y'], F['VX'], F['VY'], -x0, -y0, -vx0, -vy0)
        self.bk.xy2rt(F['X'], F['Y'], F['VX'], F['VY'], F['R'], F['T'], F['vr'], F['vt'], F['S'], F['C'])
        self.bk.shift1(F['vz'], -(-self.f_planeto*self.vz.sum()))

    def planeto2mixed(self):
        """helper_fns.planeto2mixed: planetocentric velocities to barycentric."""
        F = self.f
        self.rt2xy(F['r'], F['t'], F['vr'], F['vt'])
        f = -self.f_bary
        self.bk.shift4(F['X'], F['Y'], F['VX'], F['VY'], f*self.X.sum(), f*self.Y.sum(), f*self.VX.sum(), f*self.VY.sum())
        self.bk.xy2rt(F['X'], F['Y'], F['VX'], F['VY'], F['R'], F['T'], F['vr'], F['vt'], F['S'], F['C'])
        self.bk.shift1(F['vz'], f*self.vz.sum())

    def velocity_kick(self, dt):
        """helper_fns.velocity_kick."""
        p, F = self.p, self.f
        self.mixed2planeto()
        self.bk.sort_streamlines(self.r, self.t, self.vr, self.vt, self.z, self.vz)
        J2_15, Rp2, Rp, two_GM_Rp = self.consts
        self.bk.linear_density(F['r'], F['vr'], F['vt'], F['lam'], self.f_lambda, J2_15, Rp2, Rp, two_GM_Rp)
        if gravity_kernel_3d is not None:
            self.bk.wrap(self.z, self.zw, 0.0)
            vertical['zw'], vertical['Az'] = self.zw, None
        Ar, At = helper_fns.accelerations(self.lam, p['G_ring'], p['shear_viscosity'], p['bulk_viscosity'],
            p['c'], self.r, self.t, self.vr, self.vt, p['fast_gravity'], p['confine_inner_edge'],
            p['confine_outer_edge'])
        self.vr += Ar*dt
        self.vt += At*dt
        if vertical['Az'] is not None:
            self.vz += vertical['Az'][:, 1:-1]*dt
        self.planeto2mixed()

    def incline(self, inclination, node):
        """Initial z, vz (V2) for inclination i and node, with lambda and a from the in-plane elements."""
        p = self.p
        r, t, vr, vt = helper_fns.mixed2planeto(p['total_ring_mass'], self.r.copy(), self.t.copy(), self.vr.copy(),
                                                self.vt.copy())
        a, e, wt, M = helper_fns.coords2elem(p['J2'], p['Rp'], r, t, vr, vt)
        lam = wt + (helper_fns.Omega(p['J2'], p['Rp'], a)/helper_fns.Kappa(p['J2'], p['Rp'], a))*M
        self.z[:] = a*inclination*np.sin(lam - node)
        self.vz[:] = a*nu_a(p['J2'], p['Rp'], a)*inclination*np.cos(lam - node)
        self.vz += -self.f_bary*self.vz.sum()                 #planetocentric to barycentric, as planeto2mixed

    def snapshot(self):
        """Planetocentric z, vz, in the order that helper_fns.store_system gives r, t, vr, vt."""
        idx = self.t.argsort(axis=1)
        vz = self.vz - (-self.f_planeto*self.vz.sum())
        return np.take_along_axis(self.z, idx, axis=1), np.take_along_axis(vz, idx, axis=1)

    def monitor(self, monitor, timestep):
        """tools/patches.py's wrapped crossing check as a kernel, or else helper_fns' own check."""
        if helper_fns.monitor_streamlines.__module__ != 'patches':
            return helper_fns.monitor_streamlines(monitor, self.r, self.t, timestep)
        monitor['current_timestep'] = timestep
        if np.isnan(self.r).any() and monitor['nan_timestep'] is None:
            print('nan coordinate at timestep = ', timestep)
            monitor['nan_timestep'] = timestep
        self.bk.wrap(self.r, self.rw, 0.0)
        self.bk.wrap(self.t, self.tw, 2.0*np.pi)
        if monitor['streamline_crossing_timestep'] is None and self.bk.crossed(self.rw, self.tw):
            print('steamlines cross at timestep = ', timestep)
            monitor['streamline_crossing_timestep'] = timestep
        return monitor


def main():
    #read input parameters, as source/epi_int_lite.py does
    g = {'np': np, 'random_seed': None}
    with open('inputs.py') as fp:
        exec(fp.read(), g)
    for k in ('number_of_streamlines', 'particles_per_streamline', 'dt', 'timesteps_per_output',
              'total_number_of_outputs', 'radial_width', 'total_ring_mass'):
        print(k, '=', g[k])
    print('ring gravitation constant =', g['G_ring'])
    for k in ('fast_gravity', 'shear_viscosity', 'bulk_viscosity', 'confine_inner_edge', 'confine_outer_edge'):
        print(k, '=', g[k])
    print("Toomre's Q_ring =", g['Q_ring'])
    g.setdefault('inclination', 0.0)
    g.setdefault('node', 0.0)
    g.setdefault('epi_int_lite_gravity', False)
    for k in ('fix_lambda', 'nearpoint', 'shape_slope'):
        g.setdefault(k, not g['epi_int_lite_gravity'])
    g.setdefault('vertical_gravity', g['nearpoint'] and g['shape_slope'])
    if g['vertical_gravity'] and not (g['nearpoint'] and g['shape_slope']):
        raise ValueError('vertical_gravity requires nearpoint = True and shape_slope = True')
    for k in ('Rp', 'J2', 'initial_orbits', 'inclination', 'node', 'epi_int_lite_gravity', 'fix_lambda',
              'nearpoint', 'shape_slope', 'vertical_gravity', 'output_folder'):
        print(k, '=', g[k])
    global gravity_kernel, gravity_kernel_3d
    gravity_kernel = kernels.get_kernel(g['fix_lambda'], g['nearpoint'], g['shape_slope'])
    gravity_kernel_3d = kernels.get_kernel_3d(g['fix_lambda']) if g['vertical_gravity'] else None
    helper_fns.ring_gravity = ring_gravity
    h = helper_fns
    r, t, vr, vt, c, monitor = h.initialize_streamline(g['number_of_streamlines'], g['particles_per_streamline'],
        g['radial_width'], g['total_ring_mass'], g['G_ring'], g['fast_gravity'], g['shear_viscosity'],
        g['bulk_viscosity'], g['confine_inner_edge'], g['confine_outer_edge'], g['Q_ring'], g['J2'], g['Rp'],
        g['initial_orbits'], random_seed=g['random_seed'])
    g['c'] = c
    ring = Ring(r, t, vr, vt, g)
    if g['inclination'] != 0.0:
        ring.incline(g['inclination'], g['node'])
    dt = g['dt']

    timestep = 0
    number_of_outputs = 0
    rz, tz, vrz, vtz, timestepz = [ring.r.copy()], [ring.t.copy()], [ring.vr.copy()], [ring.vt.copy()], [timestep]
    zz, vzz = [ring.z.copy()], [ring.vz.copy()]               #as the reference, the first output is unconverted
    print('evolving system...')
    while number_of_outputs < g['total_number_of_outputs']:
        ring.velocity_kick(dt/2.0)
        for _ in range(g['timesteps_per_output']):
            ring.coordinate_kick(dt/2)
            ring.drift(dt)
            ring.coordinate_kick(dt/2)
            ring.velocity_kick(dt)
            timestep += 1
            monitor = ring.monitor(monitor, timestep)
        ring.velocity_kick(-dt/2.0)
        rz, tz, vrz, vtz, timestepz = h.store_system(rz, tz, vrz, vtz, timestepz, ring.r, ring.t, ring.vr,
                                                     ring.vt, g['total_ring_mass'], timestep)
        z, vz = ring.snapshot()
        zz.append(z)
        vzz.append(vz)
        number_of_outputs += 1
        if 20*number_of_outputs % g['total_number_of_outputs'] == 0:
            if not h.update_display(number_of_outputs, g['total_number_of_outputs'], dt, timestep, monitor):
                break

    timez = np.array(timestepz)*dt
    h.save_output(rz, tz, vrz, vtz, timez, g['modified_params'], monitor, g['output_folder'])
    np.save(g['output_folder'] + '/z.npy', np.array(zz))
    np.save(g['output_folder'] + '/vz.npy', np.array(vzz))
    print('execution time (minutes) = ', (monitor['current_time'] - monitor['start_time'])/60.0)


if __name__ == '__main__':
    main()
