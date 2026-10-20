#!/bin/zsh
# A replacement nominal ringlet for 1000 streamlines (2026-10-19): the narrow ringlet of explore/narrow
# (width 0.00125, e = 1e-5, edge taper 0.2, J2 = 0.01) with zero initial twist and viscosity and pressure
# off (as the first bending-wave runs will be). Screen to 1e3 orbits at Nt = 20 and 40, mass 3e-11
# (lambda_crit/2Delta ~ 0.09 at 1000 streamlines) and 1e-10 (~0.31). Run from explore/nominal1000.
B='"number_of_streamlines": 1000, "radial_width": 0.00125, "shear_viscosity": 0.0, "bulk_viscosity": 0.0, "Q_ring": -1.0, "timesteps_per_output": 126, "total_number_of_outputs": 100'
E='"initial_orbits": {"shape": "eccentric", "e": 1e-5, "e_prime": 0.0, "w_prime": 0.0}'
run() { rm -rf $1 && mkdir -p $1 && cp ../../tests/nominal/inputs.py $1/ && \
    (cd $1 && uv run python -u ../../../fast/epi_int_fast.py -m "{$B, $E, $2}" > output.log 2>&1)
    echo "$1: $(grep -E 'cross|terminated|execution' $1/output.log | sort -u | tr '\n' ' ')"; }
run nt20_m3e-11 '"particles_per_streamline": 20, "total_ring_mass": 3e-11'
run nt40_m3e-11 '"particles_per_streamline": 40, "total_ring_mass": 3e-11'
run nt20_m1e-10 '"particles_per_streamline": 20, "total_ring_mass": 1e-10'
run nt40_m1e-10 '"particles_per_streamline": 40, "total_ring_mass": 1e-10'
uv run python ../animate.py nt20_m3e-11   # the chosen candidate
# The chosen candidate (Nt = 20, mass 3e-11) to 5e4 orbits (629 000 steps, one output per ~100 orbits), where
# J2 alone would twist it by ~pi
L='"particles_per_streamline": 20, "total_ring_mass": 3e-11, "timesteps_per_output": 1258, "total_number_of_outputs": 500'
run long_nt20_m3e-11 "$L"
