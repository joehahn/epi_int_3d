#!/bin/zsh
# Narrow nominal ringlet (keeps J2 = 0.01): width 0.00125 (10x narrower), e = 1e-5 (q = 0.025 as
# before), pi twist, nu = 7e-13 (x2 spreading), mass bracketing 3e-11 (the m ~ J2 da^3/dw scaling
# from 3e-8). Edge taper 0.2. Run from explore/narrow.
B='"radial_width": 0.00125, "shear_viscosity": 7e-13, "bulk_viscosity": 7e-13'
PI='"initial_orbits": {"shape": "eccentric", "e": 1e-5, "e_prime": 0.0, "w_prime": 0.025148}'
run() { mkdir -p $1 && cp inputs.py $1/ && (cd $1 && uv run python ../../run.py --taper 0.2 -m "{$B, $PI, $2}" > output.log 2>&1); }
for m in 1e-11 3e-11 1e-10; do
    run nr0100_m$m '"number_of_streamlines": 100, "total_ring_mass": '$m &
done
run nr1000_m3e-11 '"number_of_streamlines": 1000, "total_ring_mass": 3e-11, "total_number_of_outputs": 40, "timesteps_per_output": 50' &
wait
uv run python ../diagnose.py nr*
