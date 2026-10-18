#!/bin/zsh
# Tapered-edge ladder (edge taper fraction 0.2; width 0.0125 keeps the middle spacing and
# surface density of the untapered W = 0.01 runs). Run from explore/taper.
B='"radial_width": 0.0125, "total_ring_mass": 3e-8, "shear_viscosity": 7e-11, "bulk_viscosity": 7e-11'
PI='"initial_orbits": {"shape": "eccentric", "e": 0.0001, "e_prime": 0.0, "w_prime": 0.025290}'    # pi twist: w_prime = pi e a_mean/radial_width
TEN='"initial_orbits": {"shape": "eccentric", "e": 0.001, "e_prime": 0.0, "w_prime": 0.014053}'    # 10 deg twist
SHORT='"total_number_of_outputs": 40, "timesteps_per_output": 10'
run() { mkdir -p $1 && cp inputs.py $1/ && (cd $1 && uv run python ../../run.py --taper 0.2 -m "{$B, $2}" > output.log 2>&1); }
run nr1000_pi       '"number_of_streamlines": 1000, '"$PI, $SHORT" &
run nr1000_10deg    '"number_of_streamlines": 1000, '"$TEN, $SHORT" &
run nr0100_pi       '"number_of_streamlines": 100, '"$PI" &
run nr0100_10deg    '"number_of_streamlines": 100, '"$TEN" &
wait
uv run python ../diagnose.py nr*

# Issue-3 diagnostics of nr1000_pi: 241 particles per streamline, and ring gravity off
# run nr1000_pi_nt241 '"number_of_streamlines": 1000, "particles_per_streamline": 241, '"$PI, $SHORT"
# run nr1000_pi_G0    '"number_of_streamlines": 1000, "G_ring": -1.0, '"$PI, $SHORT"

# Issue-4 check: nr0100_pi with the outer-edge flux corrected (compare with nr0100_pi)
# (cd nr0100_pi_fixedge && uv run python ../../run.py --taper 0.2 --fix-edge -m "{$B, \"number_of_streamlines\": 100, $PI}" > output.log 2>&1)
# (cd nr0100_pi_fixedge_own && uv run python ../../run.py --taper 0.2 --fix-edge-own -m "{$B, \"number_of_streamlines\": 100, $PI}" > output.log 2>&1)
# run nr0300_pi '"number_of_streamlines": 300, '"$PI"
