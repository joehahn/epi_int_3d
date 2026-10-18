#!/bin/zsh
# Ring pressure on: the tapered pi-twist nominal ringlet (explore/taper/nr1000_pi, nr0100_pi)
# with Toomre Q_ring = 0.5, 1, 1.2, 2. Run from explore/pressure.
B='"radial_width": 0.0125, "total_ring_mass": 3e-8, "shear_viscosity": 7e-11, "bulk_viscosity": 7e-11'
PI='"initial_orbits": {"shape": "eccentric", "e": 0.0001, "e_prime": 0.0, "w_prime": 0.025290}'    # pi twist
SHORT='"total_number_of_outputs": 40, "timesteps_per_output": 10'
run() { mkdir -p $1 && cp inputs.py $1/ && (cd $1 && uv run python ../../run.py --taper 0.2 -m "{$B, $PI, $2}" > output.log 2>&1); }
for Q in 0.5 1.0 1.2 2.0; do
    run nr1000_Q$Q '"number_of_streamlines": 1000, "Q_ring": '$Q", $SHORT" &
    run nr0100_Q$Q '"number_of_streamlines": 100, "Q_ring": '$Q &
done
wait
uv run python ../diagnose.py nr*
