#!/bin/zsh
# The nominal ringlet (tests/nominal) made 30% wider at the same surface density: width 0.01625,
# mass 3.9e-8, pi twist (w_prime = pi e a_mean/radial_width), nu = 7e-11, 100 streamlines, taper 0.2.
# Does it break into more than two co-precessing sub-ringlets? Run from explore/wide.
B='"radial_width": 0.01625, "total_ring_mass": 3.9e-8, "shear_viscosity": 7e-11, "bulk_viscosity": 7e-11'
PI='"initial_orbits": {"shape": "eccentric", "e": 0.0001, "e_prime": 0.0, "w_prime": 0.019490}'
run() { mkdir -p $1 && cp inputs.py $1/ && (cd $1 && uv run python ../../run.py --taper 0.2 -m "{$B, $2}" > output.log 2>&1); }
run nr0100_w1.3 '"number_of_streamlines": 100, '"$PI"
uv run python ../diagnose.py nr*

# Animation (figs/periapse.html): from nr0100_w1.3, make inputs_run.py (tests/nominal/inputs.py with
# radial_width 0.01625, total_ring_mass 3.9e-8, w_prime 0.019490, and the patches lines removed) and
# anim.ipynb (tests/nominal/nominal.ipynb reading inputs_run.py, make_animation = True), then
#   uv run jupyter nbconvert --to notebook --execute --inplace anim.ipynb
