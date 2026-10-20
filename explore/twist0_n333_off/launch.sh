#!/bin/zsh
# The 333-streamline zero-twist run (explore/twist0_n333, Nt = 40, crossed at step 1825, ~145 orbits) with
# ring gravity off (G0) or viscosity off (nu0), pressure off, to 1e3 orbits: which force drives the crossing?
# Run from explore/twist0_n333_off.
run() { rm -rf $1 && mkdir -p $1 && cp ../../tests/nominal/inputs.py $1/ && \
    (cd $1 && uv run python -u ../../../fast/epi_int_fast.py -m '{"number_of_streamlines": 333, "timesteps_per_output": 63, "total_number_of_outputs": 200, '$2'}' > output.log 2>&1)
    echo "$1: $(grep -E 'cross|terminated|execution' $1/output.log | sort -u | tr '\n' ' ')"; }
run G0 '"G_ring": -1.0'
run nu0 '"shear_viscosity": 0.0, "bulk_viscosity": 0.0'
run shear_only '"bulk_viscosity": 0.0'
run bulk_only '"shear_viscosity": 0.0'
# With the viscosity gating fixed (tools/patches.py fix_viscosity, issue 6): a check that both-on is unchanged
# at Nr = 100 (300 steps), then shear viscosity with its full stress, and bulk viscosity actually applied
runf() { rm -rf $1 && mkdir -p $1 && cp ${3:-inputs_fixvisc.py} $1/inputs.py && \
    (cd $1 && uv run python -u ../../../fast/epi_int_fast.py -m "{$2}" > output.log 2>&1)
    echo "$1: $(grep -E 'cross|terminated|execution|gating' $1/output.log | sort -u | tr '\n' ' ')"; }
runf check_fix '"timesteps_per_output": 30, "total_number_of_outputs": 10'
runf check_ref '"timesteps_per_output": 30, "total_number_of_outputs": 10' ../../tests/nominal/inputs.py
runf shear_only_fixed '"number_of_streamlines": 333, "timesteps_per_output": 63, "total_number_of_outputs": 200, "bulk_viscosity": 0.0'
runf bulk_only_fixed '"number_of_streamlines": 333, "timesteps_per_output": 63, "total_number_of_outputs": 200, "shear_viscosity": 0.0'
# bulk viscosity doubled, nu_b = 2 nu_s: linear theory predicts stable at R = 0.31 (needs nu_b/nu_s > 1.6)
run bulk_x2 '"bulk_viscosity": 1.4e-10'
# circular streamlines (e = 0), both viscosities: the axisymmetric analysis (notes/numerics/viscous_stability.py) predicts no crossing
run circular '"initial_orbits": {"shape": "circular"}'
