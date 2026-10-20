#!/bin/zsh
# The 333-streamline zero-twist run (explore/twist0_n333_anim/Qoff, Nt = 40, crossed at step 1825, ~145 orbits)
# with other particle counts Nt, pressure off, to 1e3 orbits (12 600 steps, one output per ~5 orbits) or a
# crossing. Does Nt change the crash? (A first pass to ~258 orbits: Nt = 20, 80, 160, 320 did not cross.)
# Run from explore/twist0_n333_Nt.
for Nt in 20 30 50 60 80 160 320; do
    d=Nt$Nt
    rm -rf $d && mkdir -p $d && cp ../../tests/nominal/inputs.py $d/ && \
    (cd $d && uv run python -u ../../../fast/epi_int_fast.py -m '{"number_of_streamlines": 333, "particles_per_streamline": '$Nt', "timesteps_per_output": 63, "total_number_of_outputs": 200}' > output.log 2>&1)
    echo "$d: $(grep -E 'cross|terminated|execution' $d/output.log | sort -u | tr '\n' ' ')"
done
