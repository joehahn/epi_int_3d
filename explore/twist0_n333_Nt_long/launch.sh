#!/bin/zsh
# The 333-streamline zero-twist runs of explore/twist0_n333_Nt that did not cross by 1e3 orbits, extended to
# 1e4 orbits (126 000 steps, one output per ~50 orbits), pressure off: do they stay saturated or cross later?
# Run from explore/twist0_n333_Nt_long.
for Nt in 30 50; do
    d=Nt$Nt
    rm -rf $d && mkdir -p $d && cp ../../tests/nominal/inputs.py $d/ && \
    (cd $d && uv run python -u ../../../fast/epi_int_fast.py -m '{"number_of_streamlines": 333, "particles_per_streamline": '$Nt', "timesteps_per_output": 630, "total_number_of_outputs": 200}' > output.log 2>&1)
    echo "$d: $(grep -E 'cross|terminated|execution' $d/output.log | sort -u | tr '\n' ' ')"
done
