#!/bin/zsh
# The zero-twist nominal ringlet at 333 streamlines (explore/twist0_n333, crossed at step 1825) with
# ring pressure on, Q_ring = 0.5, 1, 2, 10. Does pressure prevent the crossing? Run from explore/twist0_n333_Q.
for Q in 0.5 1.0 2.0 10.0; do
    mkdir -p Q$Q && cp ../../tests/nominal/inputs.py Q$Q/ && \
    (cd Q$Q && uv run python -u ../../../fast/epi_int_fast.py -m '{"number_of_streamlines": 333, "timesteps_per_output": 126, "total_number_of_outputs": 100, "Q_ring": '$Q'}' > output.log 2>&1)
done
