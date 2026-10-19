#!/bin/zsh
# The zero-twist nominal ringlet at 333 streamlines, ~1000 orbits (crossed at step 1825). Run from explore/twist0_n333.
uv run python -u ../../fast/epi_int_fast.py -m '{"number_of_streamlines": 333, "timesteps_per_output": 126, "total_number_of_outputs": 100}' > output.log 2>&1
