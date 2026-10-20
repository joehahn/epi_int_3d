#!/bin/zsh
# The 333-streamline zero-twist runs (explore/twist0_n333, twist0_n333_Q) rerun with one output per
# ~1 orbit (13 steps) for animations: pressure off and Q_ring = 0.5, 1, 2, 10. Run from explore/twist0_n333_anim.
for Q in -1.0 0.5 1.0 2.0 10.0; do
    d=Q$Q; [[ $Q == -1.0 ]] && d=Qoff
    mkdir -p $d && cp ../../tests/nominal/inputs.py $d/ && \
    (cd $d && uv run python -u ../../../fast/epi_int_fast.py -m '{"number_of_streamlines": 333, "timesteps_per_output": 13, "total_number_of_outputs": 250, "Q_ring": '$Q'}' > output.log 2>&1)
done
uv run python ../animate.py Q*
mkdir -p nr0100 && cp ../../tests/nominal/inputs.py nr0100/ && (cd nr0100 && uv run python -u ../../../fast/epi_int_fast.py -m '{"timesteps_per_output": 13, "total_number_of_outputs": 250}' > output.log 2>&1) && uv run python ../animate.py nr0100   # Nr = 100 for comparison
