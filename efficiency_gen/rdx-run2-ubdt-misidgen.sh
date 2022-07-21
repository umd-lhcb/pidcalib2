#!/usr/bin/env bash
#
# Note: Run this on glacier!

declare -A SAMPLES

for year in 16; do
    for polarity in "up" "down"; do
	for var in "Brunel_P" "Brunel_PT"; do
            for part in "Mu" "K" "Pi" "P"; do
		pidcalib2.make_eff_hists \
		    --output-dir pidcalib_output_many_precut \
		    --sample "Turbo${year}" --magnet ${polarity} \
		    --particle ${part} \
		    $(for I in $(seq 0 0.05 1); do echo "--pid-cut UBDT>${I} "; done) \
		    $(for I in $(seq 0 0.05 1); do echo "--pid-cut MC15TuneV1_ProbNNmu>${I} "; done) \
		    --cut "IsMuon==1 & DLLmu>2" \
		    --bin-var ${var}
	    done
	done
    done
done
