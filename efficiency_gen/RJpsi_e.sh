#!/usr/bin/env bash

cuts=( "MC15TuneV1_ProbNNghost>0.2" "(MC15TuneV1_ProbNNghost<0.2) & (MC15TuneV1_ProbNNe>0.1)" "MC15TuneV1_ProbNNghost<0.2 & MC15TuneV1_ProbNNe<0.1 & DLLmu>2.0" "MC15TuneV1_ProbNNghost<0.2 & MC15TuneV1_ProbNNe<0.1 & DLLmu<2.0 & MC15TuneV1_ProbNNk>0.1 & (DLLK-DLLp)>2" "MC15TuneV1_ProbNNghost<0.2 & MC15TuneV1_ProbNNe<0.1 & DLLmu<2.0 & MC15TuneV1_ProbNNk<0.1 & MC15TuneV1_ProbNNp>0.1 & (DLLp-DLLK)>2")

for year in 16; do
    for polarity in "up" "down"; do
	for part in "e_B_Jpsi"; do
            for cut in "${cuts[@]}"; do
		pidcalib2.make_eff_hists \
                    --output-dir pidcalib_output_pi_230512 \
                    --sample "Electron${year}" --magnet ${polarity} \
                    --particle ${part} \
		    --pid-cut "${cut}" \
                    --cut "InMuonAcc==1 & MuonUnbiased==1" \
                    --bin-var Brunel_P --bin-var Brunel_ETA --bin-var nTracks \
                    --binning-file "binning/customBinning${part}.json"
	    done
        done
    done
done
