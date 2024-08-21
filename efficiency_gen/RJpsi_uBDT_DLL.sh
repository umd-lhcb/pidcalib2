#!/usr/bin/env bash
#
# Note: Run this on glacier!

declare -A SAMPLES
# SAMPLES[Pi]="IsMuon==1.0 & UBDT>0.65 & DLLmu>2.0 & DLLe<-1.0 & probe_Brunel_ANNTraining_MuonNShared==1"

for year in 16 17 18; do
    for polarity in "up" "down"; do
	for var in "Brunel_P" "Brunel_PT"; do
            for part in "P" "K" "Mu" "Pi"; do
		pidcalib2.make_eff_hists \
                    --output-dir pidcalib_output_uBDT_DLL \
                    --sample "Turbo${year}" --magnet ${polarity} \
                    --particle ${part} \
		    --pid-cut "UBDT>0.65 & DLLmu>2.0 & DLLe<-1.0" \
                    --cut "probe_Brunel_ANNTraining_MuonNShared==1" \
                    --bin-var ${var} \
                    --binning-file binning/customBinningPi.json
	    done
        done
    done
done
