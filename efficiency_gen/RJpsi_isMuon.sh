#!/usr/bin/env bash
#
# Note: Run this on glacier!

declare -A SAMPLES
# SAMPLES[Pi]="IsMuon==1.0 & UBDT>0.65 & DLLmu>2.0 & DLLe<-1.0 & probe_Brunel_ANNTraining_MuonNShared==1"

for year in 16; do
    for polarity in "up" "down"; do
        for part in "Pi"; do
	    pidcalib2.make_eff_hists \
                --output-dir pidcalib_output_isMuon \
                --sample "Turbo${year}" \
		--magnet ${polarity} \
                --particle ${part} \
		--pid-cut "IsMuon==1.0 & probe_Brunel_ANNTraining_MuonNShared==0.0" \
                --cut "Brunel_InMuonAcc==1" \
                --bin-var Brunel_P \
                --binning-file binning/JpsiKBinningPi.json
        done
    done
done
