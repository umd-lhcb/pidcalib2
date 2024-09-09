#!/usr/bin/env bash
#
# Note: Run this on glacier!

declare -A SAMPLES
SAMPLES[Pi]="IsMuon==1.0 & UBDT>0.65 & DLLmu>2.0 & DLLe<-1.0 & probe_Brunel_ANNTraining_MuonNShared==1"

for year in 16; do
    for polarity in "up" "down"; do
        for part in "${!SAMPLES[@]}"; do
            pidcalib2.make_eff_hists \
                --output-dir pidcalib_output_uBDT_DLL \
                --sample "Turbo${year}" --magnet ${polarity} \
                --particle ${part} --pid-cut "${SAMPLES[${part}]}" \
                --cut "InMuonAcc==1 & MuonUnbiased==1" \
                --bin-var Brunel_P --bin-var Brunel_ETA --bin-var nTracks_Brunel \
                --binning-file binning/customBinningPi.json
        done
    done
done
