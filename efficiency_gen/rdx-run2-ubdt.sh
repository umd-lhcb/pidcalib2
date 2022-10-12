#!/usr/bin/env bash
#
# Note: Run this on glacier!

declare -A SAMPLES
SAMPLES[Mu_nopt]="DLLmu > 2.0 & DLLe < 1.0 & IsMuon == 1.0 & UBDT > 0.25"

for year in 16; do
    for polarity in "up" "down"; do
        for part in "${!SAMPLES[@]}"; do
            pidcalib2.make_eff_hists \
                --output-dir pidcalib_output \
                --sample "Turbo${year}" --magnet ${polarity} \
                --particle ${part} --pid-cut "${SAMPLES[${part}]}" \
                --bin-var Brunel_P --bin-var Brunel_ETA --bin-var nTracks_Brunel \
                --binning-file ./binning.json
        done
    done
done
