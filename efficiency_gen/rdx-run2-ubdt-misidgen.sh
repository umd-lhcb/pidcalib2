#!/usr/bin/env bash
#
# Note: Run this on glacier!
# Generates files necessary for efficiency vs. bkg efficiency plots

rm -rf pidcalib_ubdt_eff

for year in 16; do
    for polarity in "up"; do
        for var in "Brunel_P" "Brunel_PT"; do
            for part in "Mu" "K" "Pi" "P" "e_B_Jpsi"; do
                pidcalib2.make_eff_hists \
                    --output-dir pidcalib_ubdt_eff \
                    --sample "Turbo${year}" --magnet ${polarity} \
                    --particle ${part} \
                    $(for I in $(seq 0 0.05 0.3); do echo "--pid-cut UBDT>${I} "; done) \
                    $(for I in $(seq 0 0.05 0.3); do echo "--pid-cut Brunel_MC15TuneV1_ProbNNmu>${I} "; done) \
                    --cut "IsMuon==1 & MuonUnbiased==1 & DLLmu>2" \
                    --bin-var ${var}
            done
        done
    done
done
