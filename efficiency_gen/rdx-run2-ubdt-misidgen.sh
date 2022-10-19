#!/usr/bin/env bash
#
# Note: Run this on glacier!
# Generates files necessary for efficiency vs. bkg efficiency plots

OUTPUT_FOLDER=pidcalib_ubdt_eff
rm -rf ${OUTPUT_FOLDER}

declare -A PREFIX
PREFIX[K]="Turbo"
PREFIX[Pi]="Turbo"
PREFIX[P]="Turbo"
PREFIX[Mu_nopt]="Turbo"
PREFIX[e_B_Jpsi]="Electron"

for year in 16; do
    for polarity in "up"; do
        for var in "Brunel_P" "Brunel_PT"; do
            for part in "Mu_nopt" "K" "Pi" "P"; do
                pidcalib2.make_eff_hists \
                    --output-dir ${OUTPUT_FOLDER} \
                    --sample "${PREFIX[${part}]}${year}" --magnet ${polarity} \
                    --particle ${part} \
                    $(for I in $(seq 0 0.05 0.3); do echo "--pid-cut UBDT>${I} "; done) \
                    $(for I in $(seq 0 0.05 0.3); do echo "--pid-cut Brunel_MC15TuneV1_ProbNNmu>${I} "; done) \
                    --cut "IsMuon==1 & MuonUnbiased==1 & DLLmu>2" \
                    --bin-var ${var}
            done
        done
    done
done
