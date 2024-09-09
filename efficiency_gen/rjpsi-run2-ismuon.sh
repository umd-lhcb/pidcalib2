#!/usr/bin/env bash
#
# Note: Run this on glacier!

declare -A SAMPLES
# SAMPLES[Pi]="IsMuon==1.0"
SAMPLES[Pi]="IsMuon == 1.0 & probe_NShared==0.0"

declare -A POLARITY
POLARITY[up]="mu"
POLARITY[down]="md"

PARTICLE=Pi
BASE_FOLDER=pidcalib_rjpsi
rm -rf ${BASE_FOLDER}

GLOBAL_CUTS="Brunel_InMuonAcc == 1.0 & probe_Brunel_MuonUnbiased"

for year in 16; do
    for polarity in "up" "down"; do
        for part in "${!SAMPLES[@]}"; do
            folder_name="${BASE_FOLDER}/run2-rjpsi-20${year}-${polarity}-${part}-Brunel_P"
            echo "Output folder: ${folder_name}"
            pidcalib2.make_eff_hists \
                --output-dir ${folder_name} \
                --sample "Turbo${year}" --magnet ${polarity} \
                --particle ${PARTICLE} \
                --cut "${GLOBAL_CUTS}" \
                --pid-cut "${SAMPLES[${part}]}" \
                --bin-var Brunel_P \
                --binning-file ./binning/JpsiKBinningPi.json
        done
    done
done

# now rename the pkls
PKL_FOLDER=pkl-run2-rjpsi_mu_ubdt_old
rm -rf ${PKL_FOLDER}
mkdir -p ${PKL_FOLDER}

for pkl in ./${BASE_FOLDER}/*/*.pkl; do
    new_name="$(basename $(dirname ${pkl})).pkl"
    echo "Renaming $pkl to ${PKL_FOLDER}/${new_name}..."
    cp ${pkl} ${PKL_FOLDER}/${new_name}
done
