#!/usr/bin/env bash
#
# Note: Run this on glacier!

declare -A SAMPLES
SAMPLES[Mu_ubdt]="DLLmu > 2.0 & DLLe < 1.0 & IsMuon == 1.0 & UBDT > 0.25"
SAMPLES[Mu_ubdt_veto]="DLLmu > 2.0 & DLLe < 1.0 & IsMuon == 1.0 & UBDT < 0.25"

declare -A POLARITY
POLARITY[up]="mu"
POLARITY[down]="md"

PARTICLE=Mu_nopt
BASE_FOLDER=pidcalib_ubdt
rm -rf ${BASE_FOLDER}

GLOBAL_CUTS="Brunel_IPCHI2 > 45 & Brunel_TRACK_GHOSTPROB < 0.5"

for year in 16; do
    for polarity in "up" "down"; do
        for part in "${!SAMPLES[@]}"; do
            folder_name="${BASE_FOLDER}/run2-rdx-20${year}-${POLARITY[${polarity}]}-${part}-p_eta_ntracks"
            echo "Output folder: ${folder_name}"
            pidcalib2.make_eff_hists \
                --output-dir ${folder_name} \
                --sample "Turbo${year}" --magnet ${polarity} \
                --particle ${PARTICLE} \
                --cut "${GLOBAL_CUTS}" \
                --pid-cut "${SAMPLES[${part}]}" \
                --bin-var Brunel_P --bin-var Brunel_ETA --bin-var nTracks_Brunel \
                --binning-file ./binning.json
        done
    done
done

# now rename the pkls
PKL_FOLDER=pkl-run2-rdx_mu_ubdt_old
rm -rf ${PKL_FOLDER}
mkdir -p ${PKL_FOLDER}

for pkl in ./${BASE_FOLDER}/*/*.pkl; do
    new_name="$(basename $(dirname ${pkl})).pkl"
    echo "Renaming $pkl to ${PKL_FOLDER}/${new_name}..."
    cp ${pkl} ${PKL_FOLDER}/${new_name}
done
