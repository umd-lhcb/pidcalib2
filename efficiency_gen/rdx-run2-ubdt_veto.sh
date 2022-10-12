#!/usr/bin/env bash
#
# Note: Run this on glacier!

declare -A SAMPLES
SAMPLES[Mu_nopt]="DLLmu > 2.0 & DLLe < 1.0 & IsMuon == 1.0 & UBDT < 0.25"

BASE_FOLDER=pidcalib_ubdt
rm -rf ${BASE_FOLDER}

for year in 16; do
    for polarity in "up" "down"; do
        for part in "${!SAMPLES[@]}"; do
            folder_name="{BASE_FOLDER}/run2-rdx-20${year}-${POLARITY[${polarity}]}-${part}_ubdt_veto-p_eta_ntracks"
            echo "Output folder: ${folder_name}"
            pidcalib2.make_eff_hists \
                --output-dir pidcalib_output \
                --sample "Turbo${year}" --magnet ${polarity} \
                --particle ${part} --pid-cut "${SAMPLES[${part}]}" \
                --bin-var Brunel_P --bin-var Brunel_ETA --bin-var nTracks_Brunel \
                --binning-file ./binning.json
        done
    done
done


# now rename the pkls
PKL_FOLDER=pkl-run2-rdx_mu_ubdt_veto
rm -rf ${PKL_FOLDER}
mkdir -p ${PKL_FOLDER}

for pkl in ./${BASE_FOLDER}/*/*.pkl; do
    new_name="$(basename $(dirname ${pkl})).pkl"
    echo "Renaming $pkl to ${PKL_FOLDER}/${new_name}..."
    cp ${pkl} ${PKL_FOLDER}/${new_name}
done
