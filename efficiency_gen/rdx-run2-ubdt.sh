#!/usr/bin/env bash
#
# Note: Run this on glacier!

declare -A BDT_CUT
BDT_CUT[Mu_ubdt]="UBDT > 0.25"
BDT_CUT[Mu_ubdt_veto]="UBDT < 0.25"

declare -A BASE_CUT
BASE_CUT[run1rdx]="DLLmu > 2.0 & DLLe < 1.0 & IsMuon == 1.0"
BASE_CUT[run2ang]="DLLe < 1.0 & IsMuon == 1.0"

declare -A SEL_SUFFIX
SEL_SUFFIX[run1rdx]=""
SEL_SUFFIX[run2ang]="-noDLLMU"

declare -A POLARITY
POLARITY[up]="mu"
POLARITY[down]="md"

PARTICLE=Mu_nopt
BASE_FOLDER=pidcalib_ubdt
if [ -d "${BASE_FOLDER}" ]; then
    rm -r ${BASE_FOLDER}
fi

GLOBAL_CUTS="Brunel_InMuonAcc == 1.0"

for year in 16 17 18; do
    for polarity in "up" "down"; do
        for part in "${!BDT_CUT[@]}"; do
            for sel in "run1rdx" "run2ang"; do
                folder_name="${BASE_FOLDER}/run2-rdx-20${year}-${POLARITY[${polarity}]}-${part}-p_eta_ntracks${SEL_SUFFIX[${sel}]}"
                echo "Output folder: ${folder_name}"
                pidcalib2.make_eff_hists \
                    --output-dir ${folder_name} \
                    --sample "Turbo${year}" --magnet ${polarity} \
                    --particle ${PARTICLE} \
                    --ubdt-version "${sel}" \
                    --cut "${GLOBAL_CUTS}" \
                    --pid-cut "${BASE_CUT[${sel}]} & ${BDT_CUT[${part}]}" \
                    --bin-var Brunel_P --bin-var Brunel_ETA --bin-var nTracks_Brunel \
                    --binning-file ./binning.json
            done
        done
    done
done

# now rename the pkls
PKL_FOLDER=pkl-run2-rdx_mu_ubdt
if [ -d "${PKL_FOLDER}" ]; then
    rm -r ${PKL_FOLDER}
fi
mkdir -p ${PKL_FOLDER}

for pkl in ./${BASE_FOLDER}/*/*.pkl; do
    new_name="$(basename $(dirname ${pkl})).pkl"
    echo "Renaming $pkl to ${PKL_FOLDER}/${new_name}..."
    cp ${pkl} ${PKL_FOLDER}/${new_name}
done

rm -r ${BASE_FOLDER}

# Convert pkls to root
for pkl in ${PKL_FOLDER}/*.pkl; do
    echo "Converting $pkl to root..."
    pidcalib2.pklhisto2root "${pkl}"
done

# Move root files to separate directory
if [ -d "root-run2-rdx_mu_ubdt-tmp" ]; then
    rm -r root-run2-rdx_mu_ubdt-tmp
fi
mkdir -p root-run2-rdx_mu_ubdt-tmp

mv ./${PKL_FOLDER}/*.root ./root-run2-rdx_mu_ubdt-tmp/

# # Shift efficiencies
# if [ -d "root-run2-rdx_mu_ubdt-shifted-tmp" ]; then
#     rm -r root-run2-rdx_mu_ubdt-shifted-tmp
# fi
# mkdir -p root-run2-rdx_mu_ubdt-shifted-tmp

# ../../../scripts/shift_histo_efficiencies.py ./root-run2-rdx_mu_ubdt-tmp ./root-run2-rdx_mu_ubdt-shifted-tmp