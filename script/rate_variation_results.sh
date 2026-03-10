#!/usr/bin/env bash
set -euo pipefail

############################
# USER-CONFIGURABLE PATHS #
############################

if [[ $# -ge 1 ]]; then
    RESULTS_DIR="./$1"
else
    echo ""
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "!! WARNING: No results directory name was provided.  !!"
    echo "!! Defaulting to './results_rate_variation'.         !!"
    echo "!! Usage: $0 <results_dir_name>                !!"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo ""
    RESULTS_DIR="./results_rate_variation"
fi

# Template param file (never modified)
PARAM_TEMPLATE="./param_file_template.txt"

# Working param file used by param_build_and_run.sh
PARAM_WORKING="./param_file.txt"

RUN_SCRIPT="script/param_build_and_run.sh"

# Gamma distribution parameters
GAMMA_ALPHA=1.5
N_VALUES=(16 32)

################################
# ENSURE RESULTS ROOT EXISTS
################################
if [[ ! -d "${RESULTS_DIR}" ]]; then
    echo "Creating results directory: ${RESULTS_DIR}"
    mkdir -p "${RESULTS_DIR}"
fi

################################
# ITERATE OVER CATEGORY COUNTS
################################
for n in "${N_VALUES[@]}"; do
    echo "Running with Gamma n=${n} (alpha=${GAMMA_ALPHA}) ..."

    ################################
    # PREPARE RESULT DIRECTORY
    ################################
    result_dir="${RESULTS_DIR}/n_${n}"
    mkdir -p "${result_dir}"

    ################################
    # BUILD PARAM FILE
    ################################
    cp "${PARAM_TEMPLATE}" "${PARAM_WORKING}"

    # Replace rate_distribution with Gamma(n, alpha), regardless of current setting
    sed -i.bak \
        -e "s|^[#[:space:]]*rate_distribution *=.*|rate_distribution = Gamma(n=${n}, alpha=${GAMMA_ALPHA})|g" \
        "${PARAM_WORKING}"

    rm -f "${PARAM_WORKING}.bak"

    ################################
    # SAVE PARAM COPY
    ################################
    cp "${PARAM_WORKING}" "${result_dir}/params.param"

    ################################
    # RUN AND CAPTURE OUTPUT (timed)
    ################################
    start_time=$(date +%s%N)

    (
        cd "$(dirname "${RUN_SCRIPT}")/.."
        bash "${RUN_SCRIPT}"
    ) > "${result_dir}/res.log" 2>&1

    end_time=$(date +%s%N)
    elapsed_ms=$(( (end_time - start_time) / 1000000 ))
    elapsed_sec=$(echo "scale=3; ${elapsed_ms} / 1000" | bc)

    echo "  Elapsed time: ${elapsed_sec}s"
    echo "${elapsed_sec}" > "${result_dir}/time.txt"

    echo "  Results saved in: ${result_dir}"
    echo
done

echo "All runs completed successfully."
