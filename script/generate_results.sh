#!/usr/bin/env bash
set -euo pipefail

############################
# USER-CONFIGURABLE PATHS #
############################

TEST_DATA_DIR="./test_data"

if [[ $# -ge 1 ]]; then
    RESULTS_DIR="./$1"
else
    echo ""
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "!! WARNING: No results directory name was provided.  !!"
    echo "!! Defaulting to './results'.                        !!"
    echo "!! Usage: $0 <results_dir_name>                !!"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo ""
    RESULTS_DIR="./results"
fi

# Template param file (never modified)
PARAM_TEMPLATE="./param_file_template.txt"

# Working param file used by param_build_and_run.sh
PARAM_WORKING="./param_file.txt"

RUN_SCRIPT="script/param_build_and_run.sh"

################################
# ENSURE RESULTS ROOT EXISTS
################################
if [[ ! -d "${RESULTS_DIR}" ]]; then
    echo "Creating results directory: ${RESULTS_DIR}"
    mkdir -p "${RESULTS_DIR}"
fi

################################
# WALK THROUGH ALL DIRECTORIES
################################
find "${TEST_DATA_DIR}" -type d | while read -r dir; do
    # Find fasta & newick files in this directory (non-recursive)
    fasta_file=$(find "$dir" -maxdepth 1 -type f \( \
        -iname "*.fasta" -o -iname "*.fa" -o -iname "*.faa" \
    \) | head -n 1)

    tree_file=$(find "$dir" -maxdepth 1 -type f \( \
        -iname "*.newick" -o -iname "*.nwk" -o -iname "*.tree" \
    \) | head -n 1)

    # Skip directories without both files
    if [[ -z "${fasta_file}" || -z "${tree_file}" ]]; then
        continue
    fi

    echo "Processing directory: ${dir}"
    echo "  FASTA: ${fasta_file}"
    echo "  TREE:  ${tree_file}"

    ################################
    # PREPARE RESULT DIRECTORY
    ################################
    # Strip leading ./ to preserve tree structure
    rel_path="${dir#./}"
    result_dir="${RESULTS_DIR}/${rel_path}"

    mkdir -p "${result_dir}"

    ################################
    # BUILD PARAM FILE
    ################################
    cp "${PARAM_TEMPLATE}" "${PARAM_WORKING}"

    # Replace tree/data paths (even if commented)
    sed -i.bak \
        -e "s|^[#[:space:]]*_treePath *=.*|_treePath = ${tree_file}|g" \
        -e "s|^[#[:space:]]*_dataPath *=.*|_dataPath = ${fasta_file}|g" \
        "${PARAM_WORKING}"

    rm -f "${PARAM_WORKING}.bak"

    ################################
    # SAVE PARAM COPY
    ################################
    cp "${PARAM_WORKING}" "${result_dir}/params.param"

    ################################
    # RUN AND CAPTURE OUTPUT
    ################################
    (
        cd "$(dirname "${RUN_SCRIPT}")/.."
        bash "${RUN_SCRIPT}"
    ) > "${result_dir}/res.log" 2>&1

    echo "  Results saved in: ${result_dir}"
    echo
done

echo "All runs completed successfully."
