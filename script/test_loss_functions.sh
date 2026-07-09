#!/usr/bin/env bash
set -uo pipefail

############################
# USAGE
############################
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <dataset_dir> [results_dir]"
    echo "  dataset_dir  — directory containing a .fasta and a .newick file"
    echo "  results_dir  — where to write output (default: ./results/loss_func_comparison)"
    exit 1
fi

DATASET_DIR="$1"
RESULTS_DIR="${2:-./results/loss_func_comparison}"

PARAM_TEMPLATE="./param_file_template.txt"
PARAM_WORKING="./param_file.txt"
RUN_SCRIPT="script/param_build_and_run.sh"

############################
# LOCATE INPUT FILES
############################
fasta_file=$(find "$DATASET_DIR" -maxdepth 1 -type f \( \
    -iname "*.fasta" -o -iname "*.fa" -o -iname "*.faa" \
\) | while IFS= read -r f; do
    printf '%s %s\n' "$(wc -c < "$f")" "$f"
done | sort -rn | awk 'NR==1{print $2}')

tree_file=$(find "$DATASET_DIR" -maxdepth 1 -type f \( \
    -iname "*.newick" -o -iname "*.nwk" -o -iname "*.tree" \
\) | head -n 1)

if [[ -z "${fasta_file}" || -z "${tree_file}" ]]; then
    echo "ERROR: Could not find both a FASTA and a NEWICK file in '${DATASET_DIR}'."
    exit 1
fi

echo "Dataset:  ${DATASET_DIR}"
echo "FASTA:    ${fasta_file}"
echo "TREE:     ${tree_file}"
echo "Results:  ${RESULTS_DIR}"
echo

mkdir -p "${RESULTS_DIR}"

############################
# STARTING POINTS PER FUNCTION
# Each entry is a comma-separated gain init string; entries are space-separated.
############################

# Returns a newline-separated list of starting points for the given function
starts_for() {
    case "$1" in
        CONST)
            echo "0.1"
            echo "1"
            echo "5"
            ;;
        LINEAR_BD)
            echo "0.1"
            echo "1"
            echo "5"
            ;;
        LINEAR)
            echo "0.1,0"
            echo "1,0"
            echo "5,0"
            echo "0.1,0.5"
            echo "1,0.5"
            echo "1,1"
            ;;
        EXP)
            echo "1,0"
            echo "0.1,0"
            echo "5,0"
            echo "1,-0.5"
            echo "0.1,-0.5"
            ;;
        POLYNOMIAL)
            echo "1,1,0.5"
            echo "0.1,1,0.5"
            echo "1,2,1"
            echo "0.1,2,1"
            ;;
    esac
}

############################
# CSV HEADER
############################
CSV="${RESULTS_DIR}/summary.csv"
echo "lossFunc,best_start,AIC,params" > "${CSV}"

############################
# RUN EACH FUNCTION
############################
FUNCS=("CONST" "LINEAR_BD" "LINEAR" "EXP" "POLYNOMIAL")

for func in "${FUNCS[@]}"; do
    echo "--- Running lossFunc = ${func} ---"

    best_aic=""
    best_start=""
    best_log=""
    best_params=""
    run_dir="${RESULTS_DIR}/${func}_runs"
    mkdir -p "${run_dir}"

    # Iterate over starting points
    while IFS= read -r start; do
        safe_start="${start//,/_}"
        log_file="${run_dir}/${safe_start}.log"
        param_file="${run_dir}/${safe_start}.param"

        echo "  start=[${start}]"

        cp "${PARAM_TEMPLATE}" "${PARAM_WORKING}"
        sed -i.bak \
            -e "s|^[#[:space:]]*_treePath *=.*|_treePath = ${tree_file}|g" \
            -e "s|^[#[:space:]]*_dataPath *=.*|_dataPath = ${fasta_file}|g" \
            -e "s|^[#[:space:]]*_lossFunc *=.*|_lossFunc = ${func}|g" \
            -e "s|^[#[:space:]]*_loss *=.*|_loss = ${start}|g" \
            "${PARAM_WORKING}"
        rm -f "${PARAM_WORKING}.bak"
        cp "${PARAM_WORKING}" "${param_file}"

        if ! (cd "$(dirname "${RUN_SCRIPT}")/.." && bash "${RUN_SCRIPT}") > "${log_file}" 2>&1; then
            echo "    FAILED — skipping"
            continue
        fi

        aic=$(grep "^AIC score:" "${log_file}" | tail -1 | awk '{print $3}')
        if [[ -z "${aic}" ]]; then
            echo "    no AIC found — skipping"
            continue
        fi

        echo "    AIC = ${aic}"

        # Keep track of best (lowest AIC) using awk for float comparison
        if [[ -z "${best_aic}" ]] || awk "BEGIN{exit !(${aic} < ${best_aic})}"; then
            best_aic="${aic}"
            best_start="${start}"
            best_log="${log_file}"
        fi
    done < <(starts_for "$func")

    if [[ -z "${best_aic}" ]]; then
        echo "  All starts FAILED"
        echo "${func},,FAILED,\"\"" >> "${CSV}"
        echo
        continue
    fi

    echo "  Best start=[${best_start}]  AIC=${best_aic}"

    # Copy best log to the top-level result for this function
    cp "${best_log}" "${RESULTS_DIR}/${func}.log"
    cp "${run_dir}/${best_start//,/_}.param" "${RESULTS_DIR}/${func}.param"

    params=$(grep -E "^  [A-Za-z]" "${best_log}" | tail -20 \
        | awk '{gsub(/^ +/, ""); print}' \
        | paste -sd '|' -)

    echo "${func},\"${best_start}\",${best_aic},\"${params}\"" >> "${CSV}"
    echo
done

############################
# PRINT SUMMARY
############################
echo "========================================"
echo "Summary written to: ${CSV}"
echo "========================================"
column -t -s',' "${CSV}"
