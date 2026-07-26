#!/bin/bash
set -e

SIM_PARAM="${1:-./simulation_param_file.txt}"
OUTPUT_DIR="${2:-./sim_output}"
GENEVOL_PARAM="${3:-./param_file.txt}"

OUTPUT_DIR_ABS="$(realpath "$OUTPUT_DIR")"
SIM_PARAM_ABS="$(realpath "$SIM_PARAM")"
GENEVOL_PARAM_ABS="$(realpath "$GENEVOL_PARAM")"
GENEVOL_PARAM_FILENAME="$(basename "$GENEVOL_PARAM_ABS")"

# Step 1: Simulate
echo "=== Step 1: Simulating data ==="
script/simulation/docker_run_simulator.sh "$SIM_PARAM" "$OUTPUT_DIR"

# Step 2: Point param_file to simulated data and run genevol
echo ""
echo "=== Step 2: Fitting model to simulated data ==="
sed -i '' "s|_dataPath[[:space:]]*=.*|_dataPath = sim_output/simulated.fasta|" "$GENEVOL_PARAM_ABS"

GENEVOL_LOG="$(docker run --rm \
    -v "${OUTPUT_DIR_ABS}:/app/genevol/sim_output:ro" \
    -v "${GENEVOL_PARAM_ABS}:/app/genevol/${GENEVOL_PARAM_FILENAME}:ro" \
    genevol \
    ./GenEvol/genEvol "param=/app/genevol/${GENEVOL_PARAM_FILENAME}" 2>&1)"

echo "$GENEVOL_LOG"

# Step 3: Compare recovered params to ground truth
echo ""
echo "=== Step 3: Parameter comparison ==="
python3 script/simulation/compare_params.py "$SIM_PARAM_ABS" <<< "$GENEVOL_LOG"
