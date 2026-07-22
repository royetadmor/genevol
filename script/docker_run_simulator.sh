#!/bin/bash

set -e

PARAM_FILE="$(realpath "${1:-./simulation_param_file.txt}")"
OUTPUT_DIR="${2:-./sim_output}"
if [[ ! -d "${OUTPUT_DIR}" ]]; then
    mkdir -p "${OUTPUT_DIR}"
fi
OUTPUT_DIR="$(realpath "${OUTPUT_DIR}")"

PARAM_DIR="$(dirname "${PARAM_FILE}")"
PARAM_FILENAME="$(basename "${PARAM_FILE}")"

docker run --rm \
    -v "${PARAM_FILE}:/app/genevol/${PARAM_FILENAME}:ro" \
    -v "${OUTPUT_DIR}:/app/genevol/sim_output" \
    genevol \
    ./Simulator/simulator "param=/app/genevol/${PARAM_FILENAME}"
