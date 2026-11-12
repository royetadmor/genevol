set -e

if [ -z "${2}" ] && [ ! -d "./bootstrap_output" ]; then
    mkdir -p "./bootstrap_output"
fi

PARAM_FILE="$(realpath "${1:-./param_file.txt}")"
OUTPUT_DIR="$(realpath "${2:-./bootstrap_output}")"

docker build -t genevol_bootstrap -f Dockerfile.bootstrap .

docker run --rm \
    -v "${PARAM_FILE}:/app/genevol/param_file.txt:ro" \
    -v "${OUTPUT_DIR}:/app/genevol/bootstrap_output" \
    genevol_bootstrap \
    "${@:3}"
