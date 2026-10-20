set -e

PARAM_FILE="$(realpath "${1:-./param_file.txt}")"

docker build -t genevol_exec -f Dockerfile.execute .

docker run --rm \
    -v "${PARAM_FILE}:/app/genevol/param_file.txt:ro" \
    genevol_exec