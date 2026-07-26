set -e

docker build -t genevol_exec -f Dockerfile.execute . 

docker run --rm genevol_exec