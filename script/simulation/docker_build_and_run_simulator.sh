#!/bin/bash

set -e

docker build -t genevol -f Dockerfile .

script/simulation/docker_run_simulator.sh "$@"
