#!/bin/bash

set -e

docker build -t genevol -f Dockerfile .

script/infra/docker_run.sh