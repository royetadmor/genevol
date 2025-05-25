#!/bin/bash

set -e

docker build -t genevol -f Dockerfile .

script/docker_run.sh