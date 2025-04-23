#!/bin/bash


docker build -t genevol -f Dockerfile .

docker run --rm -it genevol   