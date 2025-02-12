#!/bin/bash

# Define default workdir
WORKDIR=/workspaces
if [[ "$1" == "--workdir" && -n "$2" ]]; then
    WORKDIR="$2"
fi

# Build project
cmake -DCMAKE_INSTALL_PREFIX=$WORKDIR OMP_NUM_THREADS=20 ./
make -j 3
make install