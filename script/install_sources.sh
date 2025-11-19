#!/bin/bash

# Define default workdir
WORKDIR=/workspaces
if [[ "$1" == "--workdir" && -n "$2" ]]; then
    WORKDIR="$2"
fi

# Prerequisites installation
sudo apt update
sudo apt-get install libeigen3-dev

# Set installation dir
mkdir ../sources
cd ../sources

# Install bpp-core
git clone https://github.com/BioPP/bpp-core.git
cd bpp-core
git checkout fc0695d523a6060eb8a7696a5a48076df1e349d4 # (switch to latest supported commit)
mkdir build
cd build
cmake  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$WORKDIR ..
make -j 8 # (-j 8 used 8 cores to speed the compilation).
make install
cd ../../ # (getting back to the sources directory)

# Install bpp-seq
git clone https://github.com/BioPP/bpp-seq.git
cd bpp-seq
git checkout 60edc55aa7c81ad1320589b8fbe144045ff865d9
mkdir build
cd build
cmake  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$WORKDIR ..
make -j 8
make install
cd ../../

# Install bpp-phyl
git clone https://github.com/BioPP/bpp-phyl.git
cd bpp-phyl
git checkout ea24a299fe4575f4c51357ed3ff1cb7e9e610f42
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$WORKDIR ..
make -j 4
make install
cd ../../..

