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
git checkout master # (switch to latest supported commit)
mkdir build
cd build
cmake  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$WORKDIR ..
make -j 8 # (-j 8 used 8 cores to speed the compilation).
make install
cd ../../ # (getting back to the sources directory)

# Install bpp-seq
git clone https://github.com/BioPP/bpp-seq.git
cd bpp-seq
git checkout master
mkdir build
cd build
cmake  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$WORKDIR ..
make -j 8
make install
cd ../../

# Install bpp-phyl
git clone https://github.com/BioPP/bpp-phyl.git
cd bpp-phyl
git checkout master
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$WORKDIR ..
make -j 4
make install
cd ../../..

