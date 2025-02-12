# Set installation dir
mkdir ../sources
cd ../sources

# Install bpp-core
git clone https://github.com/anatshafir1/bpp-core.git
cd bpp-core
git checkout rel_anat # (switch to rel_anat branch)
mkdir build
cd build
cmake  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=/workspaces ..
make -j 8 # (-j 8 used 8 cores to speed the compilation).
make install
cd ../../ # (getting back to the sources directory)

# Install bpp-seq
git clone https://github.com/anatshafir1/bpp-seq.git
cd bpp-seq
git checkout rel_anat
mkdir build
cd build
cmake  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=/workspaces ..
make -j 8
make install
cd ../../

# Install bpp-phyl
git clone https://github.com/anatshafir1/bpp-phyl.git
cd bpp-phyl
git checkout rel3_anat
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=/workspaces ..
make -j 1
make install
cd ../../..

