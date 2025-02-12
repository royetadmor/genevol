# GenEvol

GenEvol is a likelihood-based method for estimating gene family dynamics, including gene gain, loss, innovation, and duplication in plants. It provides a robust framework for analyzing evolutionary changes in gene families across different plant species.

## Installation

To install GenEvol, run the following script:

```sh
script/install_sources.sh
```
Once the sources are installed, execute the following commands:

```sh
cmake -DCMAKE_INSTALL_PREFIX=<working_dir> OMP_NUM_THREADS=20 ./
make -j 3
make install
```
