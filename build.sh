#!/bin/sh

module purge
module load gcc mkl hpcx/2.18-mt cmake
#module load intel/2023.1 compiler mkl/2023.1.0 hpcx/2.15-mt-icc
#module load intel/2024.0 compiler impi mkl
#module load /global/scratch/groups/nvidia/hpcx-2.15.0-umr/modulefiles/hpcx.umr

rm -rf build
mkdir build
cd build
#/global/scratch/groups/mellanox/ecmwf/fiat
#ecbuild_ROOT=$HOME/ecbuild fiat_ROOT=$HOME/FIAT-inst CC="scorep mpicc" FC="scorep mpif90" cmake -DCMAKE_INSTALL_PREFIX=$HOME/ECTRANS-scorep -DCMAKE_BUILD_TYPE=Debug .. 2>&1 | tee configure.log
ecbuild_ROOT=$HOME/ecbuild fiat_ROOT=$HOME/FIAT_papi_218 CC="mpicc" FC="mpif90" cmake -DCMAKE_INSTALL_PREFIX=$HOME/ECTRANS-main-transpose -DCMAKE_BUILD_TYPE=Release .. 2>&1 | tee configure.log
make install 2>&1 | tee make.log
