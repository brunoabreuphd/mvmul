#!/bin/bash

module load nvhpc/21.7
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/jet/home/babreu/Libraries/hdf5/pgi20.11/lib
./mvmul_oacc
