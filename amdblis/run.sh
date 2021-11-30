#!/bin/bash

module load gcc/10.2.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/jet/home/babreu/Libraries/hdf5/gcc10.2/lib
./mvmul_amdblis
