#!/bin/bash

module load intel/20.4
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/jet/home/babreu/Libraries/hdf5/intel20.4/lib
./mvmul_intelmkl
