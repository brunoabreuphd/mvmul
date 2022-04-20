#!/bin/bash

module purge
module load nvhpc/22.2
module load openmpi/4.1.2
module load hdf5/1.13.1

./mvmul_plain.exe
