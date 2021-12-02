# mvmul

This is a series of applications that perform Matrix-vector multiplication using different libraries and compilers. The goal is to have a measure of how much the final results fluctuate from build to build. For that reason, the matrix, the vector and the resulting product vector are calculated using basic Fortran data types, without calling any libraries (via straightforward implementation of loops), and the results are stored into HDF5 files, that can then be opened and read by different applications without precision loss. Makefiles and run scripts are provided, but keep in mind that they have environment variables that are set to the machines that I am using. You will need to update them.

## [plain](./plain)
This implements the matrix-vector multiplication using loops that enclose the direct definition of the operation. All of the numbers are chosen to be 64-bits floating points. The input matrix and vector elements are selected from a uniform distribution on the unitary interval using Fortran's intrinsic ```rand```, therefore having an expected value of 0.5 (in a probabilistic sense).

These three data structures (matrix and two vectors) are then stored into [HDF5](https://www.hdfgroup.org/solutions/hdf5/) datasets, enclosed in a single file. You will need to know where the HDF5 library is installed in your system to make this work (or install it, which is not too hard). 

I have intentionally and explicitly set the vector dimensions, so that a generalization to matrix-matrix multiplication is easily achieved. If you are interested in performance, I also included a stopwatch around the ```mvmul``` operation. 

## [read_h5](./read_h5)


## [intelmkl](./intelmkl)


## [amdblis](./amdblis)
