# mvmul

This is a series of applications that perform Matrix-vector multiplication using different libraries and compilers. The goal is to have a measure of how much the final results fluctuate from build to build. For that reason, the matrix, the vector and the resulting product vector are calculated using basic Fortran data types, without calling any libraries (via straightforward implementation of loops), and the results are stored into HDF5 files, that can then be opened and read by different applications without precision loss.
