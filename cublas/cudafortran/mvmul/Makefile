FC=nvfortran   # module load nvhpc/21.7 on Bridges-2

FFLAGS=-Mcuda

# HDF5 stuff
HDF5ROOT	=/jet/home/babreu/Libraries/hdf5/pgi20.11
H5INCLUDE	=-I${HDF5ROOT}/include
H5LDFLAGS	=-L${HDF5ROOT}/lib
H5LIBS		=-lhdf5 -lhdf5_fortran


exec = mvmul_cudafort

$(exec): mvmul_h5.o mvmul_cudafort.o
	$(FC) $(FFLAGS) -o $@ $^ ${H5LDFLAGS} ${H5LIBS}

mvmul_cudafort.o: mvmul_cudafort.f90
	$(FC) $(FFLAGS) ${H5INCLUDE} -c $<

mvmul_h5.o: mvmul_h5.f90
	$(FC) $(FFLAGS) ${H5INCLUDE} -c $< 

clean:
	rm *.o *.mod $(exec)
