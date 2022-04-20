FC=nvfortran   # module load nvhpc/22.2 on Delta

FFLAGS=-Mcuda

# HDF5 stuff
HDF5ROOT	=/sw/spack/delta-2022-03/apps/hdf5/1.13.1-nvhpc-22.2-vz2z7bn
H5INCLUDE	=-I${HDF5ROOT}/include
H5LDFLAGS	=-L${HDF5ROOT}/lib
H5LIBS		=-lhdf5 -lhdf5_fortran


exec = mvmul_cudafort.exe

$(exec): mvmul_h5.o mvmul_cudafort.o
	$(FC) $(FFLAGS) -o $@ $^ ${H5LDFLAGS} ${H5LIBS}

mvmul_cudafort.o: mvmul_cudafort.f90
	$(FC) $(FFLAGS) ${H5INCLUDE} -c $<

mvmul_h5.o: mvmul_h5.f90
	$(FC) $(FFLAGS) ${H5INCLUDE} -c $< 

clean:
	rm *.o *.mod $(exec)
