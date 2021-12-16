# module load nvhpc/21.7 on Bridges-2
FC=nvfortran
FFLAGS=-acc -Mcuda

# cuBLAS stuff
CUBLAS_HOME=/jet/packages/nvidia/hpc_sdk/21.7/Linux_x86_64/21.7/math_libs/11.4
LINK_CUBLAS=-L${CUBLAS_HOME}/lib64 -lcublas

# HDF5 stuff
HDF5ROOT	=/jet/home/babreu/Libraries/hdf5/pgi20.11
H5INCLUDE	=-I${HDF5ROOT}/include
H5LDFLAGS	=-L${HDF5ROOT}/lib
H5LIBS		=-lhdf5 -lhdf5_fortran


execs = mvmul_oacc

all: $(execs)

mvmul_oacc: mvmul_h5.o mvmul_oacc.o
	$(FC) $(FFLAGS) -o $@ $^ $(LINK_CUBLAS) ${H5LDFLAGS} ${H5LIBS}

mvmul_oacc.o: mvmul_oacc.f90
	$(FC) $(FFLAGS) ${H5INCLUDE} -c $<

mvmul_h5.o: mvmul_h5.f90
	$(FC) $(FFLAGS) ${H5INCLUDE} -c $<

clean:
	rm *.o *.mod $(execs)
