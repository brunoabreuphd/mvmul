FC=nvfortran

FFLAGS=

# HDF5 stuff
HDF5ROOT	=/sw/spack/delta-2022-03/apps/hdf5/1.13.1-nvhpc-22.2-vz2z7bn
H5INCLUDE	=-I${HDF5ROOT}/include
H5LDFLAGS	=-L${HDF5ROOT}/lib
H5LIBS		=-lhdf5 -lhdf5_fortran

objs = mvmul_plain.o
source = mvmul_plain.f90
exec = mvmul_plain.exe

$(exec) : $(objs)
	$(FC) -o $(exec) $(objs) ${H5LDFLAGS} ${H5LIBS}

$(objs) : $(source)
	$(FC) ${FFLAGS} ${H5INCLUDE} -c $(source) 

clean:
	rm $(objs) $(exec)
