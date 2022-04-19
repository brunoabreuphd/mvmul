FC=nvfortran

FFLAGS=
HDF5ROOT=/u/babreu/Libraries/hdf5/nvhpc
INCLUDE=-I$(HDF5ROOT)/include
LDFLAGS=-L$(HDF5ROOT)/lib
LIBS=-lhdf5 -lhdf5_fortran

objs = mvmul_plain.o
source = mvmul_plain.f90
exec = mvmul_plain

$(exec) : $(objs)
	$(FC) -o $(exec) $(objs) $(LDFLAGS) $(LIBS)

$(objs) : $(source)
	$(FC) $(FFLAGS) $(INCLUDE) -c $(source) 

clean:
	rm $(objs) $(exec)
