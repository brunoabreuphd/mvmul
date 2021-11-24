FC=gfortran

FFLAGS=-Wextra -Wall -pedantic -m64
HDF5ROOT=/jet/home/babreu/Libraries/hdf5/gcc10.2
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
