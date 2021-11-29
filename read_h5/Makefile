FC=gfortran

FFLAGS=-Wextra -Wall -pedantic -m64
HDF5ROOT=/jet/home/babreu/Libraries/hdf5/gcc10.2
INCLUDE=-I$(HDF5ROOT)/include
LDFLAGS=-L$(HDF5ROOT)/lib
LIBS=-lhdf5 -lhdf5_fortran

exec = mvmul_data

$(exec): mvmul_h5.o mvmul_data.o
	$(FC) -o $@ $^ $(LDFLAGS) $(LIBS)

mvmul_data.o: mvmul_data.f90
	$(FC) $(FFLAGS) $(INCLUDE) -c $<

mvmul_h5.o: mvmul_h5.f90
	$(FC) $(FFLAGS) $(INCLUDE) -c $< 

clean:
	rm *.o *.mod $(exec)
