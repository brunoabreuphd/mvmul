FC=ifort   # module load intel/20.4 on Bridges-2

FFLAGS=-w
HDF5ROOT=/jet/home/babreu/Libraries/hdf5/intel20.4
H5INCLUDE=-I${HDF5ROOT}/include
H5LDFLAGS=-L${HDF5ROOT}/lib
H5LIBS=-lhdf5 -lhdf5_fortran
MKLROOT=/opt/intel/compilers_and_libraries_2020.4.304/linux/mkl
MKLPATH=${MKLROOT}/lib/intel64
MKLINCLUDE=${MKLROOT}/include
MKLLDFLAGS=-mkl=parallel -L${MKLPATH} -lmkl_lapack95_lp64 -Wl,--start-group ${MKLPATH}/libmkl_intel_lp64.a ${MKLPATH}/libmkl_intel_thread.a ${MKLPATH}/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm

exec = mvmul_intelmkl

$(exec): mvmul_h5.o mvmul_intelmkl.o
	$(FC) -o $@ $^ ${H5LDFLAGS} ${H5LIBS} ${MKLLDFLAGS}

mvmul_intelmkl.o: mvmul_intelmkl.f90
	$(FC) $(FFLAGS) ${H5INCLUDE} -I${MKLINCLUDE} -I${MKLINCLUDE}/intel64/lp64 -c $<

mvmul_h5.o: mvmul_h5.f90
	$(FC) $(FFLAGS) ${H5INCLUDE} -c $< 

clean:
	rm *.o *.mod $(exec)
