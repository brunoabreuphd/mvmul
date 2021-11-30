FC=gfortran   # module load gcc/10.2.0 on Bridges-2

FFLAGS=-Wall

# HDF5 stuff
HDF5ROOT	=/jet/home/babreu/Libraries/hdf5/gcc10.2
H5INCLUDE	=-I${HDF5ROOT}/include
H5LDFLAGS	=-L${HDF5ROOT}/lib
H5LIBS		=-lhdf5 -lhdf5_fortran

# AMD BLIS stuff
BLIS_PREFIX	= /jet/packages/spack/opt/spack/linux-centos8-zen2/gcc-10.2.0/amdblis-2.2-df2bmc2zbrcgkqw5knutfcscuasntd47
BLIS_INC	= ${BLIS_PREFIX}/include/blis
BLIS_LIB	= ${BLIS_PREFIX}/lib/libblis-mt.a -lm -lpthread -fopenmp


exec = mvmul_amdblis

$(exec): mvmul_h5.o mvmul_amdblis.o
	$(FC) -o $@ $^ ${H5LDFLAGS} ${H5LIBS} ${BLIS_LIB}

mvmul_amdblis.o: mvmul_amdblis.f90
	$(FC) $(FFLAGS) ${H5INCLUDE} -I${BLIS_INC} -c $<

mvmul_h5.o: mvmul_h5.f90
	$(FC) $(FFLAGS) ${H5INCLUDE} -c $< 

clean:
	rm *.o *.mod $(exec)
