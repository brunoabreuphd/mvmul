FC=nvfortran  

FFLAGS=-Mcuda

# HDF5 stuff
HDF5ROOT	=/u/babreu/Libraries/hdf5/nvhpc
H5INCLUDE	=-I${HDF5ROOT}/include
H5LDFLAGS	=-L${HDF5ROOT}/lib
H5LIBS		=-lhdf5 -lhdf5_fortran

# MAGMA stuff
MAGMAROOT	=/sw/spack/delta-2022-03/apps/magma/2.6.1-gcc-11.2.0-bsiktjp
MAGMAINCLUDE	=-I${MAGMAROOT}/include
MAGMALDFLAGS	=-L${MAGMAROOT}/lib
MAGMALIBS	=-lmagma


exec = mvmul_magma

$(exec): mvmul_h5.o mvmul_magma.o
	$(FC) $(FFLAGS) -o $@ $^ ${H5LDFLAGS} ${H5LIBS} ${MAGMAFLAGS} ${MAGMALIBS}

mvmul_magma.o: mvmul_magma.f90
	$(FC) $(FFLAGS) ${H5INCLUDE} ${MAGMAINCLUDE} -c $<

mvmul_h5.o: mvmul_h5.f90
	$(FC) $(FFLAGS) ${H5INCLUDE} -c $< 

clean:
	rm *.o *.mod $(exec)
