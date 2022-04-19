FC=nvfortran  

FFLAGS=-Mcuda

# HDF5 vars
HDF5ROOT	=/sw/spack/delta-2022-03/apps/hdf5/1.13.1-nvhpc-22.2-vz2z7bn
H5INCLUDE	=-I${HDF5ROOT}/include
H5LDFLAGS	=-L${HDF5ROOT}/lib
H5LIBS		=-lhdf5 -lhdf5_fortran

# MAGMA vars
MAGMAROOT	=/sw/spack/delta-2022-03/apps/magma/2.6.1-nvhpc-22.2-nbjb6gh
MAGMAINCLUDE	=-I${MAGMAROOT}/include
MAGMALDFLAGS	=-L${MAGMAROOT}/lib
MAGMALIBS	=-lmagma


exec = mvmul_magma

$(exec): mvmul_h5.o mvmul_magma.o
	$(FC) $(FFLAGS) -o $@ $^ ${H5LDFLAGS} ${H5LIBS} ${MAGMALDFLAGS} ${MAGMALIBS}

mvmul_magma.o: mvmul_magma.f90
	$(FC) $(FFLAGS) ${H5INCLUDE} ${MAGMAINCLUDE} -c $<

mvmul_h5.o: mvmul_h5.f90
	$(FC) $(FFLAGS) ${H5INCLUDE} -c $< 

clean:
	rm *.o *.mod $(exec)
