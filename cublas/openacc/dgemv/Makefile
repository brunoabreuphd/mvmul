# module load nvhpc/21.7 on Bridges-2
FC=nvfortran
FFLAGS=-acc -Mcuda
CUBLAS_HOME=/jet/packages/nvidia/hpc_sdk/21.7/Linux_x86_64/21.7/math_libs/11.4
LINK_CUBLAS=-L${CUBLAS_HOME}/lib64 -lcublas

execs = dgemv_oacc

all: $(execs)

dgemv_oacc: dgemv_oacc.o
	$(FC) $(FFLAGS) -o $@ $^ $(LINK_CUBLAS)
dgemv_oacc.o: dgemv_oacc.f90
	$(FC) $(FFLAGS) -c $^

clean:
	rm *.o *.mod $(execs)
