# module load nvhpc/21.7 on Bridges-2
FC=nvfortran
FFLAGS=-acc -Mcuda
CUBLAS_HOME=/jet/packages/nvidia/hpc_sdk/21.7/Linux_x86_64/21.7/math_libs/11.4
LINK_CUBLAS=-L${CUBLAS_HOME}/lib64 -lcublas

execs = oacc_host_example

all: $(execs)

oacc_host_example: oacc_host_example.o
	$(FC) $(FFLAGS) -o $@ $^ $(LINK_CUBLAS)
oacc_host_example.o: oacc_host_example.f90
	$(FC) $(FFLAGS) -c $^

clean:
	rm *.o *.mod $(execs)
