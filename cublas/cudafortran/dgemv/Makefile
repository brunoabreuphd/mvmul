# module load nvhpc/21.7 on Bridges-2
FC=nvfortran
FFLAGS=-Mcuda

execs = dgemv

all: $(execs)

dgemv: dgemv.o
	$(FC) $(FFLAGS) -o $@ $^
dgemv.o: dgemv.f90
	$(FC) $(FFLAGS) -c $^

clean:
	rm *.o *.mod $(execs)
