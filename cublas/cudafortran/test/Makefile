# module load nvhpc/21.7 on Bridges-2
FC=nvfortran
FFLAGS=-Mcuda

execs = test

all: $(execs)

test: test.o
	$(FC) $(FFLAGS) -o $@ $^
test.o: test.f90
	$(FC) $(FFLAGS) -c $^

clean:
	rm *.o *.mod $(execs)
