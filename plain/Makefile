FC=gfortran

FFLAGS=-Wextra -Wall -pedantic -m64 
LIBS=
LDFLAGS=

objs = mvmul_plain.o
source = mvmul_plain.f90
exec = mvmul_plain

$(exec) : $(objs)
	$(FC) -o $(exec) $(objs)

$(objs) : $(source)
	$(FC) $(FFLAGS) $(LIBS) -c $(source) 

clean:
	rm $(objs) $(exec)
