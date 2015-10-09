GC=gfortran
GFLAGS=-Wall -Wextra -pedantic -std=f2008
CFLAGS=
OUT=simul


all: main.o mod_read_parameters.o mod_integrator.o mod_constants.o
	$(GC) $(CFLAGS) $^ -o $(OUT)

main.o: main.f90
	$(GC) $(CFLAGS) -c $^

mod_read_parameters.o: mod_read_parameters.f90
	$(GC) $(CFLAGS) -c $^

mod_integrator.o: mod_integrator.f90
	$(GC) $(CFLAGS) -c $^

mod_constants.o: mod_constants.f90
	$(GC) $(CFLAGS) -c $^

clean:
	\rm -rf *.o $(OUT)
