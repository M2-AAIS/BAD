GC=gfortran
GFLAGS= -Wall -Wextra -pedantic -std=f2008
CFLAGS=-llapack
OUT=simul
OUTS=courbeS


all: mod_constants.o mod_read_parameters.o mod_variables.o mod_output.o mod_timestep.o mod_integrator.o mod_S_curve.o main.o 
	$(GC) $^ $(CFLAGS) -o $(OUT)

courbeS: mod_read_parameters.o mod_constants.o mod_variables.o mod_S_curve.o main.o
	$(GC) $^ $(CFLAGS) -o $(OUTS)

main.o: main.f90
	$(GC) $(GFLAGS) -c $^

mod_read_parameters.o: mod_read_parameters.f90
	$(GC) $(GFLAGS) -c $^

mod_integrator.o: mod_integrator.f90
	$(GC) $(GFLAGS) -c $^

mod_timestep.o: mod_timestep.f90
	$(GC) $(GFLAGS) -c $^

mod_variables.o: mod_variables.f90
	$(GC) $(GFLAGS) -c $^

mod_constants.o: mod_constants.f90
	$(GC) $(GFLAGS) -c $^

mod_output.o: mod_output.f90
	$(GC) $(GFLAGS) -c $^

mod_S_curve.o: mod_S_curve.f90
	$(GC) $(GFLAGS) -c $^

clean:
	\rm -rf *.o $(OUT) *.mod

run: all
	./$(OUT)

rerun: clean all run
