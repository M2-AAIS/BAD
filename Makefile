GC=gfortran
GFLAGS=-Wall -Wextra -pedantic -std=f2008 -g -O0 -fno-second-underscore -Wno-compare-reals -fcheck=all -fimplicit-none -ffpe-trap=zero,overflow -fbacktrace -gdwarf-2 -fall-intrinsics -Wno-unused-function -fcheck=bounds
CFLAGS=-llapack -g -O0
OUT=simul
OUTS=s_curve


all: mod_constants.o mod_read_parameters.o mod_variables.o mod_output.o mod_timestep.o mod_integrator.o mod_s_curve.o mod_critical.o main.o 
	$(GC) $^ $(CFLAGS) -o $(OUT)

s_curve: mod_constants.o mod_read_parameters.o mod_variables.o mod_s_curve.o mod_critical.o main_curve.o
	$(GC) $^ $(CFLAGS) -o $(OUTS)

main.o: main.f90
	$(GC) $(GFLAGS) -c $^

main_curve.o: main_curve.f90
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

mod_s_curve.o: mod_s_curve.f90
	$(GC) $(GFLAGS) -c $^

mod_critical.o: mod_critical.f90
	$(GC) $(GFLAGS) -c $^

clean:
	\rm -rf *.o $(OUT) $(OUTS) *.mod

run: all
	./$(OUT)

run_s_curve: s_curve
	./$(OUTS)

rerun: clean all run

rebuild: clean all s_curve
