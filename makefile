FC=gfortran
FCFLAGS=-O

all: q1d_grid

q1d_grid: set_precision.o set_constants.o q1d_grid.o
	$(FC) $(FCFLAGS) set_precision.o set_constants.o q1d_grid.o -o q1d_grid

q1d_grid.o: q1d_grid.f90
	$(FC) $(FCFLAGS) -c set_precision.f90 set_constants.f90 q1d_grid.f90

set_constants.o: set_constants.f90
	$(FC) $(FCFLAGS) -c set_precision.f90 set_constants.f90

set_precision.o: set_precision.f90
	$(FC) $(FCFLAGS) -c set_precision.f90

clean:
	rm -rf *o *mod

veryclean:
	rm -rf *o *mod *grd q1d_grid