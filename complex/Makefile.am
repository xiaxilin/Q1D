
lib_LIBRARIES =

lib_LIBRARIES += libconstants.a

libconstants_a_SOURCES = \
                    ../libraries/set_precision.f90 \
                    complex_functions.f90 \
                    ../libraries/set_constants.f90 \
                    ../libraries/fluid_constants.f90 \
                    ../libraries/ddt.f90 \
                    ../libraries/set_precision_c.f90 \
                    complex_functions_c.f90 \
                    ../libraries/set_constants_c.f90 \
                    ../libraries/fluid_constants_c.f90 \
                    ../libraries/ddt_c.f90

bin_PROGRAMS =

bin_PROGRAMS += q1d q1d_grid

q1d_SOURCES = \
         ../src/initialize_grid.f90 \
         ../src/initialize_soln.f90 \
         ../src/solution_error.f90 \
         ../src/write_soln.f90 \
	 ../src/matrix_manip.f90 \
	 ../src/residual.f90 \
	 ../src/jacobians.f90 \
         ../src/lhs.f90 \
	 ../src/bc.f90 \
         ../src/solvers.f90 \
         ../src/namelist.f90 \
         ../src/initialize_grid_c.f90 \
         ../src/initialize_soln_c.f90 \
         ../src/solution_error_c.f90 \
         ../src/write_soln_c.f90 \
	 ../src/matrix_manip_c.f90 \
	 ../src/residual_c.f90 \
	 ../src/jacobians_c.f90 \
         ../src/lhs_c.f90 \
	 ../src/bc_c.f90 \
         ../src/solvers_c.f90 \
         ../src/namelist_c.f90 \
         ../src/q1d_c.f90

q1d_LDADD = libconstants.a

q1d_FCFLAGS = $(FC_MODINC)$(top_srcdir)/functions

q1d_grid_SOURCES = ../src/q1d_grid_c.f90

q1d_grid_LDADD = libconstants.a

q1d_grid_FCFLAGS = $(FC_MODINC)$(top_srcdir)/functions

DISTCLEANFILES = *.mod *.o

