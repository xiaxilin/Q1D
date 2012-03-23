
lib_LIBRARIES =

lib_LIBRARIES += libconstants.a

libconstants_a_SOURCES = \
                    ../Libraries/set_precision.f90 \
                    complex_functions.f90 \
                    ../Libraries/set_constants.f90 \
                    ../Libraries/fluid_constants.f90 \
                    ../Libraries/ddt.f90 \
                    ../Libraries/set_precision_c.f90 \
                    complex_functions_c.f90 \
                    ../Libraries/set_constants_c.f90 \
                    ../Libraries/fluid_constants_c.f90 \
                    ../Libraries/ddt_c.f90

bin_PROGRAMS =

bin_PROGRAMS += q1d q1d_grid

q1d_SOURCES = \
         ../Src/initialize_grid.f90 \
         ../Src/initialize_soln.f90 \
         ../Src/solution_error.f90 \
         ../Src/write_soln.f90 \
	 ../Src/matrix_manip.f90 \
	 ../Src/residual.f90 \
	 ../Src/jacobians.f90 \
         ../Src/lhs.f90 \
	 ../Src/bc.f90 \
         ../Src/solvers.f90 \
         ../Src/namelist.f90 \
         ../Src/initialize_grid_c.f90 \
         ../Src/initialize_soln_c.f90 \
         ../Src/solution_error_c.f90 \
         ../Src/write_soln_c.f90 \
	 ../Src/matrix_manip_c.f90 \
	 ../Src/residual_c.f90 \
	 ../Src/jacobians_c.f90 \
         ../Src/lhs_c.f90 \
	 ../Src/bc_c.f90 \
         ../Src/solvers_c.f90 \
         ../Src/namelist_c.f90 \
         ../Src/q1d_c.f90

q1d_LDADD = libconstants.a

q1d_FCFLAGS = $(FC_MODINC)$(top_srcdir)/Functions

q1d_grid_SOURCES = ../Src/q1d_grid_c.f90

q1d_grid_LDADD = libconstants.a

q1d_grid_FCFLAGS = $(FC_MODINC)$(top_srcdir)/Functions

DISTCLEANFILES = *.mod *.o

