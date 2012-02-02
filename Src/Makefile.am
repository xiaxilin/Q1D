
bin_PROGRAMS =

bin_PROGRAMS += q1d q1d_grid

q1d_SOURCES = \
         initialize_grid.f90 \
         initialize_soln.f90 \
         solution_error.f90 \
         write_soln.f90 \
	 matrix_manip.f90 \
	 jacobians.f90 \
	 bc.f90 \
         solvers.f90 \
         namelist.f90 \
         q1d.f90

q1d_LDADD = ../Libraries/libconstants.a

q1d_FCFLAGS = $(FC_MODINC)$(top_srcdir)/Functions $(FC_MODINC)../Libraries

q1d_grid_SOURCES = q1d_grid.f90

q1d_grid_LDADD = ../Libraries/libconstants.a

q1d_grid_FCFLAGS = $(FC_MODINC)$(top_srcdir)/Functions $(FC_MODINC)../Libraries

DISTCLEANFILES = *.mod *.o

