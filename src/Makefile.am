bin_PROGRAMS =

bin_PROGRAMS += q1d q1da q1d_grid

SAM_INCLUDE = -I ~/SAM/examples/Example_Tanh/build
SAM_LIB = ~/SAM/examples/Example_Tanh/build/sam.o

q1d_SOURCES = \
         initialize_grid.f90 \
         initialize_soln.f90 \
         write_soln.f90 \
	 matrix_manip.f90 \
         residual.f90 \
	 jacobians.f90 \
         lhs.f90 \
	 bc.f90 \
         solvers.f90 \
         solution_error.f90 \
         namelist.f90 \
         q1d.f90

q1d_LDADD = ../libraries/libconstants.a $(SAM_LIB)

q1d_FCFLAGS =  $(SAM_INCLUDE) $(FC_MODINC)$(top_srcdir)/functions $(FC_MODINC)../libraries

q1da_SOURCES = \
          initialize_grid.f90 \
          initialize_soln.f90 \
          solution_error.f90 \
          write_soln.f90 \
	  matrix_manip.f90 \
          residual.f90 \
	  jacobians.f90 \
          lhs.f90 \
	  bc.f90 \
          solvers.f90 \
          namelist.f90 \
          adjoint_lhs.f90 \
          adjoint_solvers.f90 \
          q1da.f90

q1da_LDADD = ../libraries/libconstants.a $(SAM_LIB)

q1da_FCFLAGS =  $(SAM_INCLUDE) $(FC_MODINC)$(top_srcdir)/functions $(FC_MODINC)../libraries

q1d_grid_SOURCES = q1d_grid.f90

q1d_grid_LDADD = ../libraries/libconstants.a

q1d_grid_FCFLAGS = $(FC_MODINC)$(top_srcdir)/functions $(FC_MODINC)../libraries

DISTCLEANFILES = *.mod *.o

