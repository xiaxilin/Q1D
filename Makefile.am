
AUTOMAKE_OPTIONS = 1.7

SUBDIRS = Libraries

bin_PROGRAMS =

bin_PROGRAMS += q1d

q1d_SOURCES = initialize_grid.f90 initialize_soln.f90 solvers.f90 namelist.f90 q1d.f90
q1d_LDADD = Libraries/libconstants.a

q1d_FCFLAGS = $(FC_MODINC)Functions $(FC_MODINC)Libraries



