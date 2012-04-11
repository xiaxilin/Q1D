! TODO: add soln restart output and visualization
! TODO: add convergence checks

program q1d_adjoint

  use namelist,        only : read_nml
  use fluid_constants, only : set_gamma_constants
  use initialize_grid, only : read_grid, cells, faces, dx, x_cc, cell_vol,     &
                              area_f, area_cc, dadx_cc, deallocate_grid
  use initialize_soln, only : allocate_soln, initial_soln, prim_cc, cons_cc,   &
                              deallocate_soln, restart
  use solvers,         only : solver
  use adjoint_solvers, only : implicit_solve
  use write_soln,      only : init_write_files

  implicit none

  continue

  print *,'*****************************************************'
  print *,'         Quasi-1D Nozzle : RELEASE 5.0, 2012         '
  print *,'*****************************************************'

! set up some files

! read inputs
  call read_nml

  solver = 'implicit'
  restart = .true.

! init files
  call init_write_files

! set derived fluid constants
  call set_gamma_constants()

! read the grid and allocated grid variables
  call read_grid

! allocate the soln variables
  call allocate_soln(3, cells)

! set up initial flow solution or read a restart
  call initial_soln(cells)

  restart = .false.

  print*, 'Beginning Implicit Solve'
  call implicit_solve(cells, faces, prim_cc, cons_cc, cell_vol,              &
                      area_f, dx, dadx_cc, x_cc)

! free memory

  print *, 'Deallocating Memory'
  call deallocate_grid
  call deallocate_soln

  print *, 'Finished'

end program q1d_adjoint
