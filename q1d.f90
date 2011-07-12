! TODO: add soln restart output and visualization
! TODO: add convergence checks

program q1d_primal

  use namelist,        only : read_nml
  use fluid_constants, only : set_gamma_constants
  use initialize_grid, only : read_grid, cells, faces, dxsi, area_f, area_cc,  &
                              dxdxsi_cc, dadx_cc, deallocate_grid, x_cc
  use initialize_soln, only : allocate_soln, initial_soln, prim_cc, cons_cc,   &
                              deallocate_soln
  use solvers,         only : explicit_solve
  use write_soln,      only : write_restart, init_write_files, write_entropy

  implicit none

  continue

  print *,'*****************************************************'
  print *,'         Quasi-1D Nozzle : RELEASE 4.0, 2011         '
  print *,'*****************************************************'

! set up some files

! read inputs
  call read_nml

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

! perform explicit solve
  call explicit_solve(cells, faces, dxsi, prim_cc, cons_cc,                    &
                      area_f, area_cc, dxdxsi_cc, dadx_cc, x_cc)

! FIXME: when implicit solver becomes live, implement this system
!  case select( solver_type )
!  case('explicit')
!    call explicit_solve
!  case('implicit')
!    call implicit_solve
!  end case select

! do solution output
  call write_restart(cells, prim_cc)

! write the entropy vars
  call write_entropy(cells, x_cc, prim_cc, cons_cc)

! free memory
  call deallocate_grid
  call deallocate_soln

end program q1d_primal
