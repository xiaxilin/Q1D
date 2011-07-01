! TODO: move initial_soln into solvers

program q1d_primal

  use namelist,        only : read_nml
  use initialize_grid, only : read_grid, cells, faces,   &
                              area_f, area_cc, dx_cc, dadx_cc
  use initialize_soln, only : allocate_soln, initial_soln, prim_cc, cons_cc
  use solvers,         only : explicit_solve

  implicit none

  continue


! read inputs
  call read_nml

! read the grid and allocated grid variables
  call read_grid

! allocate the soln variables
  call allocate_soln(3, cells)

! set up initial flow solution or read a restart
  call initial_soln(cells)

! perform explicit solve
  call explicit_solve(cells, faces, prim_cc, cons_cc,                          &
                      area_f, area_cc, dx_cc, dadx_cc)

! FIXME: when implicit solver becomes live, implement this system
!  case select( solver_type )
!  case('explicit')
!    call explicit_solve
!  case('implicit')
!    call implicit_solve
!  end case select

! do solution output
!  call write_restart
!  call write_plot_data

end program q1d_primal

! makefile layout
