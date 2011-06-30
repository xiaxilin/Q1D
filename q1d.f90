! Will hold new q1d_nozzle code

program q1d

  use namelist,        only : read_nml
  use initialize_grid, only : read_grid, cells
  use initialize_soln, only : allocate_soln, initial_soln
  use solvers,         only : explicit_solve
  use XXX, only : 
  use XXX, only : 

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
                      area_f, area_cc, dx_ccdadx_cc)

!  case select( solver_type )
!  case('explicit')
!    call explicit_solve
!  case('implicit')
!    call implicit_solve
!  end case select

! do solution output
  call write_restart
  call write_plot_data

end program q1d

! makefile layout

build all in Libraries

build initialize_grid
build initialize_soln
