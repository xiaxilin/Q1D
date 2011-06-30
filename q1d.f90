! Will hold new q1d_nozzle code

program q1d

  use XXX, only : 
  use XXX, only : 
  use XXX, only : 
  use XXX, only : 
  use XXX, only : 

  implicit none

  call read_nml

  call read_grid

  call allocate_soln(cells)

  call initialize_soln
    if ( restart ) call read_restart inside initialize_soln

!  if ( explicit ) then
    call explicit_solve
!  else
!    call implicit_solve
!  end if


end program q1d

! makefile layout

build all in Libraries

build initialize_grid
build initialize_soln
