! TODO: add soln restart output and visualization
! TODO: add convergence checks

program q1d_primal

  use set_precision,   only : dp
  use namelist,        only : read_nml
  use fluid_constants, only : set_gamma_constants
  use initialize_grid, only : read_grid, cells, faces, dx, x_cc,               &
                              area_f, area_cc, dadx_cc, deallocate_grid
  use initialize_soln, only : allocate_soln, initial_soln, prim_cc, cons_cc,   &
                              deallocate_soln
  use solvers,         only : explicit_solve, implicit_solve, solver
  use write_soln,      only : write_restart, init_write_files, write_entropy
  use solution_error,  only : calculate_exact_soln

  implicit none

  continue

  print *,'*****************************************************'
  print *,'         Quasi-1D Nozzle : RELEASE 5.0, 2012         '
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

  select case( solver )
  case('explicit')
    print*, 'Beginning Explicit Solve'
    call explicit_solve(cells, faces, prim_cc, cons_cc,                        &
                        area_f, area_cc, dx, dadx_cc, x_cc)
  case('implicit')
    print*, 'Beginning Implicit Solve'
    call implicit_solve(cells, faces, prim_cc, cons_cc,                        &
                        area_f, area_cc, dx, dadx_cc, x_cc)
  end select

! do solution output
  call write_restart(cells, prim_cc)

! write the entropy vars
  call write_entropy(cells, x_cc, prim_cc, cons_cc)

! get exact solution and plot it

  print *, 'Calculating Exact Solution'
  call calculate_exact_soln(cells, x_cc, area_cc, 0.2_dp,                      &
                            area_f(cells+1), cons_cc)

! free memory

  print *, 'Deallocating Memory'
  call deallocate_grid
  call deallocate_soln

  print *, 'Finished'

end program q1d_primal
