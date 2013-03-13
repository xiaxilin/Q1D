module adjoint_solvers

  use set_precision, only : dp

  implicit none

  private

  public :: implicit_solve

contains

!=============================== implicit_solve ==============================80
!
! Performs the dual consistent discrete adjoint solve
!
!=============================================================================80
  subroutine implicit_solve( cells, faces, prim_cc, cons_cc, cell_vol,         &
                             area_f, dx, dadx_cc, x_cc )

    use set_precision,   only : dp
    use set_constants,   only : zero, one
    use solvers,         only : iterations, itercheck, iter_out, cfl, cfl_end, &
                                set_time_step, check_convergence
    use adjoint_lhs,     only : fill_full_lhs, transpose_lhs
    use residual,        only : firstorder
    use bc,              only : subsonic_inflow_prim, set_outflow_prim
    use matrix_manip,    only : pentablocksolve
    use write_soln,      only : write_soln_line

    implicit none

    integer,                        intent(in) :: cells, faces
    real(dp), dimension(3,cells+2), intent(in) :: prim_cc, cons_cc
    real(dp), dimension(cells+2),   intent(in) :: cell_vol
    real(dp), dimension(faces),     intent(in) :: area_f
    real(dp), dimension(cells+2),   intent(in) :: dx, dadx_cc, x_cc

    integer  :: n, cell

    real(dp), dimension(3,3)         :: ident3x3
    real(dp), dimension(cells+2)     :: dt
    real(dp), dimension(3,cells+2)   :: psi, dfdq, RHS, delta_psi
    real(dp), dimension(3,3,cells+2) :: rL2, rL, rD, rU, rU2, L2, L, D, U, U2

    logical :: convergence_flag = .false.

    continue

! Initialize some stuff...

    ident3x3 = reshape( [one, zero, zero, zero, one, zero, zero, zero, one],   &
                        [3,3] )

    cfl = cfl_end

    psi = zero
    rL2 = zero
    rL  = zero
    rD  = zero
    rU  = zero
    rU2 = zero

! Get interior of LHS
    call fill_full_lhs( cells, cell_vol, area_f, dadx_cc, prim_cc,             &
                        rL2, rL, rD, rU, rU2 )

! Inflow, modify according to bc
    call subsonic_inflow_prim( firstorder+1,                                   &
                               prim_cc(:,1), prim_cc(:,2), prim_cc(:,3),       &
                               rD(:,:,1), rU(:,:,1), rU2(:,:,1), RHS(:,1) )

! Outflow
    call set_outflow_prim( firstorder+1, prim_cc(:,cells+2),                   &
                           prim_cc(:,cells+1), prim_cc(:,cells),               &
                           rD(:,:,cells+2), rL(:,:,cells+2), rL2(:,:,cells+2), &
                           RHS(:,cells+2) )

    call transpose_lhs( cells, rL2, rL, rD, rU, rU2 )

! Get the functional linearization only once and store it
    call get_dfdq( cells, dx, dfdq )

    main_loop : do n = 0, iterations

! Reinitialize the LHS destroyed by the solver
      L2 = rL2
      L  = rL
      D  = rD
      U  = rU
      U2 = rU2

! Form RHS
      call fill_rhs( cells, psi, rL2, rL, rD, rU, rU2, dfdq, RHS )

! Check residuals for convergence before they are destroyed by the direct solver
      if (mod(n,itercheck) == 0) then
        call check_convergence(cells, n, RHS, convergence_flag)
        if ( convergence_flag ) then
          write(*,*) 'Solution has converged!'
          exit
        end if
      end if

! Time step calculations
      dt = set_time_step( cells, dx, prim_cc )

! Add volume/timestep
      do cell = 2, cells+1
        D(:,:,cell) = D(:,:,cell) + ident3x3*cell_vol(cell)/dt(cell)
      end do

! Solve the system of equations
      call pentablocksolve(3, cells+2, L2, L, D, U, U2, RHS, delta_psi)

! Update the conserved variables
      psi = psi - delta_psi

      if (iter_out >= 0 .and. mod(n,iter_out) == 0) then
        call write_soln_line(n, cells, x_cc, prim_cc, cons_cc)
      end if

!      if (iter_restar >= 0 .and. mod(n,iter_restart) == 0) then
!        call write_restart(cells, prim_cc)
!      end if

    end do main_loop

    if (.not. convergence_flag ) then
      write(*,*) 'Solution failed to converge...'
      write(*,*) 'Consider continuing from q1d.rst'
    end if

    call write_soln_line(iterations, cells, x_cc, psi, cons_cc)

  end subroutine implicit_solve

!================================== get_dfdq =================================80
!
! Sets the linearized functional as the integral of pressure through the nozzle
! FIXME: add other functionals such as entropy
!
!=============================================================================80
  subroutine get_dfdq( cells, dx, dfdq )

    use set_precision,   only : dp
    use set_constants,   only : zero

    implicit none

    integer,                        intent(in)  :: cells
    real(dp), dimension(cells+2),   intent(in)  :: dx
    real(dp), dimension(3,cells+2), intent(out) :: dfdq

    integer  :: cell

    continue

    dfdq = zero

    do cell = 2,cells+1

!      dfdq(1,cell) = zero
!      dfdq(2,cell) = zero
      dfdq(3,cell) = dx(cell)

    end do

  end subroutine get_dfdq

!================================== fill_rhs =================================80
!
! Performs the matrix vector multiplication to fill the RHS and subtracts dfdq
!
!=============================================================================80
  subroutine fill_rhs( cells, psi, L2, L, D, U, U2, dfdq, RHS )

    use set_precision, only : dp

    implicit none

    integer,                          intent(in)  :: cells
    real(dp), dimension(3,cells+2),   intent(in)  :: psi, dfdq
    real(dp), dimension(3,3,cells+2), intent(in)  :: L2, L, D, U, U2
    real(dp), dimension(3,cells+2),   intent(out) :: RHS

    integer :: cell

    continue

    cell = 1
    rhs(:,cell) = matmul(D(:,:,cell),  psi(:,cell))   + &
                  matmul(U(:,:,cell),  psi(:,cell+1)) + &
                  matmul(U2(:,:,cell), psi(:,cell+2)) - dfdq(:,cell)

    cell = 2
    rhs(:,cell) = matmul(L(:,:,cell),  psi(:,cell-1)) + &
                  matmul(D(:,:,cell),  psi(:,cell))   + &
                  matmul(U(:,:,cell),  psi(:,cell+1)) + &
                  matmul(U2(:,:,cell), psi(:,cell+2)) - dfdq(:,cell)

    do cell = 3, cells
      rhs(:,cell) = matmul(L2(:,:,cell), psi(:,cell-2)) + &
                    matmul(L(:,:,cell),  psi(:,cell-1)) + &
                    matmul(D(:,:,cell),  psi(:,cell))   + &
                    matmul(U(:,:,cell),  psi(:,cell+1)) + &
                    matmul(U2(:,:,cell), psi(:,cell+2)) - dfdq(:,cell)
    end do

    cell = cells+1
    rhs(:,cell) = matmul(L2(:,:,cell), psi(:,cell-2)) + &
                  matmul(L(:,:,cell),  psi(:,cell-1)) + &
                  matmul(D(:,:,cell),  psi(:,cell))   + &
                  matmul(U(:,:,cell),  psi(:,cell+1)) - dfdq(:,cell)

    cell = cells+2
    rhs(:,cell) = matmul(L2(:,:,cell), psi(:,cell-2)) + &
                  matmul(L(:,:,cell),  psi(:,cell-1)) + &
                  matmul(D(:,:,cell),  psi(:,cell)) - dfdq(:,cell)

  end subroutine fill_rhs

end module adjoint_solvers
