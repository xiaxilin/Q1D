! Fix the Jacobian matrix so that it does the full transpose, ie,
! map rL2 to rU2 and rL to rU and vice versa.

module adjoint_solvers

  use set_precision, only : dp

  implicit none

  private

  public :: implicit_solve
  public :: check_convergence

  contains
!=============================================================================80
!
!
!
!=============================================================================80
  subroutine implicit_solve( cells, faces, prim_cc, cons_cc, cell_vol,         &
                             area_f, dx, dadx_cc, x_cc )

    use set_precision,   only : dp
    use set_constants,   only : zero, half, one, two
    use fluid_constants, only : gm1
    use solvers,         only : iterations, itercheck, iter_out, cfl, cfl_end
    use adjoint_lhs,     only : fill_full_lhs, transpose_lhs
    use residual,        only : firstorder
    use bc,              only : subsonic_inflow, set_outflow
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

    ident3x3 = reshape( (/one, zero, zero, zero, one, zero, zero, zero, one/) ,&
                        (/3,3/) )
    psi = zero
    rL2 = zero
    rL  = zero
    rD  = zero
    rU  = zero
    rU2 = zero

    call get_dfdq( cells, dadx_cc, cell_vol, prim_cc, dfdq )

    call fill_full_lhs( cells, cell_vol, area_f, dadx_cc, dt,                  &
                        cons_cc, rL2, rL, rD, rU, rU2 )

! Inflow, modify according to bc
    call subsonic_inflow( firstorder+1,                                        &
                          cons_cc(:,1), cons_cc(:,2), cons_cc(:,3),            &
                          rD(:,:,1), rU(:,:,1), rU2(:,:,1), RHS(:,1) )

! Outflow
    call set_outflow( firstorder+1,                                            &
                      cons_cc(:,cells+2), cons_cc(:,cells+1), cons_cc(:,cells),&
                      rD(:,:,cells+2), rL(:,:,cells+2), rL2(:,:,cells+2),      &
                      RHS(:,cells+2) )

    call transpose_lhs( cells, rL2, rL, rD, rU, rU2 )

    cfl = cfl_end

    main_loop : do n = 0, iterations

      L2 = rL2
      L  = rL
      D  = rD
      U  = rU
      U2 = rU2

! Time step calculations
      dt = set_time_step( cells, dx, prim_cc )

! Form RHS
      call fill_rhs( cells, psi, rL2, rL, rD, rU, rU2, RHS )

! Check residuals for convergence before they are destroyed by the direct solver
      if (n /=0 .and. mod(n,itercheck) == 0) then
        call check_convergence(cells, n, RHS, convergence_flag)
        if ( convergence_flag ) then
          write(*,*) 'Solution has converged!'
          exit
        end if
      end if

! Add volume/timestep
      do cell = 2, cells+2
        D(:,:,cell) = D(:,:,cell) + ident3x3*cell_vol(cell)/dt(cell)
      end do

      RHS = RHS - dfdq

! solve the system of equations
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

!=============================================================================80
!
!
!
!=============================================================================80
  pure function set_time_step( cells, dx, prim_cc ) result(dt)

    use set_precision, only : dp
    use set_constants, only : large
    use solvers,       only : cfl

    implicit none

    integer,                         intent(in)  :: cells
    real(dp), dimension(cells+2),    intent(in)  :: dx
    real(dp), dimension(3, cells+2), intent(in)  :: prim_cc
    real(dp), dimension(cells+2)                 :: dt

    integer  :: cell
    real(dp) :: a, dt_global

    continue

    dt_global = large

    do cell = 1, cells+2
      a = speed_of_sound( prim_cc(3,cell), prim_cc(1,cell) )
      dt(cell) = cfl*dx(cell) / ( abs(prim_cc(2,cell)) + a )
      dt_global = min(dt_global, dt(cell))
    end do

  end function set_time_step

!=============================================================================80
!
!
!
!=============================================================================80
  subroutine get_dfdq( cells, dadx_cc, cell_jac, prim_cc, dfdq )

    use set_precision,   only : dp
    use set_constants,   only : zero, half
    use fluid_constants, only : gm1

    implicit none

    integer,                        intent(in)  :: cells
    real(dp), dimension(cells+2),   intent(in)  :: dadx_cc, cell_jac
    real(dp), dimension(3,cells+2), intent(in)  :: prim_cc
    real(dp), dimension(3,cells+2), intent(out) :: dfdq

    integer  :: cell
    real(dp) :: sidewall_area

    continue

    dfdq = zero

    do cell = 2,cells+1

      sidewall_area = dadx_cc(cell)*cell_jac(cell)

      dfdq(1,cell) = half*sidewall_area*gm1*prim_cc(2,cell)**2
      dfdq(2,cell) = -gm1*prim_cc(2,cell)*sidewall_area
      dfdq(3,cell) = gm1*sidewall_area

    end do

  end subroutine get_dfdq

!=============================================================================80
!
!
!
!=============================================================================80
  subroutine fill_rhs( cells, psi, L2, L, D, U, U2, RHS )

    use set_precision, only : dp

    implicit none

    integer,                          intent(in)  :: cells
    real(dp), dimension(3,cells+2),   intent(in)  :: psi
    real(dp), dimension(3,3,cells+2), intent(in)  :: L2, L, D, U, U2
    real(dp), dimension(3,cells+2),   intent(out) :: RHS

    integer :: cell

    continue

    cell = 1
    rhs(:,cell) = matmul(D(:,:,cell),  psi(:,cell))   + &
                  matmul(U(:,:,cell),  psi(:,cell+1)) + &
                  matmul(U2(:,:,cell), psi(:,cell+2))

    cell = 2
    rhs(:,cell) = matmul(L(:,:,cell),  psi(:,cell-1)) + &
                  matmul(D(:,:,cell),  psi(:,cell))   + &
                  matmul(U(:,:,cell),  psi(:,cell+1)) + &
                  matmul(U2(:,:,cell), psi(:,cell+2))

    do cell = 3, cells
      rhs(:,cell) = matmul(L2(:,:,cell), psi(:,cell-2)) + &
                    matmul(L(:,:,cell),  psi(:,cell-1)) + &
                    matmul(D(:,:,cell),  psi(:,cell))   + &
                    matmul(U(:,:,cell),  psi(:,cell+1)) + &
                    matmul(U2(:,:,cell), psi(:,cell+2))
    end do

    cell = cells+1
    rhs(:,cell) = matmul(L2(:,:,cell), psi(:,cell-2)) + &
                  matmul(L(:,:,cell),  psi(:,cell-1)) + &
                  matmul(D(:,:,cell),  psi(:,cell))   + &
                  matmul(U(:,:,cell),  psi(:,cell+1))

    cell = cells+2
    rhs(:,cell) = matmul(L2(:,:,cell), psi(:,cell-2)) + &
                  matmul(L(:,:,cell),  psi(:,cell-1)) + &
                  matmul(D(:,:,cell),  psi(:,cell))

  end subroutine fill_rhs

!=============================================================================80
!
!
!
!=============================================================================80
  subroutine check_convergence(cells, iteration, residuals, convergence_flag)

    use set_precision,   only : dp
    use set_constants,   only : zero
    use solvers,         only : toler
    use initialize_soln, only : restart, L1_init, L2_init, Linf_init

    implicit none

    integer,                        intent(in)  :: cells, iteration
    real(dp), dimension(3,cells+2), intent(in)  :: residuals
    logical,                        intent(out) :: convergence_flag

    integer :: cell
    real(dp), dimension(3) :: L1, L2, Linf

    continue

    L1 = zero
    L2 = zero
    Linf = zero

    do cell = 2, cells+1
      L1(:) = L1(:) + abs(residuals(:,cell))
      L2(:) = L2(:) + residuals(:,cell)**2
      Linf(:) = max(Linf(:), abs(residuals(:,cell)))
    end do

    L1(:) = L1(:)/real(cells,dp)
    L2(:) = sqrt(L2(:))/real(cells,dp)

    if (iteration == 0 .and. .not. restart) then
      L1_init(:)   = L1(:)
      L2_init(:)   = L2(:)
      Linf_init(:) = Linf(:)
    end if

    L1(:)   = L1(:)/L1_init(:)
    L2(:)   = L2(:)/L2_init(:)
    Linf(:) = Linf(:)/Linf_init(:)

    write(*,300) iteration, L2(1), L2(2), L2(3)
300 format(1X,i8,2(e15.6),3(e15.6),4(e15.6))

    convergence_flag = .false.
    if (L2(1) <= toler .and. L2(2) <= toler .and. L2(3) <= toler ) then
      convergence_flag = .true.
    end if

  end subroutine check_convergence

! begin include statements

  include 'conserved_to_primitive_1D.f90'
  include 'floor_primitive_vars.f90'
  include 'primitive_to_conserved_1D.f90'
  include 'speed_of_sound.f90'

end module adjoint_solvers
