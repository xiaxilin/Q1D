!=============================================================================80
!
!
!
!=============================================================================80
module implicit

  implicit none

  private

  public :: implicit_solve

contains

  subroutine implicit_solve( cells, faces, dxsi, prim_cc, cons_cc,             &
                             area_f, area_cc, dxdxsi_cc, dadx_cc, x_cc )

    use set_precision,   only : dp
    use set_constants,   only : zero, half, one, two
    use fluid_constants, only : gm1
    use solvers,         only : iterations, itercheck, create_residual, &
                                check_convergence, set_time_step
    use bc,              only : subsonic_inflow, supersonic_outflow
    use jacobians,       only : jac_vanleer_1D
    use matrix_manip,    only : triblocksolve, mat_inv_3x3, matrix_inv
    use write_soln,      only : write_soln_line

    implicit none

    integer,                        intent(in)    :: cells, faces
    real(dp),                       intent(in)    :: dxsi
    real(dp), dimension(3,cells+2), intent(inout) :: prim_cc, cons_cc
    real(dp), dimension(faces),     intent(in)    :: area_f
    real(dp), dimension(cells+2),   intent(in)    :: area_cc, dxdxsi_cc, dadx_cc
    real(dp), dimension(cells+2),   intent(in)    :: x_cc

    integer  :: n, cell

    real(dp)                         :: divisor
    real(dp), dimension(cells+2)     :: dt
    real(dp), dimension(3,cells+2)   :: RHS, delta_cons_cc
    real(dp), dimension(3,3)         :: inv, source_jac
    real(dp), dimension(3,3)         :: left_jac_L, right_jac_L
    real(dp), dimension(3,3)         :: left_jac_C, right_jac_C
    real(dp), dimension(3,3)         :: left_jac_R, right_jac_R
    real(dp), dimension(3,3,cells+2) :: L, D, U

    logical :: convergence_flag = .false.
    real(dp), dimension(3,3) :: ident3x3

    continue

    ident3x3 = reshape( (/one, zero, zero, zero, one, zero, zero, zero, one/) ,&
                        (/3,3/) )
      print*, cons_cc(:,1)
      print*, cons_cc(:,2)
      print*, cons_cc(:,3)
      print*, cons_cc(:,4)
      print*, cons_cc(:,5)
      print*, cons_cc(:,6)
    do n = 0, iterations

      dt = set_time_step( cells, dxsi, prim_cc )

! Form RHS
      call create_residual( cells, faces, n, dxsi, prim_cc, cons_cc,           &
                            area_f, dadx_cc, dxdxsi_cc, RHS)

! Now form LHS.... should be a subroutine

! Inflow, modify according to bc
      L(:,:,1) = zero
      call subsonic_inflow( cons_cc(:,1), cons_cc(:,2), cons_cc(:,3),          &
                            D(:,:,1), U(:,:,1), L(:,:,1), RHS(:,1) )

! calculate Jacobians for inflow face

      call jac_vanleer_1D( cons_cc(:,1), cons_cc(:,2), right_jac_L, left_jac_C)

      do cell = 2, cells+1

        call jac_vanleer_1D( cons_cc(:,cell), cons_cc(:,cell+1),               &
                                  right_jac_C, left_jac_R )

        source_jac = zero
        source_jac(2,1) = half*gm1*prim_cc(2,cell)**2
        source_jac(2,2) = -gm1*prim_cc(2,cell)
        source_jac(2,3) = gm1
        source_jac = source_jac*dadx_cc(cell)

        divisor = dxdxsi_cc(cell)*dxsi*area_cc(cell)

        L(:,:,cell) = -right_jac_L*area_f(cell-1) / divisor
        D(:,:,cell) = ident3x3/dt(cell)                                        &
                    + ( right_jac_C*area_f(cell)-left_jac_C*area_f(cell-1) )   &
                    / divisor - source_jac / area_cc(cell)
        U(:,:,cell) =  left_jac_R*area_f(cell) / divisor

! shift Jacobians to avoid recalculation

        right_jac_L = right_jac_C
        left_jac_C  = left_jac_R

      end do

! Outflow, modify according to bc
      U(:,:,cells+2) = zero
      call supersonic_outflow( cons_cc(:,cells+2), cons_cc(:,cells+1),         &
                               cons_cc(:,cells),                               &
                               D(:,:,cells+2), L(:,:,cells+2), U(:,:,cells+2), &
                               RHS(:,cells+2) )

! Modify matrix to maintain block tridiagonal structure
! Needs to be made a subroutine...

! Inflow

!      call matrix_inv(3, U(:,:,2), inv)

      call mat_inv_3x3(U(:,:,2), inv)

      inv = matmul(L(:,:,1), inv)

      D(:,:,1) = D(:,:,1) - matmul(inv, L(:,:,2))
      U(:,:,1) = U(:,:,1) - matmul(inv, D(:,:,2))
      RHS(:,1) = RHS(:,1) - matmul(inv, RHS(:,2))

! Outflow
!      call matrix_inv(3, L(:,:,cells+1), inv)
      call mat_inv_3x3(L(:,:,cells+1), inv)

      inv = matmul(U(:,:,cells+2), inv)

      L(:,:,cells+2) = L(:,:,cells+2) - matmul(inv, D(:,:,cells+1))
      D(:,:,cells+2) = D(:,:,cells+2) - matmul(inv, U(:,:,cells+1))
      RHS(:,cells+2) = RHS(:,cells+2) - matmul(inv, RHS(:,cells+1))

! solve the system of equations
      call triblocksolve(3, cells+2, L, D, U, RHS, delta_cons_cc)

! Update the conserved variables
      cons_cc = cons_cc+delta_cons_cc

      print*, cons_cc(:,1)
      print*, cons_cc(:,2)
      print*, cons_cc(:,3)
      print*, cons_cc(:,cells)
      print*, cons_cc(:,cells+1)
      print*, cons_cc(:,cells+2)

      if (mod(n,itercheck) == 0) then
        call check_convergence(cells, n, RHS, convergence_flag)
      end if

!      if (iter_out >= 0 .and. mod(n,iter_out) == 0) then
!        call write_soln_line(n, cells, x_cc, prim_cc, cons_cc)
!      end if

!      if (iter_restar >= 0 .and. mod(n,iter_restart) == 0) then
!        call write_restart(cells, prim_cc)
!      end if

      if ( convergence_flag ) then
        write(*,*) 'Solution has converged!'
        exit
      end if
    end do

    if (.not. convergence_flag ) then
      write(*,*) 'Solution failed to converge...'
      write(*,*) 'Consider continuing from q1d.rst'
    end if

    call write_soln_line(iterations, cells, x_cc, prim_cc, cons_cc)

  end subroutine implicit_solve

end module implicit
