!=============================================================================80
!
!
!
!=============================================================================80

subroutine implicit_solve( cells, faces, dxsi, prim_cc, cons_cc,               &
                           area_f, area_cc, dxdxsi_cc, dadx_cc )

  use set_precision, only : dp
  use set_constants, only : two
  use bc,            only : subsonic_inflow, supersonic_outflow

  implicit none

  integer,                        intent(in)    :: cells, faces
  real(dp),                       intent(in)    :: dxsi
  real(dp), dimension(3,cells+2), intent(inout) :: prim_cc, cons_cc
  real(dp), dimension(faces),     intent(in)    :: area_f
  real(dp), dimension(cells+2),   intent(in)    :: area_cc, dxdxsi_cc, dadx_cc

  integer  :: n, cell
  real(dp) :: dt_global

  real(dp), dimension(cells+2)     :: dt
  real(dp), dimension(3,cells+2)   :: RHS, delta_cons_cc
  real(dp), dimension(3,3)         :: left_jac_L, right_jac_L
  real(dp), dimension(3,3)         :: left_jac_C, right_jac_C
  real(dp), dimension(3,3)         :: left_jac_R, right_jac_R
  real(dp), dimension(3,3,cells+2) :: L, D, U

  logical :: convergence_flag = .false.
  real(dp), dimension(3,3) :: ident3x3

  continue

  ident3x3 = reshape( (/one, zero, zero, zero, one, zero, zero, zero, one/) ,  &
                      (/3,3/) )

  do n = 0, iterations

    call set_time_step( cells, dxsi, prim_cc, dt, dt_global )

! form RHS
    call create_residual( cells, faces, n, dxsi, prim_cc, cons_cc,             &
                          area_f, dadx_cc, dxdxsi_cc, RHS)

! form LHS

! Inflow, modify according to bc
    L(:,:,1) = zero
    call subsonic_inflow( cons_cc(:,1), cons_cc(:,2), cons_cc(:,3),            &
                          D(:,:,1), U(:,:,1), DU2(:,:), RHS(:,1)

! calculate Jacobians for 1st ghost cell and 1st interior cell

    call jacobian_vanleer_1D( cons_cc(:,1), cons_cc(:,1),                      &
                              right_jac_L, left_jac_L)

    call jacobian_vanleer_1D( cons_cc(:,2), cons_cc(:,2),                      &
                              right_jac_C, left_jac_C)

    do cell = 2, cells+1

      call jacobian_vanleer_1D( cons_cc(:,cell+1), cons_cc(:,cell+1),          &
                                right_jac_R, left_jac_R )

      L(:,:,cell) = -right_jac_L/(dxdsi_cc(cell-1)*dxsi)
      D(:,:,cell) = ident/dt - source_jac                                      &
                  + (right_jac_C-left_jac_C)/(dxdsi_cc(cell)*dxsi)
      U(:,:,cell) =  left_jac_R/(dxdsi_cc(cell+1)*dxsi)

! shift Jacobians to avoid recalculation

      left_jac_L  = left_jac_C
      right_jac_L = right_jac_C

      left_jac_C  = left_jac_R
      right_jac_C = right_jac_R
    end do

! Outflow, modify according to bc
    U(:,:,cells+2) = zero
    call supersonic_outflow( cons_cc(:,cells+2), cons_cc(:,cells+1),           &
                             cons_cc(:,cells),                                 &
                             D(:,:,cells+2), L(:,:,cells+2), DL2(:,:),         &
                             RHS(:,cells+2)

! Need to perform matrix modification to preserve block tri-diagonal structure

!    bc = zero
!    bc(1,1) = one

!    call matrix_inv(3, U(:,:,2), inv)
!    call mat_inv_3x3(U(:,:,2), inv)

!    inv = matmul(bc, inv)

!    D(:,:,1) = D(:,:,1) - matmul(inv, L(:,:,2))
!    U(:,:,1) = U(:,:,1) - matmul(inv, D(:,:,2))
!    RHS(:,1) = RHS(:,1) - matmul(inv, RHS(:,2))

!    call matrix_inv(3, L(:,:,y_nodes-1), inv)
!    call mat_inv_3x3(L(:,:,y_nodes-1), inv)

!    inv = matmul(bc, inv)

!    L(:,:,y_nodes) = L(:,:,y_nodes) - matmul(inv, D(:,:,y_nodes-1))
!    D(:,:,y_nodes) = D(:,:,y_nodes) - matmul(inv, U(:,:,y_nodes-1))
!    RHS(:,y_nodes) = RHS(:,y_nodes) - matmul(inv, RHS(:,y_nodes-1))

! solve the system of equations
    call triblocksolve(3, cells+2, L, D, U, RHS, delta_cons_cc)

    cons_cc = cons_cc+deltat_cons_cc

    if (mod(n,itercheck) == 0) then
      call check_convergence(cells, n, RHS, convergence_flag)
    end if

!    if (iter_out >= 0 .and. mod(n,iter_out) == 0) then
!      call write_soln(cells, faces, prim_cc, cons_cc)
!    end if

!    if (iter_restar >= 0 .and. mod(n,iter_restart) == 0) then
!      call write_restart(cells, prim_cc)
!    end if

    if ( convergence_flag ) then
      write(*,*) 'Solution has converged!'
      exit
    end if
  end do

  if (.not. convergence_flag ) then
    write(*,*) 'Solution failed to converge...'
    write(*,*) 'Consider continuing from q1d.rst'
  end if

end subroutine implicit_solve
