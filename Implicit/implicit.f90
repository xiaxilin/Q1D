! sketch out implicit logic... how to build problem and such

! 

subroutine implicit_solve(cells, faces, ... etc)

  use set_precision, only : dp

  implicit none



  continue

! form RHS
  call create_residual(RHS)

! form LHS

! calculate Jacobians for 1st ghost cell and 1st interior cell

  call jacobian_vanleer_1D(cons_cc(:,1), cons_cc(:,1), right_jac_L, left_jac_L)

  call jacobian_vanleer_1D(cons_cc(:,2), cons_cc(:,2), right_jac_C, left_jac_C)

  do cell = 2, cells+1

    call jacobian_vanleer_1D( cons_cc(:,cell+1), cons_cc(:,cell+1),            &
                               right_jac_R, left_jac_R )

    L(:,:,cell) = -right_jac_L/(dxdsi_cc(cell-1)*dxsi)
    D(:,:,cell) = ident/dt - source_jac                                        &
                + (right_jac_C-left_jac_C)/(dxdsi_cc(cell)*dxsi)
    U(:,:,cell) =  left_jac_R/(dxdsi_cc(cell+1)*dxsi)

! shift Jacobians to avoid recalculation

    left_jac_L  = left_jac_C
    right_jac_L = right_jac_C

    left_jac_C  = left_jac_R
    right_jac_C = right_jac_R
  end do

  call triblocksolve(3, cells+2, L, D, U, RHS, cons_cc)

end subroutine implicit_solve
