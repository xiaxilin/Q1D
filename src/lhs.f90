!FIXME: make adding the time and source term its own subroutine?
!FIXME: review the full 2nd order LHS and copy it over to the adjoint
!       may want a face based indexing as well

module lhs

  implicit none

  private

  public :: fill_lhs
  public :: fill_full_lhs
  public :: lhs_order  ! Force the LHS to be 1st or 2nd order

  integer :: lhs_order

contains

!================================== fill_lhs =================================80
!
! Fills the 1st order lhs
!
!=============================================================================80
  subroutine fill_lhs( cells, cell_vol, area_f, dadx_cc, dt, prim_cc, L, D, U )

    use set_precision, only : dp
    use set_constants, only : zero
    use jacobians,     only : jac_source_q_1D, jac_vanleer_q_1D,               &
                              dconserved_dprimitive

    integer,                            intent(in)  :: cells
    real(dp), dimension(cells+1),       intent(in)  :: area_f
    real(dp), dimension(0:cells+1),     intent(in)  :: cell_vol, dadx_cc, dt
    real(dp), dimension(3,0:cells+1),   intent(in)  :: prim_cc
    real(dp), dimension(3,3,0:cells+1), intent(out) :: L, D, U

    integer                  :: cell
    real(dp), dimension(3,3) :: jac_L, jac_R, source_jac, dc_dp

    continue

    L = zero
    D = zero
    U = zero

! Inflow face
    cell = 0
    call jac_vanleer_q_1D( prim_cc(:,cell), prim_cc(:,cell+1), jac_L, jac_R )

    do cell = 1, cells
! subtract contribution from left hand face
      L(:,:,cell) = L(:,:,cell) - jac_L*area_f(cell)
      D(:,:,cell) = D(:,:,cell) - jac_R*area_f(cell)

! now add contribution from right side face
      call jac_vanleer_q_1D( prim_cc(:,cell), prim_cc(:,cell+1), jac_L, jac_R )

      D(:,:,cell) = D(:,:,cell) + jac_L*area_f(cell+1)
      U(:,:,cell) = U(:,:,cell) + jac_R*area_f(cell+1)
    end do

! Add contributions from cell volume/area and the wall pressure/source jacobian
    do cell = 1, cells
      call jac_source_q_1D( dadx_cc(cell), cell_vol(cell), source_jac )
      call dconserved_dprimitive( prim_cc(:,cell), dc_dp )
      D(:,:,cell) = D(:,:,cell) + dc_dp*cell_vol(cell)/dt(cell) - source_jac
    end do

  end subroutine fill_lhs

!=============================== fill_full_lhs ===============================80
!
! Fills the full, 2nd order LHS
!
!=============================================================================80
  subroutine fill_full_lhs( cells, cell_vol, area_f, dadx_cc, dt,             &
                            prim_cc, L2, L, D, U, U2 )

    use set_precision, only : dp
    use set_constants, only : zero, fourth, half, one
    use residual,      only : muscl_extrapolation, firstorder, kappa
    use jacobians,     only : jac_source_q_1D, jac_vanleer_q_1D,               &
                              dconserved_dprimitive

    integer,                            intent(in)  :: cells
    real(dp), dimension(cells+1),       intent(in)  :: area_f
    real(dp), dimension(0:cells+1),     intent(in)  :: cell_vol, dadx_cc, dt
    real(dp), dimension(3,0:cells+1),   intent(in)  :: prim_cc
    real(dp), dimension(3,3,0:cells+1), intent(out) :: L2, L, D, U, U2

    integer                        :: cell

    real(dp), dimension(3,3)       :: jac_L, jac_R, source_jac, dc_dp
    real(dp), dimension(3,cells+1) :: prim_L, prim_R!, limL, limR

    continue

    L2 = zero
    L  = zero
    D  = zero
    U  = zero
    U2 = zero

    call muscl_extrapolation( cells, cells+1, firstorder+1, prim_cc,           &
                              prim_L, prim_R )

! Inflow face
    cell = 0
    call jac_vanleer_q_1D( prim_L(:,cell+1), prim_R(:,cell+1), jac_L, jac_R )

! Subtract from cell to the right

! First the Jacobian on the left side of the face...
! no MUSCL extrapolation for the inflow ghost cell
    L(:,:,cell+1) = L(:,:,cell+1) - jac_L*area_f(cell+1)

! and now the right side.
    L(:,:,cell+1) = L(:,:,cell+1) - jac_R*area_f(cell+1)*fourth*(kappa+one)
    D(:,:,cell+1) = D(:,:,cell+1) - jac_R*area_f(cell+1)*(one-half*kappa)
    U(:,:,cell+1) = U(:,:,cell+1) - jac_R*area_f(cell+1)*fourth*(kappa-one)

! Take care of all interior faces
    do cell = 1, cells-1

    call jac_vanleer_q_1D( prim_L(:,cell+1), prim_R(:,cell+1), jac_L, jac_R )
! Add to cell on the left of the face

! First the Jacobian on the left side of the face...
      L(:,:,cell) = L(:,:,cell) + jac_L*area_f(cell+1)*fourth*(kappa-one)
      D(:,:,cell) = D(:,:,cell) + jac_L*area_f(cell+1)*(one-half*kappa)
      U(:,:,cell) = U(:,:,cell) + jac_L*area_f(cell+1)*fourth*(kappa+one)
! and now the right side.
      D(:,:,cell)  = D(:,:,cell)  + jac_R*area_f(cell+1)*fourth*(kappa+one)
      U(:,:,cell)  = U(:,:,cell)  + jac_R*area_f(cell+1)*(one-half*kappa)
      U2(:,:,cell) = U2(:,:,cell) + jac_R*area_f(cell+1)*fourth*(kappa-one)

! Subtract from cell to the right

! First the Jacobian on the left side of the face...
      L2(:,:,cell+1) = L2(:,:,cell+1) - jac_L*area_f(cell+1)*fourth*(kappa-one)
      L(:,:,cell+1)  = L(:,:,cell+1)  - jac_L*area_f(cell+1)*(one-half*kappa)
      D(:,:,cell+1)  = D(:,:,cell+1)  - jac_L*area_f(cell+1)*fourth*(kappa+one)

! and now the right side.
      L(:,:,cell+1) = L(:,:,cell+1) - jac_R*area_f(cell+1)*fourth*(kappa+one)
      D(:,:,cell+1) = D(:,:,cell+1) - jac_R*area_f(cell+1)*(one-half*kappa)
      U(:,:,cell+1) = U(:,:,cell+1) - jac_R*area_f(cell+1)*fourth*(kappa-one)

    end do

! Outflow face
    cell = cells
    call jac_vanleer_q_1D( prim_L(:,cell+1), prim_R(:,cell+1), jac_L, jac_R )

! Add to cell

! First the Jacobian on the left side of the face...
      L(:,:,cell) = L(:,:,cell) + jac_L*area_f(cell+1)*fourth*(kappa-one)
      D(:,:,cell) = D(:,:,cell) + jac_L*area_f(cell+1)*(one-half*kappa)
      U(:,:,cell) = U(:,:,cell) + jac_L*area_f(cell+1)*fourth*(kappa+one)
! and now the right side.
! no MUSCL extrapolation for the outflow ghost cell
      U(:,:,cell) = U(:,:,cell) + jac_R*area_f(cell+1)

! Add time term and source Jacobian
    do cell = 1, cells
      call jac_source_q_1D( dadx_cc(cell), cell_vol(cell), source_jac )
      call dconserved_dprimitive( prim_cc(:,cell), dc_dp )
      D(:,:,cell) = D(:,:,cell) + dc_dp*cell_vol(cell)/dt(cell) - source_jac
    end do

  end subroutine fill_full_lhs

!============================= modify_lhs_for_bc =============================80
!
! For 1st order LHS, changes three point BC's to two point BC's
!
!=============================================================================80
  subroutine modify_lhs_for_bc(neq, dof, lower, diag, upper, rhs)

    use set_precision, only : dp
    use matrix_manip,  only : mat_inv_3x3!, matrix_inv

    integer,                          intent(in)    :: neq, dof
    real(dp), dimension(neq,neq,dof), intent(inout) :: lower, diag, upper
    real(dp), dimension(neq,dof),     intent(inout) :: rhs

    real(dp), dimension(neq,neq) :: inv

    continue

! Inflow
!    call matrix_inv(3,upper(:,:,2),inv)
    inv = mat_inv_3x3(upper(:,:,2))
    inv = matmul(lower(:,:,1), inv)

    diag(:,:,1)  = diag(:,:,1)  - matmul(inv, lower(:,:,2))
    upper(:,:,1) = upper(:,:,1) - matmul(inv, diag(:,:,2))
    rhs(:,1)     = rhs(:,1)     - matmul(inv, rhs(:,2))

! Outflow
!    call matrix_inv(3,lower(:,:,dof-1),inv)
    inv = mat_inv_3x3(lower(:,:,dof-1))
    inv = matmul(upper(:,:,dof), inv)

    lower(:,:,dof) = lower(:,:,dof) - matmul(inv, diag(:,:,dof-1))
    diag(:,:,dof)  = diag(:,:,dof)  - matmul(inv, upper(:,:,dof-1))
    rhs(:,dof)     = rhs(:,dof)     - matmul(inv, rhs(:,dof-1))

  end subroutine modify_lhs_for_bc

  include 'conserved_to_primitive_1D.f90'
  include 'floor_primitive_vars.f90'
  include 'primitive_to_conserved_1D.f90'

end module lhs
