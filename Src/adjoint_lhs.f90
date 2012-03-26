module adjoint_lhs

  implicit none

  private

  public :: fill_full_lhs
  public :: transpose_lhs

contains

!=============================================================================80
!
!
!
!=============================================================================80
  subroutine fill_full_lhs( cells, cell_vol, area_f, dadx_cc, dt, &
                            cons_cc, L2, L, D, U, U2 )

    use set_precision, only : dp
    use set_constants, only : zero, fourth, half, one
    use residual,      only : muscl_extrapolation, firstorder, kappa
    use jacobians,     only : jac_source_1D, jac_vanleer_1D

    implicit none

    integer,                          intent(in)    :: cells
    real(dp), dimension(cells+2),     intent(in)    :: cell_vol, area_f
    real(dp), dimension(cells+2),     intent(in)    :: dadx_cc, dt
    real(dp), dimension(3,cells+2),   intent(in)    :: cons_cc
    real(dp), dimension(3,3,cells+2), intent(out)   :: L2, L, D, U, U2

    integer                        :: cell
    real(dp), dimension(3,3)       :: ident3x3, jac_L, jac_R, source_jac
    real(dp), dimension(3,cells+2) :: prim_cc, prim_L, prim_R, cons_L, cons_R

    continue

    ident3x3 = reshape( (/one, zero, zero, zero, one, zero, zero, zero, one/) ,&
                        (/3,3/) )

    L2 = zero
    L  = zero
    D  = zero
    U  = zero
    U2 = zero

    do cell = 1, cells+1
      prim_cc(:,cell) = conserved_to_primitive_1D(cons_cc(:,cell))
    end do

    call muscl_extrapolation( cells, cells+1, firstorder+1, &
                              prim_cc, prim_L, prim_R )

    do cell = 1, cells+1
      prim_L(:,cell) = floor_primitive_vars(prim_L(:,cell))
      cons_L(:,cell) = primitive_to_conserved_1D(prim_L(:,cell))
      prim_R(:,cell) = floor_primitive_vars(prim_R(:,cell))
      cons_R(:,cell) = primitive_to_conserved_1D(prim_R(:,cell))
    end do

! Inflow face
    cell = 1
    call jac_vanleer_1D( cons_L(:,cell), cons_R(:,cell), jac_L, jac_R )

! Subtract from cell to the right

! First the Jacobian on the left side of the face...
    L(:,:,cell+1)  = L(:,:,cell+1) - jac_L*area_f(cell)*(one-fourth*(kappa+one))
    D(:,:,cell+1)  = D(:,:,cell+1) - jac_L*area_f(cell)*fourth*(kappa+one)

! and now the right side.
    L(:,:,cell+1) = L(:,:,cell+1) - jac_R*area_f(cell)*fourth*(kappa+one)
    D(:,:,cell+1) = D(:,:,cell+1) - jac_R*area_f(cell)*(one-half*kappa)
    U(:,:,cell+1) = U(:,:,cell+1) - jac_R*area_f(cell)*fourth*(kappa-one)

! Take care of all interior faces
    do cell = 2, cells

      call jac_vanleer_1D( cons_L(:,cell), cons_R(:,cell), jac_L, jac_R )

! Add to cell

! First the Jacobian on the left side of the face...
      L(:,:,cell) = L(:,:,cell) + jac_L*area_f(cell)*fourth*(kappa-one)
      D(:,:,cell) = D(:,:,cell) + jac_L*area_f(cell)*(one-half*kappa)
      U(:,:,cell) = U(:,:,cell) + jac_L*area_f(cell)*fourth*(kappa+one)
! and now the right side.
      D(:,:,cell)  = D(:,:,cell)  + jac_R*area_f(cell)*fourth*(kappa+one)
      U(:,:,cell)  = U(:,:,cell)  + jac_R*area_f(cell)*(one-half*kappa)
      U2(:,:,cell) = U2(:,:,cell) + jac_R*area_f(cell)*fourth*(kappa-one)

! Subtract from cell to the right

! First the Jacobian on the left side of the face...
      L2(:,:,cell+1) = L2(:,:,cell+1) - jac_L*area_f(cell)*fourth*(kappa-one)
      L(:,:,cell+1)  = L(:,:,cell+1)  - jac_L*area_f(cell)*(one-half*kappa)
      D(:,:,cell+1)  = D(:,:,cell+1)  - jac_L*area_f(cell)*fourth*(kappa+one)

! and now the right side.
      L(:,:,cell+1) = L(:,:,cell+1) - jac_R*area_f(cell)*fourth*(kappa+one)
      D(:,:,cell+1) = D(:,:,cell+1) - jac_R*area_f(cell)*(one-half*kappa)
      U(:,:,cell+1) = U(:,:,cell+1) - jac_R*area_f(cell)*fourth*(kappa-one)

    end do

! Outflow face
    cell = cells+1
    call jac_vanleer_1D( cons_L(:,cell), cons_R(:,cell), jac_L, jac_R )

! Add to cell

! First the Jacobian on the left side of the face...
      L(:,:,cell) = L(:,:,cell) + jac_L*area_f(cell)*fourth*(kappa-one)
      D(:,:,cell) = D(:,:,cell) + jac_L*area_f(cell)*(one-half*kappa)
      U(:,:,cell) = U(:,:,cell) + jac_L*area_f(cell)*fourth*(kappa+one)
! and now the right side.
      D(:,:,cell)  = D(:,:,cell)  + jac_R*area_f(cell)*fourth*(kappa+one)
      U(:,:,cell)  = U(:,:,cell)  + jac_R*area_f(cell)*(one-fourth*kappa)

! Add time term and source Jacobian
    do cell = 2, cells+1
      call jac_source_1D( cons_cc(2,cell)/cons_cc(1,cell), dadx_cc(cell),      &
                          cell_vol(cell), source_jac )

      D(:,:,cell) = D(:,:,cell) - source_jac

    end do

  end subroutine fill_full_lhs

!=============================================================================80
!
!
!
!=============================================================================80
  subroutine transpose_lhs(cells, L2, L, D, U, U2)

    use set_precision, only : dp

    implicit none

    integer,                          intent(in)    :: cells
    real(dp), dimension(3,3,cells+2), intent(inout) :: L2, L, D, U, U2

    integer :: cell

    continue

    do cell = 1, cells+2
      L2(:,:,cell) = transpose(L2(:,:,cell))
    end do

    do cell = 1, cells+2
      L(:,:,cell) = transpose(L(:,:,cell))
    end do

    do cell = 1, cells+2
      D(:,:,cell) = transpose(D(:,:,cell))
    end do

    do cell = 1, cells+2
      U(:,:,cell) = transpose(U(:,:,cell))
    end do

    do cell = 1, cells+2
      U2(:,:,cell) = transpose(U2(:,:,cell))
    end do

  end subroutine transpose_lhs

  include 'conserved_to_primitive_1D.f90'
  include 'floor_primitive_vars.f90'
  include 'primitive_to_conserved_1D.f90'

end module adjoint_lhs
