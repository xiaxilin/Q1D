module adjoint_lhs

  implicit none

  private

  public :: fill_full_lhs
  public :: transpose_lhs

contains

!=============================== fill_full_lhs ===============================80
!
! Fills the lhs (full 2nd order linearized residual)
! FIXME: change routine so that the version on lhs.f90 can be reused here
!
!=============================================================================80
  subroutine fill_full_lhs( cells, cell_vol, area_f, dadx_cc, prim_cc,         &
                            L2, L, D, U, U2 )

    use set_precision, only : dp
    use set_constants, only : zero, fourth, third, half, one, onep5, two,      &
                              three, four, six
    use residual,      only : muscl_extrapolation, firstorder, kappa,          &
                              inflow_face, outflow_face
    use jacobians,     only : jac_source_q_1D, jac_vanleer_q_1D

    integer,                            intent(in)  :: cells
    real(dp), dimension(cells+1),       intent(in)  :: area_f
    real(dp), dimension(0:cells+1),     intent(in)  :: cell_vol, dadx_cc
    real(dp), dimension(3,0:cells+1),   intent(in)  :: prim_cc
    real(dp), dimension(3,3,0:cells+1), intent(out) :: L2, L, D, U, U2

    integer :: cell

    real(dp), dimension(3,3)       :: jac_L, jac_R, source_jac
    real(dp), dimension(3,cells+1) :: prim_L, prim_R

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
    select case( inflow_face )
    case( 0 ) ! zeroth order extrapolation
      L(:,:,cell+1) = L(:,:,cell+1) - jac_L*area_f(cell+1)
    case ( 1 ) ! zero gradient extrapolation to ghost cell
      D(:,:,cell+1) = D(:,:,cell+1) - jac_L*area_f(cell+1)*third*four
      U(:,:,cell+1) = U(:,:,cell+1) + jac_L*area_f(cell+1)*third
    case ( 2 ) ! zero curvature extrapolation to ghost cell
      D(:,:,cell+1) = D(:,:,cell+1) - jac_L*area_f(cell+1)*two
      U(:,:,cell+1) = U(:,:,cell+1) + jac_L*area_f(cell+1)
    case ( 3 ) ! third order extrapolation to ghost cell
      D(:,:,cell+1)  = D(:,:,cell+1)  - jac_L*area_f(cell+1)*three
      U(:,:,cell+1)  = U(:,:,cell+1)  + jac_L*area_f(cell+1)*three
      U2(:,:,cell+1) = U2(:,:,cell+1) - jac_L*area_f(cell+1)
    case ( -1 ) ! zero gradient extrapolation to face
      D(:,:,cell+1) = D(:,:,cell+1) - jac_L*area_f(cell+1)*7.0_dp/six
      U(:,:,cell+1) = U(:,:,cell+1) + jac_L*area_f(cell+1)/six
    case ( -2 ) ! zero curvature extrapolation to face
      D(:,:,cell+1) = D(:,:,cell+1) - jac_L*area_f(cell+1)*onep5
      U(:,:,cell+1) = U(:,:,cell+1) + jac_L*area_f(cell+1)*half
    case ( -3 ) ! third order extrapolation to face
      D(:,:,cell+1)  = D(:,:,cell+1)  - jac_L*area_f(cell+1)*11.0_dp/six
      U(:,:,cell+1)  = U(:,:,cell+1)  + jac_L*area_f(cell+1)*7.0_dp/six
      U2(:,:,cell+1) = U2(:,:,cell+1) - jac_L*area_f(cell+1)*third
    case default ! truncated MUSCL
      L(:,:,cell+1) = L(:,:,cell+1)                                            &
                    - jac_L*area_f(cell+1)*(one+fourth*(kappa+one))
      D(:,:,cell+1) = D(:,:,cell+1) - jac_L*area_f(cell+1)*fourth*(kappa+one)
    end select

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

! Add source Jacobian
    do cell = 1, cells
      call jac_source_q_1D( dadx_cc(cell), cell_vol(cell), source_jac )
      D(:,:,cell) = D(:,:,cell) - source_jac
    end do

  end subroutine fill_full_lhs

!=============================== transpose_lhs ===============================80
!
! Transposes the lhs (full 2nd order linearized residual)
!
!=============================================================================80
  subroutine transpose_lhs(cells, L2, L, D, U, U2)

    use set_precision, only : dp
    use set_constants, only : zero

    integer,                            intent(in)    :: cells
    real(dp), dimension(3,3,0:cells+1), intent(inout) :: L2, L, D, U, U2

    integer :: cell

    real(dp), dimension(3,3,0:cells+1) :: temp

    continue

    do cell = 0, cells+1
      L2(:,:,cell) = transpose(L2(:,:,cell))
    end do

    do cell = 0, cells+1
      L(:,:,cell) = transpose(L(:,:,cell))
    end do

    do cell = 0, cells+1
      D(:,:,cell) = transpose(D(:,:,cell))
    end do

    do cell = 0, cells+1
      U(:,:,cell) = transpose(U(:,:,cell))
    end do

    do cell = 0, cells+1
      U2(:,:,cell) = transpose(U2(:,:,cell))
    end do

! Now that they are transposed in place, transpose across the diagonal
    temp = L2
    L2(:,:,0) = zero
    L2(:,:,1) = zero
    do cell = 0, cells-1
      L2(:,:,cell+2) = U2(:,:,cell)
    end do
    do cell = 0, cells-1
      U2(:,:,cell) = temp(:,:,cell+2)
    end do

    temp = L
    L(:,:,0) = zero
    do cell = 0, cells
      L(:,:,cell+1) = U(:,:,cell)
    end do
    do cell = 0, cells
      U(:,:,cell) = temp(:,:,cell+1)
    end do

  end subroutine transpose_lhs

end module adjoint_lhs
