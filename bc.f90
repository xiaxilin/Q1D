! Will hold implicit BC's for Q1D Nozzle
! Subsonic inflow with Po and To specified, delta Po and delta To = 0
! Supersonic inflow with any variables specified, they get converted to cons.
! Subsonic outflow with p specified... ie delta p = 0
! Supersonic outflow... linear extrapolation from interior

module bc

  implicit none

  private

  public :: subsonic_inflow
  public :: supersonic_outflow

contains

!=============================================================================80
!
! Uses conserved variables with Po and To specified to calculate subsonic inflow
!
!=============================================================================80
  subroutine subsonic_inflow(cc_in, cc_1, cc_2, DD, DU1, DU2, RHS)

    use set_precision,   only : dp
    use set_constants,   only : zero, half, one, two
    use fluid_constants, only : gamma, gm1, gxgm1, xgm1, R, cv

    implicit none

    real(dp), dimension(3),   intent(in)  :: cc_in, cc_1, cc_2
    real(dp), dimension(3,3), intent(out) :: DD, DU1, DU2
    real(dp), dimension(3),   intent(out) :: RHS

    real(dp) :: rho, u, rhoet, p, T, m, factor
    real(dp) :: dTdrho, dTdrhou, dTdrhoet
    real(dp) :: dPdrho, dPdrhou, dPdrhoet
    real(dp) :: dM2drho, dM2drhou, dM2drhoet

    continue

    DD  = zero
    DU1 = zero
    DU2 = zero

    RHS(1) = zero
    RHS(2) = -cc_in(2)/cc_in(1) + two*cc_1(2)/cc_1(1) - cc_2(2)/cc_2(1)
    RHS(3) = zero

! Extrapolate velocity from interior
    DD(2,1) = -cc_in(2)/cc_in(1)**2
    DD(2,2) = one/cc_in(1)

    DU1(2,1) = -cc_1(2)/cc_1(1)**2
    DU1(2,2) = one/cc_1(1)

    DU1 = -two*DU1

    DU2(2,1) = -cc_2(2)/cc_2(1)**2
    DU2(2,2) = one/cc_2(1)

! Now need to account for delta Po = delta To = 0 at inflow... DD matrix

    rho   = cc_in(1)
    u     = cc_in(2)/cc_in(1)
    rhoet = cc_in(3)
    p     = gm1*(rhoet-half*rho*u**2)
    T     = p/(rho*R)

    m = u/speed_of_sound(p, rho)

    factor = one + half*gm1*m**2

    dTdrho   = -(rhoet/rho**2 - u**2/rho)/cv
    dTdrhou  = -u/(rho*cv)
    dTdrhoet = one/(rho*cv)

    dM2drho   = -m**2*(two/rho + dTdrho/T)
    dM2drhou  = two*u/(gamma*p) - dTdrhou*m**2/T
    dM2drhoet = -dTdrhoet*m**2/T

! add To equations...
    DD(1,1) = factor*dTdrho   + half*gm1*T*dM2drho
    DD(1,2) = factor*dTdrhou  + half*gm1*T*dM2drhou
    DD(1,3) = factor*dTdrhoet + half*gm1*T*dM2drhoet

    dPdrho   = half*gm1*u**2
    dPdrhou  = -gm1*u
    dPdrhoet = gm1

! add Po equations...
    DD(3,1) = factor**gxgm1*dPdrho   + half*gamma*factor**xgm1*P*dM2drho
    DD(3,2) = factor**gxgm1*dPdrhou  + half*gamma*factor**xgm1*P*dM2drhou
    DD(3,3) = factor**gxgm1*dPdrhoet + half*gamma*factor**xgm1*P*dM2drhoet
 
  end subroutine subsonic_inflow

!=============================================================================80
!
! Creates matrices for variable extrapolation
!
!=============================================================================80
  subroutine supersonic_outflow(cc_out, cc_1, cc_2, DD, DL1, DL2, RHS)

    use set_precision, only : dp
    use set_constants, only : zero, one, two

    implicit none

    real(dp), dimension(3),   intent(in)  :: cc_out, cc_1, cc_2
    real(dp), dimension(3,3), intent(out) :: DD, DL1, DL2
    real(dp), dimension(3),   intent(out) :: RHS

    real(dp), dimension(3,3) :: ident3x3

    continue

    ident3x3 = reshape( (/one, zero, zero, zero, one, zero, zero, zero, one/),&
                        (/3,3/) )
    DL2 = ident3x3
    DL1 = -two*ident3x3
    DD  = -ident3x3

    RHS = -cc_2 + two*cc_1 - cc_out

  end subroutine supersonic_outflow

  include 'speed_of_sound.f90'

end module bc
