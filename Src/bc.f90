! Will hold implicit BC's for Q1D Nozzle
! Subsonic inflow with Po and To specified, delta Po and delta To = 0
! Supersonic inflow with any variables specified, they get converted to cons.
! Subsonic outflow with p specified... ie delta p = 0
! Supersonic outflow... linear extrapolation from interior

module bc

  implicit none

  private

  public :: subsonic_inflow
  public :: set_outflow

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
    use initialize_soln, only : po, to

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

! Extrapolate velocity from interior
    DD(1,1) = -cc_in(2)/cc_in(1)**2
    DD(1,2) = one/cc_in(1)

    DU1(1,1) = -cc_1(2)/cc_1(1)**2
    DU1(1,2) = one/cc_1(1)

    DU1 = two*DU1
!    DU1 = -DU1

    DU2(1,1) = -cc_2(2)/cc_2(1)**2
    DU2(1,2) = one/cc_2(1)

    DU2 = -DU2

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
    DD(2,1) = factor*dTdrho   + half*gm1*T*dM2drho
    DD(2,2) = factor*dTdrhou  + half*gm1*T*dM2drhou
    DD(2,3) = factor*dTdrhoet + half*gm1*T*dM2drhoet

    dPdrho   = half*gm1*u**2
    dPdrhou  = -gm1*u
    dPdrhoet = gm1

! add Po equations...
    DD(3,1) = factor**gxgm1*dPdrho   + half*gamma*factor**xgm1*P*dM2drho
    DD(3,2) = factor**gxgm1*dPdrhou  + half*gamma*factor**xgm1*P*dM2drhou
    DD(3,3) = factor**gxgm1*dPdrhoet + half*gamma*factor**xgm1*P*dM2drhoet

    DD = -DD

    RHS(1) = cc_in(2)/cc_in(1) - two*cc_1(2)/cc_1(1) + cc_2(2)/cc_2(1)
    RHS(2) = T*factor - to
    RHS(3) = p*factor**gxgm1 - po


  end subroutine subsonic_inflow

!=============================================================================80
!
! Creates matrices for variable extrapolation
!
!=============================================================================80
  subroutine set_outflow(cc_out, cc_1, cc_2, DD, DL1, DL2, RHS)

    use set_precision,   only : dp
    use set_constants,   only : zero, half, one, two
    use fluid_constants, only : gm1
    use initialize_soln, only : pback

    implicit none

    real(dp), dimension(3),   intent(in)  :: cc_out, cc_1, cc_2
    real(dp), dimension(3,3), intent(out) :: DD, DL1, DL2
    real(dp), dimension(3),   intent(out) :: RHS

    real(dp) :: u_out, u_1, u_2
    real(dp) :: p_out, p_1, p_2

    continue

    u_out = cc_out(2) / cc_out(1)
    u_1   = cc_1(2)   / cc_1(1)
    u_2   = cc_2(2)   / cc_2(1)

    p_out = gm1*(cc_out(3) - half*cc_out(1)*u_out**2)
    p_1   = gm1*(cc_1(3)   - half*cc_1(1)*u_1**2)
    p_2   = gm1*(cc_2(3)   - half*cc_2(1)*u_2**2)

! Extrapolate velocity from interior
    DD(1,1) = one
    DD(2,1) = -u_out/cc_out(1)
    DD(3,1) = half*gm1*u_out**2
    DD(1,2) = zero
    DD(2,2) = one/cc_out(1)
    DD(3,2) = -gm1*u_out
    DD(1,3) = zero
    DD(2,3) = zero
    DD(3,3) = gm1

!    DD = -DD

    DL1(1,1) = one
    DL1(2,1) = -u_1/cc_1(1)
    DL1(3,1) = half*gm1*u_1**2
    DL1(1,2) = zero
    DL1(2,2) = one/cc_1(1)
    DL1(3,2) = -gm1*u_1
    DL1(1,3) = zero
    DL1(2,3) = zero
    DL1(3,3) = gm1

    DL1 = -two*DL1

    DL2(1,1) = one
    DL2(2,1) = -u_2/cc_2(1)
    DL2(3,1) = half*gm1*u_2**2
    DL2(1,2) = zero
    DL2(2,2) = one/cc_2(1)
    DL2(3,2) = -gm1*u_2
    DL2(1,3) = zero
    DL2(2,3) = zero
    DL2(3,3) = gm1

!    DL2 = -DL2

    RHS(1) = -cc_out(1) + two*cc_1(1) - cc_2(1)
    RHS(2) = -u_out     + two*u_1     - u_2
    RHS(3) = -p_out     + two*p_1     - p_2

    if (pback > zero) then
      DL1(3,:) = zero
      DL2(3,:) = zero
      RHS(3) = pback - p_out
    end if

  end subroutine set_outflow

  include 'speed_of_sound.f90'

end module bc
