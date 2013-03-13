!=============================== flux_ausm_plus ==============================80
!
! Takes left and right primitive variables.  Returns flux at interface
!
!=============================================================================80
function flux_ausm_plus(prim_L, prim_R)

  use set_precision,   only : dp
  use set_constants,   only : fourth, half, one, two
  use fluid_constants, only : xgm1

  implicit none

  real(dp), dimension(3), intent(in)  :: prim_L, prim_R
  real(dp), dimension(3)              :: flux_ausm_plus, FL, FR

  real(dp) :: aL, aR, rhoL, uL, PL, HTL, rhoR, uR, PR, HTR
  real(dp) :: a_star, a_half, m_half, p_half!, aL_tilde, aR_tilde,
  real(dp) :: ML, ML_script, PL_script, MR, MR_script, PR_script

  real(dp), parameter :: alpha = 3.0_dp/16.0_dp
  real(dp), parameter :: beta  = 1.0_dp/8.0_dp

  continue

!Calculate left (+) state
  rhoL = prim_L(1)
  uL   = prim_L(2)
  PL   = prim_L(3)
  aL   = speed_of_sound(pL,rhoL)
  HTL  = xgm1*aL**2 + half*uL**2

!Calculate right (-) state
  rhoR = prim_R(1)
  uR   = prim_R(2)
  PR   = prim_R(3)
  aR   = speed_of_sound(pR,rhoR)
  HTR  = xgm1*aR**2 + half*uR**2

!Interface speed of sound and mach numbers
  a_star = sqrt(aL*aR)

!  aL_tilde = a_star**2/max(a_star, abs(uL))
!  aR_tilde = a_star**2/max(a_star, abs(uR))

!  a_half = min(aL_tilde, aR_tilde)

  a_half = a_star

  ML = uL/a_half
  MR = uR/a_half

!Get other interface properties
!Note the 'fourth' in front of the subsonic ML/MR_script terms...
!There appears to be an error in eq. 19B of Liou, JCP 129, 1996

!Left side
  if ( abs(ML) < 1.0_dp ) then
    PL_script = fourth*(two-ML)*(ML+one)**2 + alpha*ML*(ML**2-one)**2
    ML_script = fourth*(ML+one)**2 + beta*(ML**2-one)**2
  else
    PL_script = half*(one+sign(one,ML))
    ML_script = half*(ML+abs(ML))
  end if

!Right side
  if ( abs(MR) < 1.0_dp ) then
    PR_script = fourth*(two+MR)*(MR-one)**2 - alpha*MR*(MR**2-one)**2
    MR_script = -fourth*(MR-one)**2 - beta*(MR**2-one)**2
  else
    PR_script = half*(one-sign(one,MR))
    MR_script = half*(MR-abs(MR))
  end if

  m_half = ML_script + MR_script
  p_half = PL_script*PL + PR_script*PR

  FL = [rhoL, rhoL*uL, rhoL*HTL]
  FR = [rhoR, rhoR*uR, rhoR*HTR]

  flux_ausm_plus = half*a_half*(m_half*(FL+FR)-abs(m_half)*(FR-FL))
  flux_ausm_plus(2) = flux_ausm_plus(2) + p_half

end function flux_ausm_plus
