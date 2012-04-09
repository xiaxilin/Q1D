!=============================== flux_ausm_plus ==============================80
!
! Takes left and right primitive variables.  Returns flux at interface
!
!=============================================================================80

 function flux_ausm_plus(prim_L, prim_R)

  use set_precision,   only : dp
  use set_constants,   only : fourth, half, one, two

  implicit none

  real(dp), dimension(3), intent(in)  :: prim_L, prim_R
  real(dp), dimension(3)              :: flux_ausm_plus, fL, fR, cons_L, cons_R

  real(dp) :: rhoL, aL, uL, PL, HTL, rhoR, aR, uR, PR, HTR
  real(dp) :: a_star, aL_tilde, aR_tilde, a_half, m_half, p_half
  real(dp) :: ML, ML_script, PL_script, MR, MR_script, PR_script

  real(dp), parameter :: alpha = 1.0_dp/8.0_dp
  real(dp), parameter :: beta  = 3.0_dp/16.0_dp

  continue

!Calculate left (+) state
  rhoL = prim_L(1)
  uL   = prim_L(2)
  PL   = prim_L(3)
  aL   = speed_of_sound(PL, rhoL)

  cons_L = primitive_to_conserved_1D(prim_L)

  HTL = (cons_L(3)+PL)/rhoL

!Calculate right (-) state
  rhoR = prim_R(1)
  uR   = prim_R(2)
  PR   = prim_R(3)
  aR   = speed_of_sound(PR, rhoR)

  cons_R = primitive_to_conserved_1D(prim_R)

  HTR = (cons_R(3)+PR)/rhoR

!Interface speed of sound and mach numbers
  a_star = sqrt(uL*uR)

  aL_tilde = a_star**2/max(a_star, abs(uL))
  aR_tilde = a_star**2/max(a_star, abs(uR))

  a_half = min(aL_tilde, aR_tilde)

  ML = uL/a_half
  MR = uR/a_half

!Get other interface properties

!Left side
  if (abs(Ml)<1.0_dp) then
    ML_script = half*(ML+one)**2 + beta*(ML**2-1)**2
    PL_script = fourth*(two-ML)*(ML+one)**2 + alpha*ML*(ML**2-one)**2
  else
    ML_script = half*(ML+abs(ML))
    PL_script = half*(one+sign(one,ML))
  end if

!Right side
  if (abs(Mr)<1.0_dp) then
    MR_script = -half*(MR-one)**2 - beta*(MR**2-1)**2
    PR_script = fourth*(two+MR)*(MR-one)**2 - alpha*MR*(MR**2-one)**2
  else
    MR_script = half*(MR-abs(MR))
    PR_script = half*(one-sign(one,MR))
  end if

  m_half = ML_script + MR_script
  p_half = PL_script*PL + PR_script*PR

  FL = (/rhoL, rhoL*uL, rhoL*HTL/)
  FR = (/rhoR, rhoL*uR, rhoR*HTR/)

  flux_ausm_plus = half*a_half*(m_half*(FL+FR)-abs(m_half)*(FR-FL))
  flux_ausm_plus(2) = flux_ausm_plus(2) + p_half

end function flux_ausm_plus
