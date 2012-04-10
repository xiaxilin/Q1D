!================================= flux_ausm =================================80
!
! Takes left and right primitive variables.  Returns flux
!
!=============================================================================80

 function flux_ausm(prim_L, prim_R)

  use set_precision,   only : dp
  use set_constants,   only : fourth, half, one
  use fluid_constants, only : xgm1

  implicit none

  real(dp), dimension(3), intent(in)  :: prim_L, prim_R
  real(dp), dimension(3)              :: flux_ausm

  real(dp) :: rhoL, uL, PL, aL, ML, HTL, rhoR, uR, PR, aR, MR, HTR
  real(dp), dimension(3) :: FL, FR

  continue

!Calculate left (+) state
  rhoL = prim_L(1)
  uL   = prim_L(2)
  PL   = prim_L(3)
  aL   = speed_of_sound(PL, rhoL)
  ML   = uL/aL
  HTL  = xgm1*aL**2 + half*uL**2

  if (abs(ML)<=1.0_dp) then
    PL = half*PL*(one+ML)
    ML = fourth*(ML+one)**2
  else
    PL = half*PL*(one+sign(ML,one))
    ML = half*(ML+abs(ML))
  end if

!Calculate right (-) state
  rhoR = prim_R(1)
  uR   = prim_R(2)
  PR   = prim_R(3)
  aR   = speed_of_sound(PR, rhoR)
  MR   = uR/aR
  HTR  = xgm1*aR**2 + half*uR**2

  if (abs(Mr)<=1.0_dp) then
    PR = half*PR*(one-MR)
    MR = -fourth*(MR-one)**2
  else
    PR = half*PR*(one-sign(MR,one))
    MR = half*(MR-abs(MR))
  end if

  FL = (/rhoL*aL, rhoL*aL*uL, rhoL*aL*HTL/)
  FR = (/rhoR*aR, rhoR*aR*uR, rhoR*aR*HTR/)

  flux_ausm = half*( (ML+MR)*(FL+FR) - abs(ML+MR)*(FR-FL) )
  flux_ausm(2) = flux_ausm(2) + PL + PR

end function flux_ausm
