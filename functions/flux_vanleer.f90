!=============================== flux_vanleer ================================80
!
! Takes left and right primitive variables.  Returns flux
!
!=============================================================================80
pure function flux_vanleer(qL, qR)

  use set_precision,   only : dp
  use set_constants,   only : zero, fourth, half, one, two
  use fluid_constants, only : gm1, xg, xgm1, xg2m1

  implicit none

  real(dp), dimension(3), intent(in) :: qL, qR
  real(dp), dimension(3)             :: flux_vanleer

  real(dp) :: a, M, fa, fb
  real(dp), dimension(3) :: FL, FR

  continue

!Calculate Left (+) Flux
  a = speed_of_sound(qL(3), qL(1))
  M = qL(2)/a

  if ( abs(M) < one ) then !Left subsonic flux
    fa = fourth*qL(1)*a*(M+one)**2
    fb = a*(gm1*M + two)

    FL(1) = fa
    FL(2) = fa*fb*xg
    FL(3) = half*fa*fb*fb*xg2m1
  else if ( M >= one ) then !Left supersonic flux
    FL(1) = qL(1)*a*M
    FL(2) = qL(1)*a**2*(M**2+xg)
    FL(3) = qL(1)*a**3*M*(half*M**2+xgm1)
  else
    FL = zero
  end if

!Calculate Right (-) Fluxes
  a = speed_of_sound(qR(3), qR(1))
  M = qR(2)/a

  if ( abs(M) < one ) then !Right subsonic flux
    fa = -fourth*qR(1)*a*(M-one)**2
    fb = a*(gm1*M - two)

    FR(1) = fa
    FR(2) = fa*fb*xg
    FR(3) = half*fa*fb*fb*xg2m1
  else if ( M <= -one ) then !Right supersonic flux
    FR(1) = qR(1)*a*M
    FR(2) = qR(1)*a**2*(M**2+xg)
    FR(3) = qR(1)*a**3*M*(half*M**2+xgm1)
  else
    FR = zero
  end if

  flux_vanleer = FL + FR

end function flux_vanleer
