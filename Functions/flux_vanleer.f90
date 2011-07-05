!=============================== flux_vanleer ================================80
!
! Takes left and right primitive variables.  Returns flux
!
!=============================================================================80
pure function flux_vanleer(qL, qR) result(F)

  use set_precision,   only : dp
  use set_constants,   only : zero, fourth, half, one, two
  use fluid_constants, only : gamma, gm1, xg, xgm1, xg2m1

  implicit none

  real(dp), dimension(3), intent(in)    :: qL, qR

  real(dp), dimension(3)  :: F

  real(dp) :: a, M, Mfloor, fa, fb, switch
  real(dp), dimension(3) :: Fisub, Fiss

  continue

!Calculate Left (+) Flux
  a = speed_of_sound(qL(3), qL(1))
  M = qL(2)/a

!Left sub(sonic) flux
  fa = fourth*qL(1)*a*(M+one)**2
  fb = a*(gm1*M + two)

  Fisub(1) = fa
  Fisub(2) = fa*fb*xg
  Fisub(3) = half*fa*fb*fb*xg2m1

!Floor the mach number for -supersonic flows
  Mfloor = half*(M + abs(M))
!Left supersonic flux
  Fiss(1) = qL(1)*a*Mfloor
  Fiss(2) = (Mfloor/M)*qL(1)*a**2*(Mfloor**2+xg)
  Fiss(3) = qL(1)*a**3*Mfloor*(half*Mfloor**2+xgm1)

!Combine sub and supersonic fluxes, store in final flux vector
  switch = max(zero, one-real(int(abs(M)),dp))
  F(:) = (one-switch)*Fiss(:) + switch*Fisub(:)

!Calculate Right (-) Fluxes
  a = speed_of_sound(qR(3), qR(1))
  M = qR(2)/a

!Right sub(sonic) flux
  fa = -fourth*qR(1)*a*(M-one)**2
  fb = a*(gm1*M - two)

  Fisub(1) = fa
  Fisub(2) = fa*fb*xg
  Fisub(3) = half*fa*fb*fb*xg2m1

!Floor the mach number for +supersonic flows
  Mfloor = half*(M - abs(M))

  Fiss(1) = qR(1)*a*Mfloor
  Fiss(2) = (Mfloor/M)*qR(1)*a**2*(Mfloor**2+xg)
  Fiss(3) = qR(1)*a**3*Mfloor*(half*Mfloor**2+xgm1)

!Combine sub and supersonic fluxes, finish final flux vector
  switch = max(zero, one-real(int(abs(M)),dp))
  F(:) = F(:) + (one-switch)*Fiss(:) + switch*Fisub(:)

end function flux_vanleer
