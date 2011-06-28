!=============================== flux_vanleer ================================80
!
! Takes gamma, and left and right primitive variables.  Returns flux
!
!=============================================================================80
pure function flux_vanleer(qL, qR) result(F)

  use set_precision,   only : dp
  use set_constants,   only : zero, fourth, half, one, two
  use fluid_constants, only : gamma, gm1, xg, xgm1

  implicit none

  real(dp), dimension(3), intent(in)    :: qL, qR

  real(dp), dimension(3)  :: F

  real(dp) :: a, M, Mfloor, fa, fb, switch
  real(dp), dimension(3) :: Fisub, Fiss

  continue

!Calculate Left (+) Flux
  a = sqrt(gamma*qL(3)/qL(1)) ! FIXME: make function
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
  switch = max(zero, one-int(abs(M)))
  F(:) = (one-switch)*Fiss(:) + switch*Fisub(:)

!Calculate Right (-) Fluxes
  a = sqrt(gamma*qR(3)/qR(1))
  M = qmin(2)/a

!Right sub(sonic) flux
  fa = -fourth*qR(1)*a*(M-one)**2
  fb = -a*(gm1*M + two)

  Fisub(1) = fa
  Fisub(2) = fa*fb*xg
  Fisub(3) = half*fa*fb*fb*xg2m1

!Floor the mach number for +supersonic flows
  Mfloor = half*(M - abs(M))

  Fiss(1) = qR(1)*a*M
  Fiss(2) = (Mfloor/M)*qR(1)*a**2*(Mfloor**2+xg)
  Fiss(3) = qR(1)*a**3*Mfloor*(half*Mfloor**2+xgm1)

!Combine sub and supersonic fluxes, finish final flux vector
  switch = max(zero, one-int(abs(M)))
  F(:) = F(:) + (one-switch)*Fiss(:) + switch*Fisub(:)

end function flux_vanleer
