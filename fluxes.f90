!=============================== van_leer_fvs ================================80
!
! Takes gamma, and left and right primitive variables.  Returns flux
!
!=============================================================================80
pure function van_leer_fvs(qL, qR) result(F)

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

end function van_leer_fvs

!============================ steger_warming_fvs =============================80
!
! Takes left and right primitive variables.  Returns flux
!
!=============================================================================80

pure function steger_warming_fvs(qL, qR) result(F)

  use set_precision,   only : dp
  use set_constants,   only : half, one, two, three
  use fluid_constants, only : gamma, gm1, xg

  implicit none

  real(dp), dimension(3), intent(in)  :: qL, qR

  real(dp), dimension(3)  :: F

  real(dp) :: a
  real(dp), dimension(3) :: lambda, lambdafloor, FL, FR

!Calculate Left (+) Fluxes
  a = sqrt(gamma*qL(3)/qL(1))
  lambda(1) = qL(2)
  lambda(2) = qL(2) + a
  lambda(3) = qL(2) - a

  lambdafloor(:) = half*(lambda(:)+abs(lambda(:)))
    
  FL(1) = two*gm1*lambdafloor(1) + lambdafloor(2) + lambdafloor(3)
  FL(2) = two*gm1*lambdafloor(1)*lambda(1) +                                   &
          lambdafloor(2)*lambda(2) + lambdafloor(3)*lambda(3)
  FL(3) = gm1*lambdafloor(1)*lambda(1)**2 +                                    &
          half*(lambdafloor(2)*lambda(2)**2+lambdafloor(3)*lambda(3)**2) +     &
          ((three-gamma)/(two*gm1))*(lambdafloor(2)+lambdafloor(3))*a**2
                    
  FL(:) = (half*qL(1)*xg)*FL(:)

!Calculate Right (-) Fluxes

  a = sqrt(gamma*qR(3)/qR(1))
  lambda(1) = qR(2)
  lambda(2) = qR(2) + a
  lambda(3) = qR(2) - a

  lambdafloor(:) = half*(lambda(:)-abs(lambda(:)))
    
  FR(1) = two*gm1*lambdafloor(1)+lambdafloor(2)+lambdafloor(3)
  FR(2) = two*gm1*lambdafloor(1)*lambda(1) +                                   &
          lambdafloor(2)*lambda(2) + lambdafloor(3)*lambda(3)
  FR(3) = gm1*lambdamin(1)*lambda(1)**2 +                                      &
          half*(lambdamin(2)*lambda(2)**2 + lambdamin(3)*lambda(3)**2) +       &
          ((three-gamma)/(two*gm1))*(lambdamin(2)+lambdamin(3))*a**2

  FR(:) = (half*qR(1)/gamma)*FR(:)

!Calculate Interface Flux
  F(:) = FL(:) + FR(:)

end function steger_warming_fvs

!================================== ausm_fds =================================80
!
! Takes left and right primitive variables.  Returns flux
!
!=============================================================================80


! FIXME: broken
pure function ausm_fds(gamma, qpplus, qpmin, Qplus, Qmin) result(F)

  use set_precision, only : dp
  use constants,     only : fourth, half, one

  implicit none

  real(dp), intent(in) :: gamma

  real(dp), dimension(3), intent(in)  :: qpplus, qpmin, Qplus, Qmin
  real(dp), dimension(3), intent(out) :: F

  real(dp) :: PL, PR, HTL, HTR, al, ar, Ml, Mr, Mlold, Mrold

  continue

!Calculate left (+) state
  al = sqrt(gamma*qpplus(3)/qpplus(1))
  Ml = qpplus(2)/al
  PL = qpplus(3)
  HTL = (Qplus(3)+PL)/Qplus(1)
    
  if (abs(Ml)<1.0_dp) then
    PL = PL*half*(one+Ml)
    Ml = fourth*(Ml+one)**2
  else
    Mlold = Ml
    Ml = half*(Ml+abs(Ml))
    PL = PL*Ml/Mlold	
  end if

!Calculate right (-) state
  ar = sqrt(gamma*qpmin(3)/qpmin(1))
  Mr = qpmin(2)/ar
  PR = qpmin(3)
  HTR = (Qmin(3)+PR)/Qmin(1)

  if (abs(Mr)<1.0_dp) then
    PR = PR*half*(one-Mr)
    Mr = -fourth*(Mr-one)**2
  else
    Mrold = Mr
    Mr = half*(Mr-abs(Mr))
    PR = PR*Mr/Mrold
  end if

!Combine
  F(1) = half*(Qmin(1)*ar*((Ml+Mr) - abs(Ml+Mr)) +                             &
              Qplus(1)*ar*((Ml+Mr) + abs(Ml+Mr)))
  F(2) = half*(Qmin(2)*ar*((Ml+Mr) - abs(Ml+Mr)) +                             &
              Qplus(2)*ar*((Ml+Mr) + abs(Ml+Mr))) + (PL+PR)
  F(3) = half*(Qmin(1)*HTR*ar*((Ml+Mr) - abs(Ml+Mr)) +                         &
              Qplus(1)*HTL*ar*((Ml+Mr) + abs(Ml+Mr)))
  end do

end subroutine ausm_fds

!================================== roes_fds =================================80
!
! Takes entropy fix, and left and right primitive variables. Returns flux
!
!=============================================================================80

pure function roes_fds(lambdaeps, qplus, qmin) result(F)

  use set_precision,   only : dp
  use set_constants,   only : half, one, two, four
  use fluid_constatns, only : gamma, gm1

  implicit none

  real(dp), intent(in) :: lambdaeps
  real(dp), dimension(3), intent(in)  :: qplus, qmin

  real(dp), dimension(3) :: F

  integer  :: i
  real(dp) :: rhoL, uL, pL, rhoR, uR, pR
  real(dp) :: Rhalf, RoeAvgrho, RoeAvgu, RoeAvght, RoeAvga
  real(dp), dimension(3) :: FL, FR, r1RoeAvg, r2RoeAvg, r3RoeAvg, dw, lambdaRoe

  continue

! FIXME: need rho*et term

  rhoL = qL(1)
  uL   = qL(2)
  pL   = qL(3)

  rhoR = qR(1)
  uR   = qR(2)
  pR   = qR(3)

  FL(:) = (/rhoL*uL, rhoL*uL**2 + pL, uL*(Qplus(3) + pL)/)
  FR(:) = (/rhoR*uR, rhoR*uR**2 + pR, uR*(Qmin(3)  + pR)/)

! Roe interface variable
  Rhalf = sqrt(rhoR/rhoL)

! Calculate Roe average variables
  RoeAvgrho = Rhalf*rhoL
  RoeAvgu   = (Rhalf*uR + uL) / (Rhalf+one)
  RoeAvght  = (Rhalf*(Qmin(3) + pR)/rhoR + (Qplus(3) + pL)/rhoL) / (Rhalf+one)
  RoeAvga   = sqrt(gm1*(RoeAvght-half*RoeAvgu**2))

  lambdaRoe(1) = RoeAvgu
  lambdaRoe(2) = RoeAvgu + RoeAvga
  lambdaRoe(3) = RoeAvgu - RoeAvga

!Entropy fix
  do i = 1,3
    if ( abs(lambdaRoe(i)) <= two*lambdaeps*RoeAvga ) then
      lambdaRoe(i) = (lambdaRoe(i)**2)/(four*lambdaeps*RoeAvga)                &
                   + lambdaeps*RoeAvga
    end if
  end do

  dw(1) = (rhoR-rhoL) - (pR-pL)/RoeAvga**2
  dw(2) = (uR-uL) + (pR-pL)/(RoeAvgrho*RoeAvga)
  dw(3) = (uR-uL) - (pR-pL)/(RoeAvgrho*RoeAvga)

  r1RoeAvg(:) = (/one, RoeAvgu, half*RoeAvgu**2/)
  r2RoeAvg(:) = (half*RoeAvgrho/RoeAvga) *                                     &
                (/one, RoeAvgu+RoeAvga, RoeAvght+RoeAvgu*RoeAvga/)
  r3RoeAvg(:) = (-half*RoeAvgrho/RoeAvga) *                                    &
                (/one, RoeAvgu-RoeAvga, RoeAvght-RoeAvgu*RoeAvga/)

!Calculate Interface Fluxes
  F(:) = half*((FL(:)+FR(:))                                                   &
       - (abs(lambdaRoe(1))*dw(1)*r1RoeAvg(:)                                  &
       +  abs(lambdaRoe(2))*dw(2)*r2RoeAvg(:)                                  &
       +  abs(lambdaRoe(3))*dw(3)*r3RoeAvg(:)))

end function roes_fds
