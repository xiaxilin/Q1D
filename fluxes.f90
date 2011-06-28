!=============================== van_leer_fvs =================================80
!
! Takes gamma, and left and right primitive variables.  Returns flux
!
!==============================================================================80
pure function van_leer_fvs(gamma, qL, qR) result(F)

  use set_precision, only : dp
  use constants,     only : zero, fourth, half, one, two

  implicit none

  real(dp), intent(in) :: gamma
  real(dp), dimension(3), intent(in)    :: qL, qR

  real(dp), dimension(3)  :: F

  real(dp) :: a, M, Mfloor, fa, fb, switch
  real(dp), dimension(3) :: Fisub, Fiss
  real(dp), dimension(3) :: Fiplus, Fimin

  continue

!Calculate Left (+) Flux
  a = sqrt(gamma*qL(3)/qL(1))
  M = qL(2)/a

!Left sub(sonic) flux
  fa = fourth*qL(1)*a*(M+one)**2
  fb = a*((gamma-one)*M + two)

  Fisub(1) = fa
  Fisub(2) = fa*fb / gamma
  Fisub(3) = fa*fb*fb / (two*(gamma**2-one))

!Floor the mach number for -supersonic flows
  Mfloor = half*(M + abs(M))
!Left supersonic flux
  Fiss(1) = qL(1)*a*Mfloor
  Fiss(2) = (Mfloor/M)*qL(1)*a**2*(Mfloor**2+one/gamma)
  Fiss(3) = qL(1)*a**3*Mfloor*(half*Mfloor**2+one/(gamma-one))

!Combine sub and supersonic fluxes
  switch = max(zero, one-int(abs(M)))
  Fiplus(:) = (one-switch)*Fiss(:) + switch*Fisub(:)

!Calculate Right (-) Fluxes
  a = sqrt(gamma*qR(3)/qR(1))
  M = qmin(2)/a

!Right sub(sonic) flux
  fa = -fourth*qR(1)*a*(M-one)**2
  fb = -a*((gamma-one)*M + two)

  Fisub(1) = fa
  Fisub(2) = fa*fb / gamma
  Fisub(3) = fa*fb*fb / (two*(gamma**2-one))

!Floor the mach number for +supersonic flows
  Mfloor = half*(M - abs(M))

  Fiss(1) = qR(1)*a*M
  Fiss(2) = (Mfloor/M)*qR(1)*a**2*(Mfloor**2+one/gamma)
  Fiss(3) = qR(1)*a**3*Mfloor*(half*Mfloor**2+one/(gamma-one))

!Combine sub and supersonic fluxes
  switch = max(zero, one-int(abs(M)))
  Fimin(:) = (one-switch)*Fiss(:) + switch*Fisub(:)

!Calculate Interface Flux
  F(:) = Fiplus(:) + Fimin(:)

end function van_leer_fvs

!============================ steger_warming_fvs ==============================80
!
! Takes gamma, and left and right primitive variables.  Returns flux
!
!==============================================================================80

pure function steger_warming_fvs(gamma, qplus, qmin), result(F)

  use set_precision, only : dp
  use constants,     only : half, one, two, three

  Implicit None

  real(dp), intent(in) :: gamma
  real(dp), dimension(3), intent(in)  :: qplus, qmin

  real(dp), dimension(3)  :: F

  real(dp) :: al, ar
  real(dp), dimension(3) :: lambdal,lambdar,lambdaplus,lambdamin,Fiplus,Fimin

!Calculate Left (+) Fluxes
  al = sqrt(gamma*qplus(3)/qplus(1))
  lambdal(1) = qplus(2)
  lambdal(2) = qplus(2) + al
  lambdal(3) = qplus(2) - al

  lambdaplus(:) = half*(lambdal(:)+abs(lambdal(:)))
    
  Fiplus(1) = two*(gamma-one)*lambdaplus(1) + lambdaplus(2) + lambdaplus(3)
  Fiplus(2) = two*(gamma-one)*                                                  &
              lambdaplus(1)*lambdal(1) +                                        &
              lambdaplus(2)*lambdal(2) +                                        &
              lambdaplus(3)*lambdal(3)
  Fiplus(3) = (gamma-one)*lambdaplus(1)*lambdal(1)**2 +                         &
              half*(lambdaplus(2)*lambdal(2)**2+lambdaplus(3)*lambdal(3)**2) +  &
              ((three-gamma)/(two*(gamma-one)))*                                &
              (lambdaplus(2)+lambdaplus(3))*al**2
                    
  Fiplus(:) = (half*qplus(1)/gamma)*Fiplus(:)

!Calculate Right (-) Fluxes

  ar = sqrt(gamma*qmin(3)/qmin(1))
  lambdar(1) = qmin(2)
  lambdar(2) = qmin(2) + ar
  lambdar(3) = qmin(2) - ar

  lambdamin(:) = half*(lambdar(:)-abs(lambdar(:)))
    
  Fimin(1) = two*(gamma-one)*lambdamin(1)+lambdamin(2)+lambdamin(3)
  Fimin(2) = two*(gamma-one)*                                                   &
             lambdamin(1)*lambdar(1) +                                          &
             lambdamin(2)*lambdar(2) +                                          &
             lambdamin(3)*lambdar(3)
  Fimin(3) = (gamma-one)*lambdamin(1)*lambdar(1)**2 +                           &
             half*(lambdamin(2)*lambdar(2)**2 + lambdamin(3)*lambdar(3)**2) +   &
             ((three-gamma)/(two*(gamma-one)))*                                 &
             (lambdamin(2)+lambdamin(3))*ar**2

  Fimin(:) = (half*qmin(1)/gamma)*Fimin(:)

!Calculate Interface Flux
  F(:) = Fiplus(:) + Fimin(:)

end function steger_warming_fvs

!================================== ausm_fds ==================================80
!
! Takes gamma, and left and right primitive variables.  Returns flux
!
!==============================================================================80

pure function ausm_fds(gamma, qpplus, qpmin, Qplus, Qmin), return(F)

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
  F(1) = half*(Qmin(1)*ar*((Ml+Mr) - abs(Ml+Mr)) +                              &
              Qplus(1)*ar*((Ml+Mr) + abs(Ml+Mr)))
  F(2) = half*(Qmin(2)*ar*((Ml+Mr) - abs(Ml+Mr)) +                              &
              Qplus(2)*ar*((Ml+Mr) + abs(Ml+Mr))) + (PL+PR)
  F(3) = half*(Qmin(1)*HTR*ar*((Ml+Mr) - abs(Ml+Mr)) +                          &
              Qplus(1)*HTL*ar*((Ml+Mr) + abs(Ml+Mr)))
  end do

end subroutine ausm_fds

!================================== roes_fds ==================================80
!
! Takes gamma, entropy fix, and left and right primitive variables.  Returns flux
!
!==============================================================================80

pure function roes_fds(gamma, lambdaeps, qplus, qmin), result(F)

  use set_precision, only : dp
  use constants,     only : half, one, two, four

  implicit none

  real(dp), intent(in) :: gamma, lambdaeps
  real(dp), dimension(3), intent(in)  :: qplus, qmin

  real(dp), dimension(3) :: F

  integer  :: z
  real(dp) :: rhoL, uL, PL, rhoR, uR, PR
  real(dp) ::  Rhalf, RoeAvgrho, RoeAvgu, RoeAvght, RoeAvga
  real(dp), dimension(3) :: FL, FR, r1RoeAvg, r2RoeAvg, r3RoeAvg, dw, lambdaRoe

  continue

  rhoL = qplus(1)
  uL   = qplus(2)
  PL   = qplus(3)

  rhoR = qmin(1)
  uR   = qmin(2)
  PR   = qmin(3)

  FL(:) = (/rhoL*uL, rhoL*uL**2 + PL, uL*(Qplus(3) + PL)/)
  FR(:) = (/rhoR*uR, rhoR*uR**2 + PR, uR*(Qmin(3) + PR)/)

  Rhalf = sqrt(rhoR/rhoL)                !Roe interface variable
  RoeAvgrho = Rhalf*rhoL
  RoeAvgu   = (Rhalf*uR + uL) / (Rhalf+one)
  RoeAvght  = (Rhalf*(Qmin(3) + PR)/rhoR + (Qplus(3) + PL)/rhoL) / (Rhalf+one)
  RoeAvga   = sqrt((gamma-one)*(RoeAvght-half*RoeAvgu**2))

  lambdaRoe(1) = RoeAvgu
  lambdaRoe(2) = RoeAvgu + RoeAvga
  lambdaRoe(3) = RoeAvgu - RoeAvga

!Entropy fix
  do z = 1,3
    if ( abs(lambdaRoe(z)) <= two*lambdaeps*RoeAvga ) then
      lambdaRoe(z) = (lambdaRoe(z)**2)/(four*lambdaeps*RoeAvga) +               &
           lambdaeps*RoeAvga
    end if
  end do

  dw(1) = (rhoR-rhoL) - (PR-PL)/RoeAvga**2
  dw(2) = (uR-uL) + (PR-PL)/(RoeAvgrho*RoeAvga)
  dw(3) = (uR-uL) - (PR-PL)/(RoeAvgrho*RoeAvga)

  r1RoeAvg(:) = (/one, RoeAvgu, half*RoeAvgu**2/)
  r2RoeAvg(:) = (half*RoeAvgrho/RoeAvga) *                                      &
                (/one, RoeAvgu+RoeAvga, RoeAvght+RoeAvgu*RoeAvga/)
  r3RoeAvg(:) = (-half*RoeAvgrho/RoeAvga) *                                     &
                (/one, RoeAvgu-RoeAvga, RoeAvght-RoeAvgu*RoeAvga/)

!Calculate Interface Fluxes
  F(:) = half*((FL(:)+FR(:))                                                    &
       - (abs(lambdaRoe(1))*dw(1)*r1RoeAvg(:)                                   &
       +  abs(lambdaRoe(2))*dw(2)*r2RoeAvg(:)                                   &
       +  abs(lambdaRoe(3))*dw(3)*r3RoeAvg(:)))

end function roes_fds
