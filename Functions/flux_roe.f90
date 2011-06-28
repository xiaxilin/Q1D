!================================= flux_roe ==================================80
!
! Takes entropy fix, and left and right primitive variables. Returns flux
!
!=============================================================================80

pure function flux_roe(lambdaeps, qplus, qmin) result(F)

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

end function flux_roe
