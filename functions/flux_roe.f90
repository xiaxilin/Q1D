!================================= flux_roe ==================================80
!
! Takes entropy fix, and left and right primitive variables. Returns flux
!
!=============================================================================80
pure function flux_roe(qL, qR)

  use set_precision,   only : dp
  use set_constants,   only : half, one, two, four
  use fluid_constants, only : gm1

  implicit none

  real(dp), dimension(3), intent(in) :: qL, qR
  real(dp), dimension(3)             :: flux_roe, F

  integer  :: i
  real(dp) :: rhoL, uL, pL, rhoR, uR, pR, lambdaeps
  real(dp) :: Rhalf, RoeAvgrho, RoeAvgu, RoeAvght, RoeAvga
  real(dp), dimension(3) :: FL, FR, consL, consR
  real(dp), dimension(3) :: r1RoeAvg, r2RoeAvg, r3RoeAvg, dw, lambdaRoe

  continue

! FIXME: need rho*et term
  lambdaeps = 0.1_dp

  rhoL = qL(1)
  uL   = qL(2)
  pL   = qL(3)

  rhoR = qR(1)
  uR   = qR(2)
  pR   = qR(3)

  consL = primitive_to_conserved_1D(qL)
  consR = primitive_to_conserved_1D(qR)

  FL = [rhoL*uL, rhoL*uL**2 + pL, uL*(consL(3) + pL)]
  FR = [rhoR*uR, rhoR*uR**2 + pR, uR*(consR(3) + pR)]

! Roe interface variable
  Rhalf = sqrt(rhoR/rhoL)

! Calculate Roe average variables
  RoeAvgrho = Rhalf*rhoL
  RoeAvgu   = (Rhalf*uR + uL) / (Rhalf+one)
  RoeAvght  = (Rhalf*(consR(3) + pR)/rhoR + (consL(3) + pL)/rhoL) / (Rhalf+one)
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
  dw(2) = half*( (pR-pL) + RoeAvgrho*RoeAvga*(uR-uL) ) / (RoeAvga*RoeAvga)
  dw(3) = half*( (pR-pL) - RoeAvgrho*RoeAvga*(uR-uL) ) / (RoeAvga*RoeAvga)

  r1RoeAvg = [one, RoeAvgu, half*RoeAvgu**2]
  r2RoeAvg = [one, RoeAvgu+RoeAvga, RoeAvght+RoeAvgu*RoeAvga]
  r3RoeAvg = [one, RoeAvgu-RoeAvga, RoeAvght-RoeAvgu*RoeAvga]

!Calculate Interface Fluxes
  F = half*( (FL + FR)                                                         &
    - (abs(lambdaRoe(1))*dw(1)*r1RoeAvg                                        &
    +  abs(lambdaRoe(2))*dw(2)*r2RoeAvg                                        &
    +  abs(lambdaRoe(3))*dw(3)*r3RoeAvg) )

  flux_roe = f

end function flux_roe
