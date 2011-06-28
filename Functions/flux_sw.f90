!============================ steger_warming_fvs =============================80
!
! Takes left and right primitive variables.  Returns flux
!
!=============================================================================80

pure function flux_sw(qL, qR) result(F)

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

end function flux_sw
