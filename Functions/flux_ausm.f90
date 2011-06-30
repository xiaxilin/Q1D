!================================= flux_ausm =================================80
!
! Takes left and right primitive variables.  Returns flux
!
!=============================================================================80


! FIXME: broken, comment block a lie as well
pure function flux_ausm(qpplus, qpmin, Qplus, Qmin) result(F)

  use set_precision,   only : dp
  use set_constants,   only : fourth, half, one
  use fluid_constants, only : gamma

  implicit none

  real(dp), dimension(3), intent(in)  :: qpplus, qpmin, Qplus, Qmin
  real(dp), dimension(3)              :: F

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

end function flux_ausm
