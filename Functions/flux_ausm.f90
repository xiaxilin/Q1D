!================================= flux_ausm =================================80
!
! Takes left and right primitive variables.  Returns flux
!
!=============================================================================80

 function flux_ausm(prim_L, prim_R) result(F)

  use set_precision,   only : dp
  use set_constants,   only : fourth, half, one

  implicit none

  real(dp), dimension(3), intent(in)  :: prim_L, prim_R
  real(dp), dimension(3)              :: F, cons_L, cons_R

  real(dp) :: PL, PR, HTL, HTR, al, ar, Ml, Mr, Mlold, Mrold

  continue

!Calculate left (+) state
  al = speed_of_sound(prim_L(3), prim_L(1))
  Ml = prim_L(2)/al
  PL = prim_L(3)

  cons_L = primitive_to_conserved_1D(prim_L)

  HTL = (cons_L(3)+PL)/cons_L(1)
    
  if (abs(Ml)<1.0_dp) then
    PL = PL*half*(one+Ml)
    Ml = fourth*(Ml+one)**2
  else
    Mlold = Ml
    Ml = half*(Ml+abs(Ml))
    PL = PL*Ml/Mlold
  end if

!Calculate right (-) state
  ar = speed_of_sound(prim_R(3), prim_R(1))
  Mr = prim_R(2)/ar
  PR = prim_R(3)

  cons_R = primitive_to_conserved_1D(prim_R)

  HTR = (cons_R(3)+PR)/cons_R(1)

  if (abs(Mr)<1.0_dp) then
    PR = PR*half*(one-Mr)
    Mr = -fourth*(Mr-one)**2
  else
    Mrold = Mr
    Mr = half*(Mr-abs(Mr))
    PR = PR*Mr/Mrold
  end if

!Combine FIXME: should some ar's be al's?
  F(1) = half*(cons_R(1)*ar*((Ml+Mr) - abs(Ml+Mr)) +                           &
               cons_L(1)*ar*((Ml+Mr) + abs(Ml+Mr)))
  F(2) = half*(cons_R(2)*ar*((Ml+Mr) - abs(Ml+Mr)) +                           &
               cons_L(2)*ar*((Ml+Mr) + abs(Ml+Mr))) + (PL+PR)
  F(3) = half*(cons_R(1)*HTR*ar*((Ml+Mr) - abs(Ml+Mr)) +                       &
               cons_L(1)*HTL*ar*((Ml+Mr) + abs(Ml+Mr)))

end function flux_ausm
