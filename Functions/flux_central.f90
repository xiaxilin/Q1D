!=============================== flux_central ================================80
!
! Takes left and right primitive variables.  Returns flux
!
!=============================================================================80
pure function flux_central(qL, qR)

  use set_precision,   only : dp
  use set_constants,   only : half

  implicit none

  real(dp), dimension(3), intent(in) :: qL, qR
  real(dp), dimension(3)             :: flux_central, F, consL, consR

  consL = primitive_to_conserved_1D(qL)
  consR = primitive_to_conserved_1D(qR)

  F(1) = half * ( consL(2) + consR(2) )
  F(2) = half * ( qL(1)*qL(2)**2 + qL(3) + qR(1)*qR(2)**2 + qR(3) )
  F(3) = half * ( qL(2)*(consL(3) + qL(3)) + qR(2)*(consR(3) + qR(3)) )

  flux_central = f

end function flux_central
