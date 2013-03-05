pure function conserved_to_primitive_1D(q)

  use set_precision,   only : dp
  use set_constants,   only : half
  use fluid_constants, only : gm1

  implicit none

  real(dp), dimension(3), intent(in) :: q
  real(dp), dimension(3)             :: conserved_to_primitive_1D

  conserved_to_primitive_1D(1) = q(1)
  conserved_to_primitive_1D(2) = q(2)/q(1)
  conserved_to_primitive_1D(3) = gm1*(q(3) - half*q(2)*q(2)/q(1))

end function conserved_to_primitive_1D
