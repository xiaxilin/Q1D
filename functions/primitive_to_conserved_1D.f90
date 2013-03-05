pure function primitive_to_conserved_1D(qp)

  use set_precision,   only : dp
  use set_constants,   only : half
  use fluid_constants, only : xgm1

  implicit none

  real(dp), dimension(3), intent(in) :: qp
  real(dp), dimension(3)             :: primitive_to_conserved_1D

  primitive_to_conserved_1D(1) = qp(1)
  primitive_to_conserved_1D(2) = qp(1)*qp(2)
  primitive_to_conserved_1D(3) = qp(3)*xgm1 + half*qp(1)*qp(2)*qp(2)

end function primitive_to_conserved_1D
