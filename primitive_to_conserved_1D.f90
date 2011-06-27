pure function primitive_to_conserved_1D(qp) result(q)

  use set_precision, only : dp
  use set_constants, only : half
  use fluid,         only : xgm1

  implicit none

  real(dp), dimension(3), intent(in) :: qp
  real(dp), dimension(3)             :: q

  q(1) = qp(1)
  q(2) = qp(1)*qp(2)
  q(3) = qp(3)*xgm1 + half*qp(1)*qp(2)*qp(2)

end function primitive_to_conserved_1D
