!============================= limiter_vanalbada =============================80
!
! Takes # of equations and successive gradients
!
!=============================================================================80
pure function limiter_vanalbada(neq, r)

  use set_precision,   only : dp
  use set_constants,   only : one

  implicit none

  integer,                  intent(in) :: neq
  real(dp), dimension(neq), intent(in) :: r

  integer                  :: eq
  real(dp), dimension(neq) :: limiter_vanalbada

  do eq = 1, neq
    limiter_vanalbada(eq) = (r(eq)*r(eq) + r(eq)) / (r(eq)*r(eq) + one)
    if (r(eq) < zero) limiter_vanalbada(eq) = zero
  end do

end function limiter_vanalbada
