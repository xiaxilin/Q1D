!============================= limiter_vanleer ===============================80
!
! Takes # of equations and successive gradients
!
!=============================================================================80
pure function limiter_vanleer(neq, r)

  use set_precision,   only : dp
  use set_constants,   only : one

  implicit none

  integer,                  intent(in) :: neq
  real(dp), dimension(neq), intent(in) :: r

  integer                  :: eq
  real(dp), dimension(neq) :: limiter_vanleer

  do eq = 1, neq
    limiter_vanleer(eq) = (r(eq) + abs(r(eq))) / (one + abs(r(eq)))
  end do

end function limiter_vanleer
