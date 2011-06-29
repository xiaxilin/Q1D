!============================= limiter_ospre =================================80
!
! Takes # of equations and successive gradients
!
!=============================================================================80
pure function limiter_ospre(neq, r)

  use set_precision,   only : dp
  use set_constants,   only : one, onep5

  implicit none

  integer,                  intent(in) :: neq
  real(dp), dimension(neq), intent(in) :: r

  integer                  :: eq
  real(dp), dimension(neq) :: limiter_ospre

  do eq = 1, neq
    limiter_ospre(eq) = onep5*(r(eq)*r(eq) + r(eq)) / (r(eq)*r(eq)+ r(eq) + one)
  end do

end function limiter_ospre
