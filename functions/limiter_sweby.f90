!============================= limiter_sweby =================================80
!
! Takes # of equations, successive gradients, and beta parameter
! beta = one is minmod, beta = two is superbee
!
!=============================================================================80
pure function limiter_sweby(neq, r, beta)

  use set_precision,   only : dp
  use set_constants,   only : zero, one

  implicit none

  integer,                  intent(in) :: neq
  real(dp), dimension(neq), intent(in) :: r
  real(dp),                 intent(in) :: beta

  integer                  :: eq
  real(dp), dimension(neq) :: limiter_sweby

  do eq = 1, neq
    limiter_sweby(eq) = max( zero, min(beta*r(eq), one), min(r(eq), beta) )
  end do

end function limiter_sweby
