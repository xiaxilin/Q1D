pure function floor_primitive_vars( prim ) result(floored)

  use set_precision, only : dp

  implicit none

  real(dp), dimension(3), intent(in) :: prim
  real(dp), dimension(3)             :: floored

  floored(1) = max( prim(1), 0.0001_dp )
  floored(2) = prim(2)
  floored(3) = max( prim(3), 500.0_dp )

end function floor_primitive_vars
