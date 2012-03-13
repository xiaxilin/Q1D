pure function floor_primitive_vars( prim )

  use set_precision, only : dp

  implicit none

  real(dp), dimension(3), intent(in) :: prim
  real(dp), dimension(3)             :: floor_primitive_vars

  floor_primitive_vars(1) = max( prim(1), 0.0001_dp )
  floor_primitive_vars(2) = prim(2)
  floor_primitive_vars(3) = max( prim(3), 500.0_dp )

end function floor_primitive_vars
