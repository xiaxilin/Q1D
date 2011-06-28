pure function speed_of_sound( p, rho )

  use set_precision,   only : dp
  use fluid_constants, only : gamma

  implicit none

  real(dp), intent(in) :: p, rho
  real(dp)             :: speed_of_sound

  speed_of_sound = sqrt(gamma*p/rho)

end function speed_of_sound
