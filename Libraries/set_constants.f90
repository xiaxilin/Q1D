!================================ constants ==================================80
!
! Provides useful constants
!
!=============================================================================80
module set_constants

  use set_precision, only : dp

  implicit none

  private

  public :: zero, one, onep5, two, three, four, five, six
  public :: sixth, fifth, fourth, third, half

! extremes
  public :: small, large

! transcendental numbers, etc...
  public :: pi, set_pi

  real(dp), parameter :: zero   = 0.0_dp
  real(dp), parameter :: one    = 1.0_dp
  real(dp), parameter :: onep5  = 1.5_dp
  real(dp), parameter :: two    = 2.0_dp
  real(dp), parameter :: three  = 3.0_dp
  real(dp), parameter :: four   = 4.0_dp
  real(dp), parameter :: five   = 5.0_dp
  real(dp), parameter :: six    = 6.0_dp

  real(dp), parameter :: sixth  = one/six
  real(dp), parameter :: fifth  = 0.2_dp
  real(dp), parameter :: fourth = 0.25_dp
  real(dp), parameter :: third  = one/three
  real(dp), parameter :: half   = 0.5_dp

  real(dp), parameter :: small  = tiny(1.0_dp)
  real(dp), parameter :: large  = huge(1.0_dp)

  real(dp) :: pi = 3.14159265358979_dp

contains

!==================================== set_pi =================================80
!
! Sets the value of PI by the acos intrinsic
!
!=============================================================================80
  subroutine set_pi
    use set_precision, only : dp

    implicit none

    pi = acos(-one)

  end subroutine set_pi

end module set_constants
