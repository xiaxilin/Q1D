!============================ fluid_constants ================================80
!
! Initialize constants dependent on ratio of specific heats ( gamma ) or
! the specific gas constant, R
! NOTE: x implies '1 over'
!
!=============================================================================80

module fluid_constants

  use set_precision, only : dp

  implicit none

  private

  public :: r     ! Gas constant
  public :: gamma ! Ratio of specific heats
  public :: gm1
  public :: gp1
  public :: xg
  public :: xgm1
  public :: gxgm1
  public :: gxg2m1

! Set initial values
  real(dp) :: r = 287.0_dp
  real(dp) :: gamma = 1.4_dp
  real(dp) :: gm1
  real(dp) :: gp1
  real(dp) :: xg
  real(dp) :: xgm1
  real(dp) :: gxgm1
  real(dp) :: gxg2m1

contains

!=========================== set_gamma_constants =============================80
!
!  Gamma is set by the namelist input, all other constants are defined here
!
!=============================================================================80
  subroutine set_gamma_constants()

    use set_constants, only : one

    continue

    gm1  = gamma - one
    gp1  = gamma + one
    xg   = one / gamma
    xgm1 = one / (gamma - one)
    gxgm1 = gamma / (gamma - one)
    gxg2m1 = gamma / (gamma*gamma - one)

  end subroutine set_gamma_constants

end module fluid_constants
