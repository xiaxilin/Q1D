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
  public :: xgp1
  public :: gxgm1
  public :: xg2m1
  public :: gxg2m1
  public :: gm1xgp1
  public :: gp1xgm1

  public :: set_gamma_constants

! Set initial values
  real(dp) :: r
  real(dp) :: gamma
  real(dp) :: gm1
  real(dp) :: gp1
  real(dp) :: xg
  real(dp) :: xgm1
  real(dp) :: xgp1
  real(dp) :: gxgm1
  real(dp) :: xg2m1
  real(dp) :: gxg2m1
  real(dp) :: gm1xgp1
  real(dp) :: gp1xgm1

contains

!=========================== set_gamma_constants =============================80
!
!  Gamma is set by the namelist input, all other constants are defined here
!
!=============================================================================80
  subroutine set_gamma_constants()

    use set_constants, only : one

    continue

    gm1    = gamma - one
    gp1    = gamma + one
    xg     = one / gamma
    xgm1   = one / (gamma - one)
    xgp1   = one / (gamma + one)
    xg2m1  = one / (gamma*gamma - one)
    gxgm1  = gamma / (gamma - one)
    gxg2m1 = gamma / (gamma*gamma - one)
    gm1xgp1 = gm1/gp1
    gp1xgm1 = gp1/gm1

  end subroutine set_gamma_constants

end module fluid_constants
