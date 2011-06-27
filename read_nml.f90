module namelist

  use set_precision, only : dp

  implicit none

  private

  public :: read_nml

  public :: iterations ! total number of iterations
  public :: firstorder ! if 2nd order, how many 1st order iters for stability?
  public :: itercheck  ! check for convergence every itercheck iters
  public :: rkorder    ! multistep RK order 1, 2, 3, or 4
  public :: cfl        ! CFL limit
  public :: muscl      ! MUSCL extrapolation, .true. or .false.
  public :: kappa      ! form of MUSCL
  public :: limiter    ! type of variable limiting 
  public :: toler      ! convergence tolerance
  public :: flux_type  ! flux type, '2nd', 'jst', 'sw', 'vanleer', 'roe', 'ausm'
  public :: k2         ! JST damping coefficient
  public :: k4         ! JST damping coefficient
  public :: mref       ! Initial mach number in nozzle
  public :: to         ! Inflow stag. temp
  public :: po         ! Inflow stag. pressure
  public :: pback      ! Outflow backpressure, negative for extrapolation
  public :: gamma      ! Ratio of specific heats
  public :: r          ! Gas constant

! Define namelist inputs

  integer           :: iterations, firstorder, itercheck, rkorder
  real(dp)          :: cfl, kappa, toler
  character(len=10) :: limiter
  logical           :: muscl
  namelist /code_control/ iterations, firstorder, itercheck, rkorder, cfl, &
                          limiter, muscl, kappa, toler

  character(len=10) :: flux_type
  real(dp)          :: k2, k4    ! only for flux_type = 'jst'
  namelist /flux/ flux_type, k2, k4

  real(dp) :: mref, to, po pback
  namelist /conditions/ mref, to, po, pback

  real(dp) :: gamma, r
  namelist /gas_properties/ gamma, r

contains

!================================= read_nml ==================================80
!
! Reads the q1d.nml file
!
!=============================================================================80
  subroutine read_nml

    use set_precision, only : dp

    integer :: nml_unit

    continue

    nml_unit = find_available_unit()

    open(nml_unit, file='q1d.nml', status='old')

! set defaults and read &code_control

    iterations = 100000
    firstorder = 10000
    itercheck  = 1000
    rkorder    = 1
    cfl        = 1.0
    muscl      = .false.
    kappa      = -1.0
    limiter    = 'none'
    toler      = 1.0e-13

    rewind(nml_unit)
    read(nml_unit, nml=code_control)

! set defaults and read &flux

    flux_type  = 'jst'
    k2         = 0.5
    k4         = 0.03125

    rewind(nml_unit)
    read(nml_unit, nml=flux)

! set defaults and read &conditions

    mref  = 1.5
    to    = 600.0
    po    = 300000.0
    pback = -1.0

    rewind(nml_unit)
    read(nml_unit, nml=conditions)

! set defaults and read &gas_properties

    gamma = 1.4
    r     = 287.0

    rewind(nml_unit)
    read(nml_unit, nml=gas_properties)

    close(nml_unit)

  end subroutine read_nml

end module namelist
